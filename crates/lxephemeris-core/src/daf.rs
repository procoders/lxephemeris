//! Safe, pure-Rust NAIF DAF/SPK binary parser, validator, stitcher, and Chebyshev Type 2 builder.
//!
//! Conforms to the NASA NAIF Double Precision Array File (DAF) and Spacecraft & Planet Kernel (SPK)
//! specifications without any `unsafe` code or unchecked slice indexing.

use crate::error::EphemerisError;
use crate::types::{JD_J2000, SECONDS_PER_DAY};
use serde::{Deserialize, Serialize};

const DAF_RECORD_BYTES: usize = 1024;
const DAF_RECORD_WORDS: usize = 128;
const MAX_SUMMARIES_PER_RECORD: usize = 25;
const SUMMARY_BYTES: usize = 40; // 2 f64 (16 bytes) + 6 i32 (24 bytes)
const NAME_SLOT_BYTES: usize = 40;

/// Parsed segment descriptor from a NAIF `DAF/SPK` binary file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DafSegmentDescriptor {
    pub name: String,
    pub target_id: i32,
    pub center_id: i32,
    pub frame_id: i32,
    pub spk_type: i32,
    pub start_et_s: f64,
    pub end_et_s: f64,
    pub start_jd_tdb: f64,
    pub end_jd_tdb: f64,
    /// 1-based inclusive starting `f64` word address in the DAF file.
    pub initial_word: usize,
    /// 1-based inclusive ending `f64` word address in the DAF file.
    pub final_word: usize,
}

/// Parsed structure of a NAIF `DAF/SPK` file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DafSpkFile {
    pub internal_name: String,
    pub little_endian: bool,
    pub fward: usize,
    pub bward: usize,
    pub free_word: usize,
    pub comments: String,
    pub segments: Vec<DafSegmentDescriptor>,
}

/// Input segment for building or stitching a NAIF `DAF/SPK` binary file.
#[derive(Debug, Clone)]
pub struct DafSegmentInput {
    pub name: String,
    pub target_id: i32,
    pub center_id: i32,
    pub frame_id: i32,
    pub spk_type: i32,
    pub start_et_s: f64,
    pub end_et_s: f64,
    pub data_words: Vec<f64>,
}

fn read_u32_at(bytes: &[u8], offset: usize, little_endian: bool) -> Option<u32> {
    let slice = bytes.get(offset..offset.checked_add(4)?)?;
    let arr: [u8; 4] = slice.try_into().ok()?;
    Some(if little_endian {
        u32::from_le_bytes(arr)
    } else {
        u32::from_be_bytes(arr)
    })
}

fn read_i32_at(bytes: &[u8], offset: usize, little_endian: bool) -> Option<i32> {
    read_u32_at(bytes, offset, little_endian).map(|v| v as i32)
}

fn read_f64_at(bytes: &[u8], offset: usize, little_endian: bool) -> Option<f64> {
    let slice = bytes.get(offset..offset.checked_add(8)?)?;
    let arr: [u8; 8] = slice.try_into().ok()?;
    Some(if little_endian {
        f64::from_le_bytes(arr)
    } else {
        f64::from_be_bytes(arr)
    })
}

/// Parses and validates the header, comment records, and segment descriptors of a `DAF/SPK` binary buffer.
pub fn parse_daf_spk(bytes: &[u8], path_label: &str) -> Result<DafSpkFile, EphemerisError> {
    if bytes.len() < DAF_RECORD_BYTES {
        return Err(EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: format!(
                "File size {} bytes is smaller than 1024-byte DAF record",
                bytes.len()
            ),
        });
    }

    let magic = bytes
        .get(0..8)
        .ok_or_else(|| EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: "Missing DAF magic header".to_owned(),
        })?;
    if magic != b"DAF/SPK " && magic != b"NAIF/DAF" {
        return Err(EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: format!(
                "Invalid DAF/SPK magic header: {:?}",
                String::from_utf8_lossy(magic)
            ),
        });
    }

    // Determine endianness from LOCFMT at offset 88..96 or ND==2 check at offset 8..12.
    let locfmt = bytes.get(88..96).unwrap_or(b"        ");
    let little_endian = if locfmt.starts_with(b"LTL-IEEE") {
        true
    } else if locfmt.starts_with(b"BIG-IEEE") {
        false
    } else {
        let nd_le = read_i32_at(bytes, 8, true).unwrap_or(0);
        nd_le == 2
    };

    let nd = read_i32_at(bytes, 8, little_endian).unwrap_or(0);
    let ni = read_i32_at(bytes, 12, little_endian).unwrap_or(0);
    if nd != 2 || ni != 6 {
        return Err(EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: format!("Invalid SPK summary parameters ND={nd}, NI={ni} (expected ND=2, NI=6)"),
        });
    }

    let internal_name = bytes
        .get(16..76)
        .map(|s| String::from_utf8_lossy(s).trim().to_owned())
        .unwrap_or_default();

    let fward_i32 = read_i32_at(bytes, 76, little_endian).unwrap_or(0);
    let bward_i32 = read_i32_at(bytes, 80, little_endian).unwrap_or(0);
    let free_i32 = read_i32_at(bytes, 84, little_endian).unwrap_or(0);

    if fward_i32 < 2 || bward_i32 < fward_i32 || free_i32 < 1 {
        return Err(EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: format!(
                "Invalid DAF record pointers FWARD={fward_i32}, BWARD={bward_i32}, FREE={free_i32}"
            ),
        });
    }

    let fward = fward_i32 as usize;
    let bward = bward_i32 as usize;
    let free_word = free_i32 as usize;

    // Extract DAF comment records (records 2 .. fward - 1).
    let mut comments = String::new();
    if fward > 2 {
        let mut comment_bytes = Vec::new();
        'comment_loop: for rec_idx in 2..fward {
            let rec_offset = (rec_idx - 1).saturating_mul(DAF_RECORD_BYTES);
            // Each comment record holds up to 1000 characters.
            let rec_slice =
                bytes
                    .get(rec_offset..rec_offset.saturating_add(1000))
                    .ok_or_else(|| EphemerisError::CorruptKernel {
                        path: path_label.to_owned(),
                        reason: format!("Truncated DAF comment record {rec_idx}"),
                    })?;
            for &b in rec_slice {
                if b == 4 {
                    // ASCII EOT marks end of DAF comment area.
                    break 'comment_loop;
                } else if b == 0 {
                    // Null byte represents newline in DAF comment records.
                    comment_bytes.push(b'\n');
                } else {
                    comment_bytes.push(b);
                }
            }
        }
        comments = String::from_utf8_lossy(&comment_bytes).into_owned();
    }

    // Traverse summary records starting at FWARD.
    let mut segments = Vec::new();
    let mut curr_rec = fward;
    let mut visited = 0usize;
    let max_records = bytes.len() / DAF_RECORD_BYTES + 1;

    while curr_rec != 0 {
        visited += 1;
        if visited > max_records {
            return Err(EphemerisError::CorruptKernel {
                path: path_label.to_owned(),
                reason: "Cycle detected in DAF summary record chain".to_owned(),
            });
        }

        let sum_offset = (curr_rec - 1).saturating_mul(DAF_RECORD_BYTES);
        let sum_record = bytes
            .get(sum_offset..sum_offset.saturating_add(DAF_RECORD_BYTES))
            .ok_or_else(|| EphemerisError::CorruptKernel {
                path: path_label.to_owned(),
                reason: format!("Summary record {curr_rec} extends past end of file"),
            })?;

        let name_offset = curr_rec.saturating_mul(DAF_RECORD_BYTES);
        let name_record = bytes
            .get(name_offset..name_offset.saturating_add(DAF_RECORD_BYTES))
            .unwrap_or(&[]);

        let next_f64 = read_f64_at(sum_record, 0, little_endian).unwrap_or(0.0);
        let nsum_f64 = read_f64_at(sum_record, 16, little_endian).unwrap_or(0.0);

        if !(0.0..=(MAX_SUMMARIES_PER_RECORD as f64)).contains(&nsum_f64) {
            return Err(EphemerisError::CorruptKernel {
                path: path_label.to_owned(),
                reason: format!("Invalid NSUM={nsum_f64} in summary record {curr_rec}"),
            });
        }

        let nsum = nsum_f64 as usize;
        for idx in 0..nsum {
            let entry_offset = 24 + idx.saturating_mul(SUMMARY_BYTES);
            let start_et_s =
                read_f64_at(sum_record, entry_offset, little_endian).ok_or_else(|| {
                    EphemerisError::CorruptKernel {
                        path: path_label.to_owned(),
                        reason: format!("Truncated summary entry {idx} in record {curr_rec}"),
                    }
                })?;
            let end_et_s =
                read_f64_at(sum_record, entry_offset + 8, little_endian).ok_or_else(|| {
                    EphemerisError::CorruptKernel {
                        path: path_label.to_owned(),
                        reason: format!("Truncated summary entry {idx} in record {curr_rec}"),
                    }
                })?;
            let target_id =
                read_i32_at(sum_record, entry_offset + 16, little_endian).unwrap_or(0);
            let center_id =
                read_i32_at(sum_record, entry_offset + 20, little_endian).unwrap_or(0);
            let frame_id = read_i32_at(sum_record, entry_offset + 24, little_endian).unwrap_or(0);
            let spk_type = read_i32_at(sum_record, entry_offset + 28, little_endian).unwrap_or(0);
            let initial_word_i32 =
                read_i32_at(sum_record, entry_offset + 32, little_endian).unwrap_or(0);
            let final_word_i32 =
                read_i32_at(sum_record, entry_offset + 36, little_endian).unwrap_or(0);

            if !start_et_s.is_finite()
                || !end_et_s.is_finite()
                || start_et_s >= end_et_s
                || initial_word_i32 < 1
                || final_word_i32 < initial_word_i32
            {
                return Err(EphemerisError::CorruptKernel {
                    path: path_label.to_owned(),
                    reason: format!(
                        "Invalid segment bounds for target {target_id}: et=[{start_et_s}, {end_et_s}], words=[{initial_word_i32}, {final_word_i32}]"
                    ),
                });
            }

            let initial_word = initial_word_i32 as usize;
            let final_word = final_word_i32 as usize;
            let end_byte = final_word.saturating_mul(8);
            if end_byte > bytes.len() {
                return Err(EphemerisError::CorruptKernel {
                    path: path_label.to_owned(),
                    reason: format!(
                        "Segment for target {target_id} references word {final_word} (byte {end_byte}), exceeding file size {}",
                        bytes.len()
                    ),
                });
            }

            let name_start = idx.saturating_mul(NAME_SLOT_BYTES);
            let name = name_record
                .get(name_start..name_start.saturating_add(NAME_SLOT_BYTES))
                .map(|s| String::from_utf8_lossy(s).trim().to_owned())
                .unwrap_or_default();

            let start_jd_tdb = JD_J2000 + start_et_s / SECONDS_PER_DAY;
            let end_jd_tdb = JD_J2000 + end_et_s / SECONDS_PER_DAY;

            segments.push(DafSegmentDescriptor {
                name,
                target_id,
                center_id,
                frame_id,
                spk_type,
                start_et_s,
                end_et_s,
                start_jd_tdb,
                end_jd_tdb,
                initial_word,
                final_word,
            });
        }

        curr_rec = if next_f64 <= 0.0 {
            0
        } else {
            next_f64 as usize
        };
    }

    Ok(DafSpkFile {
        internal_name,
        little_endian,
        fward,
        bward,
        free_word,
        comments,
        segments,
    })
}

/// Extracts the raw `f64` data words for a given [`DafSegmentDescriptor`] from a DAF/SPK byte slice.
pub fn extract_segment_words(
    bytes: &[u8],
    seg: &DafSegmentDescriptor,
    little_endian: bool,
    path_label: &str,
) -> Result<Vec<f64>, EphemerisError> {
    if seg.initial_word < 1 || seg.final_word < seg.initial_word {
        return Err(EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: "Invalid segment word bounds".to_owned(),
        });
    }
    let count = seg.final_word - seg.initial_word + 1;
    let start_byte = (seg.initial_word - 1).saturating_mul(8);
    let end_byte = seg.final_word.saturating_mul(8);
    let raw = bytes
        .get(start_byte..end_byte)
        .ok_or_else(|| EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: "Segment data out of bounds".to_owned(),
        })?;

    let mut words = Vec::with_capacity(count);
    for chunk in raw.chunks_exact(8) {
        let arr: [u8; 8] = chunk.try_into().map_err(|_| EphemerisError::CorruptKernel {
            path: path_label.to_owned(),
            reason: "Unaligned f64 chunk".to_owned(),
        })?;
        words.push(if little_endian {
            f64::from_le_bytes(arr)
        } else {
            f64::from_be_bytes(arr)
        });
    }
    Ok(words)
}

/// Builds a complete, valid little-endian NAIF `DAF/SPK` binary file from a list of [`DafSegmentInput`]s.
///
/// Packs up to 25 segment summaries per Summary Record, pairs each Summary Record with a Name Record,
/// assigns exact 1-based `f64` word addresses, and pads the output to 1024-byte record boundaries.
pub fn build_daf_spk_bytes(
    internal_name: &str,
    comment_text: &str,
    segments: &[DafSegmentInput],
) -> Result<Vec<u8>, EphemerisError> {
    if segments.is_empty() {
        return Err(EphemerisError::InvalidInput {
            message: "Cannot build a DAF/SPK file with zero segments".to_owned(),
        });
    }

    // Encode comment records (if any). Each 1024-byte comment record stores up to 1000 bytes.
    let mut comment_records: Vec<[u8; DAF_RECORD_BYTES]> = Vec::new();
    if !comment_text.is_empty() {
        let mut encoded = Vec::with_capacity(comment_text.len() + 2);
        for line in comment_text.lines() {
            encoded.extend_from_slice(line.as_bytes());
            encoded.push(0u8); // DAF comment newline marker
        }
        encoded.push(4u8); // DAF EOT marker

        for chunk in encoded.chunks(1000) {
            let mut rec = [b' '; DAF_RECORD_BYTES];
            if let Some(dest) = rec.get_mut(0..chunk.len()) {
                dest.copy_from_slice(chunk);
            }
            // Bytes 1000..1024 are zeroed per NAIF DAF convention.
            if let Some(tail) = rec.get_mut(1000..DAF_RECORD_BYTES) {
                tail.fill(0);
            }
            comment_records.push(rec);
        }
    }

    let fward: usize = 2 + comment_records.len();
    let num_summary_blocks = segments.len().div_ceil(MAX_SUMMARIES_PER_RECORD);
    // Each summary block takes 2 records: 1 Summary Record + 1 Name Record.
    // The last Summary Record is at `fward + 2 * (num_summary_blocks - 1)`.
    let bward: usize = fward + 2 * (num_summary_blocks - 1);
    let first_data_record: usize = fward + 2 * num_summary_blocks;
    // 1-based f64 word index where the first segment's data begins:
    let mut next_free_word: usize = (first_data_record - 1) * DAF_RECORD_WORDS + 1;

    // Precompute segment word ranges.
    let mut word_ranges: Vec<(usize, usize)> = Vec::with_capacity(segments.len());
    for seg in segments {
        if seg.data_words.is_empty() || seg.start_et_s >= seg.end_et_s {
            return Err(EphemerisError::InvalidInput {
                message: format!(
                    "Invalid segment input for target {}: words={}, et=[{}, {}]",
                    seg.target_id,
                    seg.data_words.len(),
                    seg.start_et_s,
                    seg.end_et_s
                ),
            });
        }
        let init_w = next_free_word;
        let final_w = init_w + seg.data_words.len() - 1;
        word_ranges.push((init_w, final_w));
        next_free_word = final_w + 1;
    }

    // Build Record 1: File Record.
    let mut file_record = [0u8; DAF_RECORD_BYTES];
    if let Some(s) = file_record.get_mut(0..8) {
        s.copy_from_slice(b"DAF/SPK ");
    }
    if let Some(s) = file_record.get_mut(8..12) {
        s.copy_from_slice(&2_i32.to_le_bytes()); // ND = 2
    }
    if let Some(s) = file_record.get_mut(12..16) {
        s.copy_from_slice(&6_i32.to_le_bytes()); // NI = 6
    }
    if let Some(s) = file_record.get_mut(16..76) {
        s.fill(b' ');
        let name_bytes = internal_name.as_bytes();
        let copy_len = name_bytes.len().min(60);
        if let (Some(dst), Some(src)) = (s.get_mut(0..copy_len), name_bytes.get(0..copy_len)) {
            dst.copy_from_slice(src);
        }
    }
    if let Some(s) = file_record.get_mut(76..80) {
        s.copy_from_slice(&(fward as i32).to_le_bytes());
    }
    if let Some(s) = file_record.get_mut(80..84) {
        s.copy_from_slice(&(bward as i32).to_le_bytes());
    }
    if let Some(s) = file_record.get_mut(84..88) {
        s.copy_from_slice(&(next_free_word as i32).to_le_bytes());
    }
    if let Some(s) = file_record.get_mut(88..96) {
        s.copy_from_slice(b"LTL-IEEE");
    }
    // Standard NAIF FTP validation string at offset 699..727
    let ftpstr = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP";
    if let Some(s) = file_record.get_mut(699..699 + ftpstr.len()) {
        s.copy_from_slice(ftpstr);
    }

    let mut out = Vec::new();
    out.extend_from_slice(&file_record);
    for crec in &comment_records {
        out.extend_from_slice(crec);
    }

    // Write Summary + Name record pairs.
    for block_idx in 0..num_summary_blocks {
        let curr_sum_rec = fward + 2 * block_idx;
        let prev_sum_rec = if block_idx == 0 {
            0.0_f64
        } else {
            (curr_sum_rec - 2) as f64
        };
        let next_sum_rec = if block_idx + 1 < num_summary_blocks {
            (curr_sum_rec + 2) as f64
        } else {
            0.0_f64
        };

        let start_seg = block_idx * MAX_SUMMARIES_PER_RECORD;
        let end_seg = (start_seg + MAX_SUMMARIES_PER_RECORD).min(segments.len());
        let nsum = (end_seg - start_seg) as f64;

        let mut sum_rec = [0u8; DAF_RECORD_BYTES];
        let mut name_rec = [b' '; DAF_RECORD_BYTES];

        if let Some(s) = sum_rec.get_mut(0..8) {
            s.copy_from_slice(&next_sum_rec.to_le_bytes());
        }
        if let Some(s) = sum_rec.get_mut(8..16) {
            s.copy_from_slice(&prev_sum_rec.to_le_bytes());
        }
        if let Some(s) = sum_rec.get_mut(16..24) {
            s.copy_from_slice(&nsum.to_le_bytes());
        }

        for (local_idx, seg_idx) in (start_seg..end_seg).enumerate() {
            let seg = segments
                .get(seg_idx)
                .ok_or_else(|| EphemerisError::InvalidInput {
                    message: "Segment index out of bounds".to_owned(),
                })?;
            let (init_w, final_w) =
                word_ranges
                    .get(seg_idx)
                    .copied()
                    .ok_or_else(|| EphemerisError::InvalidInput {
                        message: "Word range index out of bounds".to_owned(),
                    })?;

            let off = 24 + local_idx * SUMMARY_BYTES;
            if let Some(s) = sum_rec.get_mut(off..off + 8) {
                s.copy_from_slice(&seg.start_et_s.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 8..off + 16) {
                s.copy_from_slice(&seg.end_et_s.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 16..off + 20) {
                s.copy_from_slice(&seg.target_id.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 20..off + 24) {
                s.copy_from_slice(&seg.center_id.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 24..off + 28) {
                s.copy_from_slice(&seg.frame_id.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 28..off + 32) {
                s.copy_from_slice(&seg.spk_type.to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 32..off + 36) {
                s.copy_from_slice(&(init_w as i32).to_le_bytes());
            }
            if let Some(s) = sum_rec.get_mut(off + 36..off + 40) {
                s.copy_from_slice(&(final_w as i32).to_le_bytes());
            }

            let name_off = local_idx * NAME_SLOT_BYTES;
            let name_bytes = seg.name.as_bytes();
            let copy_len = name_bytes.len().min(NAME_SLOT_BYTES);
            if let (Some(dst), Some(src)) = (
                name_rec.get_mut(name_off..name_off + copy_len),
                name_bytes.get(0..copy_len),
            ) {
                dst.copy_from_slice(src);
            }
        }

        out.extend_from_slice(&sum_rec);
        out.extend_from_slice(&name_rec);
    }

    // Write all data words contiguously.
    for seg in segments {
        for &w in &seg.data_words {
            out.extend_from_slice(&w.to_le_bytes());
        }
    }

    // Pad to 1024-byte DAF record boundary.
    let remainder = out.len() % DAF_RECORD_BYTES;
    if remainder != 0 {
        out.resize(out.len() + (DAF_RECORD_BYTES - remainder), 0u8);
    }

    Ok(out)
}

/// Evaluates a NASA NAIF `SPK Type 21` (`Extended Modified Difference Array`, `SPKE21`) segment
/// at `et_s` (ephemeris time in seconds past J2000 TDB), returning `(pos_km, vel_km_s)` relative
/// to the segment's `center_id`.
pub fn eval_spk_type21(
    words: &[f64],
    et_s: f64,
    target_id: i32,
) -> Result<([f64; 3], [f64; 3]), EphemerisError> {
    if !et_s.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Non-finite ephemeris time et_s={et_s} in SPK Type 21 evaluation"),
        });
    }
    if words.len() < 4 {
        return Err(EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: "SPK Type 21 segment has fewer than 4 words".to_owned(),
        });
    }
    let n_rec_f64 = words.last().copied().unwrap_or(0.0);
    let maxdim_f64 = words
        .get(words.len().saturating_sub(2))
        .copied()
        .unwrap_or(0.0);
    if !n_rec_f64.is_finite()
        || !maxdim_f64.is_finite()
        || n_rec_f64.fract() != 0.0
        || maxdim_f64.fract() != 0.0
        || !(1.0..=1_000_000.0).contains(&n_rec_f64)
        || !(1.0..=100.0).contains(&maxdim_f64)
    {
        return Err(EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: format!("Invalid SPK Type 21 metadata N={n_rec_f64}, MAXDIM={maxdim_f64}"),
        });
    }
    let n_rec = n_rec_f64 as usize;
    let maxdim = maxdim_f64 as usize;
    let dlsize = 4 * maxdim + 11;
    let epochs_start = n_rec.saturating_mul(dlsize);
    let epochs_end = epochs_start.saturating_add(n_rec);
    let epochs = words
        .get(epochs_start..epochs_end)
        .ok_or_else(|| EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: "SPK Type 21 epoch table out of bounds".to_owned(),
        })?;

    // Binary search for the first epoch >= et_s.
    let idx = epochs
        .partition_point(|&ep| ep < et_s)
        .min(n_rec.saturating_sub(1));
    let rec_start = idx.saturating_mul(dlsize);
    let rec = words
        .get(rec_start..rec_start.saturating_add(dlsize))
        .ok_or_else(|| EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: "SPK Type 21 MDA record out of bounds".to_owned(),
        })?;

    let tl = rec.first().copied().unwrap_or(0.0);
    let g = rec.get(1..1 + maxdim).unwrap_or(&[]);
    let refpos = [
        rec.get(1 + maxdim).copied().unwrap_or(0.0),
        rec.get(1 + maxdim + 2).copied().unwrap_or(0.0),
        rec.get(1 + maxdim + 4).copied().unwrap_or(0.0),
    ];
    let refvel = [
        rec.get(1 + maxdim + 1).copied().unwrap_or(0.0),
        rec.get(1 + maxdim + 3).copied().unwrap_or(0.0),
        rec.get(1 + maxdim + 5).copied().unwrap_or(0.0),
    ];
    let dt0 = rec.get(1 + maxdim + 6..1 + 2 * maxdim + 6).unwrap_or(&[]);
    let dt1 = rec
        .get(1 + 2 * maxdim + 6..1 + 3 * maxdim + 6)
        .unwrap_or(&[]);
    let dt2 = rec
        .get(1 + 3 * maxdim + 6..1 + 4 * maxdim + 6)
        .unwrap_or(&[]);
    let dt = [dt0, dt1, dt2];

    let kqmax1_raw = rec.get(4 * maxdim + 7).copied().unwrap_or(0.0);
    if !kqmax1_raw.is_finite()
        || kqmax1_raw.fract() != 0.0
        || !(2.0..=(maxdim + 2) as f64).contains(&kqmax1_raw)
    {
        return Err(EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: format!("Invalid SPK Type 21 KQMAX1={kqmax1_raw} (MAXDIM={maxdim})"),
        });
    }
    let kqmax1 = kqmax1_raw as usize;

    let mut kq = [0usize; 3];
    for (axis, slot) in kq.iter_mut().enumerate() {
        let kq_raw = rec.get(4 * maxdim + 8 + axis).copied().unwrap_or(0.0);
        if !kq_raw.is_finite()
            || kq_raw.fract() != 0.0
            || !(1.0..=maxdim as f64).contains(&kq_raw)
            || (kq_raw as usize) >= kqmax1
        {
            return Err(EphemerisError::CorruptKernel {
                path: format!("NAIF:{target_id}"),
                reason: format!(
                    "Invalid SPK Type 21 KQ[{axis}]={kq_raw} (KQMAX1={kqmax1}, MAXDIM={maxdim})"
                ),
            });
        }
        *slot = kq_raw as usize;
    }
    if !tl.is_finite() {
        return Err(EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}"),
            reason: "Non-finite reference epoch TL in SPK Type 21 record".to_owned(),
        });
    }

    let delta = et_s - tl;
    let mut tp = delta;
    let mpq = kqmax1 - 2;
    let mut ks = kqmax1 - 1;

    let mut fc = vec![0.0_f64; maxdim + 4];
    let mut wc = vec![0.0_f64; maxdim + 4];
    let mut w = vec![0.0_f64; maxdim + 6];
    if let Some(f1) = fc.get_mut(1) {
        *f1 = 1.0;
    }

    for j in 1..=mpq {
        let gj = g.get(j - 1).copied().unwrap_or(0.0);
        if gj.abs() == 0.0 {
            return Err(EphemerisError::CorruptKernel {
                path: format!("NAIF:{target_id}"),
                reason: "Zero stepsize G(J) in SPK Type 21 record".to_owned(),
            });
        }
        if let Some(slot) = fc.get_mut(j + 1) {
            *slot = tp / gj;
        }
        if let Some(slot) = wc.get_mut(j) {
            *slot = delta / gj;
        }
        tp = delta + gj;
    }

    for j in 1..=kqmax1 {
        if let Some(slot) = w.get_mut(j) {
            *slot = 1.0 / (j as f64);
        }
    }

    let mut jx = 0usize;
    let mut ks1 = ks - 1;
    while ks >= 2 {
        jx += 1;
        for j in 1..=jx {
            let f_val = fc.get(j + 1).copied().unwrap_or(0.0);
            let w_ks1 = w.get(j + ks1).copied().unwrap_or(0.0);
            let wc_val = wc.get(j).copied().unwrap_or(0.0);
            let w_ks = w.get(j + ks).copied().unwrap_or(0.0);
            if let Some(slot) = w.get_mut(j + ks) {
                *slot = f_val * w_ks1 - wc_val * w_ks;
            }
        }
        ks = ks1;
        ks1 = ks1.saturating_sub(1);
    }

    let mut pos = [0.0_f64; 3];
    for i in 0..3 {
        let kqq = *kq.get(i).unwrap_or(&0);
        let dti = *dt.get(i).unwrap_or(&(&[][..]));
        let mut sum = 0.0_f64;
        for j in (1..=kqq).rev() {
            let d_val = dti.get(j - 1).copied().unwrap_or(0.0);
            let w_val = w.get(j + ks).copied().unwrap_or(0.0);
            sum += d_val * w_val;
        }
        let rp = *refpos.get(i).unwrap_or(&0.0);
        let rv = *refvel.get(i).unwrap_or(&0.0);
        if let Some(p_slot) = pos.get_mut(i) {
            *p_slot = rp + delta * (rv + delta * sum);
        }
    }

    for j in 1..=jx {
        let f_val = fc.get(j + 1).copied().unwrap_or(0.0);
        let w_ks1 = w.get(j + ks1).copied().unwrap_or(0.0);
        let wc_val = wc.get(j).copied().unwrap_or(0.0);
        let w_ks = w.get(j + ks).copied().unwrap_or(0.0);
        if let Some(slot) = w.get_mut(j + ks) {
            *slot = f_val * w_ks1 - wc_val * w_ks;
        }
    }
    ks = ks.saturating_sub(1);

    let mut vel = [0.0_f64; 3];
    for i in 0..3 {
        let kqq = *kq.get(i).unwrap_or(&0);
        let dti = *dt.get(i).unwrap_or(&(&[][..]));
        let mut sum = 0.0_f64;
        for j in (1..=kqq).rev() {
            let d_val = dti.get(j - 1).copied().unwrap_or(0.0);
            let w_val = w.get(j + ks).copied().unwrap_or(0.0);
            sum += d_val * w_val;
        }
        let rv = *refvel.get(i).unwrap_or(&0.0);
        if let Some(v_slot) = vel.get_mut(i) {
            *v_slot = rv + delta * sum;
        }
    }

    Ok((pos, vel))
}

/// Evaluates a Chebyshev polynomial of degree `degree` with coefficients `coeffs[0..=degree]`
/// at normalized coordinate `tau` in `[-1, 1]` with half-interval `radius_s` (seconds),
/// returning `(value, derivative_per_sec)`.
pub fn eval_chebyshev_1d(coeffs: &[f64], tau: f64, radius_s: f64) -> (f64, f64) {
    let degree = coeffs.len().saturating_sub(1);
    let mut w = [0.0_f64; 3];
    let mut dw = [0.0_f64; 3];

    for j in (2..=degree + 1).rev() {
        w[2] = w[1];
        w[1] = w[0];
        let cj = coeffs.get(j - 1).copied().unwrap_or(0.0);
        w[0] = cj + (2.0 * tau * w[1] - w[2]);

        dw[2] = dw[1];
        dw[1] = dw[0];
        dw[0] = w[1] * 2.0 + dw[1] * 2.0 * tau - dw[2];
    }

    let c0 = coeffs.first().copied().unwrap_or(0.0);
    let val = c0 + (tau * w[0] - w[1]);
    let deriv = if radius_s.abs() > 0.0 {
        (w[0] + tau * dw[0] - dw[1]) / radius_s
    } else {
        0.0
    };
    (val, deriv)
}

/// Fits a single continuous `SPK Type 2` (Chebyshev polynomial position) segment over `[start_et_s, end_et_s]`
/// from sorted `SPK Type 21` segments `(seg_start_et, seg_end_et, seg_words)`, verifying that the maximum
/// position residual is `< 0.5 m` and velocity residual is `< 0.01 mm/s` (`PLAN.md` §0.1).
///
/// Returns `(DafSegmentInput, max_pos_err_m, max_vel_err_mm_s)`.
pub fn fit_type21_to_type2_segment(
    name: &str,
    target_id: i32,
    center_id: i32,
    frame_id: i32,
    type21_segments: &[(f64, f64, Vec<f64>)],
) -> Result<(DafSegmentInput, f64, f64), EphemerisError> {
    if type21_segments.is_empty() {
        return Err(EphemerisError::MissingData {
            detail: format!("No Type 21 segments provided for target {target_id}"),
        });
    }
    let start_et_s = type21_segments.first().map(|s| s.0).unwrap_or(0.0);
    let end_et_s = type21_segments.last().map(|s| s.1).unwrap_or(0.0);
    let span_s = end_et_s - start_et_s;
    if span_s <= 0.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid Type 21 span for target {target_id}: [{start_et_s}, {end_et_s}]"),
        });
    }

    let eval_t21 = |et: f64| -> Result<([f64; 3], [f64; 3]), EphemerisError> {
        let clamped = et.clamp(start_et_s, end_et_s);
        let seg_idx = type21_segments
            .partition_point(|s| s.1 < clamped)
            .min(type21_segments.len() - 1);
        let seg = type21_segments
            .get(seg_idx)
            .ok_or_else(|| EphemerisError::MissingData {
                detail: format!("Missing Type 21 segment for target {target_id}"),
            })?;
        eval_spk_type21(&seg.2, clamped, target_id)
    };

    // Estimate minimum heliocentric distance over a coarse scan to select optimal interval length.
    let mut min_r_km = f64::MAX;
    let coarse_steps = 500usize;
    for k in 0..=coarse_steps {
        let et = start_et_s + span_s * (k as f64) / (coarse_steps as f64);
        if let Ok((pos, _)) = eval_t21(et) {
            let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
            if r < min_r_km {
                min_r_km = r;
            }
        }
    }
    let min_r_au = min_r_km / 149_597_870.700;

    // Candidate ladder of (interval_days, degree) from most compact to finest.
    let candidates: &[(f64, usize)] = if min_r_au >= 8.0 {
        &[(64.0, 13), (32.0, 13), (16.0, 13), (8.0, 13)]
    } else if min_r_au >= 1.7 {
        &[(32.0, 13), (16.0, 13), (8.0, 13), (4.0, 15)]
    } else if min_r_au >= 0.9 {
        &[(8.0, 13), (4.0, 13), (4.0, 15), (2.0, 15), (1.0, 15)]
    } else {
        &[(4.0, 15), (1.0, 15), (2.0, 15), (4.0, 13), (2.0, 13)]
    };

    let mut chosen = candidates.first().copied().unwrap_or((16.0, 13));
    let mut best_pos_err_m = f64::MAX;
    let mut best_vel_err_mm_s = f64::MAX;
    let mut best_score = f64::MAX;

    for &(intlen_days, degree) in candidates {
        let m = degree + 1;
        let intlen_s = intlen_days * SECONDS_PER_DAY;
        let radius_s = 0.5 * intlen_s;
        let num_records = (span_s / intlen_s).ceil() as usize;
        let stride = (num_records / 80).max(1);

        let mut nodes = Vec::with_capacity(m);
        for k in 0..m {
            let angle = (2.0 * (k as f64) + 1.0) * std::f64::consts::PI / (2.0 * (m as f64));
            nodes.push(angle.cos());
        }

        let mut max_pos_m = 0.0_f64;
        let mut max_vel_mm_s = 0.0_f64;

        for rec_idx in (0..num_records).step_by(stride) {
            let t0 = start_et_s + (rec_idx as f64) * intlen_s;
            let mid = t0 + radius_s;
            if mid + radius_s > end_et_s {
                break;
            }
            let mut samples = Vec::with_capacity(m);
            for &tau in &nodes {
                let (pos, _) = eval_t21(mid + radius_s * tau)?;
                samples.push(pos);
            }

            let mut coeffs = [vec![0.0; m], vec![0.0; m], vec![0.0; m]];
            for j in 0..m {
                let scale = if j == 0 {
                    1.0 / (m as f64)
                } else {
                    2.0 / (m as f64)
                };
                for (k, sample) in samples.iter().enumerate() {
                    let c_jk = ((j as f64)
                        * (2.0 * (k as f64) + 1.0)
                        * std::f64::consts::PI
                        / (2.0 * (m as f64)))
                    .cos();
                    for axis in 0..3 {
                        if let (Some(axis_coeffs), Some(&s_val)) =
                            (coeffs.get_mut(axis), sample.get(axis))
                        {
                            if let Some(slot) = axis_coeffs.get_mut(j) {
                                *slot += scale * s_val * c_jk;
                            }
                        }
                    }
                }
            }

            for frac in [-0.85_f64, -0.4, 0.0, 0.4, 0.85] {
                let et = mid + radius_s * frac;
                let (p_true, v_true) = eval_t21(et)?;
                let mut pos_sq = 0.0_f64;
                let mut vel_sq = 0.0_f64;
                for axis in 0..3 {
                    let (p_fit, v_fit) =
                        eval_chebyshev_1d(coeffs.get(axis).map_or(&[], Vec::as_slice), frac, radius_s);
                    let dp = p_true.get(axis).copied().unwrap_or(0.0) - p_fit;
                    let dv = v_true.get(axis).copied().unwrap_or(0.0) - v_fit;
                    pos_sq += dp * dp;
                    vel_sq += dv * dv;
                }
                let pos_err_m = pos_sq.sqrt() * 1000.0;
                let vel_err_mm_s = vel_sq.sqrt() * 1e6;
                max_pos_m = max_pos_m.max(pos_err_m);
                max_vel_mm_s = max_vel_mm_s.max(vel_err_mm_s);
            }
        }

        let score = (max_pos_m / 0.5).max(max_vel_mm_s / 0.01);
        if score < best_score {
            best_score = score;
            chosen = (intlen_days, degree);
            best_pos_err_m = max_pos_m;
            best_vel_err_mm_s = max_vel_mm_s;
        }
        if max_pos_m < 0.5 && max_vel_mm_s < 0.01 {
            chosen = (intlen_days, degree);
            best_pos_err_m = max_pos_m;
            best_vel_err_mm_s = max_vel_mm_s;
            break;
        }
    }

    // Build the full Type 2 segment using `chosen`, dividing `[start_et_s, end_et_s]` into
    // `num_records` equal intervals so the exact `[start_et_s, end_et_s]` span is preserved.
    let (intlen_days, degree) = chosen;
    let m = degree + 1;
    let nominal_intlen_s = intlen_days * SECONDS_PER_DAY;
    let num_records = ((end_et_s - start_et_s) / nominal_intlen_s).ceil().max(1.0) as usize;
    let intlen_s = (end_et_s - start_et_s) / (num_records as f64);
    let radius_s = 0.5 * intlen_s;
    let effective_end_et_s = end_et_s;
    let rsize = 2 + 3 * m;

    let mut nodes = Vec::with_capacity(m);
    let mut cos_table = vec![vec![0.0_f64; m]; m];
    for k in 0..m {
        let angle = (2.0 * (k as f64) + 1.0) * std::f64::consts::PI / (2.0 * (m as f64));
        nodes.push(angle.cos());
    }
    for (j, row) in cos_table.iter_mut().enumerate() {
        for (k, cell) in row.iter_mut().enumerate() {
            *cell = ((j as f64) * (2.0 * (k as f64) + 1.0) * std::f64::consts::PI
                / (2.0 * (m as f64)))
            .cos();
        }
    }

    let mut data_words = Vec::with_capacity(num_records * rsize + 4);
    for rec_idx in 0..num_records {
        let t0 = start_et_s + (rec_idx as f64) * intlen_s;
        let mid = t0 + radius_s;
        data_words.push(mid);
        data_words.push(radius_s);

        let mut samples = Vec::with_capacity(m);
        for &tau in &nodes {
            let (pos, _) = eval_t21((mid + radius_s * tau).clamp(start_et_s, end_et_s))?;
            samples.push(pos);
        }

        for axis in 0..3 {
            for (j, row) in cos_table.iter().enumerate() {
                let scale = if j == 0 {
                    1.0 / (m as f64)
                } else {
                    2.0 / (m as f64)
                };
                let mut c_j = 0.0_f64;
                for (k, sample) in samples.iter().enumerate() {
                    let s_val = sample.get(axis).copied().unwrap_or(0.0);
                    let c_jk = row.get(k).copied().unwrap_or(0.0);
                    c_j += s_val * c_jk;
                }
                data_words.push(scale * c_j);
            }
        }
    }

    // 4-word Type 2 directory trailer: INIT, INTLEN, RSIZE, N
    data_words.push(start_et_s);
    data_words.push(intlen_s);
    data_words.push(rsize as f64);
    data_words.push(num_records as f64);

    Ok((
        DafSegmentInput {
            name: name.to_owned(),
            target_id,
            center_id,
            frame_id,
            spk_type: 2,
            start_et_s,
            end_et_s: effective_end_et_s,
            data_words,
        },
        best_pos_err_m,
        best_vel_err_mm_s,
    ))
}

