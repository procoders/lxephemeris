# `lxephemeris` Dual Licensing Policy

Copyright (c) 2026 ProCoders (`https://procoders.tech/`) / `lxephemeris` Contributors. All rights reserved.

`lxephemeris` (`lxephemeris-core`, `lxephemeris-cli`, `lxephemeris-py`, `lxephemeris-node`, and associated build/verification tools) is **source-available** software distributed under a **Dual License** model:

1. **Personal, Academic, Research, and Non-Commercial Use (`PolyForm-Noncommercial-1.0.0`)**:
   - You may use, copy, modify, and redistribute this software free of charge for personal use or use by a non-commercial organization under the terms of the **PolyForm Noncommercial License 1.0.0**.
   - See [LICENSE-NONCOMMERCIAL.md](./LICENSE-NONCOMMERCIAL.md) for the complete, unmodified official license text.
   - **Additional Pre-Purchase Evaluation Grant**: For-profit companies and commercial developers are additionally granted a royalty-free right to evaluate and test `lxephemeris` in non-production development/staging environments for up to **30 days** solely to assess suitability for purchasing a Commercial License, provided the software is not used in production or to serve external customers.

2. **Commercial Use (`LicenseRef-Commercial`)**:
   - Any self-hosted use of `lxephemeris` in a commercial product, paid SaaS/API service, monetized mobile or web application, internal business production workflow, or any activity outside "Permitted Purpose" under the PolyForm Noncommercial License 1.0.0 requires a valid **Commercial License** granted by **ProCoders**.
   - **How to Obtain a Self-Hosted Commercial License**:
     - Submit a request via **[https://procoders.tech/contacts/](https://procoders.tech/contacts/)**
     - Or email **[lxephemeris@procoders.tech](mailto:lxephemeris@procoders.tech)**
   - **Managed Cloud Alternative ([AstroAPI — `https://astrology-api.io/`](https://astrology-api.io/))**:
     - If you prefer a turnkey cloud REST API (with built-in city geocoding, historical timezone resolution, SVG chart wheel rendering, synastry scoring, and interpretations) instead of self-hosting `lxephemeris`, commercial usage is included directly in an **[AstroAPI (`https://astrology-api.io/`)](https://astrology-api.io/)** subscription without needing a separate `lxephemeris` engine license.
   - Authorized deployments within the copyright holder's proprietary services (including **[AstroAPI](https://astrology-api.io/)** and **[Astrolium](https://astrolium.com/)**) are covered under the owner's commercial rights.

---

## Commercial License Activation (`LXEPHEMERIS_LICENSE_KEY`)

Commercial licensees receive a cryptographically signed offline license key (`LXEPH-1.<payload>.<signature>`) verified locally via `Ed25519` (RFC 8032) with **zero network calls at runtime**.

- **Environment Variable**:
  ```bash
  export LXEPHEMERIS_LICENSE_KEY="LXEPH-1...."
  ```
- **Mounted Secret File (Docker / Kubernetes)**:
  ```bash
  export LXEPHEMERIS_LICENSE_FILE="/run/secrets/lxephemeris.lic"
  ```
- **Programmatic Registration**:
  - **Rust**: `lxephemeris_core::register_commercial_license("LXEPH-1....");`
  - **Python**: `lxephemeris.register_commercial_license("LXEPH-1....")`
  - **Node.js / TypeScript**: `registerCommercialLicense("LXEPH-1....")`
- **CLI Verification**:
  ```bash
  lxephemeris-cli license-verify --key "LXEPH-1...."
  ```
- **Perpetual Fallback on Covered Releases (`updates_until`)**:
  Annual commercial subscriptions encode an `updates_until` date (`YYYY-MM-DD`) compared against the engine's build release date (`ENGINE_RELEASE_DATE`), ensuring that a running production server **never** expires or alters numerical behavior over time.

---

## Third-Party Components and Astronomical Datasets

`lxephemeris` is a clean-room implementation written in pure Rust. It depends exclusively on permissively or weak-copyleft licensed libraries (`MPL-2.0`, `MIT`, `BSD-3-Clause`, `Apache-2.0`, `IAU SOFA Software License`) and publicly available astronomical catalogs (`NASA JPL DE440s / Horizons`, `CDS VizieR I/311, I/239, V/137D`, `IERS EOP / Bulletin C`).

See [THIRD_PARTY_NOTICES.md](./THIRD_PARTY_NOTICES.md) for full provenance audit records, copyright notices, and mandatory attributions (including the 6-clause IAU SOFA Software License notice).
