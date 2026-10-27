# Архитектурная спецификация `lxephemeris` (v3.1 — Финальная спецификация после Round 2 аудита Claude Opus и Codex)

Связанные документы:
- [PLAN.md](./PLAN.md) — пошаговый план разработки, верификации и приёмки (v3.1);
- [REVIEW_AND_SPEC_ASTRA.md](./REVIEW_AND_SPEC_ASTRA.md) — декларативный реестр объектов и контрактов `astra` (68 малых тел, 23 системы домов, 12 аянамс, 3 Лилит, 4 узла, 8 уранидов + `Transpluto`, 50 звёзд);
- [REVIEW_CLAUDE.md](./REVIEW_CLAUDE.md) и [REVIEW_CLAUDE_OPUS_V3.md](./REVIEW_CLAUDE_OPUS_V3.md) — раунды 1 и 2 математического и системного аудита Claude Opus;
- [REVIEW_CODEX_V3.md](./REVIEW_CODEX_V3.md) — раунд 2 аудита совместимости и полноты контрактов от Codex (`astra`).

---

## 1. Миссия, Лицензия и Контур «Чистой комнаты»

`lxephemeris` — кроссплатформенный, неизменяемый в рантайме (`Send + Sync`), потокобезопасный эфемеридный и астрологический движок на чистом Rust (`lxephemeris-core`) с нативными бинарными пакетами для **Python (`PyO3` + `maturin`)** и **Node.js / TypeScript (`napi-rs`)**, созданный как полноценная и более точная замена C-движка Swiss Ephemeris без глобального мутабельного состояния.

### 1.1. Лицензия продукта (Dual Licensing — решено владельцем)
- **Личное / некоммерческое использование**: бесплатно по **PolyForm Noncommercial License 1.0.0** (`LicenseRef-PolyForm-Noncommercial-1.0.0`).
- **Коммерческое использование**: платная коммерческая лицензия (`LicenseRef-Commercial`).

### 1.2. Разрешённые пермиссивные источники и правила Clean-Room
1. **Разрешённые библиотеки и стандарты**:
   - **[`nyx-space/anise`](https://github.com/nyx-space/anise)** и **[`nyx-space/hifitime`](https://github.com/nyx-space/hifitime)** (**MPL-2.0**, фиксированная версия `=0.9.x` с проверкой `assert_impl_all!(Almanac: Send, Sync)`) — чтение ядер DAF/SPK и наносекундная арифметика шкал времени (`Epoch` в TAI). Доступ к ядрам в `lxephemeris-core` изолирован за внутренним трейтом `SpkSource`. Поддержка **SPK Type 21** (`ModifiedDiffType21`) проходит обязательный gate в Фазе `0.1`; при любых ограничениях `build-spk` конвертирует траектории малых тел в полиномы Чебышёва **SPK Type 2**.
   - **[`astro-xao/sofars`](https://github.com/astro-xao/sofars)** (**MIT + IAU SOFA License**) и **[`liberfa/erfa`](https://github.com/liberfa/erfa)** (**BSD-3-Clause**) — фундаментальные алгоритмы IAU 2000/2006 на чистом Rust без `unsafe` (с обязательным уведомлением IAU SOFA Board в `THIRD_PARTY_NOTICES.md`).
   - **[`jayeshmepani/jpl-ephemeris` (JME)](https://github.com/jayeshmepani/jpl-ephemeris)** (**MIT**), **[`skyfielders/python-skyfield`](https://github.com/skyfielders/python-skyfield)** (**MIT**), **[`cosinekitty/astronomy`](https://github.com/cosinekitty/astronomy)** (**MIT**), **[`TheDaniel166/moira`](https://github.com/TheDaniel166/moira)** (**MIT**) — архитектурные и алгоритмические референсы, используемые совместно с академическими первоисточниками (*Lieske et al. 1977*, *Holden 1977*, *Munkasey 1991*, *Simon et al. 1994*, *Chapront & Francou 2002*, *Neely 1980*, *Espenak & Meeus 2006*). Перед использованием материалов JME в Фазе `0.0` проводится аудит provenance коммитов и констант с записью в `THIRD_PARTY_NOTICES.md`.
2. **Строгие границы Clean-Room**:
   - **Запрещено**: открывать, читать, переводить или линковать исходный код Swiss Ephemeris (`swisseph`/`pyswisseph`), `Kerykeion`, `libephemeris`, любых GPL/AGPL-проектов, а также `Anonyfox/celestine` (содержит прямой перевод `swehouse.c`).
   - **Разрешено**: (а) читать публичную пользовательскую документацию SE исключительно для сверки внешних идентификаторов, флагов и конвенций; (б) запускать немодифицированный бинарник SE во внешнем изолированном процессе (`black-box oracle`) для снятия числовых тестовых векторов Уровня 3 (`testdata/golden/l3_se/`).

---

## 2. Принятые решения Архитектора и Доменного эксперта (Settled Calls)

### 2.1. Архитектура данных, покрытие и офлайн-рантайм
- **Строгий офлайн-режим в рантайме**: `Engine` никогда не обращается к сети при импорте или во время расчётов. Поддерживаются `Engine::from_data_dir`, `Engine::from_bsp_paths` и `Engine::from_bytes` (для встраивания и WASM). Бинарные `.bsp` поставляются отдельными пакетами данных (`lxephemeris-data-default`, `lxephemeris-data-extended`) либо через `LXEPHEMERIS_DATA_DIR`, чтобы не раздувать базовые Python wheels / npm-пакеты.
- **Два профиля `.bsp` и операционно-зависимый `safe_interval`**:
  - **Базовый профиль (`[1850-01-01T00:00:00 TT, 2100-01-01T00:00:00 TT)`)**: `de440s.bsp` (~32 МБ, `raw_interval`: `1849-12-26 .. 2150-01-22`) + единый сшитый `asteroids_68_1850_2100.bsp` (`raw_interval`: `1849-12-26 .. 2100-01-05`, все 68/68 малых тел обязательны — неполный пакет блокирует сборку).
  - **Расширенный / кастомный профиль**: загрузка полного `de440.bsp` (~114 МБ, `1549–2650`), `de441.bsp` или `asteroids_68_1600_2500.bsp` с автоматическим чтением сегментов DAF/SPK и валидацией `manifest.json` (источник, версия решения, `center_id`, `frame_id`, `spk_type`, отсутствие дыр/перекрытий, `SHA-256`).
  - **Правило `safe_interval`**: для каждой операции безопасный интервал вычисляется как $\text{raw\_interval} \setminus (\pm 2h \mp \max \tau_{\text{light}} \mp \Delta T_{\text{margin}})$, а для `Interpolated_Lilith` — с учётом окна поиска $\pm 2$ апсид Луны ($\approx \pm 60\text{ сут}$, т.е. `[1850-03-01, 2149-11-01)` на чистом `de440s.bsp`).
- **Двухуровневое хранение таблиц времени ($\Delta T$, `UT1-UTC`, `leap seconds`)**:
  - Базовый срез (исторические сплайны *Stephenson–Morrison–Hohenkerk 2016* + ряд USNO/IERS `finals2000A` `UT1-UTC` + таблица секунд координации IERS Bulletin C, ~50 КБ) вшивается в `lxephemeris-core`.
  - Команда `lxephemeris-cli data sync` скачивает обновлённые таблицы IERS во временный файл, проверяет `SHA-256` и схему, и атомарно (`rename`) сохраняет их в `LXEPHEMERIS_DATA_DIR`. При старте `Engine` выбирает таблицу с более новым `iers_release_date` (либо жёстко зафиксированную через `EngineConfig::pin_time_table_version` для 100% воспроизводимости тестов). За пределами покрытия `finals2000A` применяется параболическая экстраполяция *Morrison & Stephenson* с возвратом флага `time_extrapolated: true` и `time_table_version` в метаданных.

### 2.2. Астрометрический конвейер, скорости, режимы и метаданные
1. **Итерация светового времени и многотельная ОТО-редукция**:
   - Из `SpkSource` (`anise`) запрашиваются только сырые барицентрические состояния (`Aberration::NONE`).
   - Световое время $\tau$ находится итерацией неподвижной точки до $|\Delta\tau| < 10^{-12}\text{ сут}$ с жёстким лимитом `MAX_LT_ITER = 10` (при превышении — `Err(EphemerisError::NonConvergence)`).
   - Гравитационное отклонение света (`sofars::astro::ld`) вычисляется с ретардацией момента сближения луча для **Солнца, Юпитера и Сатурна** (плюс **Земли** при топоцентрическом наблюдении).
   - Релятивистская годичная и суточная аберрация вычисляется через `sofars::astro::ab`.
2. **Расчёт видимых скоростей (4-точечный декартов шаблон) и `Moment`**:
   - Скорость в выходной системе координат вычисляется симметричным 4-точечным шаблоном 4-го порядка (`[-2h, -h, +h, +2h]` в шкале `TT`) по 3D декартову вектору выходного положения `(x, y, z)`, после чего `(x, y, z, vx, vy, vz)` аналитически переводится в сферические `(lon, lat, r, lon_speed, lat_speed, r_speed)` с обработкой полюсов $|\beta| \to 90^\circ$.
   - Шаг $h$ фиксирован по классам тел: для Луны и лунных точек $h = 10^{-4}\text{ сут}$ ($\approx 8.64\text{ с}$), для Солнца, планет, 68 малых тел и уранидов $h = 10^{-3}\text{ сут}$ ($\approx 86.4\text{ с}$).
   - Контекст эпохи `Moment` (`engine.at(epoch, observer)`) кэширует на центральную эпоху $t_0$ шкалы `TT/TDB/UT1`, матрицу BPN (IAU 2006 / 2000A), $\varepsilon$, звёздное время GAST (строго от `Ut1`), состояние Земли/наблюдателя и офсет аянамсы, а при запросе скоростей лениво (`OnceLock`) инициализирует смещённые под-контексты $t_0 \pm h, t_0 \pm 2h$, общие для всех тел карты.
3. **Целевые тела (Target Mapping), пространства идентификаторов и режимы вывода**:
   - Раздельные типизированные пространства ID: `BodyId::Planet(PlanetId)`, `BodyId::MinorPlanet(MpcNumber)`, `BodyId::LunarPoint(LunarPointId)`, `BodyId::Hypothetical(HypotheticalId)`, `BodyId::FixedStar(HipId)`.
   - Соответствие телам `DE440s`: для `Mercury (199)`, `Venus (299)`, `Earth (399)`, `Moon (301)` используются центры масс тел (`TargetKind::BodyCenter`); для `Mars (4)`, `Jupiter (5)`, `Saturn (6)`, `Uranus (7)`, `Neptune (8)`, `Pluto (9)` — барицентры планетных систем (`TargetKind::SystemBarycenter`), `Sun (10)`, `EMB (3)`, `SSB (0)`. Метаданные ответа всегда возвращают `target_id`, `center_id` и `target_kind`.
   - Поддерживаемые оси и режимы запроса (`CalcOptions`): `Origin::{Geocentric, Topocentric(Observer), Heliocentric, Barycentric}`, `Plane::{Ecliptic, Equatorial, Horizontal}`, `Equinox::{OfDate, J2000}`, `NutationMode::{True, Mean}`, `AbMode::{Apparent, Astrometric, Geometric}`, `Zodiac::{Tropical, Sidereal(AyanamsaModel)}`.
   - Граница применимости прецессии IAU 2006 (`sofars::pnp::bp06`) — `[1000 CE, 3000 CE]`; при загрузке сверхдлинных ядер (`de441.bsp`) за пределами этого окна используется долгосрочная модель *Vondrák et al. 2011* (`eraLtp`) либо возвращается `OutOfModelValidity`.

### 2.3. Паспорта 12 аянамс (`lahiri`, `krishnamurti`, `raman`, `yukteshwar`, `true_citra`, `jn_bhasin`, `ushashashi`, `true_revati`, `true_pushya`, `ss_citra`, `suryasiddhanta`, `fagan_bradley`)
- **9 эпохо-фиксированных моделей** (заданных значением $A(t_0)$ в референсную эпоху $t_0$):
  - Распространяются двухэпоховым полиномом общей прецессии по долготе *Lieske et al. (1977, Table 5)* / *Indian Astronomical Ephemeris 1989, p. 525*:
    $$p_A(t_0 \to t) = (5029.0966'' + 2.22226'' T - 0.000042'' T^2) t + (1.11161'' - 0.000127'' T) t^2 - 0.000113'' t^3$$
    где $T = (\text{JD}_{\text{TT}}(t_0) - 2451545.0)/36525$ и $t = (\text{JD}_{\text{TT}} - \text{JD}_{\text{TT}}(t_0))/36525$ в юлианских столетиях `TT` (знаки, якоря $t_0, A(t_0)$ и скорость $\dot{A}(t)$ фиксируются в паспортах Фазы `0.3`).
  - Поскольку классические индийские таблицы фиксируют офсет аянамсы в IAU 1976, а тропические долготы планет редуцируются в IAU 2006, систематическая разница скорости прецессии ($\approx 0.3''/\text{век}$) явно документируется в `docs/SE_DIVERGENCES.md`.
- **3 звёздно-якорные модели** (`true_citra` — Spica / HIP 65474 на $180^\circ$, `true_revati` — $\zeta$ Psc A / HIP 5737 на $359^\circ 50'$, `true_pushya` — $\delta$ Cnc / HIP 42911 на $106^\circ$):
  - Базовое среднее значение: $A_{\text{mean}}(t) = \lambda_{\text{star, mean\_geo\_no\_ab}}(t) - \lambda_{\text{target}}$, где положение звезды переводится из каталога **VizieR I/311** (*van Leeuwen 2007*, с лучевыми скоростями из **XHIP V/137D**) с собственным движением и параллаксом (`pmpx`) к средней эклиптике даты **без годичной аберрации**.
  - При `NutationMode::True` к $A_{\text{mean}}(t)$ всех 12 моделей (включая 3 звёздные) прибавляется одна и та же нутация по долготе $\Delta\psi$, что гарантирует строгую инвариантность сидерической долготы к флагу нутации.

### 2.4. 23 системы домов и полярная политика (согласовано с контрактом `astra`)
- Реализуются все 23 регистрозависимых кода `astra`: `P, K, O, R, C, A, W, B, M, T, U, X, V, L, Q, I, i, F, N, S, D, Y, H` + `houses_armc` + сектора Гоклена. Для Плацидуса (`P`) уравнение полудуги решается методом Иллинойса на скобке $[-\pi, \pi]$ с лимитом 50 итераций и допуском $10^{-13}\text{ рад}$.
- **Полярная политика (`PolarPolicy`) и 3 категории `astra`**:
  - По умолчанию в SDK установлен **`PolarPolicy::Strict`** (в полном соответствии с контрактом `astra` §5); режим `PolarPolicy::FallbackToPorphyry` включается явно в `HouseOptions`. Ответ всегда содержит `requested_house_system`, `effective_house_system`, `fallback_reason` и `warnings`:
    1. `hard_fail` (`P, K, i`) при $|\phi| \ge 90^\circ - \varepsilon$ (или неразрешимости полудуги): в `Strict` возвращает `Err(EphemerisError::UndefinedHouseSystem)`; в `FallbackToPorphyry` возвращает куспиды Порфирия (`O`) с `effective_house_system: HouseSystem::Porphyry` и `HouseWarning::PolarFallbackToPorphyry`.
    2. `soft_fail` (`R, C, B, T, L, I, Y`) при $|\phi| \ge 90^\circ - \varepsilon$: вычисляет куспиды и добавляет `HouseWarning::PolarUnreliable` при любой политике (при математическом вырождении знаменателя — по правилу `PolarPolicy`).
    3. `polar_safe` (`O, A, W, M, U, X, V, Q, F, N, S, D, H`): работают за полярным кругом; в точных полюсах $|\phi| = 90^\circ$ (где не определён местный меридиан) возвращают `Err(EphemerisError::UndefinedHouseSystem)` для систем, зависящих от RAMC/ASC (кроме `N` — $0^\circ$ Овна).

### 2.5. Лунные узлы, 3 варианта Лилит, 8 Уранидов + `Transpluto` и 50 Звёзд
- **Средние лунные точки (`Mean_Node`, `Mean_South_Node`, `Mean_Lilith`)**: полиномы фундаментальных аргументов *Simon et al. (1994, A&A 282, 663)* в средней эклиптике даты; при `NutationMode::True` к долготе прибавляется $\Delta\psi$ (для совместимости с конвенцией SE и `astra`), широта $\beta = 0^\circ$, дистанция = среднее расстояние.
- **Оскулирующие лунные точки (`True_Node`, `True_South_Node`, `True_Lilith`)**: вычисляются аналитически из **геометрического** геоцентрического вектора состояния Луны $(\mathbf{r}, \mathbf{v})$ из `DE440s` в средней эклиптике даты с гравитационным параметром $\mu = GM_\oplus + GM_{\text{Moon}}$ из заголовка `DE440s` (узел — из $\mathbf{h} = \mathbf{r} \times \mathbf{v}$, апогей — из вектора Лапласа–Рунге–Ленца $-\mathbf{e}$).
- **`Interpolated_Lilith` (Natural Apogee)**: кубическая $C^2$-сплайн интерполяция развёрнутой долготы и расстояния по узлам истинных экстремумов расстояния Земля–Луна ($\dot{r} = 0$: апогеи и перигеи $+ 180^\circ$) в окне $\pm 2$ аномалистических месяцев; точный паспорт алгоритма и допуск относительно SE фиксируются в Фазе `0.3` и `docs/SE_DIVERGENCES.md`.
- **8 Уранидов (`Cupido..Poseidon`) + `Transpluto` (`Isis`)**: линейные кеплеровы элементы *Neely (1980)* для 8 гамбургских точек и опубликованные элементы *Sevin (1946) / Strubell* для `Transpluto`, с переводом из гелиоцентрических координат в геоцентрические через вектор Земли из `DE440s` (с аберрацией/световым временем согласно `AbMode`, паспорт сверяется в Фазе `0.3`).
- **50 неподвижных звёзд**: основной каталог — **VizieR I/311** (*van Leeuwen 2007*, эпоха `J1991.25`, fallback на **I/239** для кратных систем без решения в I/311) + лучевые скорости **XHIP (V/137D)**; при $\varpi \le 0$ параллакс зажимается в угловой режим с флагом `StarQualityWarning::NonPositiveParallax`.

### 2.6. Солнечные и Лунные затмения
- Контакты на фундаментальной плоскости Бесселя строятся по астрометрическим векторам (`AbMode::Astrometric` / без аберрации), а уточнение корней выполняется методом Иллинойса до $10^{-10}\text{ сут}$.
- Зафиксированные физические константы: радиусы Луны $k_{\text{mean}} = 0.2725076$ (полутень/частные фазы) и $k_{\text{umbra}} = 0.272281$ (полная/кольцеобразная фаза по IAU), сжатие Земли WGS84 $f = 1/298.257223563$, атмосферное увеличение земной тени по правилу **Danjon ($1/85$)** (с опцией традиционных $2\%$ в `EclipseOptions`). Разница времени контактов ($\sim 0.6\text{ с}$) относительно канонов с видимой аберрацией документируется в `docs/SE_DIVERGENCES.md`; допуск приёмки контактов — $\le 1.5\text{ с}$ (солнечные) и $\le 15\text{ с}$ (лунные).

---

## 3. Инженерные гарантии Rust, FFI (`Python` / `Node.js`) и бюджет производительности

- **Rust Core (`lxephemeris-core`)**: `#![forbid(unsafe_code)]`, `clippy::unwrap_used = "deny"`, `clippy::panic = "deny"`, `clippy::indexing_slicing = "warn"`; все публичные `enum` и `EphemerisError` помечены `#[non_exhaustive]` и возвращают стабильные коды (`InvalidInput`, `UnsupportedBody`, `UnsupportedModel`, `UnsupportedSpkType`, `OutOfCoverage`, `OutOfModelValidity`, `MissingData`, `CorruptKernel`, `UndefinedHouseSystem`, `NonConvergence`).
- **Параллелизм и безопасность `fork`**: `Engine` не использует глобальный пул Rayon при импорте; для `calc_batch` создаётся ленивый ограниченный пул внутри `Engine` (`EngineConfig::max_threads`), безопасный при `multiprocessing.fork` на Linux и не конфликтующий с пулом воркеров FastAPI.
- **Python (`lxephemeris-py`)**: `PyO3 0.26+` + `maturin`, `abi3-py39`, `#[pymodule(gil_used = false)]`, извлечение срезов NumPy до `py.detach`, типизированная иерархия исключений с полями `code/body/epoch/coverage`, стабы `.pyi`.
- **Node.js (`lxephemeris-node`)**: `napi-rs v3`, асинхронная загрузка `Engine::loadAsync` и `calcBatchAsync` через `AsyncTask`, безопасное владение `Float64Array`, вход времени через `SplitJd`, стабы `index.d.ts`.
- **Бюджет производительности**: $\le 50\text{ мкс}$ на одно тело (без скорости) и $\le 200\text{ мкс}$ (с 4-точечной скоростью) на прогретом `Moment`; полный расчёт натальной карты (10 светил + 23 дома + узлы/Лилит) $\le 1.5\text{ мс}$ в 1 потоке на Apple Silicon / современном x86_64.
