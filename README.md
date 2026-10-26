# HMT-S EEG dataset: preprocessing and technical validation

Code accompanying the data article *Raw and preprocessed electroencephalography dataset of verbal–visual and
active–reflective learners during a short Hagen matrix-reasoning task* (Data in Brief, under review) and the dataset on
Mendeley Data: `[dataset DOI]`.

The code covers three steps:

1. **Preprocessing.** Builds Tier 2 (8–30 Hz band-pass filtered, 1-s windows) from the raw recordings and checks
   every raw file (`eeg_preprocess/`).
2. **Timing check.** Reproduces the Tier 2 segmentation and reconstructs sample times for recordings whose
   timestamps were written in bursts (`rebuild_tier2.py`).
3. **Technical validation.** Produces the validation figure, tables and summary numbers reported in the article
   (`technical_validation.py`).

## Contents

```
README.md                  this file
requirements.txt           Python packages
config.py                  paths and all thresholds for steps 2-3 (each with its source)
eeg_preprocess/
  preprocess.py            step 1: QC report + Tier 2 from raw data
  settings.py              paths, task design and thresholds for step 1
  README.md                details of step 1
rebuild_tier2.py           step 2
technical_validation.py    step 3
build_release.py           authors only: assembles the Mendeley upload from the private working copy
LICENSE                    MIT
```

## Installation

Python 3.10 or later.

```bash
pip install -r requirements.txt
```

## Data

1. Download the dataset from Mendeley Data (`[dataset DOI]`) and unzip it. It contains:

   ```
   1. Raw EEG and Stimulus Logs/<class>/Subject N/SubjectN.csv, LS-Quiz-subjectN.xlsx
   2. Preprocessed-verbal-visual/
   3. Preprocessed-active-reflective/
   ```

2. Tell the code where the dataset is, by either:
   - setting an environment variable (recommended):
     ```bash
     export HMTS_DATASET_DIR="/path/to/unzipped/dataset"          # macOS / Linux
     set HMTS_DATASET_DIR=C:\path\to\unzipped\dataset             # Windows (cmd)
     ```
   - or editing `DATASET_DIR` in `config.py` and `eeg_preprocess/settings.py`.

Results are written inside this repository (`eeg_preprocess/output/`, `rebuilt_tier2/`, `results/`). The dataset
folder is never modified.

## Step 1 – Preprocessing: raw data → Tier 2

```bash
cd eeg_preprocess
python preprocess.py --check-only   # QC report only (fast)
python preprocess.py                # QC report + Tier 2 files
cd ..
```

For every recording, the script runs the following checks and writes them to `eeg_preprocess/output/qc_report.csv`,
with a status of OK, WARN or FAIL and the reasons:

- **Sample continuity:** samples missing from the board's sample counter.
- **Sampling rate:** effective sampling rate.
- **Timestamps:** burst-written timestamps.
- **ADC gain:** gain per channel, inferred from the quantization step.
- **Bad channels:** PREP criteria.
- **Amplitude:** signal amplitude.
- **Window outliers:** FASTER criteria.

It then builds Tier 2 in the same way as the released files:

- **Filter:** 4th-order Butterworth band-pass, 8–30 Hz, `scipy.signal.filtfilt`, applied to the whole recording.
- **Windows:** consecutive 1-s windows from `Session Started` in the answer log, each starting 1 ms after the
  previous one ends. Each window holds the first 250 samples.
- **Time axis:** with `TIME_MODEL = "auto"` (default), the board's sample counter is used for burst-written
  recordings and the timestamps for all others. This reproduces the released Tier 2 exactly.

With `EXTRA_COLUMNS = True` (default) three quality columns are added per window: `n_real_samples`, `n_interpolated`
and `faster_outlier`. Set it to `False` to obtain exactly the columns of the released files.

See `eeg_preprocess/README.md` for all options, including using the package on new recordings.

## Step 2 – Timing check (`rebuild_tier2.py`)

```bash
python rebuild_tier2.py
```

- **Input:** the raw files and the Tier 2 files of the dataset.
- **Self-check:** reapplies the timestamp-based segmentation rule to the raw data and compares the result with
  the Tier 2 file.
- **Rebuilding:** for recordings with burst-written timestamps (> 50 % identical consecutive timestamps), it rebuilds
  Tier 2 with sample times reconstructed from the sample counter. The time of sample *n* is *t*₀ + *n*/*f*ₛ:
  - *f*ₛ comes from a least-squares fit of timestamps on the counter;
  - *t*₀ is placed on the lower envelope of the timestamps.
- **Output:** `rebuilt_tier2/rebuild_report.csv` (one row per recording) and the rebuilt files.
- **On the released dataset:** the nine burst-written recordings already contain the reconstructed timing, so the
  rebuilt files are identical to the released ones. The remaining recordings reproduce exactly
  (difference < 10⁻¹⁰ µV).

## Step 3 – Technical validation (`technical_validation.py`)

```bash
python technical_validation.py
```

Writes to `results/`:

| File | Content |
|---|---|
| `figure_technical_validation.png/.pdf` | Six-panel figure (Fig. 11 of the article) |
| `summary.txt` | All numbers reported in the Technical validation section |
| `qc_windows_per_subject.csv` | (A) window completeness, sampling continuity, timestamp bursts |
| `qc_channel_quality.csv`, `qc_channel_flag_rate.csv` | (B, C) PREP statistics and flags, amplitude, inferred gain |
| `psd_per_subject.csv`, `psd_grand_average.csv` | (D) power spectra |
| `aperiodic_fooof.csv` | (E) aperiodic exponent, offset, R², number of peaks |
| `faster_window_flags.csv` | (F) percentage of outlier windows per recording |
| `behavior_per_item.csv`, `behavior_summary.csv` | Task performance (accuracy, response time, unanswered items) |

Run Step 2 first: `technical_validation.py` reads `rebuilt_tier2/rebuild_report.csv` to know which recordings have reconstructed timing (panel A) and, for panel F, uses the rebuilt files listed there. All analyses are descriptive; no
group comparison is tested.

### Criteria and sources

| Panel | Criterion | Source |
|---|---|---|
| A | Window completeness: ≥ 250 raw samples = complete, 1–249 = interpolated, 0 = zero-filled | Segmentation rule of the dataset |
| A | Effective sampling rate, timestamp gaps > 0.1 s, burst-written timestamps | Data-integrity checks defined for this dataset |
| B, C | Deviation: robust z-score of the robust SD (0.7413 × IQR) > 5 | PREP [1] |
| B, C | High-frequency noise: robust z-score of MAD(> 50 Hz) / MAD(< 50 Hz) > 5 | PREP [1] |
| B, C | Correlation: 98th percentile of \|r\| with other channels < 0.4 in > 1 % of 1-s windows | PREP [1] |
| B, C | Low SNR, NaN, constant signal | PREP [1] |
| – | RANSAC not applied: with 8 channels the 25 % predictor subset is only 2 channels | – |
| – | ADC gain from the quantization step 4.5 V / gain / (2²³ − 1); DC offset; 99.9th percentile of sample-to-sample steps | Descriptive |
| D | Welch PSD, 2-s Hann segments, 50 % overlap, shown 1–45 Hz | [3]; range as in [4] |
| E | FOOOF, fixed aperiodic mode, 2–40 Hz, peak width 1–8 Hz, max 6 peaks, min peak height 0.1, peak threshold 2.0 | [3] |
| F | Tier 2 window flagged if \|z\| > 3 for amplitude range, deviation from the channel mean or variance (epoch parameters 4–6) | FASTER [2] |

Every threshold is defined, with its source, in `config.py`.

## Citation

If you use this code, please cite the data article and the dataset: `[citation to be added]`.

## References

1. N. Bigdely-Shamlo, T. Mullen, C. Kothe, K.-M. Su, K.A. Robbins, The PREP pipeline: standardized preprocessing
   for large-scale EEG analysis, Front. Neuroinform. 9 (2015) 16. https://doi.org/10.3389/fninf.2015.00016
2. H. Nolan, R. Whelan, R.B. Reilly, FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection,
   J. Neurosci. Methods 192 (2010) 152–162. https://doi.org/10.1016/j.jneumeth.2010.07.015
3. T. Donoghue et al., Parameterizing neural power spectra into periodic and aperiodic components,
   Nat. Neurosci. 23 (2020) 1655–1665. https://doi.org/10.1038/s41593-020-00744-x
4. H. Kim, S.C. Jun, C.S. Nam, Data in Brief 68 (2026) 113064.

## License

MIT (see `LICENSE`).
