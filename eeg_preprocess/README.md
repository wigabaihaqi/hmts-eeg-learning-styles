# eeg_preprocess – raw-data QC and Tier 2 in one step

`preprocess.py` does three things for every recording:

1. **Checks the raw file.** It writes a status (OK, WARN or FAIL) and the reasons.
2. **Builds Tier 2.** Data are band-pass filtered at 8–30 Hz and cut into 1-s windows with the same rule as the released dataset.
3. **Adds per-window quality columns.** These let users of the data filter windows themselves.

On the HMT-S data the output is identical to the released Tier 2:

- **Recordings with per-sample timestamps:** difference < 10⁻⁷ µV.
- **Burst-timestamped recordings:** identical to the output of `../rebuild_tier2.py`.

## Usage

```bash
python preprocess.py              # QC + Tier 2 for all recordings
python preprocess.py --check-only # QC report only (fast)
```

**Dataset location.** Set it with the environment variable `HMTS_DATASET_DIR` or with `DATASET_DIR` in `settings.py`.

**Outputs.** Everything is written to `eeg_preprocess/output/`; the dataset folder is never modified.

**Single recording.** To process one recording only (e.g. a new or replaced participant), set in `settings.py`:

```python
ONLY = [("verbal", 7)]
EXTRA_COLUMNS = False   # exactly the columns of the released Tier 2 files
```

## Output

| File | Content |
|---|---|
| `output/<Tier 2 folder>/<class>-NN-subjectN-preprocessed.csv` | Tier 2 with the released columns, plus `n_real_samples`, `n_interpolated` and `faster_outlier` when `EXTRA_COLUMNS = True` |
| `output/qc_report.csv` | One row per recording; see the list below |
| `output/qc_channels.csv` | Per channel: robust SD, PREP flags, inferred gain |

Columns of `qc_report.csv`:

- status and reasons;
- effective sampling rate and missing samples;
- percentage of repeated timestamps;
- gain per channel;
- channels flagged by PREP;
- median robust SD;
- percentage of FASTER outlier windows;
- time model used.

## Checks

| Check | Method | Status |
|---|---|---|
| Missing samples | Jumps in the board's `Sample Index` (0–255) within the task segment | FAIL > 5 %, WARN > 0.1 % |
| Effective sampling rate | Samples / task duration | FAIL < 90 % of 250 Hz |
| Burst-written timestamps | Percentage of identical consecutive timestamps | Handled automatically (see Time axis) |
| ADC gain | Quantization step of the raw values: 4.5 V / gain / (2²³ − 1) | WARN if any channel ≠ ×24 |
| Bad channels | PREP criteria (Bigdely-Shamlo et al., 2015) | WARN if ≥ 2 channels |
| Amplitude | Median robust SD after a 1-Hz high-pass | WARN > 100 µV |
| Window outliers | FASTER epoch criteria, \|z\| > 3 (Nolan et al., 2010) | WARN > 10 % of windows |
| Empty windows | No sample in a window | FAIL |

All thresholds are in `settings.py`. Values that are not taken from a publication are marked "package choice";
report them as such.

## Time axis (`TIME_MODEL` in `settings.py`)

- **`"auto"` (default).** Sample-counter timing for burst-written recordings and the original timestamp rule for
  all others. This is how the released dataset was produced.
- **`"sample_counter"`.** The time of sample *n* is offset + *n* / rate for every recording. This is recommended for
  new datasets.
  - Missing samples are filled by linear interpolation and counted in `n_interpolated`.
  - If the counter resets or too many samples are missing, the script falls back to timestamps and records this in
    the report.
- **`"timestamps"`.** Always the original rule.

## Using the package for other recordings

Edit the **PATHS** block of `settings.py`:

- `DATASET_DIR` / `RAW_DIR`;
- `CLASSES`, i.e. class folder, Tier 2 folder and label;
- the raw and log file patterns.

If the task differs, also change `N_ITEMS`, `N_PHASES` and `WIN_PER_PHASE`. Windows start at the
`Session Started` column of the answer log.

## Recommendations for future recordings

- **Event markers.** Send stimulus markers to the EEG `Marker Channel` (e.g. via LSL or BrainFlow from the stimulus
  application), so that alignment does not depend on computer clocks.
- **Gain settings.** Set the gain of all channels to ×24 before each session and save the OpenBCI GUI settings file.
- **Immediate check.** Run `python preprocess.py --check-only` right after each recording, so that a FAIL can be
  re-recorded while the participant is still available.
