"""
Settings (settings.py) for the reusable preprocessing package (eeg_preprocess).

Edit the PATHS block for a new dataset. Every other value has a comment
saying where it comes from: the HMT-S protocol, the original preprocessing
code, a published method, or a choice made for this package (marked
"package choice", report it as such).
"""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ============================================================== PATHS
# Folder that contains one sub-folder per class, each with "Subject N" folders
# holding the raw OpenBCI csv and the answer log (xlsx).
# Set the environment variable HMTS_DATASET_DIR to the downloaded dataset folder,
# or edit the default below.
DATASET_DIR = Path(os.environ.get("HMTS_DATASET_DIR", HERE.parents[3] / "Dataset" / "2. HMT-S"))
RAW_DIR = DATASET_DIR / "1. Raw EEG and Stimulus Logs"
OUT_DIR = HERE / "output"            # Tier 2 files + QC reports are written here

# class -> (raw sub-folder, Tier 2 sub-folder, label written in the "class" column)
CLASSES = {
    "verbal":     ("3. Verbal",     "2. Preprocessed-verbal-visual",     "verbal"),
    "visual":     ("4. Visual",     "2. Preprocessed-verbal-visual",     "visual"),
    "active":     ("1. Active",     "3. Preprocessed-active-reflective", "active"),
    "reflective": ("2. Reflective", "3. Preprocessed-active-reflective", "reflective"),
}
# Process only some recordings, e.g. ONLY = [("verbal", 7)]; None = all
ONLY = None

RAW_PATTERN = r".*\.(csv|txt)"       # raw OpenBCI export in each Subject folder
EXTRA_COLUMNS = True                 # False = Tier 2 with exactly the columns of the released files
LOG_PATTERN = r"ls-quiz.*\.xlsx"     # answer log in each Subject folder
LOG_SESSION_START = "Session Started"

# ============================================================== RECORDING (HMT-S protocol)
FS = 250.0
CHANNELS = ["Fp1", "Fp2", "F3", "F4", "C3", "C4", "O1", "O2"]
EXG_COLS = [f"EXG Channel {i}" for i in range(8)]
N_ITEMS = 6
N_PHASES = 2                         # 1 = question, 2 = question + options
WIN_PER_PHASE = 15
WIN_SAMPLES = 250

# ============================================================== TIER 2 (original preprocessing code)
BAND = (8.0, 30.0)                   # Butterworth band-pass, applied with filtfilt
FILTER_ORDER = 4
WINDOW_SEC = 1.0
WINDOW_GAP_MS = 1                    # original code: next window starts 1 ms after the
                                     # previous one ends; set 0 for strictly contiguous windows

# ============================================================== TIME AXIS
# "auto"           : sample-counter timing only for burst-written recordings, the
#                    original timestamp rule otherwise (= how dataset version 2 was made)
# "sample_counter" : always reconstruct sample times from the board's counter
#                    (recommended for new datasets)
# "timestamps"     : always the original rule (version 1 behaviour)
TIME_MODEL = "auto"
BURST_PCT = 50.0                     # package choice: > 50 % identical consecutive
                                     # timestamps = burst-written
ENVELOPE_PERCENTILE = 0.5            # package choice: offset of the time model on the
                                     # lower envelope of the timestamps
MAX_FORWARD_STEP = 128               # counter steps 2..128 = lost samples; larger
                                     # (backward) steps = counter reset / reordering
TASK_MARGIN_S = 2.0

# ============================================================== QC (recording level)
# Bad channels: PREP criteria (Bigdely-Shamlo et al., 2015, Front. Neuroinform. 9:16)
PREP_HIGHPASS_HZ = 1.0
PREP_LOWPASS_HZ = 50.0
PREP_Z_DEVIATION = 5.0
PREP_Z_HF_NOISE = 5.0
PREP_CORR_THRESHOLD = 0.4
PREP_CORR_PERCENTILE = 98
PREP_CORR_BAD_FRACTION = 0.01
# Window outliers: FASTER epoch criteria (Nolan et al., 2010, J. Neurosci. Methods 192:152)
FASTER_Z = 3.0
# ADS1299 (OpenBCI Cyton) gain from the quantisation step 4.5 V / gain / (2^23 - 1)
ADS_GAINS = [1, 2, 4, 6, 8, 12, 24]
EXPECTED_GAIN = 24                   # OpenBCI default

# Status rules (package choices; adjust and report them)
FAIL_LOST_PCT = 5.0                  # FAIL if > 5 % of task samples are missing
FAIL_RATE_FRACTION = 0.9             # FAIL if effective rate < 90 % of FS
WARN_LOST_PCT = 0.1                  # WARN if > 0.1 % of task samples are missing
WARN_BAD_CHANNELS = 2                # WARN if >= 2 channels flagged by PREP
WARN_ROBUST_SD_UV = 100.0            # WARN if the median channel robust SD > 100 uV
WARN_OUTLIER_WINDOWS_PCT = 10.0      # WARN if > 10 % of windows are FASTER outliers
