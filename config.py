"""
Configuration for the HMT-S technical-validation scripts.

Edit DATASET_DIR if your dataset folder is in a different place.
Every threshold below is annotated with its source so that it can be
reported (and cited) exactly in the "Technical validation" section.

Sources
  [PREP]   Bigdely-Shamlo N, Mullen T, Kothe C, Su K-M, Robbins KA (2015).
           The PREP pipeline: standardized preprocessing for large-scale EEG
           analysis. Front. Neuroinform. 9:16. doi:10.3389/fninf.2015.00016
  [FASTER] Nolan H, Whelan R, Reilly RB (2010). FASTER: Fully Automated
           Statistical Thresholding for EEG artifact Rejection.
           J. Neurosci. Methods 192(1):152-162. doi:10.1016/j.jneumeth.2010.07.015
  [FOOOF]  Donoghue T, Haller M, Peterson EJ, et al. (2020). Parameterizing
           neural power spectra into periodic and aperiodic components.
           Nat. Neurosci. 23:1655-1665. doi:10.1038/s41593-020-00744-x
  [KIM]    Kim H, Jun SC, Nam CS (2026). Data in Brief 68:113064
           (layout of the validation figure; PSD range 1-45 Hz)
"""
import os
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
# This file lives in:  .../Semester III/Paper/SSRN Dataset/Technical Validation/
# The dataset lives in: .../Semester III/Dataset/2. HMT-S/
HERE = Path(__file__).resolve().parent
# Set the environment variable HMTS_DATASET_DIR to the downloaded dataset folder,
# or edit the default below.
DATASET_DIR = Path(os.environ.get("HMTS_DATASET_DIR", HERE.parents[2] / "Dataset" / "2. HMT-S"))
RESULTS_DIR = HERE / "results"

RAW_DIR = DATASET_DIR / "1. Raw EEG and Stimulus Logs"
PREP_DIRS = [
    DATASET_DIR / "2. Preprocessed-verbal-visual",
    DATASET_DIR / "3. Preprocessed-active-reflective",
]

# --------------------------------------------------------------------------
# Recording parameters (from the acquisition protocol)
# --------------------------------------------------------------------------
FS = 250.0
CHANNELS = ["Fp1", "Fp2", "F3", "F4", "C3", "C4", "O1", "O2"]
EXG_COLS = [f"EXG Channel {i}" for i in range(8)]
N_ITEMS = 6
N_PHASES = 2
N_WIN_PER_PHASE = 15
EXPECTED_WINDOWS = N_ITEMS * N_PHASES * N_WIN_PER_PHASE   # 180
SAMPLES_PER_WINDOW = 250

CLASS_FOLDERS = {"active": "1. Active", "reflective": "2. Reflective",
                 "verbal": "3. Verbal", "visual": "4. Visual"}
CLASS_DIMENSION = {"active": "AR", "reflective": "AR", "verbal": "VV", "visual": "VV"}
CLASS_ORDER = ["active", "reflective", "verbal", "visual"]

# --------------------------------------------------------------------------
# Panel A - window retention
# Rule taken from the dataset's own preprocessing code
# (start <= t <= end, first 250 samples, interpolate if fewer, zeros if none).
# Sampling-continuity checks below are descriptive data-integrity checks
# defined for this dataset (not taken from a published pipeline).
# --------------------------------------------------------------------------
GAP_SEC = 0.1                # jump between consecutive timestamps counted as a gap
MIN_RATE_FRACTION = 0.9      # "irregular" if effective rate < 90 % of FS
EXCLUDE_IRREGULAR_FROM_SPECTRA = True

# --------------------------------------------------------------------------
# Panels B, C - bad-channel detection  [PREP]
# Applied to the raw task segment after a temporary 1-Hz high-pass
# ("tested ... using temporary high pass filtering at 1 Hz", PREP).
# --------------------------------------------------------------------------
PREP_HIGHPASS_HZ = 1.0
PREP_LOWPASS_HZ = 50.0       # split into low / high frequency parts (PREP: 50 Hz low pass)
PREP_ROBUST_SD_FACTOR = 0.7413   # robust SD = 0.7413 x IQR  (PREP)
PREP_Z_DEVIATION = 5.0       # bad-by-deviation: robust z of robust SD > 5  (PREP)
PREP_Z_HF_NOISE = 5.0        # bad-by-HF-noise: robust z of MAD(high)/MAD(low) > 5  (PREP)
PREP_CORR_WIN_SEC = 1.0      # 1-s non-overlapping windows  (PREP)
PREP_CORR_PERCENTILE = 98    # "maximum" correlation = 98th percentile of |r|  (PREP)
PREP_CORR_THRESHOLD = 0.4    # bad-by-correlation if max |r| < 0.4 ...  (PREP)
PREP_CORR_BAD_FRACTION = 0.01  # ... in more than 1 % of windows  (PREP)
# bad-by-lowSNR = bad-by-correlation AND bad-by-HF-noise  (PREP)
# bad-by-NaN and bad-by-NoData (constant signal)  (PREP; the paper gives no
# numeric value for "constant", so a channel is treated as constant when its
# robust SD is numerically zero):
FLAT_EPS_UV = 1e-6
# bad-by-RANSAC (PREP) is NOT applied: it predicts each channel from a random
# 25 % subset of the other channels with spherical splines, which with an
# 8-channel montage means 2 predictor channels; PREP itself stops RANSAC
# when too few channels are available.

# --------------------------------------------------------------------------
# ADC gain inference (ADS1299 on the OpenBCI Cyton): the raw values are
# integer multiples of  4.5 V / gain / (2^23 - 1)  (scale factor used by the
# OpenBCI GUI), so the smallest step between values reveals the gain.
# --------------------------------------------------------------------------
ADS_GAINS = [1, 2, 4, 6, 8, 12, 24]
ADS_STEP_UV = lambda g: 4.5 / g / (2 ** 23 - 1) * 1e6

# Use Tier 2 files from ./rebuilt_tier2 (made by rebuild_tier2.py) for panel F
# when they exist. Original files are used for everything else.
USE_REBUILT_TIER2 = True

# --------------------------------------------------------------------------
# Panel F - window-level artifact screening on Tier 2  [FASTER]
# Epoch parameters 4-6 of FASTER (amplitude range, deviation from the
# channel average, variance; each averaged over channels), z-scored across
# the subject's 180 windows; a window is flagged if |z| > 3 on any parameter.
# --------------------------------------------------------------------------
FASTER_Z = 3.0

# --------------------------------------------------------------------------
# Panels D, E - spectra
# --------------------------------------------------------------------------
WELCH_SEC = 2.0              # Welch, 2-s Hann segments, 50 % overlap  [FOOOF methods]
PSD_RANGE = (1.0, 45.0)      # range shown in panel D  [KIM]
FOOOF_RANGE = (2.0, 40.0)    # fitting range  [FOOOF simulations, 2-40 Hz]
FOOOF_SETTINGS = dict(       # settings used throughout [FOOOF] ("except where noted")
    peak_width_limits=[1, 8],
    max_n_peaks=6,
    min_peak_height=0.1,
    peak_threshold=2.0,
    aperiodic_mode="fixed",  # no knee; slope = -exponent  [FOOOF eq. 3]
    verbose=False,
)
