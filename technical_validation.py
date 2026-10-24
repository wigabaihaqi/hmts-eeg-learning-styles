"""
Technical validation of the HMT-S EEG dataset.

Run from VS Code (or a terminal) inside this folder:
    python technical_validation.py

Criteria and their sources are listed in config.py ([PREP], [FASTER], [FOOOF], [KIM]).

Outputs (./results/):
    figure_technical_validation.png / .pdf   six-panel figure (A-F)
    summary.txt                              numbers ready for the manuscript
    qc_windows_per_subject.csv               (A) window retention, sampling continuity
    qc_channel_quality.csv                   (B, C) PREP statistics and flags per channel
    qc_channel_flag_rate.csv                 (C) % of recordings in which each channel is flagged
    psd_per_subject.csv / psd_grand_average.csv   (D) power spectra
    aperiodic_fooof.csv                      (E) aperiodic exponent per subject
    faster_window_flags.csv                  (F) FASTER window screening per subject
    behavior_per_item.csv / behavior_summary.csv  task performance (descriptive)
"""
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, welch

import config as C

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

try:
    from fooof import FOOOF
except ImportError:
    sys.exit("Package 'fooof' is missing. Run:  pip install -r requirements.txt")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def parse_ts(series):
    """Parse timestamps written as DD-MM-YYYY hh:mm:ss.fff or YYYY-MM-DD hh:mm:ss.fff."""
    s = series.astype(str).str.strip()
    out = pd.to_datetime(s, format="%d-%m-%Y %H:%M:%S.%f", errors="coerce")
    miss = out.isna()
    if miss.any():
        out[miss] = pd.to_datetime(s[miss], format="%Y-%m-%d %H:%M:%S.%f", errors="coerce")
    miss = out.isna()
    if miss.any():
        out[miss] = pd.to_datetime(s[miss], errors="coerce")
    return out


def find_one(folder, pattern):
    hits = [p for p in folder.iterdir() if p.is_file() and re.fullmatch(pattern, p.name, flags=re.I)]
    return sorted(hits)[0] if hits else None


def filt(x, cutoff, btype, fs=C.FS, order=4):
    sos = butter(order, cutoff, btype=btype, fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=0)


def robust_sd(x, axis=0):
    """PREP robust SD: 0.7413 x interquartile range."""
    q75, q25 = np.percentile(x, [75, 25], axis=axis)
    return C.PREP_ROBUST_SD_FACTOR * (q75 - q25)


def robust_z(v):
    """PREP robust z: (v - median) / robust SD, computed across channels."""
    v = np.asarray(v, dtype=float)
    s = robust_sd(v)
    if s == 0:
        return np.zeros_like(v)
    return (v - np.median(v)) / s


def mad(x, axis=0):
    return np.median(np.abs(x - np.median(x, axis=axis)), axis=axis)


def mean_sd(x, d=2):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return "n/a"
    return f"{x.mean():.{d}f} ± {x.std(ddof=1) if len(x) > 1 else 0:.{d}f}"


# ----------------------------------------------------------------------------
# inventory
# ----------------------------------------------------------------------------
def build_inventory():
    rows = []
    pat = re.compile(r"(active|reflective|verbal|visual)-(\d+)-subject(\d+)-preprocessed\.csv", re.I)
    for d in C.PREP_DIRS:
        if not d.exists():
            sys.exit(f"Folder not found: {d}\nCheck DATASET_DIR in config.py")
        for f in sorted(d.glob("*.csv")):
            m = pat.fullmatch(f.name)
            if not m:
                continue
            cls, num = m.group(1).lower(), int(m.group(3))
            subj_dir = C.RAW_DIR / C.CLASS_FOLDERS[cls] / f"Subject {num}"
            raw = find_one(subj_dir, r".*\.csv") if subj_dir.exists() else None
            log = find_one(subj_dir, r"ls-quiz.*\.xlsx") if subj_dir.exists() else None
            rows.append(dict(cls=cls, dim=C.CLASS_DIMENSION[cls], subject=num,
                             sid=f"{C.CLASS_DIMENSION[cls]}-S{num:02d}", prep=f, raw=raw, log=log))
    inv = pd.DataFrame(rows)
    if inv.empty:
        sys.exit("No preprocessed files found. Check DATASET_DIR in config.py")
    inv["cls"] = pd.Categorical(inv["cls"], C.CLASS_ORDER, ordered=True)
    return inv.sort_values(["dim", "cls", "subject"]).reset_index(drop=True)


# ----------------------------------------------------------------------------
# PREP bad-channel detection (deviation, correlation, HF noise, low SNR, NaN, flat)
# ----------------------------------------------------------------------------
def prep_channel_quality(x):
    nan = np.isnan(x).any(axis=0)
    x = np.nan_to_num(x)
    xh = filt(x, C.PREP_HIGHPASS_HZ, "high")

    # deviation
    rsd = robust_sd(xh)
    z_dev = robust_z(rsd)
    flat = rsd < C.FLAT_EPS_UV
    bad_dev = z_dev > C.PREP_Z_DEVIATION

    # HF noise: MAD(high) / MAD(low), low = < 50 Hz
    low = filt(xh, C.PREP_LOWPASS_HZ, "low")
    high = xh - low
    ratio = mad(high) / np.maximum(mad(low), 1e-12)
    z_hf = robust_z(ratio)
    bad_hf = z_hf > C.PREP_Z_HF_NOISE

    # correlation in 1-s windows on the < 50 Hz signal
    w = int(C.PREP_CORR_WIN_SEC * C.FS)
    nwin = low.shape[0] // w
    maxcorr = np.full((nwin, 8), np.nan)
    for k in range(nwin):
        seg = low[k * w:(k + 1) * w]
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.abs(np.corrcoef(seg.T))
        np.fill_diagonal(r, np.nan)
        maxcorr[k] = np.nanpercentile(r, C.PREP_CORR_PERCENTILE, axis=1)
    maxcorr = np.nan_to_num(maxcorr, nan=0.0)
    frac_low = (maxcorr < C.PREP_CORR_THRESHOLD).mean(axis=0)
    bad_corr = frac_low > C.PREP_CORR_BAD_FRACTION

    low_snr = bad_corr & bad_hf
    flagged = nan | flat | bad_dev | bad_hf | bad_corr
    return pd.DataFrame(dict(
        channel=C.CHANNELS, robust_sd_uv=rsd, z_deviation=z_dev, hf_noise_ratio=ratio, z_hf_noise=z_hf,
        median_max_corr=np.median(maxcorr, axis=0), pct_windows_low_corr=100 * frac_low,
        bad_nan=nan, bad_flat=flat, bad_deviation=bad_dev, bad_hf_noise=bad_hf,
        bad_correlation=bad_corr, bad_low_snr=low_snr, flagged=flagged)), xh


# ----------------------------------------------------------------------------
# FASTER epoch screening on Tier 2 windows (parameters 4, 5, 6)
# ----------------------------------------------------------------------------
def faster_windows(prep, dp_cols):
    key = ["stimulant", "session", "window"]
    prep = prep.sort_values(key + ["channel"])
    data = prep[dp_cols].to_numpy(dtype=float)
    n_ch = prep["channel"].nunique()
    data = data.reshape(-1, n_ch, data.shape[1])        # (window, channel, sample)
    ch_mean = data.mean(axis=(0, 2))                     # whole-channel mean
    p4 = (data.max(axis=2) - data.min(axis=2)).mean(axis=1)            # amplitude range
    p5 = np.abs(data.mean(axis=2) - ch_mean).mean(axis=1)               # deviation from channel average
    p6 = data.var(axis=2).mean(axis=1)                                  # variance
    z = lambda v: (v - v.mean()) / v.std() if v.std() > 0 else np.zeros_like(v)
    z4, z5, z6 = z(p4), z(p5), z(p6)
    bad = (np.abs(z4) > C.FASTER_Z) | (np.abs(z5) > C.FASTER_Z) | (np.abs(z6) > C.FASTER_Z)
    zero = np.all(data == 0, axis=(1, 2))
    return dict(n_windows=len(bad), pct_flagged=100 * bad.mean(),
                pct_range=100 * (np.abs(z4) > C.FASTER_Z).mean(),
                pct_deviation=100 * (np.abs(z5) > C.FASTER_Z).mean(),
                pct_variance=100 * (np.abs(z6) > C.FASTER_Z).mean(),
                zero_windows=int(zero.sum()))


# ----------------------------------------------------------------------------
# per-subject processing
# ----------------------------------------------------------------------------
def process_subject(r):
    out = {}
    dp_cols = [f"datapoint{i}" for i in range(1, C.SAMPLES_PER_WINDOW + 1)]

    # ---------------- Tier 2 ----------------
    prep = pd.read_csv(r.prep)
    rebuilt = C.HERE / "rebuilt_tier2" / r.prep.parent.name / r.prep.name
    out["tier2_source"] = "original"
    if C.USE_REBUILT_TIER2 and rebuilt.exists() and r.prep.name in REBUILT_FILES:
        prep_f = pd.read_csv(rebuilt)
        out["tier2_source"] = "rebuilt"
    else:
        prep_f = prep
    prep["t0"] = parse_ts(prep["timestamp_start"])
    prep["t1"] = parse_ts(prep["timestamp_end"])
    win = (prep.groupby(["stimulant", "session", "window"], as_index=False)
               .agg(t0=("t0", "first"), t1=("t1", "first"), n_ch=("channel", "nunique"))
               .sort_values("t0").reset_index(drop=True))
    out["faster"] = faster_windows(prep_f, dp_cols)

    # ---------------- Tier 1 ----------------
    raw = pd.read_csv(r.raw, skipinitialspace=True, low_memory=False)
    raw.columns = [c.strip() for c in raw.columns]
    ts = pd.to_datetime(raw["Timestamp (Formatted)"].astype(str).str.strip(), errors="coerce")
    tsv = ts.values.astype("datetime64[ns]")

    # (A) window retention, same rule as the preprocessing code
    lo = np.searchsorted(tsv, win["t0"].values.astype("datetime64[ns]"), side="left")
    hi = np.searchsorted(tsv, win["t1"].values.astype("datetime64[ns]"), side="right")
    n = hi - lo
    seg_mask = (ts >= win["t0"].min()) & (ts <= win["t1"].max())
    x = raw.loc[seg_mask, C.EXG_COLS].to_numpy(dtype=float)
    x = x - np.nanmean(x, axis=0)
    span = (win["t1"].max() - win["t0"].min()).total_seconds()
    gaps = ts[seg_mask].diff().dt.total_seconds().fillna(0)
    eff_rate = len(x) / span if span > 0 else np.nan
    out["retention"] = dict(
        windows_found=len(win), windows_expected=C.EXPECTED_WINDOWS,
        windows_complete=int((n >= C.SAMPLES_PER_WINDOW).sum()),
        windows_interpolated=int(((n > 0) & (n < C.SAMPLES_PER_WINDOW)).sum()),
        windows_empty=int((n == 0).sum()),
        median_samples_per_window=float(np.median(n)),
        channels_per_window_min=int(win["n_ch"].min()),
        task_span_s=span, task_samples=len(x), effective_rate_hz=eff_rate,
        n_gaps=int((gaps > C.GAP_SEC).sum()), max_gap_s=float(gaps.max()),
        irregular_sampling=bool(eff_rate < C.MIN_RATE_FRACTION * C.FS))

    # (B, C) PREP channel quality
    chq, xh = prep_channel_quality(x)

    # absolute amplitude descriptors (descriptive, no thresholds)
    xr = raw[C.EXG_COLS].to_numpy(dtype=float)
    gains, offs = [], []
    for c in range(8):
        v = np.unique(xr[:, c]); d = np.diff(v); d = d[d > 0]
        step = d.min() if len(d) else np.nan
        g = min(C.ADS_GAINS, key=lambda g: abs(C.ADS_STEP_UV(g) - step)) if np.isfinite(step) else np.nan
        gains.append(g)
    seg_raw = raw.loc[seg_mask, C.EXG_COLS].to_numpy(dtype=float)
    chq["detected_gain"] = gains
    chq["dc_offset_uv"] = np.median(seg_raw, axis=0)
    chq["p999_abs_sample_step_uv"] = np.percentile(np.abs(np.diff(seg_raw, axis=0)), 99.9, axis=0)
    out["chq"] = chq
    tsd = ts[seg_mask].diff().dt.total_seconds().dropna()
    out["retention"]["pct_repeated_timestamps"] = float(100 * (tsd == 0).mean())

    # (D) Welch PSD (2-s Hann, 50 % overlap)
    nper = int(C.WELCH_SEC * C.FS)
    f, pxx = welch(xh, fs=C.FS, window="hann", nperseg=nper, noverlap=nper // 2, axis=0)
    out["freqs"], out["psd"] = f, pxx
    good = ~chq["flagged"].values
    use = good if good.any() else np.ones(8, bool)
    mpsd = pxx[:, use].mean(axis=1)
    out["mean_psd"], out["n_good_channels"] = mpsd, int(good.sum())

    # (E) FOOOF, fixed aperiodic mode, 2-40 Hz
    fm = FOOOF(**C.FOOOF_SETTINGS)
    try:
        fm.fit(f, mpsd, list(C.FOOOF_RANGE))
        offset, exponent = fm.get_params("aperiodic_params")
        out["fooof"] = dict(exponent=exponent, slope=-exponent, offset=offset,
                            r_squared=fm.get_params("r_squared"), fit_error=fm.get_params("error"),
                            n_peaks=int(fm.n_peaks_))
    except Exception as e:                                   # noqa: BLE001
        print(f"    FOOOF failed for {r.sid}: {e}")
        out["fooof"] = dict(exponent=np.nan, slope=np.nan, offset=np.nan, r_squared=np.nan,
                            fit_error=np.nan, n_peaks=np.nan)

    # ---------------- behaviour ----------------
    beh = []
    if r.log is not None and not pd.isna(r.log):
        lg = pd.read_excel(r.log)
        lg.columns = [str(c).strip() for c in lg.columns]
        lg["Answered At"] = pd.to_datetime(lg["Answered At"])
        lg["logged_correct"] = lg["Correct?"].astype(str).str.strip().str.upper().isin(["TRUE", "1", "1.0"])
        lg["recomputed_correct"] = lg["Answer Index"] == lg["Correct Index"]
        opt_on = win[win["session"] == 2].groupby("stimulant")["t0"].min()
        for k, q in enumerate(sorted(lg["Question ID"].unique()), start=1):
            g_all = lg[lg["Question ID"] == q]
            g = g_all[g_all["Answer Index"].notna()].sort_values("Answered At")   # rows with an answer
            on = opt_on.get(k, pd.NaT)
            if g.empty:                                   # item shown but never answered
                beh.append(dict(sid=r.sid, cls=r.cls, dim=r.dim, item=k, question_id=q, n_selections=0,
                                final_answer=np.nan, correct_index=g_all["Correct Index"].iloc[0],
                                logged_correct=False, recomputed_correct=False, answered=False,
                                rt_first_s=np.nan, rt_final_s=np.nan, mismatch_rows=0))
                continue
            last = g.iloc[-1]
            beh.append(dict(sid=r.sid, cls=r.cls, dim=r.dim, item=k, question_id=q, n_selections=len(g), answered=True,
                            final_answer=last["Answer Index"], correct_index=last["Correct Index"],
                            logged_correct=bool(last["logged_correct"]),
                            recomputed_correct=bool(last["recomputed_correct"]),
                            rt_first_s=(g["Answered At"].iloc[0] - on).total_seconds() if pd.notna(on) else np.nan,
                            rt_final_s=(last["Answered At"] - on).total_seconds() if pd.notna(on) else np.nan,
                            mismatch_rows=int((g["logged_correct"] != g["recomputed_correct"]).sum())))
    out["beh"] = pd.DataFrame(beh)
    return out


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def add_rebuilt_retention(ret):
    """Window retention of the released Tier 2 (version 2): recordings rebuilt by
    rebuild_tier2.py have 180 complete windows with reconstructed timing."""
    ret = ret.copy()
    rep = C.HERE / "rebuilt_tier2" / "rebuild_report.csv"
    reb = set()
    if rep.exists():
        r = pd.read_csv(rep)
        reb = set(r.loc[r["action"] == "rebuilt", "sid"])
    is_reb = ret["sid"].isin(reb)
    ret["timing_reconstructed"] = is_reb
    ret["windows_reconstructed"] = np.where(is_reb, ret["windows_found"], 0)
    ret["windows_interpolated_v2"] = np.where(is_reb, 0, ret["windows_interpolated"])
    ret["windows_complete"] = np.where(is_reb, ret["windows_found"], ret["windows_complete"])
    return ret


def rebuilt_files():
    """File names that rebuild_tier2.py actually rebuilt in its latest run (from its report)."""
    rep = C.HERE / "rebuilt_tier2" / "rebuild_report.csv"
    if not rep.exists():
        return set()
    r = pd.read_csv(rep)
    return set(r.loc[r["action"] == "rebuilt", "file"])


REBUILT_FILES = set()


def main():
    global REBUILT_FILES
    REBUILT_FILES = rebuilt_files()
    C.RESULTS_DIR.mkdir(exist_ok=True)
    inv = build_inventory()
    print(f"Found {len(inv)} subjects in {C.DATASET_DIR}")
    missing = inv[inv["raw"].isna() | inv["log"].isna()]
    if len(missing):
        print("WARNING: raw or log file missing for:", ", ".join(missing["sid"]))

    ret, chq, fas, psd_rows, aper, behs, mpsds = [], [], [], [], [], [], []
    freqs = None
    for _, r in inv.iterrows():
        if r.raw is None or pd.isna(r.raw):
            continue
        print(f"  processing {r.sid} ({r.cls}) ...", flush=True)
        o = process_subject(r)
        tag = dict(sid=r.sid, cls=r.cls, dim=r.dim)
        irregular = o["retention"]["irregular_sampling"]
        ret.append({**tag, **o["retention"]})
        chq.append(o["chq"].assign(**tag))
        fas.append({**tag, **o["faster"], "tier2_source": o["tier2_source"]})
        aper.append({**tag, **o["fooof"], "n_good_channels": o["n_good_channels"], "irregular_sampling": irregular})
        freqs = o["freqs"]
        if not (irregular and C.EXCLUDE_IRREGULAR_FROM_SPECTRA):
            mpsds.append(o["mean_psd"])
        for ci, ch in enumerate(C.CHANNELS):
            psd_rows.append(pd.DataFrame(dict(sid=r.sid, channel=ch, freq_hz=freqs, psd_uv2_hz=o["psd"][:, ci])))
        behs.append(o["beh"])

    ret = pd.DataFrame(ret); chq = pd.concat(chq); fas = pd.DataFrame(fas); aper = pd.DataFrame(aper)
    beh = pd.concat(behs) if behs else pd.DataFrame()
    lp = np.log10(np.array(mpsds))
    aper_used = aper[~aper["irregular_sampling"]] if C.EXCLUDE_IRREGULAR_FROM_SPECTRA else aper

    # ---------------- tables ----------------
    R = C.RESULTS_DIR
    ret.to_csv(R / "qc_windows_per_subject.csv", index=False)
    chq.to_csv(R / "qc_channel_quality.csv", index=False)
    flag_rate = (chq.groupby("channel", sort=False)["flagged"].mean() * 100).reindex(C.CHANNELS)
    flag_rate.rename("pct_recordings_flagged").to_csv(R / "qc_channel_flag_rate.csv")
    pd.concat(psd_rows).to_csv(R / "psd_per_subject.csv", index=False)
    pd.DataFrame(dict(freq_hz=freqs, mean_log10_psd=lp.mean(0),
                      sem_log10_psd=lp.std(0, ddof=1) / np.sqrt(len(lp)) if len(lp) > 1 else 0)
                 ).to_csv(R / "psd_grand_average.csv", index=False)
    aper.to_csv(R / "aperiodic_fooof.csv", index=False)
    fas.to_csv(R / "faster_window_flags.csv", index=False)
    if len(beh):
        beh.to_csv(R / "behavior_per_item.csv", index=False)
        (beh.groupby(["cls", "item"], observed=True)
            .agg(n=("sid", "nunique"), accuracy_pct=("recomputed_correct", lambda v: 100 * v.mean()),
                 rt_final_mean_s=("rt_final_s", "mean"), rt_final_sd_s=("rt_final_s", "std"))
            .reset_index().to_csv(R / "behavior_summary.csv", index=False))

    ret_v2 = add_rebuilt_retention(ret)
    ret_v2.to_csv(R / "qc_windows_per_subject.csv", index=False)
    psd_mu = lp.mean(0)
    psd_se = lp.std(0, ddof=1) / np.sqrt(len(lp)) if len(lp) > 1 else np.zeros_like(psd_mu)
    make_figure(ret_v2, chq, flag_rate, freqs, psd_mu, psd_se, aper_used, fas)

    # ---------------- summary text ----------------
    per_subj = chq.groupby("sid")["flagged"].sum()
    L = []
    L.append(f"Subjects analysed: {len(ret)}")
    L.append("")
    L.append("[A] Window retention (rule of the dataset's preprocessing code)")
    L.append(f"  windows per subject: expected {C.EXPECTED_WINDOWS}; found {ret['windows_found'].min()}-{ret['windows_found'].max()}")
    L.append(f"  complete (>= 250 raw samples): {ret['windows_complete'].sum()} of {ret['windows_found'].sum()} "
             f"({100 * ret['windows_complete'].sum() / ret['windows_found'].sum():.2f} %)")
    L.append(f"  interpolated (1-249 samples): {ret['windows_interpolated'].sum()}; subjects affected: {(ret['windows_interpolated'] > 0).sum()}")
    L.append(f"  empty (zero-filled): {ret['windows_empty'].sum()}; subjects affected: {(ret['windows_empty'] > 0).sum()}")
    rv2 = add_rebuilt_retention(ret)
    L.append(f"  after timing reconstruction (Tier 2 version 2): {int(rv2['timing_reconstructed'].sum())} recordings rebuilt; "
             f"interpolated windows {int(rv2['windows_interpolated_v2'].sum())} in "
             f"{int((rv2['windows_interpolated_v2'] > 0).sum())} recordings")
    L.append(f"  effective sampling rate in task span (Hz): {mean_sd(ret['effective_rate_hz'], 1)} "
             f"(range {ret['effective_rate_hz'].min():.1f}-{ret['effective_rate_hz'].max():.1f})")
    L.append(f"  timestamp gaps > {C.GAP_SEC} s: {ret['n_gaps'].sum()}; subjects affected: {(ret['n_gaps'] > 0).sum()}")
    irr = ret[ret["irregular_sampling"]]
    L.append(f"  irregular sampling (< {C.MIN_RATE_FRACTION:.0%} of {C.FS:.0f} Hz): "
             + (", ".join(f"{a} ({b:.1f} Hz)" for a, b in zip(irr['sid'], irr['effective_rate_hz'])) if len(irr) else "none")
             + ("  -> excluded from D and E" if len(irr) and C.EXCLUDE_IRREGULAR_FROM_SPECTRA else ""))
    L.append("")
    L.append("[B, C] Bad channels, PREP criteria (Bigdely-Shamlo et al., 2015)")
    L.append(f"  flagged channels per 8-channel recording: {mean_sd(per_subj)} (range {per_subj.min()}-{per_subj.max()})")
    for k in ["bad_nan", "bad_flat", "bad_deviation", "bad_hf_noise", "bad_correlation", "bad_low_snr"]:
        L.append(f"    {k}: {int(chq[k].sum())} channel-recordings")
    L.append("  % of recordings in which each channel is flagged: " +
             ", ".join(f"{ch} {v:.1f}%" for ch, v in flag_rate.items()))
    L.append("")
    L.append("[Amplitude] descriptive, raw task segment")
    q1, q2, q3 = np.percentile(chq["robust_sd_uv"], [25, 50, 75])
    L.append(f"  robust SD after 1-Hz high-pass (µV): median {q2:.1f} [IQR {q1:.1f}-{q3:.1f}]")
    L.append(f"  channel-recordings with robust SD > 100 µV: {(chq['robust_sd_uv'] > 100).sum()} of {len(chq)}")
    q1, q2, q3 = np.percentile(chq["p999_abs_sample_step_uv"], [25, 50, 75])
    L.append(f"  99.9th percentile of |sample-to-sample step| (µV): median {q2:.1f} [IQR {q1:.1f}-{q3:.1f}]")
    gc = chq["detected_gain"].value_counts().sort_index()
    L.append("  ADC gain inferred from quantisation step: " + ", ".join(f"x{int(g)}: {n} channel-recordings" for g, n in gc.items()))
    mixed = chq.groupby("sid")["detected_gain"].nunique()
    L.append(f"  recordings with mixed gains across channels: {(mixed > 1).sum()} of {len(mixed)}")
    L.append(f"  recordings with > 50 % repeated timestamps (burst-written): {(ret['pct_repeated_timestamps'] > 50).sum()}: "
             + ", ".join(ret.loc[ret['pct_repeated_timestamps'] > 50, 'sid']))
    L.append("")
    L.append(f"[D] Grand-average PSD ({len(lp)} recordings, non-flagged channels)")
    ia = (freqs >= 8) & (freqs <= 13)
    L.append(f"  maximum within 8-13 Hz at {freqs[ia][np.argmax(lp.mean(0)[ia])]:.1f} Hz")
    L.append("")
    L.append(f"[E] Aperiodic component, FOOOF fixed mode {C.FOOOF_RANGE[0]:.0f}-{C.FOOOF_RANGE[1]:.0f} Hz (Donoghue et al., 2020)")
    L.append(f"  exponent: {mean_sd(aper_used['exponent'])} (range {aper_used['exponent'].min():.2f} to {aper_used['exponent'].max():.2f})")
    L.append(f"  equivalent log-log slope: {mean_sd(aper_used['slope'])}")
    L.append(f"  fit R^2: {mean_sd(aper_used['r_squared'], 3)}; peaks per spectrum: {mean_sd(aper_used['n_peaks'], 1)}")
    for cls in C.CLASS_ORDER:
        s = aper_used[aper_used["cls"] == cls]["exponent"]
        if len(s):
            L.append(f"    {cls}: exponent {mean_sd(s)}")
    L.append("")
    L.append(f"[F] Tier 2 windows flagged by FASTER epoch criteria, |z| > {C.FASTER_Z:.0f} (Nolan et al., 2010)"
             f" - Tier 2 source: {', '.join(f'{k} {v}' for k, v in fas['tier2_source'].value_counts().items())}")
    L.append(f"  per subject: {mean_sd(fas['pct_flagged'])} % (range {fas['pct_flagged'].min():.1f}-{fas['pct_flagged'].max():.1f})")
    L.append(f"    amplitude range: {mean_sd(fas['pct_range'])} %, deviation: {mean_sd(fas['pct_deviation'])} %, "
             f"variance: {mean_sd(fas['pct_variance'])} %")
    L.append(f"  all-zero windows in Tier 2: {fas['zero_windows'].sum()}")
    if len(beh):
        L.append("")
        L.append("[Behaviour] descriptive, final answer per item")
        L.append(f"  accuracy per subject: {mean_sd(beh.groupby('sid')['recomputed_correct'].mean() * 100)} %")
        L.append(f"  response time from options onset (s): {mean_sd(beh['rt_final_s'])}")
        L.append(f"  log rows where 'Correct?' disagrees with Answer Index == Correct Index: {int(beh['mismatch_rows'].sum())}")
        L.append(f"  items without any answer: {int((~beh['answered']).sum()) + len(ret) * C.N_ITEMS - len(beh)}"
                 f" ({', '.join(beh.loc[~beh['answered'], 'sid'] + ' item ' + beh.loc[~beh['answered'], 'item'].astype(str))})")
    txt = "\n".join(L)
    (R / "summary.txt").write_text(txt, encoding="utf-8")
    print("\n" + txt + f"\n\nResults written to {R}")


def make_figure(ret, chq, flag_rate, freqs, psd_mu, psd_se, aper, fas, out_stem=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"active": "#1f77b4", "reflective": "#ff7f0e", "verbal": "#2ca02c", "visual": "#d62728"}
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(2, 3, figsize=(12, 7))
    rng = np.random.default_rng(0)

    def by_class(ax, df, col, ylabel, title):
        for i, cls in enumerate(C.CLASS_ORDER):
            v = df[df["cls"] == cls][col].dropna().values
            if len(v):
                ax.boxplot([v], positions=[i], widths=0.5, showfliers=False)
                ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=12, color=colors[cls], zorder=3)
        ax.set_xticks(range(4)); ax.set_xticklabels([c.capitalize() for c in C.CLASS_ORDER])
        ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontweight="bold")

    ax = axs[0, 0]
    x = np.arange(len(ret))
    rec = ret["windows_reconstructed"] if "windows_reconstructed" in ret else np.zeros(len(ret))
    comp = ret["windows_complete"] - rec
    interp = ret["windows_interpolated_v2"] if "windows_interpolated_v2" in ret else ret["windows_interpolated"]
    col = [colors[c] for c in ret["cls"]]
    ax.bar(x, comp, color=col, label="complete")
    ax.bar(x, rec, bottom=comp, color=col, alpha=0.45, hatch="///", edgecolor="white", lw=0, label="complete (timing reconstructed)")
    ax.bar(x, interp, bottom=comp + rec, color="#bbbbbb", label="interpolated")
    ax.bar(x, ret["windows_empty"], bottom=comp + rec + interp, color="black", label="empty")
    ax.axhline(C.EXPECTED_WINDOWS, ls="--", c="k", lw=0.8)
    ax.set_ylim(0, C.EXPECTED_WINDOWS * 1.25)
    ax.set_xticks(x); ax.set_xticklabels(ret["sid"], rotation=90, fontsize=6)
    ax.set_ylabel("1-s windows per subject (Tier 2)"); ax.set_title("A  Window retention", loc="left", fontweight="bold")
    ax.legend(fontsize=6, frameon=False, loc="upper center", ncol=2)

    ps = chq.groupby(["sid", "cls"], observed=True)["flagged"].mean().mul(100).reset_index()
    by_class(axs[0, 1], ps, "flagged", "Flagged channels per recording (%)", "B  Bad-channel ratio (PREP)")

    ax = axs[0, 2]
    ax.bar(C.CHANNELS, flag_rate.values, color="#4c72b0")
    ax.set_ylabel("Recordings in which channel is flagged (%)")
    ax.set_title("C  Flag rate per channel", loc="left", fontweight="bold")

    ax = axs[1, 0]
    m = (freqs >= C.PSD_RANGE[0]) & (freqs <= C.PSD_RANGE[1])
    mu, se = np.asarray(psd_mu)[m], np.asarray(psd_se)[m]
    ax.plot(freqs[m], mu, c="k", lw=1.2)
    ax.fill_between(freqs[m], mu - se, mu + se, color="grey", alpha=0.35)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("log$_{10}$ PSD (µV²/Hz)")
    ax.set_title("D  Grand-average spectrum (task EEG)", loc="left", fontweight="bold")

    by_class(axs[1, 1], aper, "exponent", "Aperiodic exponent (FOOOF)", "E  Aperiodic exponent")
    by_class(axs[1, 2], fas, "pct_flagged", "Windows flagged (%)", "F  Tier 2 windows flagged (FASTER)")

    fig.tight_layout()
    stem = out_stem or (C.RESULTS_DIR / "figure_technical_validation")
    fig.savefig(f"{stem}.png", dpi=300)
    fig.savefig(f"{stem}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
