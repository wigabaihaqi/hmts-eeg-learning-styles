"""
Check raw OpenBCI recordings and build Tier 2 (filtered 1-s windows) in one step.

    python preprocess.py              # all recordings listed in config.py
    python preprocess.py --check-only # QC report only, no Tier 2 files

For every recording it
  1. checks the raw file: sampling continuity (board sample counter), burst-written
     timestamps, ADC gain per channel, bad channels (PREP), amplitude;
  2. builds the time axis (config.TIME_MODEL) and cuts the 1-s windows exactly as
     the original HMT-S code did (filter whole recording, 8-30 Hz, then cut);
  3. adds per-window quality columns to Tier 2:
        n_real_samples   raw samples found in the window (before any interpolation)
        n_interpolated   samples in the window that are interpolated
        faster_outlier   window is an outlier by the FASTER epoch criteria
  4. writes output/qc_report.csv with a status (OK / WARN / FAIL) and reasons.

Outputs (config.OUT_DIR):
    <Tier 2 sub-folder>/<class>-NN-subjectN-preprocessed.csv
    qc_report.csv, qc_channels.csv
"""
import argparse
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, sosfiltfilt

import settings as C

warnings.filterwarnings("ignore")
DP = [f"datapoint{i}" for i in range(1, C.WIN_SAMPLES + 1)]


# ----------------------------------------------------------------- input
def find(folder, pattern):
    hits = sorted(p for p in folder.iterdir() if p.is_file() and re.fullmatch(pattern, p.name, flags=re.I))
    return hits[0] if hits else None


def recordings():
    out = []
    for cls, (raw_sub, tier2_sub, label) in C.CLASSES.items():
        base = C.RAW_DIR / raw_sub
        if not base.exists():
            continue
        for sd in base.iterdir():
            m = re.fullmatch(r"subject\s*(\d+)", sd.name, flags=re.I)
            if not (sd.is_dir() and m):
                continue
            num = int(m.group(1))
            if C.ONLY and (cls, num) not in C.ONLY:
                continue
            out.append(dict(cls=cls, num=num, label=label, tier2_sub=tier2_sub,
                            raw=find(sd, C.RAW_PATTERN), log=find(sd, C.LOG_PATTERN)))
    return sorted(out, key=lambda r: (r["tier2_sub"], list(C.CLASSES).index(r["cls"]), r["num"]))


def load_raw(path):
    r = pd.read_csv(path, skipinitialspace=True, low_memory=False)
    r.columns = [c.strip() for c in r.columns]
    t = pd.to_datetime(r["Timestamp (Formatted)"].astype(str).str.strip(), errors="coerce")
    t = (t.values.astype("datetime64[ns]").astype("int64") // 1_000_000).astype(float)   # ms, exact
    x = r[C.EXG_COLS].to_numpy(dtype=float)
    counter = r["Sample Index"].to_numpy(dtype=float).astype(int)
    return x, t, counter


def window_bounds(session_start):
    """Original rule: windows of WINDOW_SEC from the session start, next window WINDOW_GAP_MS later.
    All times in integer milliseconds so that boundary samples are included exactly as before."""
    rows, t = [], session_start
    for item in range(1, C.N_ITEMS + 1):
        for phase in range(1, C.N_PHASES + 1):
            for w in range(1, C.WIN_PER_PHASE + 1):
                rows.append((item, phase, w, t, t + C.WINDOW_SEC * 1000))
                t = t + C.WINDOW_SEC * 1000 + C.WINDOW_GAP_MS
    return pd.DataFrame(rows, columns=["stimulant", "session", "window", "t0", "t1"])


# ----------------------------------------------------------------- time axis
def counter_steps(counter):
    step = np.diff(counter) % 256
    fwd = (step >= 2) & (step <= C.MAX_FORWARD_STEP)
    back = step > C.MAX_FORWARD_STEP
    return step, fwd, back


def sample_time_axis(x, t, counter, task):
    """Insert missing samples (linear interpolation) and return (x_full, t_full, is_filled)."""
    step, fwd, back = counter_steps(counter)
    step_n = np.where(back, 1, step)
    n = np.concatenate([[0], np.cumsum(step_n)]).astype(float)
    full = np.arange(n[-1] + 1)
    filled = ~np.isin(full, n)
    if filled.any():
        x = np.column_stack([np.interp(full, n, x[:, c]) for c in range(x.shape[1])])
        t = np.interp(full, n, t)
    fit = (t >= task[0]) & (t <= task[1])
    slope, _ = np.polyfit(full[fit], t[fit], 1)
    offset = np.percentile(t[fit] - slope * full[fit], C.ENVELOPE_PERCENTILE)
    return x, offset + slope * full, filled, 1000 / slope


def bandpass(x):
    b, a = butter(C.FILTER_ORDER, list(C.BAND), btype="band", fs=C.FS)
    return filtfilt(b, a, x, axis=0)


def cut(xf, t, filled, t0, t1):
    lo = np.searchsorted(t, t0, side="left")
    hi = np.searchsorted(t, t1, side="right")
    seg = xf[lo:hi][:C.WIN_SAMPLES]
    fil = filled[lo:hi][:C.WIN_SAMPLES]
    n = len(seg)
    n_real = int(n - fil.sum())
    if n == C.WIN_SAMPLES:
        return seg, n_real, int(fil.sum())
    if n == 0:
        return np.zeros((C.WIN_SAMPLES, xf.shape[1])), 0, C.WIN_SAMPLES
    xi = np.linspace(0, n - 1, C.WIN_SAMPLES)
    seg = np.column_stack([np.interp(xi, np.arange(n), seg[:, c]) for c in range(seg.shape[1])])
    return seg, n_real, C.WIN_SAMPLES - n_real


# ----------------------------------------------------------------- QC
def robust_sd(v, axis=0):
    q75, q25 = np.percentile(v, [75, 25], axis=axis)
    return 0.7413 * (q75 - q25)


def robust_z(v):
    s = robust_sd(v)
    return np.zeros_like(v) if s == 0 else (v - np.median(v)) / s


def prep_flags(x):
    """PREP bad-channel criteria (deviation, HF noise, correlation, low SNR, NaN, constant)."""
    nan = np.isnan(x).any(axis=0)
    x = np.nan_to_num(x)
    xh = sosfiltfilt(butter(4, C.PREP_HIGHPASS_HZ, "high", fs=C.FS, output="sos"), x, axis=0)
    rsd = robust_sd(xh)
    low = sosfiltfilt(butter(4, C.PREP_LOWPASS_HZ, "low", fs=C.FS, output="sos"), xh, axis=0)
    mad = lambda v: np.median(np.abs(v - np.median(v, axis=0)), axis=0)
    hf = mad(xh - low) / np.maximum(mad(low), 1e-12)
    w = int(C.FS)
    k = low.shape[0] // w
    mc = np.zeros((k, x.shape[1]))
    for i in range(k):
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.abs(np.corrcoef(low[i * w:(i + 1) * w].T))
        np.fill_diagonal(r, np.nan)
        mc[i] = np.nan_to_num(np.nanpercentile(r, C.PREP_CORR_PERCENTILE, axis=1))
    bad_corr = (mc < C.PREP_CORR_THRESHOLD).mean(axis=0) > C.PREP_CORR_BAD_FRACTION
    bad_dev = robust_z(rsd) > C.PREP_Z_DEVIATION
    bad_hf = robust_z(hf) > C.PREP_Z_HF_NOISE
    flat = rsd < 1e-6
    df = pd.DataFrame(dict(channel=C.CHANNELS, robust_sd_uv=rsd, bad_deviation=bad_dev, bad_hf_noise=bad_hf,
                           bad_correlation=bad_corr, bad_low_snr=bad_corr & bad_hf, bad_nan=nan, bad_flat=flat))
    df["flagged"] = bad_dev | bad_hf | bad_corr | nan | flat
    return df


def gains(x_raw):
    out = []
    for c in range(x_raw.shape[1]):
        v = np.unique(x_raw[:, c]); d = np.diff(v); d = d[d > 0]
        step = d.min() if len(d) else np.nan
        out.append(min(C.ADS_GAINS, key=lambda g: abs(4.5 / g / (2 ** 23 - 1) * 1e6 - step)) if np.isfinite(step) else np.nan)
    return out


def faster_outliers(win):
    """win: (windows, channels, samples). FASTER epoch parameters 4-6, |z| > FASTER_Z."""
    chm = win.mean(axis=(0, 2))
    p = [(win.max(2) - win.min(2)).mean(1), np.abs(win.mean(2) - chm).mean(1), win.var(2).mean(1)]
    z = [(v - v.mean()) / v.std() if v.std() > 0 else np.zeros_like(v) for v in p]
    return np.any([np.abs(v) > C.FASTER_Z for v in z], axis=0)


# ----------------------------------------------------------------- one recording
def process(rec, write=True):
    sid = f"{rec['cls']}-S{rec['num']:02d}"
    q = dict(recording=sid, raw_file=rec["raw"].name if rec["raw"] else None,
             log_file=rec["log"].name if rec["log"] else None)
    reasons_fail, reasons_warn = [], []
    if rec["raw"] is None or rec["log"] is None:
        q.update(status="FAIL", reasons="raw or log file missing")
        return q, None
    x_raw, t, counter = load_raw(rec["raw"])
    log = pd.read_excel(rec["log"])
    log.columns = [str(c).strip() for c in log.columns]
    s0 = float(pd.to_datetime(log[C.LOG_SESSION_START].iloc[0]).value // 1_000_000)   # ms
    wb = window_bounds(s0)
    task = (wb["t0"].min() - C.TASK_MARGIN_S * 1000, wb["t1"].max() + C.TASK_MARGIN_S * 1000)
    in_task = (t >= task[0]) & (t <= task[1])
    if in_task.sum() < 10:
        q.update(status="FAIL", reasons="raw file does not cover the session in the log")
        return q, None

    # ---- continuity and timestamps
    step, fwd, back = counter_steps(counter)
    st_task = in_task[1:] & in_task[:-1]
    lost = int((step[fwd & st_task] - 1).sum())
    lost_pct = 100 * lost / max(st_task.sum() + lost, 1)
    span = (wb["t1"].max() - wb["t0"].min()) / 1000
    eff_rate = ((t >= wb["t0"].min()) & (t <= wb["t1"].max())).sum() / span
    rep_pct = 100 * float((np.diff(t)[st_task] == 0).mean())
    burst = rep_pct > C.BURST_PCT
    q.update(n_samples=len(t), effective_rate_hz=round(eff_rate, 2), lost_samples_task=lost,
             lost_pct_task=round(lost_pct, 3), backward_counter_steps_task=int(back[st_task].sum()),
             pct_repeated_timestamps=round(rep_pct, 2), burst_written=burst)
    if lost_pct > C.FAIL_LOST_PCT: reasons_fail.append(f"{lost_pct:.1f}% samples missing")
    elif lost_pct > C.WARN_LOST_PCT: reasons_warn.append(f"{lost_pct:.2f}% samples missing (interpolated)")
    if eff_rate < C.FAIL_RATE_FRACTION * C.FS: reasons_fail.append(f"effective rate {eff_rate:.1f} Hz")

    # ---- gain, channels, amplitude
    g = gains(x_raw)
    q.update(gains=" ".join(f"{ch}:x{int(v)}" for ch, v in zip(C.CHANNELS, g)))
    if any(v != C.EXPECTED_GAIN for v in g):
        reasons_warn.append("gain " + ",".join(f"{ch} x{int(v)}" for ch, v in zip(C.CHANNELS, g) if v != C.EXPECTED_GAIN))
    seg = x_raw[(t >= wb["t0"].min()) & (t <= wb["t1"].max())]
    chq = prep_flags(seg - seg.mean(axis=0))
    chq.insert(0, "recording", sid); chq["gain"] = g
    bad = chq.loc[chq["flagged"], "channel"].tolist()
    q.update(prep_flagged=" ".join(bad), n_prep_flagged=len(bad),
             median_robust_sd_uv=round(float(chq["robust_sd_uv"].median()), 1))
    if len(bad) >= C.WARN_BAD_CHANNELS: reasons_warn.append(f"PREP flagged {', '.join(bad)}")
    if chq["robust_sd_uv"].median() > C.WARN_ROBUST_SD_UV:
        reasons_warn.append(f"median robust SD {chq['robust_sd_uv'].median():.0f} uV")

    # ---- time axis and windows
    model = C.TIME_MODEL
    if model == "auto":
        model = "sample_counter" if burst else "timestamps"
    if model == "sample_counter" and (lost_pct > C.FAIL_LOST_PCT or back[st_task].any()):
        model = "timestamps"
        reasons_warn.append("sample-counter timing not possible (sample loss / counter reset); timestamps used")
    if model == "sample_counter":
        xs, ts, filled, rate = sample_time_axis(x_raw, t, counter, task)
        q["fitted_rate_hz"] = round(rate, 3)
    else:
        xs, ts, filled = x_raw, t, np.zeros(len(t), bool)
    q["time_model"] = model
    xf = bandpass(xs)
    wins, n_real, n_int = [], [], []
    for w in wb.itertuples():
        s, nr, ni = cut(xf, ts, filled, w.t0, w.t1)
        wins.append(s); n_real.append(nr); n_int.append(ni)
    wins = np.stack(wins).transpose(0, 2, 1)             # (window, channel, sample)
    n_int = np.array(n_int)
    out_flag = faster_outliers(wins)
    q.update(windows=len(wb), windows_interpolated=int((n_int > 0).sum()),
             pct_faster_outliers=round(100 * out_flag.mean(), 2))
    if (n_int == C.WIN_SAMPLES).any(): reasons_fail.append(f"{int((n_int == C.WIN_SAMPLES).sum())} empty windows")
    if 100 * out_flag.mean() > C.WARN_OUTLIER_WINDOWS_PCT:
        reasons_warn.append(f"{100 * out_flag.mean():.1f}% FASTER outlier windows")

    q["status"] = "FAIL" if reasons_fail else ("WARN" if reasons_warn else "OK")
    q["reasons"] = "; ".join(reasons_fail + reasons_warn)

    # ---- Tier 2 table (same columns as the released files + quality columns)
    if not write:
        return q, chq
    rows = []
    fmt = lambda v: pd.Timestamp(int(round(v)), unit="ms").strftime("%d-%m-%Y %H:%M:%S.%f")[:-3]   # as released
    for k, w in enumerate(wb.itertuples()):
        trial = int((w.stimulant - 1) * C.N_PHASES + w.session)
        nomor = int((rec["num"] - 1) * C.N_ITEMS * C.N_PHASES + trial)
        for c in range(len(C.CHANNELS)):
            rows.append([trial, f"subject{rec['num']}", fmt(w.t0), fmt(w.t1), int(w.stimulant), int(w.session), int(w.window),
                         c, nomor, *wins[k, c], rec["label"], n_real[k], int(n_int[k]), bool(out_flag[k])])
    cols = ["trial", "subject", "timestamp_start", "timestamp_end", "stimulant", "session", "window", "channel",
            "nomor", *DP, "class", "n_real_samples", "n_interpolated", "faster_outlier"]
    df = pd.DataFrame(rows, columns=cols)
    out = C.OUT_DIR / rec["tier2_sub"] / f"{rec['cls']}-{rec['num']:02d}-subject{rec['num']}-preprocessed.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not C.EXTRA_COLUMNS:
        df = df.drop(columns=["n_real_samples", "n_interpolated", "faster_outlier"])
    df.to_csv(out, index=False)
    q["tier2_file"] = str(out.relative_to(C.OUT_DIR))
    return q, chq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true", help="write the QC report only")
    a = ap.parse_args()
    recs = recordings()
    if not recs:
        sys.exit(f"No recordings found under {C.RAW_DIR}. Check config.py")
    C.OUT_DIR.mkdir(parents=True, exist_ok=True)
    qs, chs = [], []
    for rec in recs:
        print(f"  {rec['cls']} subject {rec['num']} ...", flush=True)
        q, ch = process(rec, write=not a.check_only)
        qs.append(q)
        if ch is not None: chs.append(ch)
    rep = pd.DataFrame(qs)
    rep.to_csv(C.OUT_DIR / "qc_report.csv", index=False)
    if chs: pd.concat(chs).to_csv(C.OUT_DIR / "qc_channels.csv", index=False)
    pd.set_option("display.width", 200, "display.max_colwidth", 90)
    print("\n" + rep[["recording", "status", "time_model", "reasons"]].to_string(index=False))
    print(f"\nOK {sum(rep.status == 'OK')}, WARN {sum(rep.status == 'WARN')}, FAIL {sum(rep.status == 'FAIL')}"
          f"  ->  {C.OUT_DIR / 'qc_report.csv'}")


if __name__ == "__main__":
    main()
