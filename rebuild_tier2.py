"""
Rebuild Tier 2 (preprocessed 1-s windows) using sample-based timing.

Why
  In some HMT-S recordings the OpenBCI GUI wrote timestamps in bursts
  (one timestamp for a whole packet, every ~0.44 s). The original Tier 2
  code cut windows by timestamp, so those windows received fewer than
  250 samples and were interpolated, although no samples were missing.

What this script does (nothing in the dataset folder is changed)
  1. Self-check: rebuilds one recording with regular timestamps using the
     original rule and confirms that the result is identical to the
     published Tier 2 file.
  2. For every recording, unwraps the board's Sample Index counter
     (0-255) and checks for lost samples.
  3. Recordings whose original Tier 2 contains interpolated windows and
     whose sample counter shows no (or negligible) loss are rebuilt:
       - time of sample n = offset + n / rate, where rate is the
         least-squares slope of timestamp vs. n and offset places the
         line on the lower envelope of the timestamps (a timestamp can
         only be later than the true sampling time, never earlier);
       - the full channel is filtered exactly as in the original code
         (4th-order Butterworth 8-30 Hz, scipy filtfilt);
       - the same window boundaries as in the published file are used
         (timestamp_start <= t <= timestamp_end, first 250 samples).
     Only the datapoint1..datapoint250 columns change; all other columns
     are copied from the published file.
  4. Recordings with real sample loss are NOT rebuilt and are listed.

Output: ./rebuilt_tier2/<same sub-folder>/<same file name>
        ./rebuilt_tier2/rebuild_report.csv
Run:    python rebuild_tier2.py
"""
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

import config as C

warnings.filterwarnings("ignore")

OUT_DIR = C.HERE / "rebuilt_tier2"
BAND = (8.0, 30.0)            # as in the original preprocessing code
ORDER = 4                     # as in the original preprocessing code
MAX_LOST_FRACTION = 0.001     # rebuild only if < 0.1 % of the samples inside the task
                              # segment are missing (missing samples are linearly
                              # interpolated before filtering; gaps outside the task
                              # segment, e.g. at the start of the file, are tolerated)
TASK_MARGIN_S = 2.0
BURST_PCT = 50.0              # recording is "burst-written" if > 50 % of consecutive
                              # timestamps in the task segment are identical; only
                              # these recordings are rebuilt (in the HMT-S data the
                              # values are either near 0 % or above 60 %)
MAX_FORWARD_STEP = 128        # counter steps of 2..128 = lost samples; larger steps
                              # (= going backwards modulo 256) are counter resets or
                              # re-ordered packets, not loss           # task segment = first window start - 2 s .. last window end + 2 s
ENVELOPE_PERCENTILE = 0.5     # lower envelope of timestamp residuals (robust minimum)
DP = [f"datapoint{i}" for i in range(1, C.SAMPLES_PER_WINDOW + 1)]


def parse_ts(s):
    s = s.astype(str).str.strip()
    out = pd.to_datetime(s, format="%d-%m-%Y %H:%M:%S.%f", errors="coerce")
    m = out.isna()
    if m.any():
        out[m] = pd.to_datetime(s[m], format="%Y-%m-%d %H:%M:%S.%f", errors="coerce")
    m = out.isna()
    if m.any():
        out[m] = pd.to_datetime(s[m], errors="coerce")
    return out


def load_raw(path):
    r = pd.read_csv(path, skipinitialspace=True, low_memory=False)
    r.columns = [c.strip() for c in r.columns]
    ts = pd.to_datetime(r["Timestamp (Formatted)"].astype(str).str.strip(), errors="coerce")
    return r, ts.values.astype("datetime64[ns]").astype("int64") / 1e9   # seconds


def bandpass_full(x):
    b, a = butter(ORDER, list(BAND), btype="band", fs=C.FS)
    return filtfilt(b, a, x, axis=0)


def cut(xf, t, t0, t1):
    """Original rule: start <= t <= end, first 250 samples; interpolate if fewer; zeros if none."""
    lo = np.searchsorted(t, t0, side="left")
    hi = np.searchsorted(t, t1, side="right")
    seg = xf[lo:hi][:C.SAMPLES_PER_WINDOW]
    n = len(seg)
    if n == C.SAMPLES_PER_WINDOW:
        return seg, n
    if n == 0:
        return np.zeros((C.SAMPLES_PER_WINDOW, xf.shape[1])), 0
    xi = np.linspace(0, n - 1, C.SAMPLES_PER_WINDOW)
    return np.column_stack([np.interp(xi, np.arange(n), seg[:, c]) for c in range(seg.shape[1])]), n


def rebuild(prep, xf, t):
    """Return a copy of the Tier 2 frame with datapoints recomputed on time axis t."""
    out = prep.copy()
    t0 = parse_ts(prep["timestamp_start"]).values.astype("datetime64[ns]").astype("int64") / 1e9
    t1 = parse_ts(prep["timestamp_end"]).values.astype("datetime64[ns]").astype("int64") / 1e9
    vals = np.empty((len(prep), C.SAMPLES_PER_WINDOW))
    nsamp = np.empty(len(prep), int)
    cache = {}
    for i, (a, b, ch) in enumerate(zip(t0, t1, prep["channel"].astype(int))):
        if (a, b) not in cache:
            cache[(a, b)] = cut(xf, t, a, b)
        seg, n = cache[(a, b)]
        vals[i] = seg[:, ch]
        nsamp[i] = n
    out[DP] = vals
    return out, nsamp


def sample_counter(r):
    si = r["Sample Index"].to_numpy(dtype=float).astype(int)
    step = np.diff(si) % 256
    return step


def inventory():
    pat = re.compile(r"(active|reflective|verbal|visual)-(\d+)-subject(\d+)-preprocessed\.csv", re.I)
    rows = []
    for d in C.PREP_DIRS:
        for f in sorted(d.glob("*.csv")):
            m = pat.fullmatch(f.name)
            if not m:
                continue
            cls, num = m.group(1).lower(), int(m.group(3))
            sd = C.RAW_DIR / C.CLASS_FOLDERS[cls] / f"Subject {num}"
            raws = [p for p in sd.glob("*") if p.suffix.lower() == ".csv"] if sd.exists() else []
            rows.append(dict(sid=f"{C.CLASS_DIMENSION[cls]}-S{num:02d}", cls=cls, prep=f,
                             raw=sorted(raws)[0] if raws else None))
    return pd.DataFrame(rows)


def main():
    inv = inventory()
    if inv.empty:
        sys.exit("No Tier 2 files found. Check DATASET_DIR in config.py")
    # start from an empty output folder so that files from an earlier run do not remain
    if OUT_DIR.exists():
        for f in OUT_DIR.rglob("*.csv"):
            f.unlink()
    OUT_DIR.mkdir(exist_ok=True)
    report, selfcheck_done = [], False
    print(f"Found {len(inv)} recordings")

    for _, row in inv.iterrows():
        print(f"  {row.sid} ...", flush=True)
        prep = pd.read_csv(row.prep)
        r, t_orig = load_raw(row.raw)
        x = r[C.EXG_COLS].to_numpy(dtype=float)

        # original rule on original timestamps -> must equal the published file
        xf_orig = bandpass_full(x)
        repro, n_orig = rebuild(prep, xf_orig, t_orig)
        repro_err = float(np.abs(repro[DP].to_numpy() - prep[DP].to_numpy()).max())
        n_interp_orig = int((n_orig < C.SAMPLES_PER_WINDOW).sum() // 8)

        step = sample_counter(r)
        fwd = (step >= 2) & (step <= MAX_FORWARD_STEP)
        back = step > MAX_FORWARD_STEP
        lost = int((step[fwd] - 1).sum())
        ts0 = parse_ts(prep["timestamp_start"]).min().value / 1e9 - TASK_MARGIN_S
        ts1 = parse_ts(prep["timestamp_end"]).max().value / 1e9 + TASK_MARGIN_S
        in_task = (t_orig[1:] >= ts0) & (t_orig[:-1] <= ts1)
        st_task = step[in_task]
        lost_task = int((st_task[fwd[in_task]] - 1).sum())
        dup_task = int((st_task == 0).sum())
        back_task = int(back[in_task].sum())
        pct_rep = 100 * float((np.diff(t_orig)[in_task] == 0).mean())
        lost_frac = lost_task / max(in_task.sum() + lost_task, 1)
        rec = dict(sid=row.sid, file=row.prep.name, n_samples=len(r),
                   lost_samples_file=lost, lost_samples_task=lost_task,
                   duplicate_counter_steps_task=dup_task, backward_counter_steps_task=back_task,
                   lost_fraction_task=lost_frac, pct_repeated_timestamps_task=pct_rep, original_rule_reproduced_max_abs_diff_uv=repro_err,
                   windows_interpolated_original=n_interp_orig)

        if not selfcheck_done and n_interp_orig == 0:
            rec["note"] = "self-check: original Tier 2 reproduced" if repro_err < 1e-6 else "self-check FAILED"
            print(f"    self-check on {row.sid}: max |difference| = {repro_err:.2e} µV")
            selfcheck_done = True

        if n_interp_orig == 0:
            rec.update(action="kept (no interpolated windows)")
            report.append(rec)
            continue
        if pct_rep <= BURST_PCT:
            rec.update(action="kept (timestamps per sample; few interpolated windows)")
            report.append(rec)
            continue
        if lost_frac > MAX_LOST_FRACTION or dup_task > 0 or back_task > 0:
            rec.update(action="NOT rebuilt: sample loss in the recording")
            report.append(rec)
            continue

        # sample-based time axis
        step_n = np.where(back, 1, step)          # resets outside the task: no loss assumed
        n = np.concatenate([[0], np.cumsum(step_n)]).astype(float)
        if (step_n > 1).any():
            full_n = np.arange(n[-1] + 1)
            x = np.column_stack([np.interp(full_n, n, x[:, c]) for c in range(8)])
            tt = np.interp(full_n, n, t_orig)
            n, t_obs = full_n, tt
        else:
            t_obs = t_orig
        fit = (t_obs >= ts0) & (t_obs <= ts1)        # time model fitted on the task segment only
        slope, _ = np.polyfit(n[fit], t_obs[fit], 1)
        offset = np.percentile(t_obs[fit] - slope * n[fit], ENVELOPE_PERCENTILE)
        t_new = offset + slope * n

        xf = bandpass_full(x)
        new, n_new = rebuild(prep, xf, t_new)
        out_path = OUT_DIR / row.prep.parent.name / row.prep.name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        new.to_csv(out_path, index=False)
        r_corr = np.corrcoef(new[DP].to_numpy().ravel(), prep[DP].to_numpy().ravel())[0, 1]
        rec.update(action="rebuilt", fitted_rate_hz=1 / slope,
                   windows_interpolated_rebuilt=int((n_new < C.SAMPLES_PER_WINDOW).sum() // 8),
                   corr_with_original=r_corr)
        report.append(rec)

    rep = pd.DataFrame(report)
    rep.to_csv(OUT_DIR / "rebuild_report.csv", index=False)
    print("\n" + rep[["sid", "action", "windows_interpolated_original",
                      *[c for c in ["windows_interpolated_rebuilt", "fitted_rate_hz", "lost_samples_task"] if c in rep]]]
          .to_string(index=False))
    print(f"\nRebuilt files and report written to {OUT_DIR}")


if __name__ == "__main__":
    main()
