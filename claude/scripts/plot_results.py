#!/usr/bin/env python3
"""Plot metrics from run_experiment.sh CSV files.

  ./plot_results.py results.csv -x samples -y max_rss_kb --step build -o rss_vs_samples.png
  ./plot_results.py results.csv -x chunks -y wall_s --step merge

Columns come from run_experiment.sh: exp,step,name,samples,partitions,threads,
bf_size,chunks,batch,z,wall_s,max_rss_kb,max_fds,out_bytes,tmp_files_max,status
"""

import argparse
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path, step):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if step and r["step"] != step:
                continue
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-x", required=True)
    ap.add_argument("-y", required=True)
    ap.add_argument("--step", default="")
    ap.add_argument("--group", default="", help="column used to split series")
    ap.add_argument("--logx", action="store_true")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("-o", "--output", default="")
    a = ap.parse_args()

    rows = load(a.csv, a.step)
    if not rows:
        raise SystemExit("no rows")
    series = defaultdict(list)
    for r in rows:
        try:
            x, y = float(r[a.x]), float(r[a.y])
        except (KeyError, ValueError):
            continue
        series[r.get(a.group, "") if a.group else "all"].append((x, y))

    fig, ax = plt.subplots(figsize=(6, 4))
    for name, pts in sorted(series.items()):
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o",
                label=f"{a.group}={name}" if a.group else None)
    ax.set_xlabel(a.x)
    ax.set_ylabel(a.y)
    ax.set_title(f"{a.y} vs {a.x}" + (f" ({a.step})" if a.step else ""))
    if a.logx:
        ax.set_xscale("log")
    if a.logy:
        ax.set_yscale("log")
    if a.group:
        ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = a.output or f"{a.y}_vs_{a.x}{'_' + a.step if a.step else ''}.png"
    fig.savefig(out, dpi=120)
    print(out)


if __name__ == "__main__":
    main()
