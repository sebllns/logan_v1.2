# claude/

Performance review and experiment kit for indexing the 2025 SRA (~40M samples) with
kmindex on TACC Vista. The index is built in chunks (one kmtricks build each) and
combined with `kmindex merge`; this folder gathers the models, the experiment plan and
the scripts needed to choose the chunking strategy.

Produced with Claude Code (see `method.md` for versions, sources and prompts).

## Disclaimer

This folder is reserved for automated work by Claude and must not be edited by hand.
Its content is not validated: formulas, citations and conclusions may be wrong. The
intended workflow is to use it as a source of information, review each item, and
integrate what holds (with corrections if needed) into the official documentation.

## Contents

| File | Purpose |
| --- | --- |
| `performances_review.md` | Main report. Verifies every formula of `performances.md` against kmindex 0.6.1 and kmtricks 1.5/1.6 (`file:line` citations), extracts the parameters driving open files, peak RAM, build/merge/query time and storage, defines experiments E0-E7 and confronts the results with the Logan Pareto merge calculator |
| `plan.md` | Approved plan the report was written from (structure, deliverables, checks) |
| `method.md` | Methodology: tool versions, code commits reviewed, process, original prompts |
| `scripts/probe_limits.sh` | E0: prints host limits as `key=value` (open files, `vm.max_map_count`, RAM, disks, quotas). `PROBE_DISK=1` adds a 1 GiB write test |
| `scripts/estimate.py` | Stand-alone resource model (no kmhelpers dependency): per-stage RAM, open files, temp disk, storage and time for one configuration, or `--sweep` over the number of chunks |
| `scripts/make_fofs.py` | Splits a sample list or directory into `chunk_<i>.fof` files in kmtricks format |
| `scripts/run_experiment.sh` | E1-E6 runner: builds chunks with `--from`, merges (flat or two-level), queries, and records wall time, peak RSS, peak fds, output size and temp-disk peak per step in `results.csv` |
| `scripts/plot_results.py` | Plots one metric against one parameter from `results.csv` |

`performances.md`, the document under review, lives outside this folder.

## Requirements

- `kmindex` >= 0.6.1 and `kmtricks` >= 1.6.0 in `$PATH` (runner only)
- Python 3, `matplotlib` for `plot_results.py`
- `/usr/bin/time` (GNU time) for the runner

`estimate.py` and `make_fofs.py` need only the Python standard library.

## Usage

Probe the limits of a compute node:

```bash
sbatch -p gg -N 1 -t 00:05:00 --wrap "bash scripts/probe_limits.sh" -o probe_%j.out
```

Estimate resources for a configuration, or sweep the number of chunks:

```bash
python3 scripts/estimate.py --kmers 5e8 --samples 5000 --total-samples 100000 \
    --threads 64 --ulimit 1048576 --ram 237e9 --fpr 0.25
python3 scripts/estimate.py --kmers 5e8 --total-samples 100000 --threads 64 \
    --ulimit 1048576 --ram 237e9 --sweep 1,2,4,8,16,32 --csv
```

Prepare chunk FOFs, run a build/merge/query experiment, plot:

```bash
python3 scripts/make_fofs.py /data/unitigs -c 8 -o fofs
sbatch -p gg -N 1 -t 08:00:00 -o e4_%j.out --wrap \
    "ulimit -S -n \$(ulimit -H -n); bash scripts/run_experiment.sh -f fofs -o e4_c8 -c 8 -t 64 -p 64 -b 300000000 -q query.fa -B 10 -T 8"
python3 scripts/plot_results.py e4_c8/results.csv -x chunks -y wall_s --step merge
```

Each script prints its options with `-h` (Python) or in its header comment (bash).

## Status

- `run_experiment.sh` is syntax-checked only: no kmindex/kmtricks binary was available
  on the review machine.
- `estimate.py` reproduces the kmhelpers `kmparams.py` values for the tested inputs.
- Formulas come from reading the code, not from measurements; the experiments in
  section 10 of the report turn them into constants. Open points are listed in section 12.
