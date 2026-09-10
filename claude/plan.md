# Plan: verification report of `seb/performances.md` + experiment design for the SRA-2025 index

## Context

`seb/performances.md` collects heuristic formulas (kmhelpers `auto-params`, vendored
`pykmhelpers/vendor/kmparams.py`) that estimate peak RAM, open files, storage and time for
`kmindex build` / `query`. The team must index ~40M SRA samples on TACC Vista (GG nodes: 144
cores, 237 GB RAM, 286 GB local `/tmp`, 48 h queue limit, no documented ulimit). A single
kmtricks build cannot hold all samples (open-file hard limit in the kmtricks merge stage), so
the index will be built as chunks then combined with `kmindex merge`. The goal of the new
report is (1) to check every formula/assertion against the code, docs and papers, with
`file:line` citations, (2) to state clearly which parameters drive each metric, (3) to define
the experiments that decide the chunking strategy, and (4) to confront this with the Logan
Pareto merge calculator.

Exploration is done (kmindex agent report + inline reading of kmtricks, GATB, kmhelpers
`dev/simplify-options`, sibling repos, docs, Vista docs, calculator repo). No code changes are
needed in kmindex; the deliverable is a report plus helper scripts.

## Deliverables

1. `seb/performances_review.md` (English, same notation as `performances.md`).
2. `seb/scripts/probe_limits.sh` - prints the system limits to measure on a Vista compute node
   (`ulimit -H -n`, `ulimit -S -n`, `/proc/sys/fs/nr_open`, `file-max`, `vm.max_map_count`,
   nproc, `free -g`, `df` of `/tmp`, `$SCRATCH`, `$WORK`, inode quota, `uname -m`). Meant to be
   run via `sbatch -p gg --wrap`.
3. `seb/scripts/estimate.py` - stand-alone model (no kmhelpers dependency) that, given
   `K, S_chunk, S_total, P, T, bf_size, fpr, n_chunks, ulimit, ram`, prints per-stage
   estimates: RAM (count, kmtricks merge, kmindex merge, query), open files per stage,
   temp-file count, storage, and a merge/build time model with user-supplied disk throughput.
   Includes a `--sweep n_chunks` mode producing a CSV/plot for the chunk trade-off.
4. `seb/scripts/run_experiment.sh` - generic runner: builds N chunks with `kmindex build`
   (first chunk, then `--from`), merges them, queries; every step wrapped in
   `/usr/bin/time -v`, plus an fd sampler (`ls /proc/PID/fd | wc -l` every second) and a
   temp-file counter; appends one CSV line per step (params, wall, max RSS, max fds, size).
5. `seb/scripts/plot_results.py` - reads the CSV from (4) and plots metric vs parameter.

Scripts stay small and sober (no emojis, no em-dashes, concise comments).

## Report structure (each section starts with the questions it answers)

### 0. Scope, sources, versions
- kmindex 0.6.1 (`CMakeLists.txt:2`), requires kmtricks binary >= 1.3.0
  (`lib/include/kmindex/index/index.hpp:12`). The submodule `thirdparty/kmtricks` is pinned at
  v1.5.0; the user targets kmtricks 1.6.0. Verified: `git diff v1.5.0 v1.6.0` only adds BAM
  input support (BankBam, filtering options); the count, merge, memory and file-handling code
  cited below is identical, so citations are given against v1.5.0 lines and stated valid for
  1.6.0. Flag: the Vista install must use the same version pair (0.6.1 / 1.6.0).
- kmindex 0.6.x additions relevant to the report: `kmindex query2` (several global index
  paths queried together; sub-indexes are still processed sequentially, `app/kmindex/
  query2.cpp:187-245`), `kmindex compress` (zstd block-compressed matrices, `--fast` ignored
  for compressed indexes, `app/kmindex/query.cpp:324`), official ARM support (changelog
  v0.6.0, kmtricks PR #45 sse2neon). `--static-report` (changelog v0.6.1) not found in
  `build.cpp`; kmhelpers passes `--static-repart`: check the real flag name during writing.
- Papers: kmindex (Lemane et al., Nat. Comput. Sci. 4, 2024), kmtricks (Lemane et al.,
  Bioinf. Adv. 2022), findere (Robidou and Peterlongo 2021). Docs: tlemane.github.io/kmindex,
  Vista guide, calculator repo.
- Vista GG nodes are aarch64 (NVIDIA Grace): ARM is supported since kmindex 0.6.0, but the
  conda `tlemane` channel packages must be checked for aarch64 (else build from source with
  `conda/kmindex/build_local.sh`).

### 1. Pipeline recap and what `kmindex build` actually controls
- `kmindex build` forwards exactly: `--file --run-dir --kmer-size --hard-min --bloom-size
  --minimizer-size --nb-partitions --threads --mode hash:bf:bin` (+`--repart-from`, `--cpr`)
  (`app/kmindex/build.cpp:159-190`). `--max-memory` is NOT a kmtricks CLI option: it is a
  hard-coded 8000 MB (`thirdparty/kmtricks/include/kmtricks/cmd/all.hpp:60`) only used by GATB
  to auto-compute partitions; `--focus` default 0.5 (`src/cli.cpp:318-320`), not reachable
  from kmindex. Correct `performances.md` accordingly: M is not a kmtricks knob, it is a
  budget the user must enforce through P and T.
- Pipeline in the code: config -> repart -> superk+count interleaved (`task_scheduler.hpp:
  exec_superk_count`, `max_running = floor(T*focus)` at :265) -> merge -> (format only for
  BFT, `task_scheduler.hpp:530`). kmindex uses `hash:bf:bin`, so the kmtricks merge output
  (`write_as_bf`, `merge.hpp:575-597`) IS the index (`matrices/matrix_p.cmbf`).
- Auto partitions (`--nb-partitions 0`): GATB `ConfigurationAlgorithm.cpp:319-425`
  (`volume = kmersNb*8/MB`, `*0.6`, `P = volume/8000 + 1`, then rounded and capped by
  `getrlimit(RLIMIT_NOFILE)/2/3`), min 4 (`task.hpp:75-76`); bloom window rounded to 64 bits
  per partition (`hash.hpp:31-38`). Recommend always passing P explicitly.

### 2. Build: open files (hard-limited metric)
Questions: what is the exact per-stage fd count, what is the binding stage, what is the hard
limit on Vista.
- superk: each SuperK task opens P files (`io/superk_storage.hpp:174-185`), at most
  `floor(T*focus)` tasks in flight; count tasks open 1 reader + 1 writer. So
  `superk = floor(T*f)*P + 2*(T - floor(T*f))` (kmhelpers' `T*P + n_w` is a safe upper bound).
- kmtricks merge: one `HashReader` per sample + 1 output per partition task
  (`merge.hpp:404-409`, `task.hpp:695-701`), P tasks over a T-thread pool
  (`task_scheduler.hpp:exec_merge`) -> `min(T,P)*(S+1)`. Pierre's formula is the exact one;
  kmhelpers `T*(S+1)` is the safe bound (equal when T <= P). Confirm this is the binding stage.
- kmindex merge: `T_eff*(n_chunks+1)` with `T_eff = min(T, cores, P)` (`lib/src/index/
  merge.cpp:145-159, 239-252`; `lib/src/threadpool.cpp:5-12`). No rlimit handling anywhere.
- GATB itself reads the SOFT limit (`FileSystemLinux.cpp:71-79`); kmhelpers `get_max_open_files`
  also uses the soft limit (`resources.py`), `maximize_nofile()` raises soft to hard. Vista:
  the docs say nothing; local reference values (this laptop): hard 1073741816, soft 1048576,
  systemd default 524288. Must be measured with `probe_limits.sh` (login node vs compute node).
- Extra hard-limited resource missed by `performances.md`: temp file COUNT: superk + count
  files = 2*S*P per chunk (plus deletion) and inode quotas ($HOME 500k files); also
  `vm.max_map_count` for kmindex merge/query mmaps.

### 3. Build: peak RAM
Questions: which stage peaks, formula per stage, role of K, P, T, S.
- count: `nbk*8 + 8192` bytes per task (`utils.hpp:125-128`, `task.hpp:355-368`) with nbk =
  k-mers (with multiplicity, from super-k-mers) of one sample in one partition. Peak ~
  `T * max_{s,p} nbk * 8` = the `K*8*1.05*T/P` formula; verified, with caveats: K is
  occurrences not distinct k-mers (fine for Logan unitigs, wrong for raw reads), partition
  imbalance (GATB assumes 1.2; kmhelpers uses 1.05), superk tasks add GATB partition cache
  (to check in `gatb/fill_partitions.hpp`) and P*32 KB buffers.
- kmtricks merge stage (missing from `performances.md`): each HashReader holds ~160 KB of
  buffers (`hash_file.hpp:222-228` + `io_common.hpp:229`, buf_size 32768) so
  `RAM_merge ~ min(T,P) * S_chunk * 0.16 MB` (+lz4 buffers with `--cpr`). For S=10k, T=32:
  ~50 GB. This makes RAM depend on S and can dominate the count stage for big chunks.
- kmindex merge: anonymous RAM negligible (`P * ceil(S_total/8)` bytes, `merge.cpp:199`);
  RSS is page cache from `T_eff * n_chunks` mmaps.
- Update the build arrow table: `nb_samples` -> `↑` on peak RAM (merge stage), `nb_partitions`
  `↓` (count) and `·` (merge stage), threads `↑↑` both.

### 4. Build time
Questions: which volumes drive each stage, how to extrapolate from small runs.
- Model: `t_build ~ a*input_bytes (superk, parse) + b*total_kmers (count+merge CPU) +
  c*S*bf_size/8 (merge writes the full dense matrix: `write_as_bf` writes empty rows too)`.
  Merge output volume equals index size, independent of k-mer content. Threads: superk/count
  scale with T up to focus limits; merge scales with min(T,P).
- Reference points to confront: Logan/Genouest campaign 0.13328 s/sample at 40 threads (k=25,
  1.58x speedup 16->40 threads), Azure Logan runs (`building_logan_search`: P=256, 16-32
  threads, `ulimit -n 98304`), calculator estimate ~11 days per full-corpus configuration.

### 5. `kmindex merge` cost (the chunk trade-off)
Questions: is merge IO- or CPU-bound, does cost depend on n_chunks, multi-level merge cost.
- Algorithm (`lib/src/index/merge.cpp:143-237`): per partition, mmap all chunk files, memcpy
  the first chunk's row, then bit-by-bit copy for every other sample
  (`BITCHECK/BITSET`). Cost = read(total size) + write(total size) + CPU
  `bf_size * (S_total - S_first)` bit ops. Hence merge time is ~independent of n_chunks
  (except the first chunk) and roughly `2 * index_size / disk_throughput` when IO-bound;
  every extra merge level re-reads and re-writes the whole index: prefer one flat merge,
  bounded by `T_eff*(n_chunks+1) <= ulimit`.
- Output size <= sum of inputs (row padding coalesced) exactly
  `P*49 + bf_size*ceil(S*bw/8)` (verified on `tests/data/indexes`). Disk peak = inputs +
  output unless `-d` (per-partition deletion, `merge.cpp:236`).
- Constraints: identical sha1 over bloom size, P, k, minim size, bw and the full repartition
  table (`index_infos.cpp:253-281`) -> build chunks with `--from`; unique sample ids.

### 6. Query: peak RAM and time
Questions: does RAM depend on index size, on S, on P, on n sub-indexes, on threads/batch.
- Index is mmapped, never read (`kindex.cpp:10-27`), so RAM is not the index size.
- Per sequence response buffer `(len - s + 1) * ceil(S*bw/8)` bytes (`query.cpp:6-10`),
  s-mer list 24 B per s-mer, `-b 0` (default) keeps the whole query file's responses live
  (`app/kmindex/query.cpp:329`); RAM ~ `T * b * L * S/8`. With S=40M one 1 kb sequence needs
  ~5 GB: the merged sub-index size S is the dominant query-RAM driver, so "less chunks =
  bigger index = more RAM" is TRUE for query, and `-b` must be set. `z` adds CPU only
  (`query_results.cpp:44-54`), P is negligible for RAM.
- Sub-indexes are processed sequentially, each re-reads the fastx and re-maps partitions per
  batch (`query.cpp:276-346`, `kindex.hpp:58-71`); `--fast` keeps P files mapped. Time ~
  `sum over sub-indexes [fixed cost + partitions touched * (mmap + page faults) +
  n_smers * ceil(S/8) bytes]`. This matches the calculator's score
  `sum_g min(P_g, M(L))` and the docs' "optimal query times" after merge.
- Table of query drivers (S, n_smers, b, T, P, z, format json_vec, registered sub-indexes
  via repartition tables `2^(2m)*2 B`).

### 7. Storage
- Verify the bits/k-mer table and `f = -ln(p)/ln(2)^2` against kmhelpers `bloom_filter.py`
  (`f_value`, `get_bf_size` rounding to 8, `bf_max_kmers`) and the 64-bit per-partition
  rounding in kmtricks; note 1 hash function (docs). Exact size formula above; P adds only
  49 B per partition. Index size is the same whether chunked or merged (minus padding).
- Mention `kmindex compress` (0.6.0) as a post-merge storage lever with its query-side cost
  (no `--fast`, block decompression per lookup); out of scope for the formulas, listed as an
  optional experiment (E7) only.

### 8. Verification table of `performances.md`
One row per formula/assertion: status (exact / safe bound / wrong / unverifiable), correction,
source. Includes the `get_best_params` procedure (checked against
`pykmhelpers/core/build_params.py` on `dev/simplify-options`) and the note that
`safety_margin` scales the soft ulimit, not the hard one.

### 9. System limits on Vista
Facts from the Vista guide (nodes, RAM, /tmp, queues, 48 h, SU cost, VAST scratch purge,
no striping) and the list of values to measure with `probe_limits.sh`. Flag the 48 h wall
time as a hard limit on chunk size (build time per chunk) and the local 286 GB `/tmp` as a
possible temp-file target.

### 10. Experiment plan (one subsection per question, each with a parameter table)
Columns: parameter, values, fixed/varied, known/unknown, and expected driver. Small synthetic
or Logan subsets (kmhelpers `test create-db`, or a few thousand real unitig files) suffice
because the models are linear; each experiment states the minimum size needed.
- E0 probe limits (no data).
- E1 open files vs (T, P, S_chunk): `ulimit -n` lowered artificially, find failure threshold;
  validates `min(T,P)*(S+1)`.
- E2 build peak RAM vs T, P, K, S_chunk (isolate count stage vs kmtricks merge stage using
  `--until`-free runs with fd/RSS sampler over time).
- E3 build time vs S, bf_size, T; extract a, b, c coefficients for the time model.
- E4 kmindex merge: time and RSS vs n_chunks at fixed S_total, vs T, one-level vs two-level.
- E5 query: RAM and time vs S (merged vs not merged), n sub-indexes, P, T, b, z, `--fast`.
- E6 end-to-end strategy comparison: for a fixed dataset, n_chunks in {1, 2, 4, 8, 16}:
  total build time (parallel chunks on separate nodes counted as max, not sum), merge time,
  query time, query RAM, storage.
- Decision matrix: how each result feeds `n_chunks`, `T`, `P`, `-b` choices, plus the
  cluster-level parallelism argument (chunks on separate nodes, SU cost).

### 11. Confrontation with the Logan Pareto merge calculator
- What it models (query cost `sum_g min(P_g, M(L))`, P from storage `clamp(pow2(1 + bytes/4GB),
  4, 256)`, 16 measured configs, Spearman 0.84/0.73/0.53) and what it does not (build,
  merge, RAM, open files). Note that its partition formula is the old kmhelpers `compose`
  rule (now commented out in `pipeline/composer.py:229-249`) and differs from the RAM-driven
  `auto_params` partitions; agreement with sections 5-6; gaps to fill with E4-E6.

### 12. Open points / flags
Vista ulimit unknown; kmtricks version on Vista; K distinct vs occurrences; GATB superk cache
size; aarch64 build; `--cpr` interaction (docs say indexes must be built without `--cpr`,
`kmindex build` still forwards it); `index-infos` mergeable dead code
(`lib/src/index/index.cpp:124`); abundance-mode query buffers never freed (not used here).

## Files
- Create: `seb/performances_review.md`, `seb/scripts/{probe_limits.sh, estimate.py,
  run_experiment.sh, plot_results.py}`.
- Reuse (reference only, no edits): `pykmhelpers/vendor/kmparams.py` and
  `pykmhelpers/core/build_params.py` (`dev/simplify-options`) for the formulas;
  `pykmhelpers/core/bloom_filter.py` for bf sizing; `kmhelpers test create-db` for synthetic
  data; kmindex `tests/data/indexes` for size checks.
- Do not touch `seb/performances.md` (corrections go in the new report).

## Verification
- Every citation re-checked with `sed -n` on the referenced lines before writing.
- Numeric checks: size formula against `tests/data/indexes/*/matrices` (`ls -l`);
  `estimate.py` reproduces kmhelpers `auto-params` outputs for a few (K, S, ulimit, threads)
  inputs (`python -m pykmhelpers ... test auto-params` on `dev/simplify-options`, or the
  vendored `kmparams.py`).
- `bash -n` / `python -m py_compile` on the scripts; a dry run of `run_experiment.sh` on
  `tests/data` or a `kmhelpers test create-db` dataset if a kmindex/kmtricks binary is
  available locally (`which kmindex kmtricks`), otherwise documented as untested.
- Suggested commit message at the end (per user instructions).


The kmindex agent finished, and I've replaced the two killed subagents (kmtricks and kmhelpers) with inline exploration, now nearly complete—no subagents remain running. I'm finalizing a few kmtricks details (superk cache memory, hash reader buffers, format stage relevance to kmindex) before writing the plan.

Exploration is complete and the plan is written. Submitting it for approval.

Thanks. I'll re-check the versions now and re-verify the key code paths against kmindex 0.6.1 and kmtricks 1.6.0 before updating the plan.

kmindex is confirmed at 0.6.1 with kmtricks pinned to v1.5.0; the cited code paths still hold, just with shifted line numbers. Now checking the new 0.6.x features (query2, compressed indexes, --static-report) and the 1.5.0/1.6.0 diff before updating the plan.

All verifications are done, and local checks confirm no kmindex/kmtricks binaries are installed here (so the runner script will be syntax-checked only). Writing the scripts first, then the report.

Finishing up: the default kmtricks count type is uint32 (MAX_C=4294967295), so the hash-file record is 12 bytes per k-mer, not 9. Fixing that in the report and the model, then final checks.

Done. The review report and four helper scripts are in `seb/`.

**Deliverables**
- `seb/performances_review.md`: 13 sections, each opening with the questions it answers. Every formula of `performances.md` gets a verdict (`EXACT`, `BOUND`, `APPROX`, `WRONG`, `OPEN`) with `file:line` citations against kmindex 0.6.1 and kmtricks v1.5.0 (checked identical to 1.6.0 for the cited code, which only adds BAM input). Sections 10 and 11 hold the experiment plan (E0 to E7, one parameter table each) and the confrontation with the Logan Pareto calculator.
- `seb/scripts/probe_limits.sh`: prints all limits to measure on a Vista compute node (tested locally).
- `seb/scripts/estimate.py`: stand-alone resource model, single configuration or `--sweep` over the chunk count; reproduces the `kmparams.py` values.
- `seb/scripts/run_experiment.sh`: builds chunks with `--from`, merges (flat or two-level), queries, and records wall time, peak RSS, peak open fds, output size and temp-file peak per step in a CSV. Syntax-checked only, since no kmindex/kmtricks binary exists on this machine.
- `seb/scripts/plot_results.py`: metric vs parameter plots from that CSV (tested on a fake CSV).

**Main corrections to `performances.md`**
- `max_memory` is not a kmtricks option. It is a hard-coded 8000 MB used only for the automatic partition count, so $M$ is a user budget enforced through $P$ and $T$.
- The kmtricks merge stage allocates about 160 KB per open sample per merge task, so build peak RAM also grows with $S$: roughly $\min(T,P) \cdot S \cdot 160$ KB, which reaches a full Vista node around 20k samples at 64 threads. This term is absent from kmhelpers `auto-params`.
- Open files in the merge stage are exactly $\min(T,P)(S+1)$; the `T(S+1)` formula is a safe bound. Super-k-mer files are $\lfloor Tf \rfloor P + 2(T - \lfloor Tf \rfloor)$.
- Query peak RAM is linear in the sub-index sample count: $\min(N_{seq}, Tb) \cdot n \cdot \lceil S/8 \rceil$, and `-b 0` keeps the whole query file live. A merged 40M-sample index needs about 5 GB per 1 kb sequence.
- `kmindex merge` cost does not depend on the number of chunks (read plus write of the whole index, IO-bound by the model); a two-level merge doubles it.
- Temporary disk per chunk is about 12 bytes per k-mer occurrence (uncompressed hash files), far beyond the 286 GB local `/tmp` for large spans.
- The "p = 0.25" sizing uses the optimal-hash Bloom formula while kmindex uses one hash function; the real rate is about 0.29.

**Open points flagged in the report**: Vista `ulimit -n` and `vm.max_map_count` (probe first), the merge-stage RAM constant and partition imbalance (E2), `--cpr` ratio and queryability, build-time coefficients (E3), the 48 h wall-time limit per chunk, aarch64 binaries.

Two memory notes were saved for future sessions (project strategy, calculator reference). Nothing was committed.

Suggested commit message:

```
Add performance review report and experiment scripts for the SRA-2025 index

Verify the formulas of seb/performances.md against kmindex 0.6.1 and
kmtricks 1.5/1.6 with file:line citations, add the missing merge-stage
RAM and temp-disk terms, define the chunk/merge experiment plan for TACC
Vista, and add probe_limits.sh, estimate.py, run_experiment.sh and
plot_results.py under seb/scripts.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01MtW6mFPMKgYpCKyruXVZ8o
```
