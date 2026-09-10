# Defects found while writing `performances_review.md`

Versions: kmindex 0.6.1 (`d326bbfc`), kmtricks v1.5.0 (`561554f3`, identical code in
v1.6.0), kmhelpers `dev/simplify-options` (`f5dc0bd2`). Paths as in the review: kmtricks
paths relative to `thirdparty/kmtricks/`, kmindex paths to the repo root. None of these
affects the index content or the query results; they affect tooling, documentation and
memory.

| # | Component | Where | What | Effect | Fix |
| --- | --- | --- | --- | --- | --- |
| 1 | kmindex | `lib/src/index/index.cpp:199` | `index::mergeable()` tests `vec.size() > 1` (the output vector) instead of `v.size() > 1` (the group of indexes sharing a SHA-1) | `kmindex index-infos` never lists mergeable sets, the output vector stays empty | replace `vec.size()` by `v.size()` |
| 2 | kmindex | `app/kmindex/query.cpp:104` (help), `:375` (loop) | help of `-b/--batch-size` says `0≈nb_seq/nb_thread`; the worker loop only solves a batch when `nq == batch_size` or at end of file, so `0` means "all the sequences a worker reads, in one batch" | with the default `-b 0` every response buffer of the file is live at once, peak RAM $\approx N_{seq} \cdot n \cdot \lceil S \cdot bw/8 \rceil$ (review 6.1, 6.4); on a 40M-sample index one 1 kb sequence is 4.9 GB | fix the help text, or implement the documented default |
| 3 | kmindex docs | `docs/changelogs/v0.6.1.md:3` | changelog names the new flag `--static-report`; the flag is `--static-repart` (`app/kmindex/build.cpp:84, 186`) | user passes a flag that does not exist | fix the changelog |
| 4 | kmindex docs | `docs/kmindex/docs/construction.md:8` | says indexes built by kmtricks are usable "without `--cpr`", while `kmindex build` forwards `--cpr` itself (`build.cpp:80`, changelog v0.3.0). In kmtricks v1.5.0 `pipeline --cpr` compresses intermediate files only: the Bloom matrices are always written with `compressed = false` (`include/kmtricks/task.hpp:768-774`). Only a manual `kmtricks merge --cpr` (`src/cli.cpp:626`) compresses matrices | the warning is stale for the pipeline and hides that `--cpr` is a safe temp-disk lever (review 7.2); it still holds for a manual `kmtricks merge --cpr` | reword: "not built with `kmtricks merge --cpr`" |
| 5 | kmindex | `lib/src/query/query_results.cpp:114-200` | `compute_abs()` and `compute_abs_pos()` never call `m_qr->free()`, unlike `compute_ratios()` (`:70`) and `compute_ratios_pos()` (`:111`) | for abundance indexes (`--bitw` > 1) the response buffers of a batch stay allocated until the batch is destroyed, so the batch peak adds the responses to the results; not hit by presence/absence indexes | add `m_qr->free()` at the end of both functions |
| 6 | kmtricks (GATB layer) | `include/kmtricks/gatb/fill_partitions.hpp:37-56` | the `cache_items` constructor argument of `KmFillPartitions` is never used | no per-partition k-mer cache exists, so the super-k-mer stage RAM is $\lfloor Tf \rfloor (P \cdot 32\,\text{KB} + \text{PartiInfo})$ (review 3.1) and GATB's cache term in its memory estimate is void | remove the argument or the estimate term; informational |
| 7 | kmhelpers | `core/build_params.py:92, 104`; `pipeline/index_ops.py:134`; `cli/build.py:152`; `cli/shared.py:129` | `safety_margin` has four defaults: `auto_params` signature 1.0, its docstring 0.9, `IndexOpsConfig` 0.9, `kmhelpers build` hard-codes 0.75, the CLI option 1.0 | the same dataset gets different chunk sizes and thread counts depending on the entry point | pick one default, expose it once |
| 8 | kmhelpers | `core/resources.py:34-37` | `get_max_open_files()` reads the soft `RLIMIT_NOFILE`; `safety_margin` scales that soft value | `maximize_nofile()` (soft raised to hard) runs only just before the kmindex command (`core/kmindex_wrapper.py:253`), after `auto_params` (`pipeline/index_ops.py:580`) computed the plan: on a node where the soft limit is far below the hard one, the chunk size comes from the soft limit while the run could use the hard one (review 2.2) | call `maximize_nofile()` before `auto_params`, or read the hard limit in `get_max_open_files()` |

Not defects, but easy to mistake for one:

- `--fast` on a compressed index is ignored with a warning (`app/kmindex/query.cpp:322-325`), not silently.
- `kmindex merge` has no partition-range option (`app/kmindex/merge.cpp:19-57`); a merge cannot be split across nodes by design (review 5.1).
- kmtricks `--max-memory` is not exposed by `kmindex build`; it is a constant used only by the automatic partition count (review section 1).

Errors in `performances.md` (formulas and drivers, not code) are in the verification
table, section 8 of the review.
