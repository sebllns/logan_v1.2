# Methodology used to produce `performances_review.md`

Date: 2026-09-09. Not committed at the time of writing.

## Tooling

| Item | Value |
| --- | --- |
| Assistant | Claude Fable 5.1 (model id `claude-fable-5-1`) in Claude Code (terminal), effort level `high` |
| Mode | plan mode first (read-only exploration, plan written to a plan file, approved by the user), then execution mode;  |

## Code versions

| Repository | Reference | Commit | Note |
| --- | --- | --- | --- |
| kmindex | `main`, `v0.6.1-5-gd326bbf` | `d326bbfc5b3b64b085bf9a3e048e1dfc468bfc6c` | working copy, all kmindex citations |
| kmtricks (submodule `thirdparty/kmtricks`) | `v1.5.0` | `561554f34e08a448e4f4a85ff3b0188ae938c9ae` | all kmtricks and GATB citations |
| kmtricks (target version) | `v1.6.0` | `9bccc774304b1303049399cbf1ffb562d27bd9ea` | target version; `git diff v1.5.0 v1.6.0` only adds BAM input, cited code identical |
| kmhelpers | `dev/simplify-options` | `f5dc0bd24a68a48fb72f367da6684106c7951ee3` | read with `git show`, same as `origin/dev/simplify-options` |

## Process

1. Exploration (plan mode). One `Explore` subagent read the kmindex sources (build, merge,
   query, storage layout, docs) and returned a cited report. Two other subagents
   (kmtricks, kmhelpers and sibling repos) were stopped by accident and replaced by direct
   reading with `grep`, `sed -n` and `git show` in the main session. Every line cited in
   the report was re-read after the version update to kmindex 0.6.1.
2. External sources fetched: kmindex documentation (`construction/`, `merge/`, `query/`),
   the TACC Vista user guide, the Logan Pareto merge calculator page, its README,
   `docs/index.html`, `code/enumerate_designs.py`, `data/training_campaign/training16.csv`
   and the `genouest_training_campaign` bundle (`README.md`, `run_one_config.sh`,
   `config.env`, `TUTORIAL.md`). GATB code was read in
   `thirdparty/kmtricks/thirdparty/gatb-core-stripped`.
3. Local reference values (open-file limits, `vm.max_map_count`) were measured on the
   laptop with the same commands as `claude/scripts/probe_limits.sh`.
4. Plan written and approved, then: scripts written and tested (`bash -n`,
   `py_compile`, sample runs of `estimate.py` cross-checked against the vendored
   `kmparams.py` of kmhelpers, `plot_results.py` on a fake CSV, `make_fofs.py` on dummy
   files), report written, then the hash-record size corrected after checking the default
   `MAX_C` of kmtricks.
5. Second round (this file): run commands added to each experiment, calculator section
   reworded, runner extended (RSS timeline, temp-disk sampling, `--cpr` flag).
6. Third round: question and answer loop on the review ("what to conclude" blocks, `-b`
   sizing, terms and notation in section 0), and `bugs.md` listing the defects found.

Limits of the method: formulas come from reading buffer declarations and loops, not from
measurements; the experiment plan is what turns them into constants. Line numbers refer
to kmindex 0.6.1 and kmtricks v1.5.0 (identical to v1.6.0 for the cited code).

## Prompts (user, typos fixed, content unchanged)

Prompt 1 (task definition):

> Here is the kmindex cloned repo (documentation link: https://tlemane.github.io/kmindex/).
>
> We started a report about performance vs parameters in `claude/performances.md`. I would like
> you to generate another markdown report in the same folder, `claude/`, where you verify that
> the formulas and assertions there are exact. You would source your assertions by
> referencing code (file:line) or URLs (documentation) or paper references. An important
> point is that we do not always need exact formulas; consistent estimations are good. We
> mainly need to extract the main parameters that impact each important metric. You'll see
> that we focus on 2 steps: build and query. One of the most important metrics is usually
> peak RAM usage. The others are how to estimate build time and query time, and which
> parameters impact them. Also note the presence of open files in the metrics: this metric
> is capped by a hard limit (see `ulimit -H -n`). We need to identify the system hard
> limits. At the start of each section of your document, don't hesitate to summarize the
> main questions we try to answer in this section.
>
> The main objective:
>
> We are going to index a huge dataset (complete 2025 SRA, ~40 million samples). We are
> going to do it on this cluster (https://docs.tacc.utexas.edu/hpc/vista/). We need to
> perform a bunch of experiments, especially to know how we should design the index.
> Indeed, we'll need to use the `kmindex merge` command, because building a sub-index in
> one go may sometimes not be possible due to a hard limit. Thus, we split it into chunks
> and then merge them. But the merge also has a cost (in time) that we need to estimate.
> So we would like to find the best strategy:
> - how many chunks?
> - positive impact of the number of chunks? (allows more threads; each chunk can also be
>   built in parallel on separate nodes; does this have an impact on query peak RAM too?
>   fewer chunks = bigger index = more RAM?)
> - negative impact of the number of chunks on build time (due to the merge step at the
>   end)?
>
> Still, it is important that you expose clearly each question we want to answer, then the
> tests/solutions.
>
> For this, we may need to use kmhelpers
> (https://github.com/sebllns/kmhelpers/tree/dev/simplify-options), especially the
> `auto-params` command in `kmhelpers/pykmhelpers/cli/test.py` (which implements the
> formulas in `performances.md`).
>
> Note that the experiments must be consistent, but we do not need to use millions of
> samples. Often, a small number of samples to be indexed should be enough for most
> experiments (still to be confirmed by your expertise, case by case).
>
> This is why it could be interesting that each experiment subsection exposes a table with
> parameters (clearly defined or not; you may not always have an answer).
>
> Bonus: another work has been done to evaluate the impact of index design on
> https://rebeccasalles.github.io/logan-pareto-merge-calculator/ (code here:
> https://github.com/RebeccaSalles/logan-pareto-merge-calculator) that you can also
> analyze and/or confront with what you found.
>
> PS: if useful, don't hesitate to generate bash or python scripts (or put the code in the
> markdown report) that we can use either to run actual experiments, or simulations, or
> plot generation.
>
> PS2: when you can't conclude or really hesitate on something, don't hesitate to flag it;
> you must be as consistent and concise as possible.

