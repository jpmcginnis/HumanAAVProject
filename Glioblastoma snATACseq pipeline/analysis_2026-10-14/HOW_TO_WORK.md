# How to work on analysis-heavy / data-pulling projects

A playbook distilled from running a cross-cohort snATAC-seq atlas rebuild
(8 cohorts, 1.5M cells, ~16 hours of overnight pipeline on EC2). Written to be
reusable for any similar project: long-running jobs, remote compute, large
intermediates, non-trivial failure modes.

Audience: Claude (or anyone else) picking up a similar project with my style.

---

## Part 1 — Status update format

Every substantive progress message to the user follows this shape:

```
**Working:** <one line: what's actually running right now>

**Issues hit + handled:** <one bullet per issue, with root cause + fix>
  (omit entirely if nothing new since last update)

**ETA current step:** <concrete number + what "done" looks like>

**ETA whole project:** <sum of remaining steps, biggest unknown called out>

**Running cost:** <only every ~6 hr or when crossing a round threshold>

**Next check:** <when I'll surface next + what will trigger it>
```

### Rules

- **Lead with action, not status.** "macs3 sample 7/14, X008 skipped" beats "pipeline running."
- **No fake numbers.** If ETA is a guess, say so ("±30 min, X008 is the unknown"). Never report a progress bar you don't believe.
- **"Handled" means verified fixed.** Not "I think it should work." Verify with a byte count or a real output file before claiming resolution.
- **Pivots count as issues.** Mid-stream strategy change: log it as `old-approach → new-approach — reason`.
- **"Routine" is a valid update.** When nothing changed, one word beats repackaging the same facts.
- **Flag problems immediately.** Don't wait to look confident. Users prefer "we hit X, fixing it as y" over polished "all good" that was wrong.
- **Running cost only occasionally.** It barely moves in 20-min ticks; mentioning it every tick creates noise.
- **"Next check" is a commitment.** If I say "when sample 7 lands," the next update begins at that trigger, not drift.

---

## Part 2 — Troubleshooting discipline

### Silence is not success

The single biggest failure mode: assuming no error = success. Verify with a
positive signal (file byte count, row count, known-good value), never with
"nothing crashed." My "7/14 OK" count was wrong in this project because my
regex matched 0-byte output files — if I'd trusted it, I'd have reported
hours of fake progress.

**Rule:** before reporting "step done," check an output you can measure.

### Diagnose before switching strategies

When something fails, spend 2 minutes probing before pivoting:
- `/proc/<pid>/io` — bytes read/written, lets you distinguish IO-bound from CPU-bound
- `/proc/<pid>/status` — process state (S=sleeping, D=disk-wait, R=running)
- `df -h` on all mounts — disk-full is a silent killer
- `free -g` — memory pressure
- `tail -N <log>` — most recent output

In this project, reading `State: D` on the greenwald process + `read_bytes` growing
at 5 MB/s pointed to HDF5 chunk-level compression as the bottleneck. That led
to a 10-line fix (sequential col-sums over data+indices) instead of a 100-line
rearchitecture.

**Rule:** every failure has a diagnosis window of ~2 min. Spend it before
switching libraries.

### Multiple paths to the goal

Each failure has at least three reasonable workarounds. Pick cheapest first;
only escalate if cheapest fails:

| Failure mode | Cheap | Medium | Expensive |
|---|---|---|---|
| Rust panic on obs read | catch BaseException, skip sample | split per-chromosome | rebuild h5ad clean |
| Disk fills | redirect TMPDIR | delete redundant tarball | mount new EBS volume |
| Slow HDF5 read | increase chunk cache | sequential bulk reads | re-materialize uncompressed |
| Snakatac2 won't read a sample | fall back to matrix-mode projection | manual h5py read | regenerate h5ad |

### Catch BaseException, not just Exception

snapatac2 raises `pyo3_runtime.PanicException` which does not inherit from
Exception. A vanilla `except Exception` lets the whole run die on one bad
sample. For any Rust-backed library (snapatac2, polars, pyarrow in some paths):

```python
try:
    a = snap.read(path)
except BaseException as e:
    log(f"ERROR {type(e).__name__}: {e}")
    continue  # skip this sample, keep going
```

### Pre-flight check on one unit before scaling

Before running 8 cohorts × hours each, test the per-unit operation on one unit
for seconds. In this project: `snap.pp.make_peak_matrix` on one sundaram sample
succeeded in 35 s. If I'd just scheduled it on all 14 without testing, I'd have
hit the API-signature mismatch at hour 2 of a 6-hour run.

### Don't accept absurd ETAs

8 hours for a projection of 1M cells × 438k peaks is wrong — the computation
is cheap (sparse matmul). If a reasonable operation is slow, you're doing it
inefficiently. Dig until you understand *why*, then fix the *why*. "I'll just
wait 8 hours" is almost always a bug you haven't found.

### Save intermediates aggressively

Every stage writes a checkpoint parquet/BED/log. Downstream stages skip steps
whose output already exists (`if out.exists(): continue`). If the chain dies at
stage 6/8, restart doesn't redo stages 1-5. This saved ~3 hours in this project
after the first greenwald attempt failed.

---

## Part 3 — Launch-in-background pattern

For any step that takes >2 min:

1. Write the script to a file (not inline), so it has a stable `__main__`.
2. Run as `nohup bash -c "python x.py && python y.py" > wrapper.log 2>&1 </dev/null & disown`.
3. Arm a `Monitor` that polls every 3-10 min for state changes.
4. Keep working on the next step's script while this runs.
5. On monitor event, update user with structured status.

**Monitor pattern:**
```
# Emit a one-line state summary; dedup consecutive identical lines; stall alert at 10 min
while true; do
  state=$(check_everything_in_one_ssh_call)
  if [ "$state" != "$prev" ]; then echo "$state"; prev=$state; last_change=now; fi
  if [ $((now - last_change)) -gt 600 ]; then echo "STALL — investigating"; last_change=now; fi
  if echo "$state" | grep -qE "done:1"; then echo "DONE"; break; fi
  sleep 300
done
```

Monitor emits only on *change* (noise suppression) and on *stall* (so you catch
hangs). Each notification goes to the user as a brief status.

---

## Part 4 — Precompute likely objections

Before the user raises a concern you can anticipate, write the answer. In this
project the user asked "what about the union-collapse case where 3 near-identical
peaks get merged into one wider interval?" — I had already written the
per-peak provenance table + the `pileup_region.py` helper for exactly this
concern.

**Rule:** when designing a result that collapses data (unions, averages,
merges), ship a traceback tool alongside it.

---

## Part 5 — Don't ask the user what you can decide

If the answer is in the brief, the code, or the data, decide. "Should I use
projection or recount?" was my mistake earlier in this project — the brief
clearly said "recount." Bouncing it back to the user wastes their attention.

Save escalation for:
- Genuine architectural tradeoffs with no technical "best" answer
- Scope boundaries the brief didn't anticipate
- Destructive actions (`rm -rf`, force-push, overwrite)
- Spending > $N when a budget was set

---

## Part 6 — Specific technical traps I'd want next-person to avoid

Not analysis-specific, project-wide. The ones that cost real time:

1. **macs3 writes temp files to $TMPDIR, defaulting to /tmp.** If /tmp is on a
   small root partition, it fills silently and macs3 emits zero-byte peak files.
   Set `export TMPDIR=/data/project/tmp` before any macs3 call.

2. **`bedtools merge -d 0` chains 500 bp peaks into chromosome-wide intervals**
   when peaks are densely packed. Use `-d -100` or do a midpoint clustering
   pass first.

3. **snapatac2 h5ads can hit a Rust panic on obs read** (`Expecting a dataframe
   but found: Mapping`) if obs has been written through certain code paths.
   Catch BaseException and skip; the panic is reproducible on specific h5ads.

4. **Trust `pgrep -c` returns, not `grep -c log_pattern`.** The latter is a lie
   if the log doesn't record the event you grepped for. Count output files
   instead.

5. **Cohort barcodes often have `-<sample_idx>` suffixes** that per-sample h5ads
   don't carry. If every cell shows "unassigned" after annotation join, suspect
   this. Auto-detect with a short routine that tries suffixes -0 through -20
   against a known cohort barcode set.

6. **anndata backed-mode row-slice reads are slow on gzip-compressed HDF5.**
   For large matrices, read `X.data` + `X.indices` sequentially with a 4 GB
   `rdcc_nbytes` chunk cache and accumulate column sums via `np.add.at`.
   Beats row-slicing by ~30×.

7. **Delete redundant tarballs before extraction, not after.** The GDC
   `full_cancer_scatacseq.tar.gz` was 134 GB sitting next to its own extracted
   contents — found while debugging a disk-full event. Checking `du -sh */` at
   project start would have surfaced this immediately.

8. **Pyo3 PanicException propagates past `except Exception`.** Always
   `except BaseException` around Rust-library calls in production code.

9. **"ok" in a status log is not success.** Verify with:
   `find reports -name "*.parquet" -size +0 | wc -l`

10. **Store things with canonical names from day one.** Legacy slugs (`gse276177_khan_astro`
    turned out to be misattributed to Khan when it was Sojka et al.) propagate
    everywhere — into filenames, column names, documentation, Zenodo metadata.
    Changing them later means a cross-tree rename that's easy to miss. Use
    `<first-author>_<year>` from the start.

---

## Part 7 — Project packaging discipline

Every deliverable project ends with:

- **CHANGELOG.md** with semantic versioning. Each version has: headline,
  what's new in data, what's new in method, what's removed, caveats.
- **README.md** describing the structure, the primary deliverables, how to
  interrogate specific peaks, known limitations.
- **Preserved regeneration scripts.** If I had to regenerate the atlas, I
  should be able to run `scripts/build_v5_reports.py` and have it produce
  the current `reports/` content without extra hand-steps.
- **Tests that use the actual data.** `tests/test_schema.py` reads the live
  parquet files and asserts column presence. If someone renames a column,
  the test fails loudly.
- **Trash dir for the previous version** (`_v3_trash_delete_me/`), not immediate
  `rm`. Lets the user audit before confirming. Delete only after sign-off.
- **Attribution.** `cohort_attribution.csv` with proper authorship. Correcting
  misattributions is worth doing explicitly.

---

## Part 8 — When to call it

A step is done when:
1. The expected output file exists
2. Its size / row count is in the expected range
3. A sanity check on one row matches biology (e.g., GABA_neuron top peak = GAD2
   promoter)
4. The schema test passes
5. The user has what they need in a human-readable form (HTML + CSV + BED)

If any of these are failing silently, you're not done.

---

## Final meta-point

Users remember the texture of how you worked more than any single answer.
Being persistent, diagnosing before pivoting, verifying claims, catching
problems early, preserving regeneration, and keeping status updates honest
and tight — these compound. One session of this discipline buys trust for
the next session. Don't ship vibes; ship checked work.
