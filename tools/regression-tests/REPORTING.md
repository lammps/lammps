# Regression test results: what the output files contain

This describes the information that `run_tests.py` and `merge_results.py` produce beyond
"passed" and "failed", and how a dashboard can turn it into something a developer can act
on.  See the `README` in this folder for how to run the tests.

## Why the extra information is needed

The example inputs under `examples/` were not written to be tests.  They are there to show
how to use LAMMPS, and many of them are old contributions from a time when submissions were
reviewed less strictly.  Using them as an additional source of validation is worthwhile
anyway: LAMMPS has far more functionality than the `unittest/` tree can cover, and running
the examples catches side effects of a change that nothing else would notice.  But it means
the results need interpretation, and a plain pass/fail number is misleading:

- A classical MD trajectory is chaotic.  The tiniest difference in the computed forces --
  from another machine, compiler, level of optimization, or number of MPI processes --
  grows exponentially.  After roughly 1000 steps it is visible in the thermo output, and a
  few thousand steps later it has reached the printed precision of every quantity.  A test
  that "fails" that way says nothing about the correctness of the code.
- A real difference in what is computed shows up in the first steps, or in subtle cases
  after one or two hundred.
- Some inputs cannot agree with their reference log file at all, because the input itself
  is not reproducible across a different number of MPI processes.  Those need a fix in the
  repository, not an investigation of the code.
- Some inputs run a production number of steps and use most of the machine time of a test
  run, or hit the timeout and produce no result at all.  Those need to be shortened in the
  repository.

So the goal of the reporting is to separate *look at this* from *expected* and from
*fix the example*, and to make the third group actionable.

## Where the information is

| file | written by | purpose |
|---|---|---|
| `progress.yaml` | `run_tests.py` | one line per input, YAML mapping keyed by the input script name |
| `output.xml` | `run_tests.py` | JUnit XML, one `<testcase>` per input |
| `run.json` | `merge_results.py` | merged results of all shards, for archiving and for the dashboard |
| `costs.json` | `merge_results.py` | measured run times, fed back into the next `--analyze` run |
| `summary.md` | `merge_results.py` | Markdown summary, already grouped by the categories below |

`run.json` is the file a dashboard should read.  Its `tests` mapping is keyed by
`<folder>/<input script>` relative to the top-level examples folder, so results from
different machines and configurations line up.

## Per-test fields in `run.json`

| field | type | meaning |
|---|---|---|
| `status` | string | `passed`, `failed`, `error`, `skipped` |
| `time` | float | measured wall-clock time of the run in seconds, including a run that crashed or timed out |
| `message` | string | one-line status, e.g. `completed, 6 abs diff and 6 rel diff checks failed` |
| `diverged` | int | number of thermo quantities that leave their tolerance somewhere in the run |
| `diverged_at` | int | the timestep at which the earliest of them does so; absent if the thermo output has no `Step` column |
| `diverged_row` | int | how many thermo outputs into the trajectory that is; `0` means the very first one |
| `attention` | string | a problem with the input script itself that a developer has to fix; see below |
| `details` | string | the full per-quantity report; only present when `merge_results.py --max-details` was used, always present in the merged JUnit XML |

The same values appear as attributes on the JUnit `<testcase>` element (`diverged`,
`diverged-at`, `diverged-row`, `attention`) and in the `failed_checks` mapping of
`progress.yaml` (`diverged`, `diverged_at`, `diverged_row`), with `attention` as a
top-level key of the progress entry.

## Classifying a failure

Use `diverged_row` and `diverged_at` together.  `diverged_row == 0` means the run already
differs in its very first thermo output, before any trajectory divergence is possible.
Otherwise `diverged_at` says after how many steps the deviation appears:

| condition | reading | what to do |
|---|---|---|
| `attention` is set | the input cannot match its reference | fix the example; check this first, it explains an early deviation on its own |
| `diverged_row == 0` | differs before the trajectory can diverge | a difference in the setup or in the computation; investigate |
| `diverged_at <= 200` | too early for a rounding difference to have grown | a difference in what is computed until proven otherwise; investigate |
| `200 < diverged_at <= 1000` | early, but not impossible | look at it when nothing more urgent is open |
| `diverged_at > 1000` | consistent with chaotic divergence | expected; shorten the run or widen the tolerance |

Check `attention` first: an input whose initial velocities depend on the number of MPI
processes deviates at the second line of output for a reason that has nothing to do with
the code, and in a configuration that uses fewer MPI processes than the reference log
files there are enough of those to bury everything else.  The Markdown summary therefore
leaves them out of the "worth investigating" group and lists them separately.

One caveat the dashboard should carry: a deviation cannot be seen before the first thermo
output after it starts.  An input with thermo output every 5000 steps cannot distinguish a
bug from chaos with this rule; `run_tests.py` says so in the detailed report when it
applies, and the dashboard can spot it as `diverged_row <= 1` with a large `diverged_at`.

## The `attention` field: what a developer has to fix

This is the field to build a work list from.  It is set independently of the verdict, so a
test that passes can still carry one.  There are currently four kinds, and an entry can
carry more than one, separated by `; `.

**Initial velocities that depend on the number of MPI processes.**

    velocity create with the default "loop all" and atoms from create_atoms: cannot match
    the reference log file, which was written with 4 MPI processes instead of 2

The `velocity create` command only assigns reproducible velocities when `loop geom` is
used, which seeds the random numbers from the coordinates of each atom.  The default
`loop all` assigns them by atom ID, which is only reproducible when the atoms were read
from a data file -- with `create_atoms` the IDs depend on the domain decomposition -- and
`loop local` is never reproducible.  Such an input diverges from its reference log file at
the second line of output whenever the number of MPI processes differs from the one the
reference was made with, which is exactly what happens in the configurations that use
fewer MPI processes and more threads.  The fix is in the repository: add `loop geom` to
the `velocity` command and recreate the reference log files.  This is only reported when
the run actually used a different number of processes than the reference log.

**Production sized runs.**

    the run needs 149 s: a production sized run that should be shortened in the repository
    (its run commands cannot be shortened automatically)

Reported when a run takes more than a quarter of the configured timeout or hits it.  These
are the legacy inputs that were never trimmed for testing.  They dominate the machine time
of a full test run, and the ones that hit the timeout produce no result at all.  The fix is
to reduce the number of steps in the input and recreate its reference log files.

Inputs without a reference log file are run with their `run` commands reduced to
`smoke_steps` steps (see the config file), since without a reference there is nothing to
compare and the only thing left to check is that the input does not crash.  Two things can
go wrong with that, and each adds a parenthesis that says what kind of work is needed:

    ... (shortening its run commands makes it fail, so it has to be done by hand)

The shortened input stops with an error, which happens when a variable refers to a fix
whose output is only produced after more steps than the shortened run has, for example
`fix ave/correlate` with a long correlation length.  The input is run unchanged instead, so
it costs what it always did, and shortening it needs a look at what the input actually does.

    ... (it does not even finish with its run commands shortened to 100 steps, so the
    cost is in the setup)

The shortened input still hits the timeout, so the cost is not in the number of steps at
all but in the system size, the potential, or a minimization.  Shortening the runs will not
help; the example needs a smaller system.

**A style that exists in no package.**

    names a style that does not exist in any package, so this is a typo or an outdated
    name that has to be fixed in the input

When a style is not found, LAMMPS names the package it belongs to
(`utils::check_packages_for_style()`), so `Unrecognized fix style 'qeqr/reaxff'` without
that part means no package provides it under that name at all.  That is a typo or a name
that was changed, not a build that is missing something, and the earlier reporting hid it
by calling every unrecognized style a missing package.

**A reference log file that belongs to no input script.**

    no reference log file matches this input, while log.29Aug24.tracker.g++.1 in this
    folder match no input script at all: check the names

A reference log file is only found when its name follows
`log.{date}.{basename}.{compiler}.{nprocs}` with the basename of an input script in the
same folder.  A typo in either name, or renaming one of the two without the other, turns
the input into one that cannot be checked against anything -- silently, because a missing
reference log file is not an error.  When an input has no usable reference log file and
its folder contains log files that belong to no input at all, the two facts are almost
always the same problem, and the message names the candidates.

## Statuses that are not verdicts

Several statuses mean "this input was not really tested".  A dashboard should count them
separately rather than folding them into "skipped", since each implies different work:

| status text contains | meaning | what it implies |
|---|---|---|
| `needs a multi-partition run` | the input uses `neb`, `prd`, `tad`, `temper`, `fix pimd`, a `universe` variable, ... | needs its own test configuration with `-partition NxM` and a matching number of MPI processes; the thermo output goes to one log file per partition |
| `no reference log file, only checked that ... does not crash` | no reference log exists for this number of processes, so the input was run shortened as a crash test | a reference log file could be added with `--gen-ref` |
| `numerical checks skipped due to missing the reference log file` | same, but the input could not be shortened | see the `attention` field |
| `numerical checks skipped, unsupported log file format` | the log file could not be converted into thermo data | the input probably uses `thermo_style multi` or prints something the parser does not handle |
| `package not installed` | the style is not in this build | not a problem with the input |

## Suggested dashboard sections

1. **Input scripts that need to be fixed** -- every test with an `attention` field, grouped
   by kind.  This is a work list against the repository, not against the code, and it is
   stable across runs.
2. **Worth investigating** -- failed tests with `diverged_row == 0` or
   `diverged_at <= 200`, most recent first.
3. **Expected divergence** -- failed tests with `diverged_at > 1000`, collapsed by default.
4. **Not tested** -- the statuses above, with counts per kind.
5. **Most expensive runs** -- sorted by `time`; the sum of the top entries is the lower
   bound for the walltime of the whole test run.

Comparing consecutive runs is what makes the difference between a new problem and a known
one; `merge_results.py --compare` already reports new failures, fixed tests, new tests and
removed tests, and the same comparison applied to the `attention` list shows whether the
work list is shrinking.
