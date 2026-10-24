Working with AI coding agents
-----------------------------

.. contents::
   :local:

------------

This page collects practical suggestions for using AI coding agents when
working on the LAMMPS source code.  They complement the general remarks
on the :doc:`Information for Developers <Developer>` page and reflect
the experience of the LAMMPS developers.  Since this is a rapidly
evolving field, the suggestions will be revised as needed.

Instruction files for coding agents
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The LAMMPS git repository contains instruction files that coding agents
read automatically.  They are organized so that only a small core is
always read and more specific information is added when it is needed:

- ``.github/copilot-instructions.md`` is the always-read core: how to
  build and test LAMMPS, general conventions, the development workflow,
  and an index of the more specific guides.
- ``.github/instructions/*.instructions.md`` contain rules for specific
  parts of the source tree, e.g. the C++ sources, the documentation, the
  KOKKOS package, or the unit tests.  The ``applyTo:`` line at the top
  of each file lists the paths it applies to.
- ``.github/dev-docs/`` contains longer background documents (e.g. on
  porting styles to the KOKKOS or OPENMP package, or on testing and
  verification strategies) that agents are told to read before starting
  the corresponding kind of work.

GitHub Copilot reads these files directly.  Claude Code reads
``.claude/CLAUDE.md``, which imports the core file, and the small files
in ``.claude/rules/`` tell it which of the specific guides to read when
it works on matching files.  Most other coding agents read the
``AGENTS.md`` file in the top-level folder, which points to the same
files.

These files are not included in the LAMMPS release tarballs.  To work
with a coding agent, use a git clone of the `LAMMPS GitHub repository
<https://github.com/lammps/lammps>`_ instead of an unpacked tarball.
This is also the recommended way to prepare a contribution, since pull
requests are submitted against the *develop* branch.

Information that only applies to your own setup, e.g. the location of
your build folder, the names of your git remotes, or additional
checkouts created with ``git worktree``, should not go into these shared
files.  Keep it in a local
file that is ignored by git, e.g. ``CLAUDE.local.md`` in the top-level
folder for Claude Code.

Verify, do not trust
^^^^^^^^^^^^^^^^^^^^

Coding agents are very good at producing code that looks plausible and
at explaining convincingly why it is correct.  Neither is proof that it
works.  The following habits have caught many errors in agent-written
code:

- Ask for tests that demonstrate a bug fix: the new or modified test
  must fail with the unmodified code and pass with the fix.
- Run changes to parallel code with several MPI processes and, where
  applicable, several OpenMP threads.  Changes to communication or
  threading that look correct by inspection often fail in actual
  parallel runs.
- For changes that are not supposed to change results (refactoring,
  code cleanup, porting), compare the output with that of the unmodified
  code, ideally to the last digit.
- Check that the unit tests actually ran: CTest reports a test as passed
  even when all of its test cases were skipped, e.g. because a required
  package was not included in the build.
- Build and test the code yourself (or have the agent do it) before
  accepting a change suggested in an automated code review; a suggestion
  can be plausible and still be wrong.

Scope and organization of changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Many styles exist in several variants (e.g. for the OPENMP, KOKKOS,
  GPU, or OPT package) that were created by copying and adapting code.
  After a bug is found, ask the agent to check the base style, all its
  accelerated variants, and similar styles for the same problem.
- Agents tend to "improve" code next to what they were asked to change.
  Ask them to report such findings instead of fixing them right away, so
  you can decide whether they belong into the same pull request.
- Several small, related fixes can be combined into one pull request;
  unrelated changes are easier to review in separate pull requests.

Git and licensing
^^^^^^^^^^^^^^^^^

- Do not let an agent rewrite history that was already pushed (no
  ``git rebase``, ``git commit --amend``, or ``git push --force`` on
  published branches); bring branches up to date by merging instead.
- Disclose the use of AI tools only in the corresponding section of the
  pull request template.  Do not add AI attribution lines (like
  ``Co-Authored-By:``) to commit messages or pull request descriptions.
- Code added to LAMMPS must be compatible with its license (GNU GPL,
  version 2).  Agents may reproduce code from other projects they have
  seen during training; ask them to name the sources they used and do
  not accept code copied from software with an incompatible license.

Keeping the instructions useful
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a session reveals a lesson that applies beyond the current task
(a pitfall, a convention, a debugging technique), consider adding it to
the instruction files as part of a pull request.  Keep the always-read
core file short, put rules for specific parts of the source tree into
the matching ``.github/instructions`` file (with a stub in
``.claude/rules`` for new files), and put longer background material
into ``.github/dev-docs``.  Information that an agent can easily obtain
by reading the source code does not need to be repeated there.
