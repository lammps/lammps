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

.. _ai-contributions-policy:

Current policy for AI-assisted contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the time of writing (fall 2026), LAMMPS does not have an official
policy for contributions that were created with the help of AI tools.
Formulating such a policy is a pending task for the LAMMPS `Technical
Steering Committee
<https://github.com/lammps/lammps/blob/develop/.github/GOVERNANCE.md>`_.
Until then, the following rules apply:

- Contributions created with the help of AI tools must meet the same
  :doc:`requirements <Modify_requirements>` as all other contributions,
  and the LAMMPS developers apply these requirements more strictly to
  them (see the remarks on the :doc:`Information for Developers
  <Developer>` page).
- Disclose the use of AI tools in the "Artificial Intelligence (AI)
  Tools Usage" section of the pull request template: name the tools and
  models that were used and describe which parts of the changes they
  generated.  Using AI tools only to check spelling and syntax does not
  need to be reported.
- This section of the pull request is the only place for AI
  attribution.  Do not add AI attribution lines (like
  ``Co-Authored-By:``) to commit messages or pull request descriptions,
  and do not list an AI tool or model as contributing author in the
  header of source files.

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
  you can decide whether they belong in the same pull request.
- Several small, related fixes can be combined into one pull request;
  unrelated changes are easier to review in separate pull requests.

Git and licensing
^^^^^^^^^^^^^^^^^

- Do not let an agent rewrite history that was already pushed (no
  ``git rebase``, ``git commit --amend``, or ``git push --force`` on
  published branches); bring branches up to date by merging instead.
- Disclose the use of AI tools as described in the :ref:`current policy
  <ai-contributions-policy>` above.
- Code added to LAMMPS must be compatible with its license (GNU GPL,
  version 2).  Agents may reproduce code from other projects they have
  seen during training; ask them to name the sources they used and do
  not accept code copied from software with an incompatible license.

Keeping the instructions useful
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The instruction files are only as good as the experience that went into
them, and the LAMMPS developers welcome pull requests that correct,
clarify, or extend them.

We encourage everyone who works on LAMMPS with AI support to regularly
review the lessons learned in their sessions and to integrate those
that apply beyond the current task (a pitfall, a convention, a
debugging technique, a mistake the agent keeps repeating) into the
instruction files.  The current files are the result of such reviews,
which the LAMMPS developers do roughly once a month.  Coding agents can
help with this: ask them to summarize what they learned and to propose
changes to the instruction files.

When changing the instruction files, keep the always-read core file
short, put rules for specific parts of the source tree into the
matching ``.github/instructions`` file (with a stub in
``.claude/rules`` for new files), and put longer background material
into ``.github/dev-docs``.  Information that an agent can easily obtain
by reading the source code does not need to be repeated there, and
information that only applies to your own setup belongs in a local
file as explained above.
