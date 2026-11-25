# Instructions for AI Coding Agents

The LAMMPS instructions for coding agents are shared between tools and live under
`.github/`.  Before starting a task:

1. Read `.github/copilot-instructions.md`: the always-applicable core (build, tests,
   conventions, development workflow) with an index of the task-specific guides.
2. Before working on files, read the matching guides in `.github/instructions/`; the
   `applyTo:` header of each file lists the paths it covers (e.g. `src/**` for the C++
   rules, `doc/**` for the documentation rules).
3. For larger tasks, read the deep-dive documents in `.github/dev-docs/` listed in the
   index at the end of the core file.

Machine- or user-specific notes (local build directory, git remotes) belong in an
untracked, gitignored file such as `CLAUDE.local.md`, never in these shared files.

These files are excluded from release tarballs; work in a git clone of
https://github.com/lammps/lammps.
