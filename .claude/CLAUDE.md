# LAMMPS Project Instructions (Claude Code entry point)

@../.github/copilot-instructions.md

The imported file above is the always-loaded core shared with GitHub Copilot.  Its
"Task-Specific Guides" index lists the detailed guides in `.github/instructions/` and
`.github/dev-docs/`; read the matching guide BEFORE starting that kind of work.  The
files in `.claude/rules/` load the path-scoped guides automatically when matching
source files are read.

Machine- and user-specific facts (local build directory, git remotes, worktrees) live
in the untracked, gitignored `CLAUDE.local.md` at the repository root, which Claude
Code loads automatically when present.
