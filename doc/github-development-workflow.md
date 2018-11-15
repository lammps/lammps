# Outline of the GitHub Development Workflow

This purpose of this document is to provide a point of reference for the
core LAMMPS developers and other LAMMPS contibutors to understand the
choices the LAMMPS developers have agreed on. Git and GitHub provide the
tools, but do not set policies, so it is up to the developers to come to
an agreement as to how to define and interpret policies. This document
is likely to change as our experiences and needs change and we try to
adapt accordingly. Last change 2018-11-15.

## Table of Contents

  * [GitHub Merge Management](#github-merge-management)
  * [Pull Requests](#pull-requests)
    * [Pull Request Assignments](#pull-request-assignments)
    * [Pull Request Reviews](#pull-request-reviews)
    * [Pull Request Discussions](#pull-request-discussions)
    * [Checklist for Pull Requests](#checklist-for-pull-requests)
  * [GitHub Issues](#github-issues)
    * [Feature Request Issues](#feature-request-issues)
    * [Bug Report Issues](#bug-report-issues)
  * [Milestones and Release Planning](#milestones-and-release-planning)

## GitHub Merge Management

In the interest of consistency, ONLY ONE of the core LAMMPS developers
should doing the merging itself.  This is currently @akohlmey (Axel
Kohlmeyer).  If this assignment needs to be changed, it shall be done
right after a stable release.

## Pull Requests

ALL changes to the LAMMPS code and documentation, however trivial, MUST
be submitted as a pull request to GitHub. All changes to the "master"
branch must be made exclusively through merging pull requests. The
"unstable" and "stable" branches are only to be updated upon patch or
stable releases with fast-forward merges based on the associated
tags. Pull requests may also be submitted to (long-running) feature
branches created by LAMMPS developers inside the LAMMPS project, if
needed. Those are not subject to the merge and review restrictions
discussed in this document, though.

### Pull Request Assignments

Pull requests can be "chaperoned" by one of the LAMMPS core developers.
This is indicated by who the pull request is assigned to. LAMMPS core
developers can self-assign or they can decide to assign a pull request
to a different LAMMPS developer. Being assigned to a pull request means,
that this pull request may need some work and the assignee is tasked to
determine what this might be needed or not, and may either implement the
required changes or ask the submitter of the pull request to implement
them.  Even though, all LAMMPS developers may have write access to pull
requests (if enabled by the submitter, which is the default), only the
submitter or the assignee of a pull request may do so.  During this
period the "work_in_progress" label shall be applied to the pull
request.  The assignee gets to decide what happens to the pull request
next, e.g. whether it should be assigned to a different developer for
additional checks and changes, or is recommended to be merged.  Removing
the "work_in_progress" label and assigning the pull request to the
developer tasked with merging signals that a pull request is ready to be
merged.

### Pull Request Reviews

People can be assigned to review a pull request in two ways:
  * They can be assigned manually to review a pull request
    by the submitter or a LAMMPS developer
  * They can be automatically assigned, because 
developer

### Pull Request Discussions

### Checklist for Pull Requests

## GitHub Issues

### Feature Request Issues

### Bug Report Issues

## Milestones and Release Planning

