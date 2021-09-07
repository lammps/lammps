# Outline of the GitHub Development Workflow

This purpose of this document is to provide a point of reference for the
core LAMMPS developers and other LAMMPS contributors to understand the
choices the LAMMPS developers have agreed on. Git and GitHub provide the
tools, but do not set policies, so it is up to the developers to come to
an agreement as to how to define and interpret policies. This document
is likely to change as our experiences and needs change and we try to
adapt accordingly. Last change 2021-09-02.

## Table of Contents

  * [GitHub Merge Management](#github-merge-management)
  * [Pull Requests](#pull-requests)
    * [Pull Request Assignments](#pull-request-assignments)
    * [Pull Request Reviews](#pull-request-reviews)
    * [Pull Request Discussions](#pull-request-discussions)
    * [Checklist for Pull Requests](#checklist-for-pull-requests)
  * [GitHub Issues](#github-issues)
  * [Milestones and Release Planning](#milestones-and-release-planning)

## GitHub Merge Management

In the interest of consistency, ONLY ONE of the core LAMMPS developers
should doing the merging itself.  This is currently
[@akohlmey](https://github.com/akohlmey) (Axel Kohlmeyer).  If this
assignment needs to be changed, it shall be done right after a stable
release.  If the currently assigned developer cannot merge outstanding
pull requests in a timely manner, or in other extenuating circumstances,
other core LAMMPS developers with merge rights can merge pull requests,
when necessary.

## Pull Requests

ALL changes to the LAMMPS code and documentation, however trivial, MUST
be submitted as a pull request to GitHub. All changes to the "master"
branch must be made exclusively through merging pull requests. The
"unstable" and "stable" branches, respectively are only to be updated
upon patch or stable releases with fast-forward merges based on the
associated tags. Pull requests may also be submitted to (long-running)
feature branches created by LAMMPS developers inside the LAMMPS project,
if needed. Those are not subject to the merge and review restrictions
discussed in this document, though, but get managed as needed on a
case-by-case basis.

### Pull Request Assignments

Pull requests can be "chaperoned" by one of the LAMMPS core developers.
This is indicated by who the pull request is assigned to. LAMMPS core
developers can self-assign or they can decide to assign a pull request
to a different LAMMPS developer. Being assigned to a pull request means,
that this pull request may need some work and the assignee is tasked to
determine whether this might be needed or not, and may either implement
the required changes or ask the submitter of the pull request to implement
them.  Even though, all LAMMPS developers may have write access to pull
requests (if enabled by the submitter, which is the default), only the
submitter or the assignee of a pull request may do so.  During this
period the `work_in_progress` label may be applied to the pull
request.  The assignee gets to decide what happens to the pull request
next, e.g. whether it should be assigned to a different developer for
additional checks and changes, or is recommended to be merged.  Removing
the `work_in_progress` label and assigning the pull request to the
developer tasked with merging signals that a pull request is ready to be
merged. In addition, a `ready_for_merge` label may also be assigned
to signal urgency to merge this pull request quickly.

### Pull Request Reviews

People can be assigned to review a pull request in two ways:

  * They can be assigned manually to review a pull request
    by the submitter or a LAMMPS developer
  * They can be automatically assigned, because a developers matches
    a file pattern in the `.github/CODEOWNERS` file, which associates
    developers with the code they contributed and maintain.

Reviewers are requested to state their appraisal of the proposed changes
and either approve or request changes. People may unassign themselves
from review, if they feel not competent about the changes proposed. At
least two approvals from LAMMPS developers with write access are required
before merging in addition to the automated compilation tests.
Merging counts as implicit approval, so does submission of a pull request
(by a LAMMPS developer). So the person doing the merge may not also submit
an approving review. The feature, that reviews from code owners are "hard"
reviews (i.e. they must all be approved before merging is allowed), is
currently disabled and it is in the discretion of the merge maintainer to
assess when a sufficient degree of approval, especially from external
contributors, has been reached in these cases.  Reviews may be
(automatically) dismissed, when the reviewed code has been changed,
and then approval is required a second time.

### Pull Request Discussions

All discussions about a pull request should be kept as much as possible
on the pull request discussion page on GitHub, so that other developers
can later review the entire discussion after the fact and understand the
rationale behind choices made.  Exceptions to this policy are technical
discussions, that are centered on tools or policies themselves
(git, GitHub, c++) rather than on the content of the pull request.

## GitHub Issues

The GitHub issue tracker is the location where the LAMMPS developers
and other contributors or LAMMPS users can report issues or bugs with
the LAMMPS code or request new features to be added. Bug reports have
a `[Bug]` marker in the subject line; suggestions for changes or
adding new functionality are indicated by a `[Feature Request]`
marker in the subject. This is automatically done when using the
corresponding template for submitting an issue.  Issues may be assigned
to one or more developers, if they are working on this feature or
working to resolve an issue.  Issues that have nobody working
on them at the moment or in the near future, have the label
`volunteer needed` attached.

When an issue, say `#125` is resolved by a specific pull request,
the comment for the pull request shall contain the text `closes #125`
or `fixes #125`, so that the issue is automatically deleted when
the pull request is merged.  The template for pull requests includes
a header where connections between pull requests and issues can be listed
and thus were this comment should be placed.

## Milestones and Release Planning

LAMMPS uses a continuous release development model with incremental
changes, i.e. significant effort is made - including automated pre-merge
testing - that the code in the branch "master" does not get easily
broken.  These tests are run after every update to a pull request.  More
extensive and time consuming tests (including regression testing) are
performed after code is merged to the "master" branch. There are patch
releases of LAMMPS every 3-5 weeks at a point, when the LAMMPS
developers feel, that a sufficient amount of changes have happened, and
the post-merge testing has been successful. These patch releases are
marked with a `patch_<version date>` tag and the "unstable" branch
follows only these versions (and thus is always supposed to be of
production quality, unlike "master", which may be temporary broken, in
the case of larger change sets or unexpected incompatibilities or side
effects.

About 1-2 times each year, there are going to be "stable" releases of
LAMMPS.  These have seen additional, manual testing and review of
results from testing with instrumented code and static code analysis.
Also, the last 1-3 patch releases before a stable release are "release
candidate" versions which only contain bugfixes and documentation
updates.  For release planning and the information of code contributors,
issues and pull requests being actively worked on are assigned a
"milestone", which corresponds to the next stable release or the stable
release after that, with a tentative release date.
