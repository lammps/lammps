# LAMMPS Release Steps

The following notes chronicle the current steps for preparing and publishing LAMMPS releases.
For definition of what LAMMPS versions and the different kinds of releases mean, please
refer to [the corresponding section in the LAMMPS manual](https://docs.lammps.org/Manual_version.html).

## LAMMPS Feature Release

A LAMMPS feature release is currently prepared after about 500 to 750 commits to the
'develop' branch or after a period of four weeks up to two months.

### Preparing a 'next\_release' branch

Create a 'next\_release' branch off 'develop' and make the following changes:
- set the LAMMPS\_VERSION define to the planned release date in src/version.h in the format "D Mmm YYYY" or "DD Mmm YYYY"
- remove the LAMMPS\_UPDATE define in src/version.h
- update the release date in doc/lammps.1
- update all TBD arguments for ..versionadded::, ..versionchanged:: ..deprecated:: to the
  planned release date in the format "DMmmYYYY" or "DDMmmYYYY"

Submit this pull request, rebase if needed. This is the last pull request merged for the release
and should not contain any other changes. (Exceptions: this document, last minute trivial(!) changes).
This PR shall not be merged before **all** pending tests have completed and cleared. If needed, a
bugfix pull request should be created and merged to clear all tests.

## LAMMPS Stable Release

A LAMMPS stable release is prepared about once per year in the months July, August, or September.
One (or two, if needed) feature releases before the stable release shall contain only bug fixes
or minor feature updates in optional packages.  Also substantial changes to the core of the code
shall be applied rather toward the beginning of a development cycle between two stable releases
than toward the end.  The intention is to stablilize significant change to the core and have
outside users and developers try them out during the development cycle; the sooner the changes
are included, the better chances for spotting peripheral bugs and issues.

### Prerequesites

Before making a stable release all remaining backported bugfixes shall be released as a (final)
stable update release (see below).

A LAMMPS stable release process starts like a feature release (see above), only that this feature
release is called a "Stable Release Candidate" and no assets are uploaded to GitHub.

### Synchronize 'maintenance' branch with 'release'

The state of the 'release' branch is then transferred to the 'maintenance' branch (which will
have diverged significantly from 'release' due to the selectively backported bug fixes).

### Fast-forward merge of 'maintenance' into 'stable' and apply tag

At this point it should be possible to do a fast-forward merge of 'maintenance' to 'stable'
and then apply the stable\_DMmmYYYY tag.

### Push branches and tags



## LAMMPS Stable Update Release
