# Security Policy

[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/15199/badge)](https://www.bestpractices.dev/projects/15199)

LAMMPS is designed as a user-level application to conduct computer
simulations for research using classical mechanics.  As such LAMMPS
depends to some degrees on users providing correctly formatted input and
LAMMPS needs to read and write files based on uncontrolled user input.
As a parallel application for use in high-performance computing
environments, performance critical steps are also done without checking
data.

LAMMPS also is interfaced to a number of external libraries, including
libraries with experimental research software, that are not validated
and tested by the LAMMPS developers, so it is easy to import bad
behavior from calling functions in one of those libraries.  The section
[Security considerations](https://docs.lammps.org/latest/Build_security.html)
of the LAMMPS manual lists the external projects that may be used when
compiling LAMMPS, explains how files that are downloaded during the
build are checked, and what applies to LAMMPS provided by other projects
and to add-on packages for LAMMPS that are maintained elsewhere.

Thus it is quite easy to crash LAMMPS through malicious input and do all
kinds of file system manipulations.  A LAMMPS input can also run other
programs through shell commands, load plugins, run Python code, and
access the network.  And because of that LAMMPS should
**NEVER** be compiled or **run** as superuser, either from a "root" or
"administrator" account directly or indirectly via "sudo" or "su".
The LAMMPS executable prints a warning when it is started that way, and
so does CMake when configuring LAMMPS.

Therefore what could be seen as a security vulnerability is usually
either a user mistake or a bug in the code.  Bugs can be reported in the
LAMMPS project [issue tracker on
GitHub](https://github.com/lammps/lammps/issues).

To mitigate issues with using homoglyphs or bidirectional reordering in
unicode, which have been demonstrated as a vector to obfuscate and hide
malicious changes to the source code, all LAMMPS submissions are checked
for unicode characters and only all-ASCII source code is accepted.

# Reporting a Vulnerability

If you have found a problem that could be used to harm other LAMMPS
users or the LAMMPS project, please do **not** report it in a public
issue.  Examples for such problems are ways to get malicious changes
into the LAMMPS source code, into the libraries that are downloaded
when compiling LAMMPS, or into the pre-compiled LAMMPS packages, and
also passwords or access keys that have become public by accident.
Please report those problems privately in one of these two ways:

- with the "Report a vulnerability" button on the [Security
  page](https://github.com/lammps/lammps/security) of the LAMMPS
  repository on GitHub
- with an email to developers@lammps.org which is forwarded to the
  LAMMPS core developers

Please include which version of LAMMPS is affected and how the problem
can be reproduced.  You can expect a first response within one week.
The LAMMPS developers will then work with you to confirm and correct the
problem and will agree with you on when and how it is made public.

# Version Updates

LAMMPS follows a continuous release development model.  We aim to keep
the development version (`develop` branch) always fully functional and
employ a variety of automatic testing procedures to detect failures of
existing functionality from adding or modifying features.  Most of those
tests are run on pull requests and must be passed *before* merging to
the `develop` branch.  The `develop` branch is protected, so all changes
*must* be submitted as a pull request and thus cannot avoid the
automated tests.

Additional tests are run *after* merging.  Before releases are made
*all* tests must have cleared.  Then a release tag is applied and the
`release` branch is fast-forwarded to that tag.  This is referred to
as a "feature release".  Bug fixes and updates are applied first to the
`develop` branch.  Later, they appear in the `release` branch when the
next patch release occurs.  For stable releases, backported bug fixes
and infrastructure updates are first applied to the `maintenance` branch
and then merged to `stable` and published as "updates".  For a new
stable release the `stable` branch is updated to the corresponding state
of the `release` branch and a new stable tag is applied in addition to
the release tag.

# Supported Versions

Corrections are always applied to the development version first and
thus are part of the next feature release.  For the most recent stable
release, selected corrections are back-ported and published as stable
update releases.  Older versions are in general not updated, so please
upgrade to a current version.

# Integrity of Downloaded Archives

For *all* files that can be downloaded from the "lammps.org" web server
we provide SHA-256 checksum data in files named SHA256SUMS or
SHA256SUM.  These checksums can be used to validate the integrity of
the downloaded archives.  Please note that we also use symbolic links
to point to the latest or stable releases and the checksums for those
files *will* change (and so their checksums) because the symbolic links
will be updated for new releases.

Starting with the first release after the stable release of 30 Sep 2026,
the releases published on GitHub also contain a file `SHA256SUMS` with
the SHA-256 checksums of all files of that release, and a file
`SHA256SUMS.asc` with a digital signature for this list.  After
downloading these two files into the same folder as the other downloaded
files, the downloads can be checked on Linux with:

```
gpg --verify SHA256SUMS.asc SHA256SUMS
sha256sum --ignore-missing -c SHA256SUMS
```

The first command confirms that the list of checksums was created by
the LAMMPS developer who published the release (see below for the key),
and the second command confirms that the downloaded files are identical
to the published files.  On macOS use `shasum -a 256` instead of
`sha256sum`.

# Signed Release Tags

The tags for LAMMPS releases in the git repository are digitally signed.
In a clone of the LAMMPS repository a tag can be checked with:

```
git verify-tag stable_30Sep2026
```

The signatures for tags and checksum lists are currently made by Axel
Kohlmeyer with the key that has the fingerprint

```
EEA1 0376 4C6C 633E DC8A  C428 D9B4 4E93 BF0C 375A
```

The public part of this key is needed for the checks and can be imported
with:

```
curl -sL https://github.com/akohlmey.gpg | gpg --import
```

# Immutable GitHub Releases

Starting with LAMMPS version 10 Sep 2025 the LAMMPS releases published
on GitHub are configured as `immutable`.  This means that after the
release is published the release tag cannot be changed or any of the
uploaded assets, i.e. the source tarball, the static Linux executable
tarball and the pre-compiled packages of LAMMPS with LAMMPS-GUI included.
GitHub will generate a release attestation JSON file which can be
used to verify the integrity of the files provided with the release.
With the [GitHub command line program](https://cli.github.com/) this
can be done for a release as a whole, or for a single downloaded file,
with:

```
gh release verify --repo lammps/lammps stable_30Sep2026
gh release verify-asset --repo lammps/lammps stable_30Sep2026 lammps-src-30Sep2026.tar.gz
```
