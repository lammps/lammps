# Security Policy

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
behavior from calling functions in one of those libraries.

Thus is is quite easy to crash LAMMPS through malicious input and do all
kinds of file system manipulations.  And because of that LAMMPS should
**NEVER** be compiled or **run** as superuser, either from a "root" or
"administrator" account directly or indirectly via "sudo" or "su".

Therefore what could be seen as a security vulnerability is usually
either a user mistake or a bug in the code.  Bugs can be reported in the
LAMMPS project [issue tracker on
GitHub](https://github.com/lammps/lammps/issues).

To mitigate issues with using homoglyphs or bidirectional reordering in
unicode, which have been demonstrated as a vector to obfuscate and hide
malicious changes to the source code, all LAMMPS submissions are checked
for unicode characters and only all-ASCII source code is accepted.

# Version Updates

LAMMPS follows a continuous release development model.  We aim to keep
the development version (`develop` branch) always fully functional and
employ a variety of automatic testing procedures to detect failures
of existing functionality from adding or modifying features.  Most of
those tests are run on pull requests *before* merging to the `develop`
branch.  The `develop` branch is protected, so all changes *must* be
submitted as a pull request and thus cannot avoid the automated tests.

Additional tests are run *after* merging.  Before releases are made
*all* tests must have cleared.  Then a release tag is applied and the
`release` branch is fast-forwarded to that tag.  This is often referred
to as a patch release. Bug fixes and updates are
applied first to the `develop` branch.  Later, they appear in the `release`
branch when the next patch release occurs.
For stable releases, selected bug fixes, updates, and new functionality
are pushed to the `stable` branch and a new stable tag is applied.
