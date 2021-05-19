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
kinds of filesystem manipulations.  And because of that LAMMPS should
**NEVER** be compiled or **run** as superuser, either from a "root" or
"administrator" account directly or indirectly via "sudo" or "su".

Therefore what could be seen as a security vulnerability is usually
either a user mistake or a bug in the code.  Bugs can be reported in
the LAMMPS project
[issue tracker on GitHub](https://github.com/lammps/lammps/issues).


# Version Updates

LAMMPS follows continuous release development model.  We aim to keep all
release versions (stable or patch) fully functional and employ a variety
of automatic testing procedures to detect failures of existing
functionality from adding new features before releases are made.  Thus
bugfixes and updates are only integrated into the current development
branch and thus the next (patch) release and users are recommended to
update regularly.
