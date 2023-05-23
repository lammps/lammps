What does a LAMMPS version mean
-------------------------------

The LAMMPS "version" is the date when it was released, such as 1 May
2014.  LAMMPS is updated continuously, and we aim to keep it working
correctly and reliably at all times.  You can follow its development
in a public `git repository on GitHub <https://github.com/lammps/lammps>`_.

Modifications of the LAMMPS source code (like bug fixes, code refactors,
updates to existing features, or addition of new features) are organized
into pull requests.  Pull requests will be merged into the *develop*
branch of the git repository after they pass automated testing and code
review by the LAMMPS developers.  When a sufficient number of changes
have accumulated *and* the *develop* branch version passes an extended
set of automated tests, we release it as a *feature release*, which are
currently made every 4 to 8 weeks.  The *release* branch of the git
repository is updated with every such release.  A summary of the most
important changes of the patch releases are on `this website page
<https://www.lammps.org/bug.html>`_.  More detailed release notes are
`available on GitHub <https://github.com/lammps/lammps/releases/>`_.

Once or twice a year, we have a "stabilization period" where we apply
only bug fixes and small, non-intrusive changes to the *develop*
branch.  At the same time, the code is subjected to more detailed and
thorough manual testing than the default automated testing.  Also,
several variants of static code analysis are run to improve the overall
code quality, consistency, and compliance with programming standards,
best practices and style conventions.

The release after such a stabilization period is called a *stable*
version and both, the *release* and the *stable* branches are updated
with it.  Between stable releases, we collect back-ported bug fixes and
updates from the *develop* branch in the *maintenance* branch.  From the
*maintenance* branch we make occasional update releases and update the
*stable* branch accordingly.

Each version of LAMMPS contains all the documented *features* up to and
including its version date.  For recently added features, we add markers
to the documentation at which specific LAMMPS version a feature or
keyword was added or significantly changed.

The version date is printed to the screen and log file every time you run
LAMMPS.  It is also in the file src/version.h and in the LAMMPS
directory name created when you unpack a tarball.  And it is on the
first page of the :doc:`manual <Manual>`.

* If you browse the HTML pages of the online version of the LAMMPS
  manual, they will by default describe the most current feature release
  version of LAMMPS.  In the navigation bar on the bottom left, there is
  the option to view instead the documentation for the most recent
  *stable* version or the documentation corresponding to the state of
  the development branch.
* If you browse the HTML pages included in your downloaded tarball, they
  describe the version you have, which may be older than the online
  version.
