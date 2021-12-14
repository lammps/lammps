What does a LAMMPS version mean
-------------------------------

The LAMMPS "version" is the date when it was released, such as 1 May
2014.  LAMMPS is updated continuously and we aim to keep it working
correctly and reliably at all times.  You can follow its development
in a public `git repository on GitHub <https://github.com/lammps/lammps>`_.

Whenever we fix a bug or update or add a feature, it will be merged into
the *develop* branch of the git repository.  When a sufficient number of
changes have accumulated *and* the software passes a set of automated
tests, we release it in the next *patch* release, which are made every
few weeks.  The *release* branch of the git repository is updated with
every such release.  Info on patch releases are on `this website page
<https://www.lammps.org/bug.html>`_.

Once or twice a year, we apply only bug fixes and small, non-intrusive
changes to the *develop* branch and the code is subjected to more detailed
and thorough testing than the default automated testing.  The latest
patch release after such a period is then also labeled as a *stable* version
and the *stable* branch is updated with it.  Between stable releases
we occasionally release some updates to the stable release containing
only bug fixes and updates back-ported from *develop* but no new features
and update the *stable* branch accordingly.

Each version of LAMMPS contains all the documented features up to and
including its version date.

The version date is printed to the screen and logfile every time you
run LAMMPS. It is also in the file src/version.h and in the LAMMPS
directory name created when you unpack a tarball.  And it is on the
first page of the :doc:`manual <Manual>`.

* If you browse the HTML pages on the LAMMPS WWW site, they will by
  default describe the most current patch release version of LAMMPS.
  In the navigation bar on the bottom left, there is the option to
  view instead the documentation for the most recent *stable* version
  or the latest version from the current development branch.
* If you browse the HTML pages included in your tarball, they
  describe the version you have, which may be older.
