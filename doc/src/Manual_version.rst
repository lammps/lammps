What does a LAMMPS version mean
-------------------------------

The LAMMPS "version" is the date when it was released, such as 1 May
2014.  LAMMPS is updated continuously and we aim to keep it working
correctly and reliably at all times.  You can follow its development
in a public `git repository on GitHub <https://github.com/lammps/lammps>`_.

Whenever we fix a bug or update or add a feature, it will be merged into
the `master` branch of the git repository.  When a sufficient number of
changes have accumulated *and* the software passes a set of automated
tests, we release it in the next *patch* release, which are made every
few weeks.  Info on patch releases are on `this website page
<https://www.lammps.org/bug.html>`_.

Once or twice a year, only bug fixes and small, non-intrusive changes are
included for a period of time, and the code is subjected to more detailed
and thorough testing than the default automated testing.  The latest
patch release after such a period is then labeled as a *stable* version.

Each version of LAMMPS contains all the features and bug-fixes up to
and including its version date.

The version date is printed to the screen and logfile every time you
run LAMMPS. It is also in the file src/version.h and in the LAMMPS
directory name created when you unpack a tarball.  And it is on the
first page of the :doc:`manual <Manual>`.

* If you browse the HTML pages on the LAMMPS WWW site, they always
  describe the most current patch release of LAMMPS.
* If you browse the HTML pages included in your tarball, they
  describe the version you have, which may be older.
