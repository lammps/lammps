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
review by the LAMMPS developers.

_Feature_ V.S _Stable_ release channels
'''''''''''''''''''''''''''''''''''''''

LAMMPS has 2 release channels - _Feature_ and _Stable_. By default commits
reach first the `develop` branch, and every certain amount of time we
create a GitHub release, along with a Git tag for either of the channels.
In parallel to the Git tags, we also maintain a Git branch that is
identical of the latest Git tag of the release (either a _Feature_ or
_Stable_). The table below lists the names

Particularly for the _Stable_ channel, we use a branch named `maintenance`,
to backport commits from the `develop` branch, that we want the stable
channel to include. Hence the `maintenance` branch is a mediator between
the main `develop` branch and the `stable` branch that is parallel to the
latest tag of the _Stable_ channel.

+--------------+----------------+-----------------------------------+----------------------------------------------+------------------+
| Release name | Git tag prefix | Git Branch parallel to latest tag | Git branches workflow                        | Release schedule |
+--------------+----------------+-----------------------------------+----------------------------------------------+------------------+
| Feature      | ``patch_``     | ``release``                       | ``develop`` -> ``release``                   | Every 4-8 weeks  |
+--------------+----------------+-----------------------------------+----------------------------------------------+------------------+
| Stable       | ``stable_``    | ``stable``                        | ``develop`` -> ``maintenance`` -> ``stable`` | 1-2 times a year |
+--------------+----------------+-----------------------------------+----------------------------------------------+------------------+

A summary of the most important changes of the _Stable_ channel releases
are on `this website page <https://www.lammps.org/bug.html>`_.  More
detailed release notes are `available on GitHub
<https://github.com/lammps/lammps/releases/>`_.

The code in the _Stable_ release channel is subjected to more detailed and
thorough manual testing than the default automated testing. Also, several
variants of static code analysis are run to improve the overall code
quality, consistency, and compliance with programming standards, best
practices and style conventions.

Each version of LAMMPS contains all the documented *features* up to and
including its version date.  For recently added features, we add markers to
the documentation at which specific LAMMPS version a feature or keyword was
added or significantly changed.

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
