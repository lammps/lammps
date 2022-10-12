Reporting bugs
==============

If you are confident that you have found a bug in LAMMPS, please follow
the steps outlined below:

 * Check the `New features and bug fixes
   <https://www.lammps.org/bug.html>`_ section of the `LAMMPS WWW site
   <https://www.lammps.org>`_ or the
   `GitHub Releases page <https://github.com/lammps/lammps/releases>`_ to
   see if the bug has already been addressed in a patch release.
 * Check that your issue can be reproduced with the latest development
   version of LAMMPS.
 * Check the manual carefully to verify that the unexpected behavior you
   are observing is indeed in conflict with the documentation
 * Check the `GitHub Issue page <https://github.com/lammps/lammps/issues>`_
   if your issue has already been reported and if it is still open.
 * Check the `GitHub Pull Requests page <https://github.com/lammps/lammps/pulls>`_
   to see if there is already a fix for your bug pending.
 * Check the `LAMMPS forum at MatSci <https://matsci.org/lammps/>`_
   to see if the issue has been discussed before.

If none of these steps yields any useful information, please file a new
bug report on the `GitHub Issue page <https://github.com/lammps/lammps/issues>`_.
The website will offer you to select a suitable template with explanations
and then you should replace those explanations with the information that
you can provide to reproduce your issue.

The most useful thing you can do to help us verify and fix a bug is to
isolate the problem.  Run it on the smallest number of atoms and fewest
number of processors with the simplest input script that reproduces the
bug.  Try to identify what command or combination of commands is causing
the problem and upload the complete input deck as a tar or zip archive.
Please avoid using binary restart files unless the issue requires it.
In the latter case you should also include an input deck to quickly
generate this restart from a data file or a simple additional input.
This input deck can be used with tools like a debugger or `valgrind
<https://valgrind.org>`_ to further :doc:`debug the crash <Errors_debug>`.

You may also post a message in the `development category of the LAMMPS
forum at MatSci <https://matsci.org/c/lammps/lammps-development/>`_
describing the problem with the same kind of information.  The forum can
provide a faster response, especially if the bug reported is actually
expected behavior or other LAMMPS users have come across it before.

