Download source via Git
=======================

All LAMMPS development is coordinated through the "LAMMPS GitHub
site".  If you clone the LAMMPS repository onto your local machine, it
has several advantages:

* You can stay current with changes to LAMMPS with a single git
  command.
* You can create your own development branches to add code to LAMMPS.
* You can submit your new features back to GitHub for inclusion in
  LAMMPS.

You must have `Git <git_>`_ installed on your system to communicate with
the public Git server for LAMMPS.

.. warning::

   As of Oct 2016, the official home of public LAMMPS
   development is on GitHub.  The previously advertised LAMMPS git
   repositories on git.lammps.org and bitbucket.org are now deprecated,
   may not be up-to-date, and may go away at any time.

.. _git: http://git-scm.com



You can follow LAMMPS development on 3 different Git branches:

* **stable**   :  this branch is updated with every stable release
* **unstable** :  this branch is updated with every patch release
* **master**   :  this branch continuously follows ongoing development

To access the Git repositories on your box, use the clone command to
create a local copy of the LAMMPS repository with a command like:


.. parsed-literal::

   git clone -b unstable https://github.com/lammps/lammps.git mylammps

where "mylammps" is the name of the directory you wish to create on
your machine and "unstable" is one of the 3 branches listed above.
(Note that you actually download all 3 branches; you can switch
between them at any time using "git checkout <branch name>".)

Once the command completes, your directory will contain the same files
as if you unpacked a current LAMMPS tarball, with the exception, that
the HTML documentation files are not included.  They can be fetched
from the LAMMPS website by typing "make fetch" in the doc directory.
Or they can be generated from the content provided in doc/src by
typing "make html" from the doc directory.

After initial cloning, as bug fixes and new features are added to
LAMMPS, as listed on :doc:`this page <Errors_bugs>`, you can stay
up-to-date by typing the following Git commands from within the
"mylammps" directory:


.. parsed-literal::

   git checkout unstable      # not needed if you always stay in this branch
   git checkout stable        # use one of the 3 checkout commands
   git checkout master
   git pull

Doing a "pull" will not change any files you have added to the LAMMPS
directory structure.  It will also not change any existing LAMMPS
files you have edited, unless those files have changed in the
repository.  In that case, Git will attempt to merge the new
repository file with your version of the file and tell you if there
are any conflicts.  See the Git documentation for details.

If you want to access a particular previous release version of LAMMPS,
you can instead "checkout" any version with a published tag. See the
output of "git tag -l" for the list of tags.  The Git command to do
this is as follows.


.. parsed-literal::

   git checkout tagID

Stable versions and what tagID to use for a particular stable version
are discussed on :doc:`this page <Errors_bugs>`.  Note that this command
will print some warnings, because in order to get back to the latest
revision and to be able to update with "git pull" again, you first
will need to first type "git checkout unstable" (or check out any
other desired branch).

Once you have updated your local files with a "git pull" (or "git
checkout"), you still need to re-build LAMMPS if any source files have
changed.  To do this, you should cd to the src directory and type:


.. parsed-literal::

   make purge             # remove any deprecated src files
   make package-update    # sync package files with src files
   make foo               # re-build for your machine (mpi, serial, etc)

just as described on the :doc:`Install patch <Install_patch>` doc page,
after a patch has been installed.

.. warning::

   If you wish to edit/change a src file that is from a
   package, you should edit the version of the file inside the package
   sub-directory with src, then re-install the package.  The version in
   the src dir is merely a copy and will be wiped out if you type "make
   package-update".

.. warning::

   The GitHub servers support both the "git://" and
   "https://" access protocols for anonymous read-only access.  If you
   have a correspondingly configured GitHub account, you may also use SSH
   with "git@github.com:/lammps/lammps.git".

The LAMMPS GitHub project is managed by Christoph Junghans (LANL,
junghans at lanl.gov), Axel Kohlmeyer (Temple U, akohlmey at
gmail.com) and Richard Berger (Temple U, richard.berger at
temple.edu).


