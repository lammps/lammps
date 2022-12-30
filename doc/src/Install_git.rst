Download the LAMMPS source with git
-----------------------------------

All LAMMPS development is coordinated through the "LAMMPS GitHub
site".  If you clone the LAMMPS repository onto your local machine, it
has several advantages:

* You can stay current with changes to LAMMPS with a single git
  command.
* You can create your own development branches to add code to LAMMPS.
* You can submit your new features back to GitHub for inclusion in
  LAMMPS.  For that you should first create your own :doc:`fork on
  GitHub <Howto_github>`.

You must have `git <git_>`_ installed on your system to use the
commands explained below to communicate with the git servers on
GitHub.  For people still using subversion (svn), GitHub also
provides `limited support for subversion clients <svn_>`_.

.. note::

   As of October 2016, the official home of public LAMMPS development is
   on GitHub.  The previously advertised LAMMPS git repositories on
   git.lammps.org and bitbucket.org are now offline or deprecated.

.. _git: https://git-scm.com
.. _svn: https://help.github.com/en/github/importing-your-projects-to-github/working-with-subversion-on-github

You can follow the LAMMPS development on 3 different git branches:

* **stable**   :  this branch is updated from the *release* branch with
  every stable release version and also has selected bug fixes and updates
  back-ported from the *develop* branch
* **release**  :  this branch is updated with every patch release;
  updates are always "fast forward" merges from *develop*
* **develop**  :  this branch follows the ongoing development and
  is updated with every merge commit of a pull request

To access the git repositories on your box, use the clone command to
create a local copy of the LAMMPS repository with a command like:

.. code-block:: bash

   git clone -b release https://github.com/lammps/lammps.git mylammps

where "mylammps" is the name of the directory you wish to create on
your machine and "release" is one of the 3 branches listed above.
(Note that you actually download all 3 branches; you can switch
between them at any time using "git checkout <branch name>".)

.. admonition:: Saving time and disk space when using ``git clone``

   The complete git history of the LAMMPS project is quite large because
   it contains the entire commit history of the project since fall 2006,
   which includes the time when LAMMPS was managed with subversion.
   This includes a few commits that have added and removed some large
   files (mostly by accident).  If you do not need access to the entire
   commit history (most people don't), you can speed up the "cloning"
   process and reduce local disk space requirements by using the
   *--depth* git command line flag.  That will create a "shallow clone"
   of the repository containing only a subset of the git history.  Using
   a depth of 1000 is usually sufficient to include the head commits of
   the *develop* and the *release* branches.  To include the head commit
   of the *stable* branch you may need a depth of up to 10000.  If you
   later need more of the git history, you can always convert the
   shallow clone into a "full clone".

Once the command completes, your directory will contain the same files
as if you unpacked a current LAMMPS tarball, with the exception, that
the HTML documentation files are not included.  They can be fetched
from the LAMMPS website by typing ``make fetch`` in the doc directory.
Or they can be generated from the content provided in ``doc/src`` by
typing ``make html`` from the ``doc`` directory.

After initial cloning, as bug fixes and new features are added to
LAMMPS you can stay up-to-date by typing the following git commands
from within the "mylammps" directory:

.. code-block:: bash

   git checkout release      # not needed if you always stay in this branch
   git checkout stable       # use one of these 3 checkout commands
   git checkout develop      # to choose the branch to follow
   git pull

Doing a "pull" will not change any files you have added to the LAMMPS
directory structure.  It will also not change any existing LAMMPS
files you have edited, unless those files have changed in the
repository.  In that case, git will attempt to merge the new
repository file with your version of the file and tell you if there
are any conflicts.  See the git documentation for details.

If you want to access a particular previous release version of LAMMPS,
you can instead "check out" any version with a published tag. See the
output of ``git tag -l`` for the list of tags.  The git command to do
this is as follows.

.. code-block:: bash

   git checkout tagID

Stable versions and what tagID to use for a particular stable version
are discussed on `this page <https://www.lammps.org/bug.html#version>`_.
Note that this command will print some warnings, because in order to get
back to the latest revision and to be able to update with ``git pull``
again, you will need to do ``git checkout release`` (or
check out any other desired branch) first.

Once you have updated your local files with a ``git pull`` (or ``git
checkout``), you still need to re-build LAMMPS if any source files have
changed.  How to do this depends on the build system you are using.

.. tabs::

   .. tab:: CMake build

      Change to your build folder and type:

      .. code-block:: bash

         cmake . --build

      CMake should auto-detect whether it needs to re-run the CMake
      configuration step and otherwise redo the build for all files
      that have been changed or files that depend on changed files.
      In case some build options have been changed or renamed, you
      may have to update those by running:

      .. code-block:: bash

         cmake .

      and then rebuild.

   .. tab:: Traditional make

      Switch to the src directory and type:

      .. code-block:: bash

         make purge             # remove any deprecated src files
         make package-update    # sync package files with src files
         make foo               # re-build for your machine (mpi, serial, etc)

      to enforce consistency of the source between the src folder
      and package directories.  This is OK to do even if you don't
      use any packages. The ``make purge`` command removes any deprecated
      src files if they were removed by the patch from a package
      sub-directory.

      .. warning::

         If you wish to edit/change a src file that is from a package,
         you should edit the version of the file inside the package
         sub-directory with src, then re-install the package.  The
         version in the source directory is merely a copy and will be
         wiped out if you type "make package-update".

.. admonition:: Git protocols
   :class: note

   The servers at github.com support the "https://" access protocol for
   anonymous, read-only access.  If you have a suitably configured
   GitHub account, you may also use SSH protocol with the URL
   "git@github.com:lammps/lammps.git".

The LAMMPS GitHub project is currently managed by Axel Kohlmeyer
(Temple U, akohlmey at gmail.com).
