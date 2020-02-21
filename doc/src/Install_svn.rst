Download source via SVN
=======================

.. warning::

   As of Oct 2016, SVN support is now implemented via a
   git-to-subversion interface service on GitHub and no longer through a
   mirror of the internal SVN repository at Sandia.

You must have the `Subversion (SVN) client software <svn_>`_ installed on
your system to communicate with the Git server in this mode.

.. _svn: http://subversion.apache.org



You can follow LAMMPS development on 3 different SVN branches:

* **stable**   :  this branch is updated with every stable release
* **unstable** :  this branch is updated with every patch release
* **master**   :  this branch continuously follows ongoing development

The corresponding command lines to do an initial checkout are as
follows.  (Note that unlike Git, you must perform a separate checkout
into a unique directory for each of the 3 branches.)


.. parsed-literal::

   svn checkout https://github.com/lammps/lammps.git/branches/unstable mylammps
   svn checkout https://github.com/lammps/lammps.git/branches/stable mylammps
   svn checkout https://github.com/lammps/lammps.git/trunk mylammps

where "mylammps" is the name of the directory you wish to create on
your machine.

Once the command completes, your directory will contain the same files
as if you unpacked a current LAMMPS tarball, with the exception, that
the HTML documentation files are not included.  They can be fetched
from the LAMMPS website by typing "make fetch" in the doc directory.
Or they can be generated from the content provided in doc/src by
typing "make html" from the doc directory.

After initial checkout, as bug fixes and new features are added to
LAMMPS, as listed on :doc:`this page <Errors_bugs>`, you can stay
up-to-date by typing the following SVN commands from within the
"mylammps" directory:


.. parsed-literal::

   svn update

You can also check if there are any updates by typing:


.. parsed-literal::

   svn -qu status

Doing an "update" will not change any files you have added to the
LAMMPS directory structure.  It will also not change any existing
LAMMPS files you have edited, unless those files have changed in the
repository.  In that case, SVN will attempt to merge the new
repository file with your version of the file and tell you if there
are any conflicts.  See the SVN documentation for details.

Please refer to the `subversion client support help pages on GitHub <https://help.github.com/articles/support-for-subversion-clients>`_
if you want to use advanced features like accessing particular
previous release versions via tags.

Once you have updated your local files with an "svn update" (or "svn
co"), you still need to re-build LAMMPS if any source files have
changed.  To do this, you should cd to the src directory and type:


.. parsed-literal::

   make purge             # remove any deprecated src files
   make package-update    # sync package files with src files
   make foo               # re-build for your machine (mpi, serial, etc)

just as described on the :doc:`Install patch <Install_patch>` doc page,
after a patch has been installed.

.. warning::

   If you wish to edit/change a source file that is from a package, you
   should edit the version of the file inside the package sub-directory
   with src, then re-install the package.  The version in the src
   directory is merely a copy and will be wiped out if you type "make
   package-update".

The LAMMPS GitHub project is managed by Christoph Junghans (LANL,
junghans at lanl.gov), Axel Kohlmeyer (Temple U, akohlmey at
gmail.com) and Richard Berger (Temple U, richard.berger at
temple.edu).
