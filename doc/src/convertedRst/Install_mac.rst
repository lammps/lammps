Download an executable for Mac
==============================

LAMMPS can be downloaded, built, and configured for OS X on a Mac with
`Homebrew <homebrew_>`_.  Only four of the LAMMPS packages are unavailable
at this time because of additional needs not yet met: KIM, GPU,
USER-INTEL, USER-ATC.

After installing Homebrew, you can install LAMMPS on your system with
the following commands:


.. parsed-literal::

   % brew tap homebrew/science
   % brew install lammps              # serial version
   % brew install lammps --with-mpi   # mpi support

This will install the executable "lammps", a python module named
"lammps", and additional resources with all the standard packages.  To
get the location of the additional resources type this:


.. parsed-literal::

   % brew info lammps

This command also tells you additional installation options available.
The user-packages are available as options, just install them like
this example for the USER-OMP package:


.. parsed-literal::

   % brew install lammps --enable-user-omp

It is usually best to install LAMMPS with the most up to date source
files, which can be done with the "--HEAD" option:


.. parsed-literal::

   % brew install lammps --HEAD

To re-install the LAMMPS HEAD, run this command occasionally (make sure
to use the desired options).


.. parsed-literal::

   % brew install --force lammps --HEAD ${options}

Once LAMMPS is installed, you can test the installation with the
Lennard-Jones benchmark file:


.. parsed-literal::

   % brew test lammps -v

If you have problems with the installation you can post issues to
`this link <homebrew_>`_.

.. _homebrew: https://github.com/Homebrew/homebrew-science/issues

Thanks to Derek Thomas (derekt at cello.t.u-tokyo.ac.jp) for setting
up the Homebrew capability.



.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
