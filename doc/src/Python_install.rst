Installing LAMMPS in Python
===========================

For Python to invoke LAMMPS, there are 2 files it needs to know about:

* python/lammps.py
* liblammps.so or liblammps.dylib

The python source code in lammps.py is the Python wrapper on the
LAMMPS library interface. The liblammps.so or liblammps.dylib file
is the shared LAMMPS library that Python loads dynamically.

You can achieve that Python can find these files in one of two ways:

* set two environment variables pointing to the location in the source tree
* run "make install-python" or run the python/install.py script explicitly

When calling "make install-python" LAMMPS will try to install the
python module and the shared library into the python site-packages folders;
either the system-wide ones, or the local users ones (in case of insufficient
permissions for the global install). Python will then find the module
and shared library file automatically. The exact location of these folders
depends on your python version and your operating system.

If you set the paths to these files as environment variables, you only
have to do it once.  For the csh or tcsh shells, add something like
this to your ~/.cshrc file, one line for each of the two files:


.. parsed-literal::

   setenv PYTHONPATH ${PYTHONPATH}:/home/sjplimp/lammps/python
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home/sjplimp/lammps/src

On MacOSX you may also need to set DYLD\_LIBRARY\_PATH accordingly.
For Bourne/Korn shells accordingly into the corresponding files using
the "export" shell builtin.

If you use "make install-python" or the python/install.py script, you need
to invoke it every time you rebuild LAMMPS (as a shared library) or
make changes to the python/lammps.py file, so that the site-packages
files are updated with the new version.

If the default settings of "make install-python" are not what you want,
you can invoke install.py from the python directory manually as


.. parsed-literal::

   % python install.py -m \<python module\> -l <shared library> -v <version.h file> [-d \<pydir\>]

* The -m flag points to the lammps.py python module file to be installed,
* the -l flag points to the LAMMPS shared library file to be installed,
* the -v flag points to the version.h file in the LAMMPS source
* and the optional -d flag to a custom (legacy) installation folder

If you use a legacy installation folder, you will need to set your
PYTHONPATH and LD\_LIBRARY\_PATH (and/or DYLD\_LIBRARY\_PATH) environment
variables accordingly, as described above.

Note that if you want Python to be able to load different versions of
the LAMMPS shared library (see :doc:`this section <Python_shlib>`), you will
need to manually copy files like liblammps\_g++.so into the appropriate
system directory.  This is not needed if you set the LD\_LIBRARY\_PATH
environment variable as described above.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
