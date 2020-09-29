Installation
************

The LAMMPS Python module enables calling the :ref:`LAMMPS C library API <lammps_c_api>`
from Python by dynamically loading functions in the LAMMPS shared library through the
Python ``ctypes`` module. Because of the dynamic loading, it is required that
LAMMPS is compiled in *shared mode*.

Two files are necessary for Python to be able to invoke LAMMPS code:

* LAMMPS Python Module (``python/lammps.py``)
* LAMMPS Shared Library (e.g., ``liblammps.so``)


.. _python_virtualenv: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment
.. _python_venv: https://docs.python.org/3.8/library/venv.html
.. _python_pep405: https://www.python.org/dev/peps/pep-0405

.. _python_install_guides:

Installing the LAMMPS Python Module and Shared Library
======================================================

Making LAMMPS usable within Python and vice versa requires putting the LAMMPS
Python module into a location that the Python interpreter can find and
installing the LAMMPS shared library into a folder that the dynamic loader
searches. For some potentials LAMMPS also needs to know where it can find the
necessary potential files.

Both CMake and traditional make build options offer ways to automate these tasks.

.. tabs::

   .. tab:: CMake (local user)

      LAMMPS can be configured and compiled as shared library with CMake by enabling the ``BUILD_SHARED_LIBS`` option.
      The file name of the shared library depends on the platform (Unix/Linux, MacOS, Windows) and the build configuration
      being used. See :ref:`Build the LAMMPS executable and library <library>` for more details and how the name is
      determined.

      After compilation, the generated binaries, shared library, Python module,
      and other files can be installed to a custom location defined by the
      ``CMAKE_INSTALL_PREFIX`` setting.  By default, this is set to the current
      user's ``$HOME/.local`` directory. This leads to an installation to the following locations:

      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                  | Notes                                                       |
      +========================+===========================================================+=============================================================+
      | LAMMPS Python Module   | * ``$HOME/.local/lib/pythonX.Y/site-packages/`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$HOME/.local/lib64/pythonX.Y/site-packages/`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$HOME/.local/lib/`` (32bit)                           |                                                             |
      |                        | * ``$HOME/.local/lib64/`` (64bit)                         |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | ``$HOME/.local/share/lammps/potentials/``                 |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+

      The following is a minimal working example:

      1. Install LAMMPS Shared Library and Python module using CMake

         .. code-block:: bash

            # create and change into build directory
            mkdir build
            cd build

            # configure LAMMPS compilation
            # compile with shared library, PYTHON package, and C++ exceptions
            # TODO: add more options to customize your LAMMPS installation
            cmake -C ../cmake/presets/minimal.cmake    \
                  -D BUILD_SHARED_LIBS=on              \
                  -D PKG_PYTHON=on                     \
                  -D LAMMPS_EXCEPTIONS=on              \
                  ../cmake

             # compile LAMMPS (in parallel for faster builds)
             cmake --build . --parallel

             # install LAMMPS into myvenv
             cmake --install .

      2. Configure Environment Variables

         To use this installation you have to ensure that the folder containing
         the LAMMPS shared library is part of the ``LD_LIBRARY_PATH`` environment variable (or
         ``DYLD_LIBRARY_PATH`` on MacOS). This allows the dynamic library loader of your system
         to find the LAMMPS shared library when needed.

         .. code-block:: bash

            # Unix/Linux
            export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH

            # MacOS
            export DYLD_LIBRARY_PATH=$HOME/.local/lib:$DYLD_LIBRARY_PATH


         LAMMPS will also need to know the location of the folder
         containing its potential files.  This can be set with the ``LAMMPS_POTENTIALS``
         environment variable:

         .. code-block::

            export LAMMPS_POTENTIALS=$HOME/.local/share/lammps/potentials

         To set these environment variables for each new shell, add the above
         ``export`` commands at the end of the ``$HOME/.bashrc`` file.

      3. Verify if LAMMPS can be successfully started from Python

         .. code-block:: bash

            $ python
            Python 3.8.5 (default, Sep  5 2020, 10:50:12)
            [GCC 10.2.0] on linux
            Type "help", "copyright", "credits" or "license" for more information.
            >>> import lammps
            >>> lmp = lammps.lammps()
            LAMMPS (18 Sep 2020)
              using 1 OpenMP thread(s) per MPI task
            >>>

      .. note::

         If you recompile LAMMPS, you will have to also rerun the install step to
         ensure the latest Python module and shared library are installed.

   .. tab:: CMake (system-wide)

      A system-wide installation allows all users to run Python with LAMMPS
      included.  Note that during the installation step you will need to either be
      root or use ``sudo`` to elevate your write privileges. The compilation steps are identical
      to the local user installation, with the only difference that
      ``CMAKE_INSTALL_PREFIX`` is set to system folder such as ``/usr``. This leads to
      the following installation locations:

      +------------------------+---------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                          | Notes                                                       |
      +========================+===================================================+=============================================================+
      | LAMMPS Python Module   | * ``/usr/lib/pythonX.Y/site-packages/`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``/usr/lib64/pythonX.Y/site-packages/`` (64bit) |                                                             |
      +------------------------+---------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``/usr/lib/`` (32bit)                           |                                                             |
      |                        | * ``/usr/lib64/`` (64bit)                         |                                                             |
      +------------------------+---------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | ``/usr/share/lammps/potentials/``                 |                                                             |
      +------------------------+---------------------------------------------------+-------------------------------------------------------------+

      The following is a minimal working example:

      1. Install LAMMPS shared library and Python module into system folder

         .. code-block:: bash

            # configure LAMMPS compilation
            # compile with shared library, PYTHON package, and C++ exceptions
            # TODO: add more options to customize your LAMMPS installation
            cmake -C ../cmake/presets/minimal.cmake    \
                  -D BUILD_SHARED_LIBS=on              \
                  -D PKG_PYTHON=on                     \
                  -D LAMMPS_EXCEPTIONS=on              \
                  -D CMAKE_INSTALL_PREFIX=/usr         \
                  ../cmake

             # compile LAMMPS (in parallel for faster builds)
             cmake --build . --parallel

             # install LAMMPS into /usr (requires write access)
             sudo cmake --install .

         Unlike the local user installation, no additional environment
         variables need to be set.  The system locations such as ``/usr/lib`` and
         ``/usr/lib64`` are already part of the search path of the dynamic library
         loader.  Therefore ``LD_LIBRARY_PATH`` or ``DYLD_LIBRARY_PATH`` on MacOS do not
         have be set.

         All other environment variables will be automatically set when
         launching a new shell.  This is due to files installed in system folders
         ``/etc/profile.d/``, such as ``/etc/profile.d/lammps.sh``, that are loaded when a
         login shell is started.

      2. Open a new shell

         Close the current shell and open a new one or use ``source /etc/profile`` to
         update your environment

         .. note::

            On some systems you might also need to log out your current user and log back in.

      3. Verify if LAMMPS can be successfully started from Python

         Open a new terminal and test if LAMMPS can be started from within Python:

         .. code-block:: bash

            $ python
            Python 3.8.5 (default, Sep  5 2020, 10:50:12)
            [GCC 10.2.0] on linux
            Type "help", "copyright", "credits" or "license" for more information.
            >>> import lammps
            >>> lmp = lammps.lammps()
            LAMMPS (18 Sep 2020)
              using 1 OpenMP thread(s) per MPI task
            >>>

      .. note::

         If you recompile LAMMPS, you will have to also rerun the install step to
         ensure the latest Python module and shared library are installed.

   .. tab:: CMake (virtual environment)

      LAMMPS and its Python module can be installed together into a
      Python virtual environment.

      A virtual environment is a minimalistic Python installation inside of a
      folder. It allows isolating and customizing a Python environment that is
      independent from a user or system installation. This gives you the flexibility
      to install (newer) versions of Python packages that would potentially conflict
      with already installed system packages. It also does not requite any superuser
      privileges. See `PEP 405: Python Virtual Environments <python_pep405>`_
      for more information.

      To install into the virtual environment, it is first activated and the
      ``CMAKE_INSTALL_PREFIX`` is set to value of the ``$VIRTUAL_ENV`` environment
      variable. This leads to the following installation locations:

      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                  | Notes                                                       |
      +========================+===========================================================+=============================================================+
      | LAMMPS Python Module   | * ``$VIRTUAL_ENV/lib/pythonX.Y/site-packages/`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$VIRTUAL_ENV/lib64/pythonX.Y/site-packages/`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$VIRTUAL_ENV/lib/`` (32bit)                           |                                                             |
      |                        | * ``$VIRTUAL_ENV/lib64/`` (64bit)                         |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | ``$VIRTUAL_ENV/share/lammps/potentials/``                 |                                                             |
      +------------------------+-----------------------------------------------------------+-------------------------------------------------------------+

      The following is a minimal working example using CMake:

      1. Create a virtual environment

         Use the `venv <python_venv>`_ module to create a new environment
         inside of the folder ``$HOME/myenv``. For Python versions prior 3.3,
         you can use `virtualenv <python_virtualenv>`_ instead.

         .. code-block:: bash

            # create virtual environment in folder $HOME/myenv
            python3 -m venv $HOME/myenv

      2. Modify the ``$HOME/myenv/bin/activate`` script

         The ``activate`` script initializes the environment for use. For convienience,
         add two additional lines at the end of this script:

         * To allow the dynamic library loader to find the LAMMPS shared library, add
           the folder where it will be installed to ``LD_LIBRARY_PATH`` environment
           variable (``DYLD_LIBRARY_PATH`` on MacOS). When installing LAMMPS into a
           virtual environment this location will be ``$VIRTUAL_ENV/lib``.
           Run the following command to add the necessary line to the ``activate`` script:

           .. code-block:: bash

              # Unix/Linux
              echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> $HOME/myenv/bin/activate

              # MacOS
              echo 'export DYLD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> $HOME/myenv/bin/activate

         * Any LAMMPS installation will need to know the location of the folder containing its potential files.
           This can be set with the ``LAMMPS_POTENTIALS`` environment variable. When installing LAMMPS into a
           virtual environment this location will be ``$VIRTUAL_ENV/share/lammps/potentials``.
           Run the following command to add the change in the ``activate`` script:

           .. code-block:: bash

              echo 'export LAMMPS_POTENTIALS=$VIRTUAL_ENV/share/lammps/potentials' >> $HOME/myenv/bin/activate

      3. Compile LAMMPS and install it into virtual environment

         .. code-block:: bash

            # create and change into build directory
            mkdir build
            cd build

            # activate environment, this sets VIRTUAL_ENV and other environment variables
            source $HOME/myenv/bin/activate

            # configure LAMMPS compilation
            # compile with shared library, PYTHON package, and C++ exceptions
            # and install into virtual environment folder (VIRTUAL_ENV)
            # TODO: add more options to customize your LAMMPS installation
            (myenv)$ cmake -C ../cmake/presets/minimal.cmake    \
                           -D BUILD_SHARED_LIBS=on              \
                           -D PKG_PYTHON=on                     \
                           -D LAMMPS_EXCEPTIONS=on              \
                           -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
                           ../cmake

             # compile LAMMPS (in parallel for faster builds)
             (myenv)$ cmake --build . --parallel

             # install LAMMPS into myenv
             (myenv)$ cmake --install .

      4. Verify if LAMMPS can be successfully started from Python

         .. code-block:: bash

            (myenv)$ python
            Python 3.8.5 (default, Sep  5 2020, 10:50:12)
            [GCC 10.2.0] on linux
            Type "help", "copyright", "credits" or "license" for more information.
            >>> import lammps
            >>> lmp = lammps.lammps()
            LAMMPS (18 Sep 2020)
              using 1 OpenMP thread(s) per MPI task
            >>>

      .. note::

         If you recompile LAMMPS, you will have to also rerun the install step to
         ensure the virtual environment contains the latest Python module and shared
         library.


   .. tab:: Traditional make

      Instructions on how to build LAMMPS as a shared library are given on
      the :doc:`Build_basics <Build_basics>` doc page.  A shared library is
      one that is dynamically loadable, which is what Python requires to
      wrap LAMMPS.  On Linux this is a library file that ends in ``.so``, not
      ``.a``.

      From the src directory, type

      .. code-block:: bash

         make foo mode=shared

      where ``foo`` is the machine target name, such as ``mpi`` or ``serial``.
      This should create the file ``liblammps_foo.so`` in the ``src`` directory, as
      well as a soft link ``liblammps.so``, which is what the Python wrapper will
      load by default.  Note that if you are building multiple machine
      versions of the shared library, the soft link is always set to the
      most recently built version.

      .. note::

         If you are building LAMMPS with an MPI or FFT library or other
         auxiliary libraries (used by various packages), then all of these
         extra libraries must also be shared libraries.  If the LAMMPS
         shared-library build fails with an error complaining about this, see
         the :doc:`Build_basics <Build_basics>` doc page.

      You can achieve that Python can find these files in one of two ways:

      * set two environment variables pointing to the location in the source tree
      * run ``make install-python`` or run the ``python/install.py`` script explicitly

      When calling ``make install-python`` LAMMPS will try to install the
      python module and the shared library into the python site-packages folders;
      either the system-wide ones, or the local users ones (in case of insufficient
      permissions for the global install). Python will then find the module
      and shared library file automatically. The exact location of these folders
      depends on your python version and your operating system.

      You can override the python version to version x.y when calling
      ``make`` with ``PYTHON=pythonX.Y``.

      If you set the paths to these files as environment variables, you only
      have to do it once.  For the csh or tcsh shells, add something like
      this to your ~/.cshrc file, one line for each of the two files:

      .. code-block:: csh

         setenv PYTHONPATH ${PYTHONPATH}:/home/sjplimp/lammps/python
         setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home/sjplimp/lammps/src

      On MacOS you may also need to set ``DYLD_LIBRARY_PATH`` accordingly.
      For Bourne/Korn shells accordingly into the corresponding files using
      the ``export`` shell builtin.

      If you use ``make install-python`` or the ``python/install.py`` script, you need
      to invoke it every time you rebuild LAMMPS (as a shared library) or
      make changes to the ``python/lammps.py`` file, so that the site-packages
      files are updated with the new version.

      If the default settings of ``make install-python`` are not what you want,
      you can invoke ``install.py`` from the python directory manually as

      .. code-block:: bash

         $ python install.py -m <python module> -l <shared library> -v <version.h file> [-d <pydir>]

      * The ``-m`` flag points to the ``lammps.py`` python module file to be installed,
      * the ``-l`` flag points to the LAMMPS shared library file to be installed,
      * the ``-v`` flag points to the ``version.h`` file in the LAMMPS source
      * and the optional ``-d`` flag to a custom (legacy) installation folder

      If you use a legacy installation folder, you will need to set your
      ``PYTHONPATH`` and ``LD_LIBRARY_PATH`` (and/or ``DYLD_LIBRARY_PATH``) environment
      variables accordingly, as described above.

      Note that if you want Python to be able to load different versions of
      the LAMMPS shared library (see :ref:`python_create_lammps`), you will
      need to manually copy files like ``liblammps_mpi.so`` into the appropriate
      system directory.  This is not needed if you set the ``LD_LIBRARY_PATH``
      environment variable as described above.


Extending Python to run in parallel
===================================

If you wish to run LAMMPS in parallel from Python, you need to extend
your Python with an interface to MPI.  This also allows you to
make MPI calls directly from Python in your script, if you desire.

We have tested this with mpi4py and pypar:

* `MPI for Python <https://mpi4py.readthedocs.io/>`_
* `pypar <https://github.com/daleroberts/pypar>`_

We recommend the use of mpi4py as it is the more complete MPI interface,
and as of version 2.0.0 mpi4py allows passing a custom MPI communicator
to the LAMMPS constructor, which means one can easily run one or more
LAMMPS instances on subsets of the total MPI ranks.

To install mpi4py (version 3.0.3 as of Sep 2020),

.. tabs::

   .. tab:: local user

      .. code-block:: bash

         pip install --user mpi4py

   .. tab:: system-wide

      .. code-block:: bash

         sudo pip install mpi4py

   .. tab:: virtual environment

      .. code-block:: bash

         $ source $HOME/myenv/activate
         (myenv)$ pip install mpi4py

.. _mpi4py_install: https://mpi4py.readthedocs.io/en/stable/install.html

For more detailed installation instructions, please see the `mpi4py installation <mpi4py_install>`_ page.

If you have successfully installed mpi4py, you should be able to run
Python and type

.. code-block:: python

   from mpi4py import MPI

without error.  You should also be able to run Python in parallel
on a simple test script

.. code-block:: bash

   $ mpirun -np 4 python test.py

where ``test.py`` contains the lines

.. code-block:: python

   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())

and see one line of output for each processor you run on.

.. note::

   To use mpi4py and LAMMPS in parallel from Python, you must
   insure both are using the same version of MPI.  If you only have one
   MPI installed on your system, this is not an issue, but it can be if
   you have multiple MPIs.  Your LAMMPS build is explicit about which MPI
   it is using, since you specify the details in your low-level
   src/MAKE/Makefile.foo file.  mpi4py uses the "mpicc" command to find
   information about the MPI it uses to build against.  And it tries to
   load "libmpi.so" from the ``LD_LIBRARY_PATH``.  This may or may not find
   the MPI library that LAMMPS is using.  If you have problems running
   both mpi4py and LAMMPS together, this is an issue you may need to
   address, e.g. by moving other MPI installations so that mpi4py finds
   the right one.

