Installation
************

The LAMMPS Python module enables calling the :ref:`LAMMPS C library API
<lammps_c_api>` from Python by dynamically loading functions in the
LAMMPS shared library through the Python `ctypes <ctypes_>`_
module.  Because of the dynamic loading, it is required that LAMMPS is
compiled in :ref:`"shared" mode <exe>`.  It is also recommended to
compile LAMMPS with :ref:`C++ exceptions <exceptions>` enabled.

Two components are necessary for Python to be able to invoke LAMMPS code:

* The LAMMPS Python Package (``lammps``) from the ``python`` folder
* The LAMMPS Shared Library (``liblammps.so``, ``liblammps.dylib`` or
  ``liblammps.dll``) from the folder where you compiled LAMMPS.

.. _ctypes: https://docs.python.org/3/library/ctypes.html
.. _python_virtualenv: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment
.. _python_venv: https://docs.python.org/3/library/venv.html
.. _python_pep405: https://www.python.org/dev/peps/pep-0405

.. _python_install_guides:

Installing the LAMMPS Python Module and Shared Library
======================================================

Making LAMMPS usable within Python and vice versa requires putting the
LAMMPS Python package (``lammps``) into a location where the
Python interpreter can find it and installing the LAMMPS shared library
into a folder that the dynamic loader searches or inside of the installed
``lammps`` package folder.  There are multiple ways to achieve
this.

#. Do a full LAMMPS installation of libraries, executables, selected
   headers, documentation (if enabled), and supporting files (only
   available via CMake), which can also be either system-wide or into
   user specific folders.

#. Install both components into a Python ``site-packages`` folder, either
   system-wide or in the corresponding user-specific folder. This way no
   additional environment variables need to be set, but the shared
   library is otherwise not accessible.

#. Do an installation into a virtual environment. This can either be an
   installation of the Python package only or a full installation of LAMMPS.

#. Leave the files where they are in the source/development tree and
   adjust some environment variables.

.. tabs::

   .. tab:: Full install (CMake-only)

      :ref:`Build the LAMMPS executable and library <library>` with
      ``-DBUILD_SHARED_LIBS=on``, ``-DLAMMPS_EXCEPTIONS=on`` and
      ``-DPKG_PYTHON=on`` (The first option is required, the other two
      are optional by recommended).  The exact file name of the shared
      library depends on the platform (Unix/Linux, MacOS, Windows) and
      the build configuration being used.  The installation base folder
      is already set by default to the ``$HOME/.local`` directory, but
      it can be changed to a custom location defined by the
      ``CMAKE_INSTALL_PREFIX`` CMake variable.  This uses a folder
      called ``build`` to store files generated during compilation.

      .. code-block:: bash

         # create build folder
         mkdir build
         cd build

         # configure LAMMPS compilation
         cmake -C ../cmake/presets/basic.cmake -D BUILD_SHARED_LIBS=on \
               -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on ../cmake

         # compile LAMMPS
         cmake --build .

         # install LAMMPS into $HOME/.local
         cmake --install .


      This leads to an installation to the following locations:

      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                        | Notes                                                       |
      +========================+=================================================================+=============================================================+
      | LAMMPS Python package  | * ``$HOME/.local/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$HOME/.local/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$HOME/.local/lib/`` (32bit)                                 | Set shared loader environment variable to this path         |
      |                        | * ``$HOME/.local/lib64/`` (64bit)                               | (see below for more info on this)                           |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS executable      | * ``$HOME/.local/bin/``                                         |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | * ``$HOME/.local/share/lammps/potentials/``                     | Set ``LAMMPS_POTENTIALS`` environment variable to this path |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+

      For a system-wide installation you need to set
      ``CMAKE_INSTALL_PREFIX`` to a system folder like ``/usr`` (or
      ``/usr/local``); the default is ``${HOME}/.local``.  The
      installation step for a system folder installation (**not** the
      configuration/compilation) needs to be done with superuser
      privilege, e.g. by using ``sudo cmake --install .``.  The
      installation folders will then be changed to (assuming ``/usr`` as
      prefix):

      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                | Notes                                                       |
      +========================+=========================================================+=============================================================+
      | LAMMPS Python package  | * ``/usr/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``/usr/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``/usr/lib/`` (32bit)                                 |                                                             |
      |                        | * ``/usr/lib64/`` (64bit)                               |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS executable      | * ``/usr/bin/``                                         |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | * ``/usr/share/lammps/potentials/``                     |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+

      To be able to use the "user" installation you have to ensure that
      the folder containing the LAMMPS shared library is either included
      in a path searched by the shared linker (e.g. like
      ``/usr/lib64/``) or part of the ``LD_LIBRARY_PATH`` environment
      variable (or ``DYLD_LIBRARY_PATH`` on MacOS).  Otherwise you will
      get an error when trying to create a LAMMPS object through the
      Python module.

      .. code-block:: bash

         # Unix/Linux
         export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH

         # MacOS
         export DYLD_LIBRARY_PATH=$HOME/.local/lib:$DYLD_LIBRARY_PATH

      If you plan to use the LAMMPS executable (e.g., ``lmp``), you may
      also need to adjust the ``PATH`` environment variable (but many
      newer Linux distributions already have ``$HOME/.local/bin``
      included). Example:

      .. code-block:: bash

         export PATH=$HOME/.local/bin:$PATH

      To make those changes permanent, you can add the commands to your
      ``$HOME/.bashrc`` file.  For a system-wide installation is is not
      necessary due to files installed in system folders that are loaded
      automatically when a login shell is started.

   .. tab:: Python package only

      Compile LAMMPS with either :doc:`CMake <Build_cmake>` or the
      :doc:`traditional make <Build_make>` procedure in :ref:`shared
      mode <exe>`.  After compilation has finished type (in the
      compilation folder):

      .. code-block:: bash

         make install-python

      This will try to install (only) the shared library and the Python
      package into a system folder and if that fails (due to missing
      write permissions) will instead do the installation to a user
      folder under ``$HOME/.local``.  For a system-wide installation you
      would have to gain superuser privilege, e.g. though ``sudo``

      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                        | Notes                                                       |
      +========================+=================================================================+=============================================================+
      | LAMMPS Python package  | * ``$HOME/.local/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$HOME/.local/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$HOME/.local/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$HOME/.local/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+

      For a system-wide installation those folders would then become.

      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                | Notes                                                       |
      +========================+=========================================================+=============================================================+
      | LAMMPS Python package  | * ``/usr/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``/usr/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``/usr/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``/usr/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+---------------------------------------------------------+-------------------------------------------------------------+

      No environment variables need to be set for those, as those
      folders are searched by default by Python or the LAMMPS Python
      package.

      For the traditional make process you can override the python
      version to version x.y when calling ``make`` with
      ``PYTHON=pythonX.Y``.  For a CMake based compilation this choice
      has to be made during the CMake configuration step.

      If the default settings of ``make install-python`` are not what you want,
      you can invoke ``install.py`` from the python directory manually as

      .. code-block:: bash

         $ python install.py -p <python package> -l <shared library> -v <version.h file> [-d <pydir>]

      * The ``-p`` flag points to the ``lammps`` Python package folder to be installed,
      * the ``-l`` flag points to the LAMMPS shared library file to be installed,
      * the ``-v`` flag points to the ``version.h`` file in the LAMMPS source
      * and the optional ``-d`` flag to a custom (legacy) installation folder

      If you use a legacy installation folder, you will need to set your
      ``PYTHONPATH`` and ``LD_LIBRARY_PATH`` (and/or ``DYLD_LIBRARY_PATH``) environment
      variables accordingly as explained in the description for "In place use".

   .. tab:: Virtual environment

      A virtual environment is a minimal Python installation inside of a
      folder.  It allows isolating and customizing a Python environment
      that is mostly independent from a user or system installation.
      For the core Python environment, it uses symbolic links to the
      system installation and thus it can be set up quickly and will not
      take up much disk space.  This gives you the flexibility to
      install (newer/different) versions of Python packages that would
      potentially conflict with already installed system packages.  It
      also does not requite any superuser privileges. See `PEP 405:
      Python Virtual Environments <python_pep405>`_ for more
      information.

      To create a virtual environment in the folder ``$HOME/myenv``,
      use the `venv <python_venv>`_ module as follows.

      .. code-block:: bash

         # create virtual environment in folder $HOME/myenv
         python3 -m venv $HOME/myenv

      For Python versions prior 3.3 you can use `virtualenv
      <python_virtualenv>`_ command instead of "python3 -m venv".  This
      step has to be done only once.

      To activate the virtual environment type:

      .. code-block:: bash

         source $HOME/myenv/bin/activate

      This has to be done every time you log in or open a new terminal
      window and after you turn off the virtual environment with the
      ``deactivate`` command.

      When using CMake to build LAMMPS, you need to set
      ``CMAKE_INSTALL_PREFIX`` to the value of the ``$VIRTUAL_ENV``
      environment variable during the configuration step. For the
      traditional make procedure, no additional steps are needed.
      After compiling LAMMPS you can do a "Python package only"
      installation with ``make install-python`` and the LAMMPS Python
      package and the shared library file are installed into the
      following locations:

      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                        | Notes                                                       |
      +========================+=================================================================+=============================================================+
      | LAMMPS Python Module   | * ``$VIRTUAL_ENV/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$VIRTUAL_ENV/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$VIRTUAL_ENV/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$VIRTUAL_ENV/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+

      If you do a full installation (CMake only) with "install", this
      leads to the following installation locations:

      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | File                   | Location                                                        | Notes                                                       |
      +========================+=================================================================+=============================================================+
      | LAMMPS Python Module   | * ``$VIRTUAL_ENV/lib/pythonX.Y/site-packages/lammps`` (32bit)   | ``X.Y`` depends on the installed Python version             |
      |                        | * ``$VIRTUAL_ENV/lib64/pythonX.Y/site-packages/lammps`` (64bit) |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS shared library  | * ``$VIRTUAL_ENV/lib/`` (32bit)                                 | Set shared loader environment variable to this path         |
      |                        | * ``$VIRTUAL_ENV/lib64/`` (64bit)                               | (see below for more info on this)                           |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS executable      | * ``$VIRTUAL_ENV/bin/``                                         |                                                             |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+
      | LAMMPS potential files | * ``$VIRTUAL_ENV/share/lammps/potentials/``                     | Set ``LAMMPS_POTENTIALS`` environment variable to this path |
      +------------------------+-----------------------------------------------------------------+-------------------------------------------------------------+

      In that case you need to modify the ``$HOME/myenv/bin/activate``
      script in a similar fashion you need to update your
      ``$HOME/.bashrc`` file to include the shared library and
      executable locations in ``LD_LIBRARY_PATH`` (or
      ``DYLD_LIBRARY_PATH`` on MacOS) and ``PATH``, respectively.

      For example with:

      .. code-block:: bash

         # Unix/Linux
         echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> $HOME/myenv/bin/activate

         # MacOS
         echo 'export DYLD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$DYLD_LIBRARY_PATH' >> $HOME/myenv/bin/activate

   .. tab:: In place usage

      You can also :doc:`compile LAMMPS <Build>` as usual in
      :ref:`"shared" mode <exe>` leave the shared library and Python
      package inside the source/compilation folders. Instead of
      copying the files where they can be found, you need to set the environment
      variables ``PYTHONPATH`` (for the Python package) and
      ``LD_LIBRARY_PATH`` (or ``DYLD_LIBRARY_PATH`` on MacOS

      For Bourne shells (bash, ksh and similar) the commands are:

      .. code-block:: bash

         export PYTHONPATH=${PYTHONPATH}:${HOME}/lammps/python
         export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/lammps/src

      For the C-shells like csh or tcsh the commands are:

      .. code-block:: csh

         setenv PYTHONPATH ${PYTHONPATH}:${HOME}/lammps/python
         setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HOME}/lammps/src

      On MacOS you may also need to set ``DYLD_LIBRARY_PATH`` accordingly.
      You can make those changes permanent by editing your ``$HOME/.bashrc``
      or ``$HOME/.login`` files, respectively.

      .. note::

         The ``PYTHONPATH`` needs to point to the parent folder that contains the ``lammps`` package!


To verify if LAMMPS can be successfully started from Python, start the
Python interpreter, load the ``lammps`` Python module and create a
LAMMPS instance.  This should not generate an error message and produce
output similar to the following:

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

   Unless you opted for "In place use", you will have to rerun the installation
   any time you recompile LAMMPS to ensure the latest Python package and shared
   library are installed and used.

.. note::

   If you want Python to be able to load different versions of the
   LAMMPS shared library with different settings, you will need to
   manually copy the files under different names
   (e.g. ``liblammps_mpi.so`` or ``liblammps_gpu.so``) into the
   appropriate folder as indicated above. You can then select the
   desired library through the *name* argument of the LAMMPS object
   constructor (see :ref:`python_create_lammps`).

.. _python_install_mpi4py:

Extending Python to run in parallel
===================================

If you wish to run LAMMPS in parallel from Python, you need to extend
your Python with an interface to MPI.  This also allows you to
make MPI calls directly from Python in your script, if you desire.

We have tested this with `MPI for Python <https://mpi4py.readthedocs.io/>`_
(aka mpi4py) and you will find installation instruction for it below.

Installation of mpi4py (version 3.0.3 as of Sep 2020) can be done as
follows:

- Via ``pip`` into a local user folder with:

  .. code-block:: bash

     pip install --user mpi4py

- Via ``dnf`` into a system folder for RedHat/Fedora systems:

  .. code-block:: bash

     # for use with OpenMPI
     sudo dnf install python3-mpi4py-openmpi
     # for use with MPICH
     sudo dnf install python3-mpi4py-openmpi

- Via ``pip`` into a virtual environment (see above):

  .. code-block:: bash

     $ source $HOME/myenv/activate
     (myenv)$ pip install mpi4py

- Via ``pip`` into a system folder (not recommended):

  .. code-block:: bash

     sudo pip install mpi4py

.. _mpi4py_install: https://mpi4py.readthedocs.io/en/stable/install.html

For more detailed installation instructions and additional options,
please see the `mpi4py installation <mpi4py_install>`_ page.


To use ``mpi4py`` and LAMMPS in parallel from Python, you **must** make
certain that **both** are using the **same** implementation and version
of MPI library.  If you only have one MPI library installed on your
system this is not an issue, but it can be if you have multiple MPI
installations (e.g. on an HPC cluster to be selected through environment
modules).  Your LAMMPS build is explicit about which MPI it is using,
since it is either detected during CMake configuration or in the
traditional make build system you specify the details in your low-level
``src/MAKE/Makefile.foo`` file. The installation process of ``mpi4py``
uses the ``mpicc`` command to find information about the MPI it uses to
build against.  And it tries to load "libmpi.so" from the
``LD_LIBRARY_PATH``.  This may or may not find the MPI library that
LAMMPS is using.  If you have problems running both mpi4py and LAMMPS
together, this is an issue you may need to address, e.g. by loading the
module for different MPI installation so that mpi4py finds the right
one.

If you have successfully installed mpi4py, you should be able to run
Python and type

.. code-block:: python

   from mpi4py import MPI

without error.  You should also be able to run Python in parallel
on a simple test script

.. code-block:: bash

   $ mpirun -np 4 python3 test.py

where ``test.py`` contains the lines

.. code-block:: python

   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   print("Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size()))

and see one line of output for each processor you run on.

.. code-block:: bash

   # NOTE: the line order is not deterministic
   $ mpirun -np 4 python3 test.py
   Proc 0 out of 4 procs
   Proc 1 out of 4 procs
   Proc 2 out of 4 procs
   Proc 3 out of 4 procs

