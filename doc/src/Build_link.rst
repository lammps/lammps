Link LAMMPS as a library to another code
========================================

LAMMPS is designed as a library of C++ objects that can be integrated
into other applications, including Python scripts.  The files
``src/library.cpp`` and ``src/library.h`` define a C-style API for using
LAMMPS as a library.  See the :doc:`Howto_library` page for a
description of the interface and how to use it for your needs.

The :doc:`Build_basics` page explains how to build LAMMPS as either a
shared or static library.  This results in a file in the compilation
folder called ``liblammps.a`` or ``liblammps_<name>.a`` in case of
building a static library.  In case of a shared library, the name is the
same only that the suffix is going to be either ``.so`` or ``.dylib`` or
``.dll`` instead of ``.a`` depending on the OS.  In some cases, the
``.so`` file may be a symbolic link to a file with the suffix ``.so.0``
(or some other number).

.. note::

   Care should be taken to use the same MPI library for the calling code
   and the LAMMPS library, unless LAMMPS is to be compiled without (real)
   MPI support using the included STUBS MPI library.

Link with LAMMPS as a static library
------------------------------------

The calling application can link to LAMMPS as a static library with
compilation and link commands, as in the examples shown below.  These
are examples for a code written in C in the file ``caller.c``.
The benefit of linking to a static library is, that the resulting
executable is independent of that library since all required
executable code from the library is copied into the calling executable.

CMake installs the LAMMPS library, its header files, and the
``liblammps.pc`` file for ``pkg-config`` only when LAMMPS is built as a
shared library.  To link with a static library, you need to use the
header files in the ``src`` folder and the library file in the build
folder directly, for example:

.. code-block:: bash

   mpicc -c -O -I${HOME}/lammps/src caller.c
   mpicxx -o caller caller.o -L${HOME}/lammps/build -llammps

However, it is only as simple as shown above for the case of a plain
LAMMPS library without any optional packages that depend on libraries
(bundled or external).  Otherwise, you need to include all flags,
libraries, and paths for the coupled executable, that are also required
to link the LAMMPS executable.  Additional libraries with sources in the
``lib`` folder are built as separate static libraries in the build
folder and are not included in ``liblammps.a``, so you have to determine
the necessary link flags manually.  It is therefore recommended to build
and link with a shared library instead.

Link with LAMMPS as a shared library
------------------------------------

When linking to LAMMPS built as a shared library, the situation becomes
much simpler, as all dependent libraries and objects are either included
in the shared library or registered as a dependent library in the shared
library file.  Thus, those libraries need not be specified when linking
the calling executable.  This assumes that LAMMPS has been configured
without setting a ``LAMMPS_MACHINE`` name, installed with ``cmake
--install build``, and the ``PKG_CONFIG_PATH`` environment variable has
been updated to include the folder with the ``liblammps.pc`` file in the
configured destination folder.  The commands to compile and link a
coupled executable are then:

.. code-block:: bash

   mpicc -c -O $(pkg-config --cflags liblammps) caller.c
   mpicxx -o caller caller.o $(pkg-config --libs liblammps)

Locating liblammps.so at runtime
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike with a static link, now the ``liblammps.so`` file is required at
runtime and needs to be in a folder, where the shared linker program of
the operating system can find it.  This would be either a folder like
``/usr/local/lib64`` or ``${HOME}/.local/lib64`` or a folder pointed to
by the ``LD_LIBRARY_PATH`` environment variable. You can type

.. code-block:: bash

   printenv LD_LIBRARY_PATH

to see what directories are in that list.

Or you can add the folder you performed the CMake build in to your
``LD_LIBRARY_PATH`` environment variable,
so that the current version of the shared library is always available to
programs that use it.

For the Bourne or Korn shells (/bin/sh, /bin/ksh, /bin/bash etc.), you
would add something like this to your ``${HOME}/.profile`` file:

.. code-block:: bash

   LD_LIBRARY_PATH=${LD_LIBRARY_PATH-/usr/lib64}:${HOME}/lammps/build
   export LD_LIBRARY_PATH

For the csh or tcsh shells, you would equivalently add something like this
to your ``${HOME}/.cshrc`` file:

.. code-block:: csh

   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HOME}/lammps/build

You can verify whether all required shared libraries are found with the
``ldd`` tool.  Example:

.. code-block:: bash

   LD_LIBRARY_PATH=/home/user/lammps/build ldd caller
        linux-vdso.so.1 (0x00007ffe729e0000)
        liblammps.so => /home/user/lammps/build/liblammps.so (0x00007fc91bb9e000)
        libstdc++.so.6 => /lib64/libstdc++.so.6 (0x00007fc91b984000)
        libm.so.6 => /lib64/libm.so.6 (0x00007fc91b83e000)
        libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007fc91b824000)
        libc.so.6 => /lib64/libc.so.6 (0x00007fc91b65b000)
        /lib64/ld-linux-x86-64.so.2 (0x00007fc91c094000)

If a required library is missing, you would get a 'not found' entry:

.. code-block:: bash

   ldd caller
        linux-vdso.so.1 (0x00007ffd672fe000)
        liblammps.so => not found
        libstdc++.so.6 => /usr/lib64/libstdc++.so.6 (0x00007fb7c7e86000)
        libm.so.6 => /usr/lib64/libm.so.6 (0x00007fb7c7d40000)
        libgcc_s.so.1 => /usr/lib64/libgcc_s.so.1 (0x00007fb7c7d26000)
        libc.so.6 => /usr/lib64/libc.so.6 (0x00007fb7c7b5d000)
        /lib64/ld-linux-x86-64.so.2 (0x00007fb7c80a2000)
