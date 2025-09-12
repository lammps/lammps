Platform abstraction functions
------------------------------

The ``platform`` sub-namespace inside the ``LAMMPS_NS`` namespace
provides a collection of wrapper and convenience functions and utilities
that perform common tasks for which platform specific code would be
required or for which a more high-level abstraction would be convenient
and reduce duplicated code.  This reduces redundant implementations and
encourages consistent behavior and thus has some overlap with the
:doc:`"utils" sub-namespace <Developer_utils>`.

Time functions
^^^^^^^^^^^^^^

.. doxygenfunction:: cputime
   :project: progguide

.. doxygenfunction:: walltime
   :project: progguide

.. doxygenfunction:: usleep
   :project: progguide

Platform information functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: os_info
   :project: progguide

.. doxygenfunction:: compiler_info
   :project: progguide

.. doxygenfunction:: cxx_standard
   :project: progguide

.. doxygenfunction:: openmp_standard
   :project: progguide

.. doxygenfunction:: mpi_vendor
   :project: progguide

.. doxygenfunction:: mpi_info
   :project: progguide

.. doxygenfunction:: compress_info
   :project: progguide


File and path functions and global constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since we are requiring C++17 to compile LAMMPS, you can also make use of
the functionality of the `C++ filesystem library
<https://en.cppreference.com/w/cpp/filesystem.html>`_.  The following
functions are in part convenience functions or emulate the behavior of
similar Python functions or Unix shell commands.  Please note that the
you need to use the ``string()`` member function of the
`std::filesystem::path class
<https://en.cppreference.com/w/cpp/filesystem/path.html>`_ to get access
to the path as a C++ string class instance.

.. doxygenvariable:: filepathsep
   :project: progguide

.. doxygenvariable:: pathvarsep
   :project: progguide

.. doxygenfunction:: guesspath
   :project: progguide

.. doxygenfunction:: path_basename
   :project: progguide

.. doxygenfunction:: path_dirname
   :project: progguide

.. doxygenfunction:: path_join
   :project: progguide

.. doxygenfunction:: file_is_readable
   :project: progguide

.. doxygenfunction:: is_console
   :project: progguide

.. doxygenfunction:: disk_free
   :project: progguide

.. doxygenfunction:: list_directory
   :project: progguide

.. doxygenfunction:: chdir
   :project: progguide

.. doxygenfunction:: mkdir
   :project: progguide

.. doxygenfunction:: rmdir
   :project: progguide

.. doxygenfunction:: unlink
   :project: progguide

Standard I/O function wrappers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenvariable:: END_OF_FILE
   :project: progguide

.. doxygenfunction:: ftell
   :project: progguide

.. doxygenfunction:: fseek
   :project: progguide

.. doxygenfunction:: ftruncate
   :project: progguide

.. doxygenfunction:: popen
   :project: progguide

.. doxygenfunction:: pclose
   :project: progguide

Environment variable functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: putenv
   :project: progguide

.. doxygenfunction:: unsetenv
   :project: progguide

.. doxygenfunction:: list_pathenv
   :project: progguide

.. doxygenfunction:: find_exe_path
   :project: progguide

Dynamically loaded object or library functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: dlopen
   :project: progguide

.. doxygenfunction:: dlclose
   :project: progguide

.. doxygenfunction:: dlsym
   :project: progguide

.. doxygenfunction:: dlerror
   :project: progguide

Compressed file I/O functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: has_compress_extension
   :project: progguide

.. doxygenfunction:: compressed_read
   :project: progguide

.. doxygenfunction:: compressed_write
   :project: progguide
