Retrieving LAMMPS configuration information
===========================================

This section documents the following functions:

- :cpp:func:`lammps_config_has_mpi_support`
- :cpp:func:`lammps_config_has_gzip_support`
- :cpp:func:`lammps_config_has_png_support`
- :cpp:func:`lammps_config_has_jpeg_support`
- :cpp:func:`lammps_config_has_ffmpeg_support`
- :cpp:func:`lammps_config_has_exceptions`
- :cpp:func:`lammps_config_has_package`
- :cpp:func:`lammps_config_package_count`
- :cpp:func:`lammps_config_package_name`
- :cpp:func:`lammps_has_style`
- :cpp:func:`lammps_style_count`
- :cpp:func:`lammps_style_name`

--------------------

The following library functions can be used to query the LAMMPS library
about compile time settings and included packages and styles.  This
enables programs that use the library interface to run LAMMPS
simulations to determine, whether the linked LAMMPS library is compatible
with the requirements of the application without crashing during the
LAMMPS functions (e.g. due to missing pair styles from packages) or to
choose between different options (e.g. whether to use ``lj/cut``,
``lj/cut/opt``, ``lj/cut/omp`` or ``lj/cut/intel``).  Most of the
functions can be called directly without first creating a LAMMPS
instance.  While crashes within LAMMPS may be recovered from through
enabling :ref:`exceptions <exceptions>`, avoiding them proactively is
a safer approach.

.. code-block:: C
   :caption: Example for using configuration settings functions

   #include "library.h"
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;

     handle = lammps_open_no_mpi(0, NULL, NULL);
     lammps_file(handle, "in.missing");
     if (lammps_has_error(handle)) {
       char errmsg[256];
       int errtype;
       errtype = lammps_get_last_error_message(handle, errmsg, 256);
       fprintf(stderr, "LAMMPS failed with error: %s\n", errmsg);
       return 1;
     }
     /* write compressed dump file depending on available of options */
     if (lammps_has_style(handle, "dump", "atom/zstd")) {
       lammps_command(handle, "dump d1 all atom/zstd 100 dump.zst");
     } else if (lammps_has_style(handle, "dump", "atom/gz")) {
       lammps_command(handle, "dump d1 all atom/gz 100 dump.gz");
     } else if (lammps_config_has_gzip_support()) {
       lammps_command(handle, "dump d1 all atom 100 dump.gz");
     } else {
       lammps_command(handle, "dump d1 all atom 100 dump");
     }
     lammps_close(handle);
     return 0;
   }

-----------------------

.. doxygenfunction:: lammps_config_has_mpi_support
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_gzip_support
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_png_support
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_jpeg_support
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_ffmpeg_support
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_exceptions
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_has_package
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_package_count
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_config_package_name
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_has_style
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_style_count
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_style_name
   :project: progguide

