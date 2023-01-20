System properties
=================

This section documents the following functions:

- :cpp:func:`lammps_get_natoms`
- :cpp:func:`lammps_get_thermo`
- :cpp:func:`lammps_extract_box`
- :cpp:func:`lammps_reset_box`
- :cpp:func:`lammps_memory_usage`
- :cpp:func:`lammps_get_mpi_comm`
- :cpp:func:`lammps_extract_setting`
- :cpp:func:`lammps_extract_global_datatype`
- :cpp:func:`lammps_extract_global`

--------------------

The library interface allows the extraction of different kinds of
information about the active simulation instance and also - in some
cases - to apply modifications to it.  This enables combining of a
LAMMPS simulation with other processing and simulation methods computed
by the calling code, or by another code that is coupled to LAMMPS via
the library interface.  In some cases the data returned is direct
reference to the original data inside LAMMPS, cast to a void pointer.
In that case the data needs to be cast to a suitable pointer for the
calling program to access it, and you may need to know the correct
dimensions and lengths.  This also means you can directly change those
value(s) from the calling program (e.g., to modify atom positions).  Of
course, changing values should be done with care.  When accessing per-atom
data, please note that these data are the per-processor **local** data and are
indexed accordingly. Per-atom data can change sizes and ordering at
every neighbor list rebuild or atom sort event as atoms migrate between
sub-domains and processors.

.. code-block:: c

   #include "library.h"
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int i;

     handle = lammps_open_no_mpi(0, NULL, NULL);
     lammps_file(handle,"in.sysinit");
     printf("Running a simulation with %g atoms.\n",
            lammps_get_natoms(handle));

     printf(" %d local and %d ghost atoms. %d atom types\n",
            lammps_extract_setting(handle,"nlocal"),
            lammps_extract_setting(handle,"nghost"),
            lammps_extract_setting(handle,"ntypes"));

     double  *dt = (double *)lammps_extract_global(handle,"dt");
     printf("Changing timestep from %g to 0.5\n", *dt);
     *dt = 0.5;

     lammps_command(handle,"run 1000 post no");

     for (i=0; i < 10; ++i) {
       lammps_command(handle,"run 100 pre no post no");
       printf("PE = %g\nKE = %g\n",
              lammps_get_thermo(handle,"pe"),
              lammps_get_thermo(handle,"ke"));
     }
     lammps_close(handle);
     return 0;
   }


-----------------------

.. doxygenfunction:: lammps_get_natoms
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_thermo
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_box
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_reset_box
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_memory_usage
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_mpi_comm
   :project: progguide

-------------------

.. doxygenfunction:: lammps_extract_setting
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global_datatype
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global
   :project: progguide

