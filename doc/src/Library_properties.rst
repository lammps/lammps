Retrieving or setting LAMMPS system properties
==============================================

The library interface allows to extract different kinds of information
about the active simulation instance and also to modify some of it.
This enables combining MD simulation steps with other processing and
simulation methods computed in the calling code or another code that is
coupled to LAMMPS via the library interface.  In some cases the data
returned is direct reference to the original data inside LAMMPS cast
to a void pointer.  In that case the data needs to be cast to a suitable
pointer to be able to access it, and you need to know the correct dimensions
and lengths.  When accessing per-atom data, please note that this data
is the per-processor **local** data and indexed accordingly. These arrays
can change sizes and order at every neighbor list rebuild and atom sort
event as atoms are migrating between sub-domains.

.. code-block:: C

   #include "library.h"
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int i;

     handle = lammps_open_no_mpi(0, NULL, NULL);
     lammps_file(handle,"in.sysinit");
     printf("running simulation with %g atoms\n",
            lammps_get_natoms(handle));

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

.. doxygenfunction:: lammps_version
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_memory_usage
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_mpi_comm
   :project: progguide

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

-------------------

.. doxygenfunction:: lammps_extract_setting
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global_datatype
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_atom_datatype
   :project: progguide


-----------------------

.. doxygenfunction:: lammps_extract_atom
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_create_atoms(void *handle, int n, int *id, int *type, double *x, double *v, int *image, int bexpand)
   :project: progguide


