Using the C++ API directly
**************************

Using the C++ classes of the LAMMPS library is lacking some of the
convenience of the C library API, but it allows a more direct access to
simulation data and thus more low-level manipulations and tighter
integration of LAMMPS into another code.  While for the complete C
library API is provided in the ``library.h`` header file, for using
the C++ API it is required to include the individual header files
defining the individual classes in use.  Typically the name of the
class and the name of the header follow some simple rule.  Examples
are given below.


Creating or deleting a LAMMPS object
*************************************

When using the LAMMPS library interfaces, the core task is to create an
instance of the :cpp:class:`LAMMPS_NS::LAMMPS` class.  In C++ this can
be done directly through the ``new`` operator.  All further operations
are then initiated through calling member functions of some of the
components of the LAMMPS class or accessing their data members.  The
destruction of the LAMMPS instance is correspondingly initiated by using
the ``delete`` operator.  Here is a simple example:

.. code-block:: c++

   #include "lammps.h"

   #include <mpi.h>
   #include <iostream>

   int main(int argc, char **argv)
   {
       LAMMPS_NS::LAMMPS *lmp;
       // custom argument vector for LAMMPS library
       const char *lmpargv[] {"liblammps", "-log", "none"};
       int lmpargc = sizeof(lmpargv)/sizeof(const char *);

       // explicitly initialize MPI
       MPI_Init(&argc, &argv);

       // create LAMMPS instance
       lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
       // output numerical version string
       std::cout << "LAMMPS version ID: " << lmp->num_ver << std::endl;
       // delete LAMMPS instance
       delete lmp;

       // stop MPI environment
       MPI_Finalize();
       return 0;
   }

This minimal example only requires to include the ``lammps.h`` header
file since it only accesses a non-pointer member of the LAMMPS class.


Executing LAMMPS commands
*************************

Once a LAMMPS instance is created by your C++ code, you need to set up a
simulation and that is most conveniently done by "driving" it through
issuing commands like you would do when running a LAMMPS simulation from
an input script. Processing of input in LAMMPS is handled by the
:cpp:class:`Input <LAMMPS_NS::Input>` class an instance of which is a
member of the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class.  You have
two options: reading commands from a file, or executing a single
command from a string. See below for a small example:

.. code-block:: c++

   #include "lammps.h"
   #include "input.h"
   #include <mpi.h>

   using namespace LAMMPS_NS;

   int main(int argc, char **argv)
   {
       const char *lmpargv[] {"liblammps", "-log", "none"};
       int lmpargc = sizeof(lmpargv)/sizeof(const char *);

       MPI_Init(&argc, &argv);
       LAMMPS *lmp = new LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
       lmp->input->file("in.melt");
       lmp->input->one("run 100 post no");
       delete lmp;
       return 0;
   }
