Executing LAMMPS commands
=========================

Once a LAMMPS instance is created, there are multiple ways to "drive" a
simulation.  In most cases it is easiest to process single or multiple
LAMMPS commands like in an input file.  This can be done through reading
a file or passing single commands or lists of commands or blocks of
commands with the following functions.

Via these functions, the calling code can have the LAMMPS instance act
on a series of :doc:`input file commands <Commands_all>` that are either
read from a file or passed as strings.  This for, for example, allows to
setup a problem from a template file and then run it in stages while
performing other operations in between or concurrently.  The caller can
interleave the LAMMPS function calls with operations it performs, calls
to extract information from or set information within LAMMPS, or calls
to another code's library.

Also equivalent to regular :doc:`input script parsing <Commands_parse>`
is the handling of comments and expansion of variables with ``${name}``
or ``$(expression)`` syntax before the commands are parsed and
executed. Below is a short example using some of these functions.

.. code-block:: C

   #include <library.h>
   #include <mpi.h>
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int i;

     MPI_Init(&argc, &argv);
     lammps_open(0, NULL, MPI_COMM_WORLD, &handle);
     lammps_file(handle,"in.sysinit");
     lammps_command(handle,"run 1000 post no");

     for (i=0; i < 100; ++i) {
       lammps_commands_string(handle,"run 100 pre no post no\n"
                                     "print 'PE = $(pe)'\n"
                                     "print 'KE = $(ke)'\n")
     }
     lammps_close(handle);
     MPI_Finalize();
     return 0;
   }

-----------------------

.. doxygenfunction:: lammps_file
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_command
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_commands_list
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_commands_string
   :project: progguide

