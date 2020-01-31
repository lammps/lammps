Input script command style
==========================

New commands can be added to LAMMPS input scripts by adding new
classes that have a "command" method.  For example, the create\_atoms,
read\_data, velocity, and run commands are all implemented in this
fashion.  When such a command is encountered in the LAMMPS input
script, LAMMPS simply creates a class with the corresponding name,
invokes the "command" method of the class, and passes it the arguments
from the input script.  The command method can perform whatever
operations it wishes on LAMMPS data structures.

The single method your new class must define is as follows:

+---------+-----------------------------------------+
| command | operations performed by the new command |
+---------+-----------------------------------------+

Of course, the new class can define other methods and variables as
needed.


