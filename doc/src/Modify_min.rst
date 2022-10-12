Minimization styles
===================

Classes that perform energy minimization derived from the Min class.
New styles can be created to add new minimization algorithms to
LAMMPS.

Min_cg.cpp is an example of conjugate gradient minimization.

Here is a brief description of methods you define in your new derived
class.  See min.h for details.

+---------------+------------------------------------------+
| init          | initialize the minimization before a run |
+---------------+------------------------------------------+
| run           | perform the minimization                 |
+---------------+------------------------------------------+
| memory_usage  | tally of memory usage                    |
+---------------+------------------------------------------+
