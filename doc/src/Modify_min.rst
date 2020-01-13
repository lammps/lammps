Minimization styles
===================

Classes that perform energy minimization derived from the Min class.
New styles can be created to add new minimization algorithms to
LAMMPS.

Min\_cg.cpp is an example of conjugate gradient minimization.

Here is a brief description of methods you define in your new derived
class.  See min.h for details.

+---------------+------------------------------------------+
| init          | initialize the minimization before a run |
+---------------+------------------------------------------+
| run           | perform the minimization                 |
+---------------+------------------------------------------+
| memory\_usage | tally of memory usage                    |
+---------------+------------------------------------------+


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
