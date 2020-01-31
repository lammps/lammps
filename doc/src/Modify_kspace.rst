Kspace styles
=============

Classes that compute long-range Coulombic interactions via K-space
representations (Ewald, PPPM) are derived from the KSpace class.  New
styles can be created to add new K-space options to LAMMPS.

Ewald.cpp is an example of computing K-space interactions.

Here is a brief description of methods you define in your new derived
class.  See kspace.h for details.

+---------------+----------------------------------------------+
| init          | initialize the calculation before a run      |
+---------------+----------------------------------------------+
| setup         | computation before the 1st timestep of a run |
+---------------+----------------------------------------------+
| compute       | every-timestep computation                   |
+---------------+----------------------------------------------+
| memory\_usage | tally of memory usage                        |
+---------------+----------------------------------------------+


