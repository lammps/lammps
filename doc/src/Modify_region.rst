Region styles
=============

Classes that define geometric regions are derived from the Region
class.  Regions are used elsewhere in LAMMPS to group atoms, delete
atoms to create a void, insert atoms in a specified region, etc.  New
styles can be created to add new region shapes to LAMMPS.

Region\_sphere.cpp is an example of a spherical region.

Here is a brief description of methods you define in your new derived
class.  See region.h for details.

+-------------------+---------------------------------------------------------------------+
| inside            | determine whether a point is in the region                          |
+-------------------+---------------------------------------------------------------------+
| surface\_interior | determine if a point is within a cutoff distance inside of surface  |
+-------------------+---------------------------------------------------------------------+
| surface\_exterior | determine if a point is within a cutoff distance outside of surface |
+-------------------+---------------------------------------------------------------------+
| shape\_update     | change region shape if set by time-dependent variable               |
+-------------------+---------------------------------------------------------------------+


