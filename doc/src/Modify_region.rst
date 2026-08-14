Region styles
=============

Classes that define geometric regions are derived from the Region
class.  Regions are used elsewhere in LAMMPS to group atoms, delete
atoms to create a void, insert atoms in a specified region, etc.  New
styles can be created to add new region shapes to LAMMPS.

The files ``region_sphere.cpp`` and ``region_block.cpp`` and their
corresponding headers are representative examples of region style
implementations.

Here is a brief description of methods you define in your new derived
class.  See ``src/region.h`` for details.

+---------------------------+---------------------------------------------------------------------+
| Required                  | "pure" methods that *must* be overridden in a derived class         |
+===========================+=====================================================================+
| inside                    | determine whether a point is in the region                          |
+---------------------------+---------------------------------------------------------------------+
| surface_interior          | determine if a point is within a cutoff distance inside of surface  |
+---------------------------+---------------------------------------------------------------------+
| surface_exterior          | determine if a point is within a cutoff distance outside of surface |
+---------------------------+---------------------------------------------------------------------+

+---------------------------+---------------------------------------------------------------------+
| Optional                  | methods that have a default or empty implementation                 |
+===========================+=====================================================================+
| shape_update              | change region shape if set by time-dependent variable               |
+---------------------------+---------------------------------------------------------------------+
| bbox_update               | update the bounding box when the shape is variable                  |
+---------------------------+---------------------------------------------------------------------+
| init                      | check style-specific conditions before a run                        |
+---------------------------+---------------------------------------------------------------------+
| set_velocity              | compute the velocity of a moving region surface                     |
+---------------------------+---------------------------------------------------------------------+
| set_velocity_shape        | shape-specific part of set_velocity                                 |
+---------------------------+---------------------------------------------------------------------+
| velocity_contact_shape    | return velocity at a contact point for the specific shape           |
+---------------------------+---------------------------------------------------------------------+
| write_restart             | write region state to a restart file                                |
+---------------------------+---------------------------------------------------------------------+
| restart                   | read region state from a restart file                               |
+---------------------------+---------------------------------------------------------------------+
| length_restart_string     | return the size of the restart data                                 |
+---------------------------+---------------------------------------------------------------------+
| reset_vel                 | reset velocity when the region is reset                             |
+---------------------------+---------------------------------------------------------------------+
