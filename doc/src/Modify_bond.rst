Bond, angle, dihedral, improper styles
======================================

Classes that compute molecular interactions are derived from the Bond,
Angle, Dihedral, and Improper classes.  New styles can be created to
add new potentials to LAMMPS.

Bond_harmonic.cpp is the simplest example of a bond style.  Ditto for
the harmonic forms of the angle, dihedral, and improper style
commands.

Here is a brief description of common methods you define in your
new derived class.  See bond.h, angle.h, dihedral.h, and improper.h
for details and specific additional methods.

+-----------------------+---------------------------------------------------------------------+
| Required              | "pure" methods that *must* be overridden in a derived class         |
+=======================+=====================================================================+
| compute               | compute the molecular interactions for all listed items             |
+-----------------------+---------------------------------------------------------------------+
| coeff                 | set coefficients for one type                                       |
+-----------------------+---------------------------------------------------------------------+
| equilibrium_distance  | length of bond, used by SHAKE (bond styles only)                    |
+-----------------------+---------------------------------------------------------------------+
| equilibrium_angle     | opening of angle, used by SHAKE (angle styles only)                 |
+-----------------------+---------------------------------------------------------------------+
| write & read_restart  | writes/reads coeffs to restart files                                |
+-----------------------+---------------------------------------------------------------------+
| single                | force/r (bond styles only) and energy of a single bond or angle     |
+-----------------------+---------------------------------------------------------------------+


+--------------------------------+----------------------------------------------------------------------+
| Optional                       | methods that have a default or dummy implementation                  |
+================================+======================================================================+
| init                           | check if all coefficients are set, calls init_style()                |
+--------------------------------+----------------------------------------------------------------------+
| init_style                     | check if style specific conditions are met                           |
+--------------------------------+----------------------------------------------------------------------+
| settings                       | apply global settings for all types                                  |
+--------------------------------+----------------------------------------------------------------------+
| write & read_restart_settings  | writes/reads global style settings to restart files                  |
+--------------------------------+----------------------------------------------------------------------+
| write_data                     | write corresponding Coeffs section(s) in data file                   |
+--------------------------------+----------------------------------------------------------------------+
| memory_usage                   | tally memory allocated by the style                                  |
+--------------------------------+----------------------------------------------------------------------+
| extract                        | provide access to internal data  (bond or angle styles only)         |
+--------------------------------+----------------------------------------------------------------------+
| reinit                         | reset all type-based parameters, called by fix adapt (bonds only)    |
+--------------------------------+----------------------------------------------------------------------+
| pack & unpack_forward_comm     | copy data to and from buffer in forward communication (bonds only)   |
+--------------------------------+----------------------------------------------------------------------+
| pack & unpack_reverse_comm     | copy data to and from buffer in reverse communication (bonds only)   |
+--------------------------------+----------------------------------------------------------------------+

Here is a list of flags or settings that should be set in the
constructor of the derived class when they differ from the default
setting.

+---------------------------------+------------------------------------------------------------------------------+---------+
| Name of flag                    | Description                                                                  | default |
+=================================+==============================================================================+=========+
| writedata                       | 1 if write_data() is implemented                                             | 1       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| single_extra                    | number of extra single values calculated (bond styles only)                  | 0       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| partial_flag                    | 1 if bond type can be set to 0 and deleted (bond styles only)                | 0       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| reinitflag                      | 1 if style has reinit() and is compatible with fix adapt                     | 1       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| comm_forward                    | size of buffer (in doubles) for forward communication  (bond styles only)    | 0       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| comm_reverse                    | size of buffer (in doubles) for reverse communication  (bond styles only)    | 0       |
+---------------------------------+------------------------------------------------------------------------------+---------+
| comm_reverse_off                | size of buffer for reverse communication with newton off (bond styles only)  | 0       |
+---------------------------------+------------------------------------------------------------------------------+---------+
