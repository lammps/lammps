Compute styles
==============

Classes that compute scalar and vector quantities like temperature
and the pressure tensor, as well as classes that compute per-atom
quantities like kinetic energy and the centro-symmetry parameter
are derived from the Compute class.  New styles can be created
to add new calculations to LAMMPS.

Compute_temp.cpp is a simple example of computing a scalar
temperature.  Compute_ke_atom.cpp is a simple example of computing
per-atom kinetic energy.

Here is a brief description of methods you define in your new derived
class.  See compute.h for details.

+-----------------------+------------------------------------------------------------------+
| init                  | perform one time setup (required)                                |
+-----------------------+------------------------------------------------------------------+
| init_list             | neighbor list setup, if needed (optional)                        |
+-----------------------+------------------------------------------------------------------+
| compute_scalar        | compute a scalar quantity (optional)                             |
+-----------------------+------------------------------------------------------------------+
| compute_vector        | compute a vector of quantities (optional)                        |
+-----------------------+------------------------------------------------------------------+
| compute_peratom       | compute one or more quantities per atom (optional)               |
+-----------------------+------------------------------------------------------------------+
| compute_local         | compute one or more quantities per processor (optional)          |
+-----------------------+------------------------------------------------------------------+
| pack_comm             | pack a buffer with items to communicate (optional)               |
+-----------------------+------------------------------------------------------------------+
| unpack_comm           | unpack the buffer (optional)                                     |
+-----------------------+------------------------------------------------------------------+
| pack_reverse          | pack a buffer with items to reverse communicate (optional)       |
+-----------------------+------------------------------------------------------------------+
| unpack_reverse        | unpack the buffer (optional)                                     |
+-----------------------+------------------------------------------------------------------+
| remove_bias           | remove velocity bias from one atom (optional)                    |
+-----------------------+------------------------------------------------------------------+
| remove_bias_all       | remove velocity bias from all atoms in group (optional)          |
+-----------------------+------------------------------------------------------------------+
| restore_bias          | restore velocity bias for one atom after remove_bias (optional)  |
+-----------------------+------------------------------------------------------------------+
| restore_bias_all      | same as before, but for all atoms in group (optional)            |
+-----------------------+------------------------------------------------------------------+
| pair_tally_callback   | callback function for *tally*\ -style computes (optional).       |
+-----------------------+------------------------------------------------------------------+
| memory_usage          | tally memory usage (optional)                                    |
+-----------------------+------------------------------------------------------------------+

Tally-style computes are a special case, as their computation is done
in two stages: the callback function is registered with the pair style
and then called from the Pair::ev_tally() function, which is called for
each pair after force and energy has been computed for this pair. Then
the tallied values are retrieved with the standard compute_scalar or
compute_vector or compute_peratom methods. The USER-TALLY package
provides *examples*\ _compute_tally.html for utilizing this mechanism.
