Compute styles
==============

Classes that compute scalar and vector quantities like temperature
and the pressure tensor, as well as classes that compute per-atom
quantities like kinetic energy and the centro-symmetry parameter
are derived from the Compute class.  New styles can be created
to add new calculations to LAMMPS.

Compute\_temp.cpp is a simple example of computing a scalar
temperature.  Compute\_ke\_atom.cpp is a simple example of computing
per-atom kinetic energy.

Here is a brief description of methods you define in your new derived
class.  See compute.h for details.

+-----------------------+------------------------------------------------------------------+
| init                  | perform one time setup (required)                                |
+-----------------------+------------------------------------------------------------------+
| init\_list            | neighbor list setup, if needed (optional)                        |
+-----------------------+------------------------------------------------------------------+
| compute\_scalar       | compute a scalar quantity (optional)                             |
+-----------------------+------------------------------------------------------------------+
| compute\_vector       | compute a vector of quantities (optional)                        |
+-----------------------+------------------------------------------------------------------+
| compute\_peratom      | compute one or more quantities per atom (optional)               |
+-----------------------+------------------------------------------------------------------+
| compute\_local        | compute one or more quantities per processor (optional)          |
+-----------------------+------------------------------------------------------------------+
| pack\_comm            | pack a buffer with items to communicate (optional)               |
+-----------------------+------------------------------------------------------------------+
| unpack\_comm          | unpack the buffer (optional)                                     |
+-----------------------+------------------------------------------------------------------+
| pack\_reverse         | pack a buffer with items to reverse communicate (optional)       |
+-----------------------+------------------------------------------------------------------+
| unpack\_reverse       | unpack the buffer (optional)                                     |
+-----------------------+------------------------------------------------------------------+
| remove\_bias          | remove velocity bias from one atom (optional)                    |
+-----------------------+------------------------------------------------------------------+
| remove\_bias\_all     | remove velocity bias from all atoms in group (optional)          |
+-----------------------+------------------------------------------------------------------+
| restore\_bias         | restore velocity bias for one atom after remove\_bias (optional) |
+-----------------------+------------------------------------------------------------------+
| restore\_bias\_all    | same as before, but for all atoms in group (optional)            |
+-----------------------+------------------------------------------------------------------+
| pair\_tally\_callback | callback function for *tally*\ -style computes (optional).       |
+-----------------------+------------------------------------------------------------------+
| memory\_usage         | tally memory usage (optional)                                    |
+-----------------------+------------------------------------------------------------------+

Tally-style computes are a special case, as their computation is done
in two stages: the callback function is registered with the pair style
and then called from the Pair::ev\_tally() function, which is called for
each pair after force and energy has been computed for this pair. Then
the tallied values are retrieved with the standard compute\_scalar or
compute\_vector or compute\_peratom methods. The USER-TALLY package
provides *examples*\ \_compute\_tally.html for utilizing this mechanism.
