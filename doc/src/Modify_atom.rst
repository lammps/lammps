Atom styles
===========

Classes that define an :doc:`atom style <atom_style>` are derived from
the AtomVec class and managed by the Atom class.  The atom style
determines what attributes are associated with an atom.  A new atom
style can be created if one of the existing atom styles does not
define all the attributes you need to store and communicate with
atoms.

Atom_vec_atomic.cpp is a simple example of an atom style.

Here is a brief description of methods you define in your new derived
class.  See atom_vec.h for details.

+-------------------------+--------------------------------------------------------------------------------+
| init                    | one time setup (optional)                                                      |
+-------------------------+--------------------------------------------------------------------------------+
| grow                    | re-allocate atom arrays to longer lengths (required)                           |
+-------------------------+--------------------------------------------------------------------------------+
| grow_reset              | make array pointers in Atom and AtomVec classes consistent (required)          |
+-------------------------+--------------------------------------------------------------------------------+
| copy                    | copy info for one atom to another atom's array locations (required)            |
+-------------------------+--------------------------------------------------------------------------------+
| pack_comm               | store an atom's info in a buffer communicated every timestep (required)        |
+-------------------------+--------------------------------------------------------------------------------+
| pack_comm_vel           | add velocity info to communication buffer (required)                           |
+-------------------------+--------------------------------------------------------------------------------+
| pack_comm_hybrid        | store extra info unique to this atom style (optional)                          |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_comm             | retrieve an atom's info from the buffer (required)                             |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_comm_vel         | also retrieve velocity info (required)                                         |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_comm_hybrid      | retrieve extra info unique to this atom style (optional)                       |
+-------------------------+--------------------------------------------------------------------------------+
| pack_reverse            | store an atom's info in a buffer communicating partial forces  (required)      |
+-------------------------+--------------------------------------------------------------------------------+
| pack_reverse_hybrid     | store extra info unique to this atom style (optional)                          |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_reverse          | retrieve an atom's info from the buffer (required)                             |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_reverse_hybrid   | retrieve extra info unique to this atom style (optional)                       |
+-------------------------+--------------------------------------------------------------------------------+
| pack_border             | store an atom's info in a buffer communicated on neighbor re-builds (required) |
+-------------------------+--------------------------------------------------------------------------------+
| pack_border_vel         | add velocity info to buffer (required)                                         |
+-------------------------+--------------------------------------------------------------------------------+
| pack_border_hybrid      | store extra info unique to this atom style (optional)                          |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_border           | retrieve an atom's info from the buffer (required)                             |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_border_vel       | also retrieve velocity info (required)                                         |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_border_hybrid    | retrieve extra info unique to this atom style (optional)                       |
+-------------------------+--------------------------------------------------------------------------------+
| pack_exchange           | store all an atom's info to migrate to another processor (required)            |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_exchange         | retrieve an atom's info from the buffer (required)                             |
+-------------------------+--------------------------------------------------------------------------------+
| size_restart            | number of restart quantities associated with proc's atoms (required)           |
+-------------------------+--------------------------------------------------------------------------------+
| pack_restart            | pack atom quantities into a buffer (required)                                  |
+-------------------------+--------------------------------------------------------------------------------+
| unpack_restart          | unpack atom quantities from a buffer (required)                                |
+-------------------------+--------------------------------------------------------------------------------+
| create_atom             | create an individual atom of this style (required)                             |
+-------------------------+--------------------------------------------------------------------------------+
| data_atom               | parse an atom line from the data file (required)                               |
+-------------------------+--------------------------------------------------------------------------------+
| data_atom_hybrid        | parse additional atom info unique to this atom style (optional)                |
+-------------------------+--------------------------------------------------------------------------------+
| data_vel                | parse one line of velocity information from data file (optional)               |
+-------------------------+--------------------------------------------------------------------------------+
| data_vel_hybrid         | parse additional velocity data unique to this atom style (optional)            |
+-------------------------+--------------------------------------------------------------------------------+
| memory_usage            | tally memory allocated by atom arrays (required)                               |
+-------------------------+--------------------------------------------------------------------------------+

The constructor of the derived class sets values for several variables
that you must set when defining a new atom style, which are documented
in atom_vec.h.  New atom arrays are defined in atom.cpp.  Search for
the word "customize" and you will find locations you will need to
modify.

.. note::

   It is possible to add some attributes, such as a molecule ID, to
   atom styles that do not have them via the :doc:`fix property/atom <fix_property_atom>` command.  This command also
   allows new custom attributes consisting of extra integer or
   floating-point values to be added to atoms.  See the :doc:`fix property/atom <fix_property_atom>` doc page for examples of cases
   where this is useful and details on how to initialize, access, and
   output the custom values.

New :doc:`pair styles <pair_style>`, :doc:`fixes <fix>`, or
:doc:`computes <compute>` can be added to LAMMPS, as discussed below.
The code for these classes can use the per-atom properties defined by
fix property/atom.  The Atom class has a find_custom() method that is
useful in this context:

.. code-block:: c++

   int index = atom->find_custom(char *name, int &flag);

The "name" of a custom attribute, as specified in the :doc:`fix property/atom <fix_property_atom>` command, is checked to verify
that it exists and its index is returned.  The method also sets flag =
0/1 depending on whether it is an integer or floating-point attribute.
The vector of values associated with the attribute can then be
accessed using the returned index as

.. code-block:: c++

   int *ivector = atom->ivector[index];
   double *dvector = atom->dvector[index];

Ivector or dvector are vectors of length Nlocal = # of owned atoms,
which store the attributes of individual atoms.
