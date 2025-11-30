Atom Class
**********

The Atom class provides access to atom style related global settings and
per-atom data that is stored with atoms and migrates with them from
sub-domain to sub-domain as atoms move around.  This includes topology
data, which is stored with either one specific atom or all atoms
involved depending on the settings of the :doc:`newton command
<newton>`.  It also contains and manages the list of molecule templates
that are usually created by the :doc:`molecule command <molecule>`.

The actual per-atom data is allocated and managed by one of the various
classes derived from the AtomVec class as determined by the
:doc:`atom_style command <atom_style>`.  The pointers in the Atom class
are updated by the AtomVec class as needed.

.. doxygenclass:: LAMMPS_NS::Atom
   :project: progguide
   :members:
