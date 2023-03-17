.. index:: comm_modify

comm_modify command
===================

Syntax
""""""

.. code-block:: LAMMPS

   comm_modify keyword value ...

* one or more keyword/value pairs may be appended
* keyword = *mode* or *cutoff* or *cutoff/multi* or *group* or *reduce/multi* or *vel*

  .. parsed-literal::

       *mode* value = *single*, *multi*, or *multi/old* = communicate atoms within a single or multiple distances
       *cutoff* value = Rcut (distance units) = communicate atoms from this far away
       *cutoff/multi* collection value
          collection = atom collection or collection range (supports asterisk notation)
          value = Rcut (distance units) = communicate atoms for selected types from this far away
       *reduce/multi* arg = none = reduce number of communicated ghost atoms for multi style
       *cutoff/multi/old* type value
          type = atom type or type range (supports asterisk notation)
          value = Rcut (distance units) = communicate atoms for selected types from this far away
       *group* value = group-ID = only communicate atoms in the group
       *vel* value = *yes* or *no* = do or do not communicate velocity info with ghost atoms

Examples
""""""""

.. code-block:: LAMMPS

   comm_modify mode multi reduce/multi
   comm_modify mode multi group solvent
   comm_modify mode multi cutoff/multi 1 10.0 cutoff/multi 2*4 15.0
   comm_modify vel yes
   comm_modify mode single cutoff 5.0 vel yes
   comm_modify cutoff/multi * 0.0

Description
"""""""""""

This command sets parameters that affect the inter-processor
communication of atom information that occurs each timestep as
coordinates and other properties are exchanged between neighboring
processors and stored as properties of ghost atoms.

.. note::

   These options apply to the currently defined comm style.  When
   you specify a :doc:`comm_style <comm_style>` or
   :doc:`read_restart <read_restart>` command, all communication settings
   are restored to their default or stored values, including those
   previously reset by a comm_modify command.  Thus if your input script
   specifies a comm_style or read_restart command, you should use the
   comm_modify command after it.

The *mode* keyword determines whether a single or multiple cutoff
distances are used to determine which atoms to communicate.

The default mode is *single* which means each processor acquires
information for ghost atoms that are within a single distance from its
subdomain.  The distance is by default the maximum of the neighbor
cutoff across all atom type pairs.

For many systems this is an efficient algorithm, but for systems with
widely varying cutoffs for different type pairs, the *multi* or *multi/old* mode can
be faster.  In *multi*, each atom is assigned to a collection which should
correspond to a set of atoms with similar interaction cutoffs.
In this case, each atom collection is assigned its own distance
cutoff for communication purposes, and fewer atoms will be
communicated. in *multi/old*, a similar technique is used but atoms
are grouped by atom type. See the :doc:`neighbor multi <neighbor>`  and
:doc:`neighbor multi/old <neighbor>` commands for
neighbor list construction options that may also be beneficial for
simulations of this kind. The *multi* communication mode is only compatible
with the *multi* neighbor style. The *multi/old* communication mode is comparable
with both the *multi* and *multi/old* neighbor styles.

The *cutoff* keyword allows you to extend the ghost cutoff distance
for communication mode *single*, which is the distance from the borders
of a processor's subdomain at which ghost atoms are acquired from other
processors.  By default the ghost cutoff = neighbor cutoff = pairwise
force cutoff + neighbor skin.  See the :doc:`neighbor <neighbor>` command
for more information about the skin distance.  If the specified Rcut is
greater than the neighbor cutoff, then extra ghost atoms will be acquired.
If the provided cutoff is smaller, the provided value will be ignored,
the ghost cutoff is set to the neighbor cutoff and a warning will be
printed. Specifying a cutoff value of 0.0 will reset any previous value
to the default. If bonded interactions exist and equilibrium bond length
information is available, then also a heuristic based on that bond length
is computed. It is used as communication cutoff, if there is no pair
style present and no *comm_modify cutoff* command used. Otherwise a
warning is printed, if this bond based estimate is larger than the
communication cutoff used.

The *cutoff/multi* option is equivalent to *cutoff*, but applies to
communication mode *multi* instead. Since the communication cutoffs are
determined per atom collections, a collection specifier is needed and
cutoff for one or multiple collections can be extended. Also ranges of
collections using the usual asterisk notation can be given.  Collections
are indexed from 1 to N where N is the total number of collections.
Note that the arguments for *cutoff/multi* are parsed right before each
simulation to account for potential changes in the number of
collections.  Custom cutoffs are preserved between runs but if
collections are redefined, one may want to re-specify the communication
cutoffs.  For granular pair styles,the default cutoff is set to the sum
of the current maximum atomic radii for each collection.  The
*cutoff/multi/old* option is similar to *cutoff/multi* except it
operates on atom types as opposed to collections.

The *reduce/multi* option applies to *multi* and sets the communication
cutoff for a particle equal to the maximum interaction distance between particles
in the same collection. This reduces the number of
ghost atoms that need to be communicated. This method is only compatible with the
*multi* neighbor style and requires a half neighbor list and Newton on.
See the :doc:`neighbor multi <neighbor>` command for more information.

These are simulation scenarios in which it may be useful or even
necessary to set a ghost cutoff > neighbor cutoff:

* a single polymer chain with bond interactions, but no pairwise interactions
* bonded interactions (e.g. dihedrals) extend further than the pairwise cutoff
* ghost atoms beyond the pairwise cutoff are needed for some computation

In the first scenario, a pairwise potential is not defined.  Thus the
pairwise neighbor cutoff will be 0.0.  But ghost atoms are still
needed for computing bond, angle, etc interactions between atoms on
different processors, or when the interaction straddles a periodic
boundary.

The appropriate ghost cutoff depends on the :doc:`newton bond <newton>`
setting.  For newton bond *off*, the distance needs to be the furthest
distance between any two atoms in the bond, angle, etc.  E.g. the
distance between 1-4 atoms in a dihedral.  For newton bond *on*, the
distance between the central atom in the bond, angle, etc and any
other atom is sufficient.  E.g. the distance between 2-4 atoms in a
dihedral.

In the second scenario, a pairwise potential is defined, but its
neighbor cutoff is not sufficiently long enough to enable bond, angle,
etc terms to be computed.  As in the previous scenario, an appropriate
ghost cutoff should be set.

In the last scenario, a :doc:`fix <fix>` or :doc:`compute <compute>` or
:doc:`pairwise potential <pair_style>` needs to calculate with ghost
atoms beyond the normal pairwise cutoff for some computation it
performs (e.g. locate neighbors of ghost atoms in a manybody pair
potential).  Setting the ghost cutoff appropriately can ensure it will
find the needed atoms.

.. note::

   In these scenarios, if you do not set the ghost cutoff long
   enough, and if there is only one processor in a periodic dimension
   (e.g. you are running in serial), then LAMMPS may "find" the atom it
   is looking for (e.g. the partner atom in a bond), that is on the far
   side of the simulation box, across a periodic boundary.  This will
   typically lead to bad dynamics (i.e. the bond length is now the
   simulation box length).  To detect if this is happening, see the
   :doc:`neigh_modify cluster <neigh_modify>` command.

The *group* keyword will limit communication to atoms in the specified
group.  This can be useful for models where no ghost atoms are needed
for some kinds of particles.  All atoms (not just those in the
specified group) will still migrate to new processors as they move.
The group specified with this option must also be specified via the
:doc:`atom_modify first <atom_modify>` command.

The *vel* keyword enables velocity information to be communicated with
ghost particles.  Depending on the :doc:`atom_style <atom_style>`,
velocity info includes the translational velocity, angular velocity,
and angular momentum of a particle.  If the *vel* option is set to
*yes*, then ghost atoms store these quantities; if *no* then they do
not.  The *yes* setting is needed by some pair styles which require
the velocity state of both the I and J particles to compute a pairwise
I,J interaction, as well as by some compute and fix commands.

Note that if the :doc:`fix deform <fix_deform>` command is being used
with its "remap v" option enabled, then the velocities for ghost atoms
(in the fix deform group) mirrored across a periodic boundary will
also include components due to any velocity shift that occurs across
that boundary (e.g. due to dilation or shear).

Restrictions
""""""""""""

Communication mode *multi* is currently only available for
:doc:`comm_style <comm_style>` *brick*\ .

Related commands
""""""""""""""""

:doc:`comm_style <comm_style>`, :doc:`neighbor <neighbor>`

Default
"""""""

The option defaults are mode = single, group = all, cutoff = 0.0, vel =
no.  The cutoff default of 0.0 means that ghost cutoff = neighbor
cutoff = pairwise force cutoff + neighbor skin.
