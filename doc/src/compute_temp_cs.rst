.. index:: compute temp/cs

compute temp/cs command
=======================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID temp/cs group1 group2

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/cs = style name of this compute command
* group1 = group-ID of either cores or shells
* group2 = group-ID of either shells or cores

Examples
""""""""


.. parsed-literal::

   compute oxygen_c-s all temp/cs O_core O_shell
   compute core_shells all temp/cs cores shells

Description
"""""""""""

Define a computation that calculates the temperature of a system based
on the center-of-mass velocity of atom pairs that are bonded to each
other.  This compute is designed to be used with the adiabatic
core/shell model of :ref:`(Mitchell and Finchham) <MitchellFinchham1>`.  See
the :doc:`Howto coreshell <Howto_coreshell>` doc page for an overview of
the model as implemented in LAMMPS.  Specifically, this compute
enables correct temperature calculation and thermostatting of
core/shell pairs where it is desirable for the internal degrees of
freedom of the core/shell pairs to not be influenced by a thermostat.
A compute of this style can be used by any command that computes a
temperature via :doc:`fix_modify <fix_modify>` e.g. :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix npt <fix_nh>`, etc.

Note that this compute does not require all ions to be polarized,
hence defined as core/shell pairs.  One can mix core/shell pairs and
ions without a satellite particle if desired. The compute will
consider the non-polarized ions according to the physical system.

For this compute, core and shell particles are specified by two
respective group IDs, which can be defined using the
:doc:`group <group>` command.  The number of atoms in the two groups
must be the same and there should be one bond defined between a pair
of atoms in the two groups.  Non-polarized ions which might also be
included in the treated system should not be included into either of
these groups, they are taken into account by the *group-ID* (2nd
argument) of the compute.

The temperature is calculated by the formula KE = dim/2 N k T, where
KE = total kinetic energy of the group of atoms (sum of 1/2 m v\^2),
dim = 2 or 3 = dimensionality of the simulation, N = number of atoms
in the group, k = Boltzmann constant, and T = temperature.  Note that
the velocity of each core or shell atom used in the KE calculation is
the velocity of the center-of-mass (COM) of the core/shell pair the
atom is part of.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that v\^2 is replaced by vx\*vy for the xy
component, etc.  The 6 components of the vector are ordered xx, yy,
zz, xy, xz, yz.  In contrast to the temperature, the velocity of
each core or shell atom is taken individually.

The change this fix makes to core/shell atom velocities is essentially
computing the temperature after a "bias" has been removed from the
velocity of the atoms.  This "bias" is the velocity of the atom
relative to the COM velocity of the core/shell pair.  If this compute
is used with a fix command that performs thermostatting then this bias
will be subtracted from each atom, thermostatting of the remaining COM
velocity will be performed, and the bias will be added back in.  This
means the thermostatting will effectively be performed on the
core/shell pairs, instead of on the individual core and shell atoms.
Thermostatting fixes that work in this way include :doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix temp/berendsen <fix_temp_berendsen>`, and :doc:`fix langevin <fix_langevin>`.

The internal energy of core/shell pairs can be calculated by the
:doc:`compute temp/chunk <compute_temp_chunk>` command, if chunks are
defined as core/shell pairs.  See the :doc:`Howto coreshell <Howto_coreshell>` doc page doc page for more discussion
on how to do this.

**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1-6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""


The number of core/shell pairs contributing to the temperature is
assumed to be constant for the duration of the run.  No fixes should
be used which generate new molecules or atoms during a simulation.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute temp/chunk <compute_temp_chunk>`

**Default:** none


----------


.. _MitchellFinchham1:



**(Mitchell and Finchham)** Mitchell, Finchham, J Phys Condensed Matter,
5, 1031-1038 (1993).
