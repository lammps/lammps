.. index:: compute bpm/peri/damage/atom

compute bpm/peri/damage/atom command
====================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID bpm/peri/damage/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* bpm/peri/damage/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all bpm/peri/damage/atom

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the volume-weighted bond damage of
each peridynamic node for the :doc:`bond_style bpm/peri <bond_bpm_peri>`
model.  The damage of node *i* is

.. math::

   d_i = 1 - \frac{\sum_{j\,\mathrm{(intact)}} V_j}{\sum_{j\,\mathrm{(initial)}} V_j}

where the numerator sums the nodal volumes :math:`V_j` of the bonds that
are still intact and the denominator is the reference interaction volume
(the same sum over all bonds present in the initial configuration).
Damage is 0 for a fully intact node and approaches 1 as a node loses
bonds, so it is a convenient measure for visualizing cracks and fracture
surfaces.  It is the BPM-framework equivalent of the legacy :doc:`compute
damage/atom <compute_damage_atom>`.

The reference interaction volume is recorded by :doc:`bond_style bpm/peri
<bond_bpm_peri>` when the bonds are first set up, so that compute style
must be defined for this compute to be used.  Damage is 0 for nodes not
in the specified compute group.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input.  See the
:doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS
output options.  The per-atom values are unitless, between 0 and 1.

Restrictions
""""""""""""

This compute is part of the BPM package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This compute requires :doc:`bond_style bpm/peri <bond_bpm_peri>`, which
supplies the reference interaction volume, and a per-atom *vfrac*
property.

Related commands
""""""""""""""""

:doc:`bond_style bpm/peri <bond_bpm_peri>`, :doc:`compute damage/atom
<compute_damage_atom>`, :doc:`compute nbond/atom <compute_nbond_atom>`

Default
"""""""

none
