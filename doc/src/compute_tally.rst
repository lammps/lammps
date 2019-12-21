.. index:: compute force/tally

compute force/tally command
===========================

compute heat/flux/tally command
===============================

compute pe/tally command
========================

compute pe/mol/tally command
============================

compute stress/tally command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID style group2-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *force/tally* or *pe/tally* or *pe/mol/tally* or *stress/tally*
* group2-ID = group ID of second (or same) group

Examples
""""""""


.. parsed-literal::

   compute 1 lower force/tally upper
   compute 1 left pe/tally right
   compute 1 lower stress/tally lower

Description
"""""""""""

Define a computation that calculates properties between two groups of
atoms by accumulating them from pairwise non-bonded computations.  The
two groups can be the same. This is similar to :doc:`compute group/group <compute_group_group>` only that the data is
accumulated directly during the non-bonded force computation. The
computes *force/tally*\ , *pe/tally*\ , *stress/tally*\ , and
*heat/flux/tally* are primarily provided as example how to program
additional, more sophisticated computes using the tally callback
mechanism. Compute *pe/mol/tally* is one such style, that can
- through using this mechanism - separately tally intermolecular
and intramolecular energies. Something that would otherwise be
impossible without integrating this as a core functionality into
the based classes of LAMMPS.


----------


The pairwise contributions are computing via a callback that the
compute registers with the non-bonded pairwise force computation.
This limits the use to systems that have no bonds, no Kspace, and no
many-body interactions. On the other hand, the computation does not
have to compute forces or energies a second time and thus can be much
more efficient. The callback mechanism allows to write more complex
pairwise property computations.


----------


**Output info:**

Compute *pe/tally* calculates a global scalar (the energy) and a per
atom scalar (the contributions of the single atom to the global
scalar). Compute *pe/mol/tally* calculates a global 4-element vector
containing (in this order): *evdwl* and *ecoul* for intramolecular pairs
and *evdwl* and *ecoul* for intermolecular pairs. Since molecules are
identified by their molecule IDs, the partitioning does not have to be
related to molecules, but the energies are tallied into the respective
slots depending on whether the molecule IDs of a pair are the same or
different. Compute *force/tally* calculates a global scalar (the force
magnitude) and a per atom 3-element vector (force contribution from
each atom).  Compute *stress/tally* calculates a global scalar
(average of the diagonal elements of the stress tensor) and a per atom
vector (the 6 elements of stress tensor contributions from the
individual atom).

Both the scalar and vector values calculated by this compute are
"extensive".

Restrictions
""""""""""""


This compute is part of the USER-TALLY package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Not all pair styles can be evaluated in a pairwise mode as required by
this compute.  For example, 3-body and other many-body potentials,
such as :doc:`Tersoff <pair_tersoff>` and
:doc:`Stillinger-Weber <pair_sw>` cannot be used.  :doc:`EAM <pair_eam>`
potentials only include the pair potential portion of the EAM
interaction when used by this compute, not the embedding term.  Also
bonded or Kspace interactions do not contribute to this compute.

The computes in this package are not compatible with dynamic groups.

Related commands
""""""""""""""""

*compute group/group*\ \_compute\_group\_group.html, *compute
heat/flux*\ \_compute\_heat\_flux.html

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
