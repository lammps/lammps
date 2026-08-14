.. index:: pair_style template/offset

pair_style template/offset command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style template/offset

   pair_coeff * * mol-ID offset

* mol-ID = ID of a molecule template (from the :doc:`molecule <molecule>` command)
* offset = energy added to the system for each distinct match of that template (energy units)

Examples
""""""""

.. code-block:: LAMMPS

   molecule h3o h3o.mol
   pair_style hybrid/overlay lj/cut/coul/long 10.0 template/offset
   pair_coeff * * lj/cut/coul/long ...
   pair_coeff * * template/offset h3o 8.8115

Description
"""""""""""

Pair style *template/offset* adds a constant energy *offset* to the system for
every distinct match of a molecule template.  It is intended to give a chemical
species a fixed diabatic reference energy: the internal energy of a molecule is
arbitrary in classical mechanics, so shifting each species by a constant lets
reacting species be placed on potential energy surfaces of comparable magnitude.

For each registered template the style counts the number of distinct, non
overlapping matches of the template in the current bonded topology and adds
``offset * count`` to the potential energy.  Matching uses a graph-isomorphism
kernel on the 1-2 bond topology (types and connectivity must both match); the
kernel is shared with :doc:`fix msevb <fix_msevb>`.  Only connected
(single-molecule) templates are supported.  A one-atom template counts atoms of
its type.

Because the offset is a per-configuration constant, this style contributes
energy only. The forces and the per-atom virial it produces are zero.  It
is therefore normally used through :doc:`pair_style hybrid/overlay
<pair_hybrid>` alongside the real force field, and one ``pair_coeff * *`` line is
given per template (templates accumulate).  The style applies to all atom type
pairs.

This pair style is independent of :doc:`fix msevb <fix_msevb>`, but pairs
naturally with it: because fix msevb evaluates every EVB diabatic state through
the normal pair pipeline, adding this style automatically includes the
appropriate species offset in each state's energy, with no bookkeeping inside
the fix.

Restrictions
""""""""""""

This pair style is part of the MSEVB package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>` page.

It requires a molecular :doc:`atom style <atom_style>` (bond topology is needed
for matching) and at least one molecule template defined with the
:doc:`molecule <molecule>` command.

Related commands
""""""""""""""""

:doc:`fix msevb <fix_msevb>`, :doc:`molecule <molecule>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`

Default
"""""""

none
