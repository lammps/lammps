.. index:: compute cluster/atom
.. index:: compute fragment/atom
.. index:: compute aggregate/atom

compute cluster/atom command
============================

compute fragment/atom command
=============================

compute aggregate/atom command
==============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID cluster/atom cutoff
   compute ID group-ID fragment/atom keyword value ...
   compute ID group-ID aggregate/atom cutoff

* ID, group-ID are documented in :doc:`compute <compute>` command
* *cluster/atom* or *fragment/atom* or *aggregate/atom* = style name of this compute command
* cutoff = distance within which to label atoms as part of same cluster (distance units)
* zero or more keyword/value pairs may be appended to *fragment/atom*
* keyword = *single*

    .. parsed-literal::

       *single* value = *yes* or *no* to treat single atoms (no bonds) as fragments

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all cluster/atom 3.5
   compute 1 all fragment/atom
   compute 1 all fragment/atom single no
   compute 1 all aggregate/atom 3.5

Description
"""""""""""

Define a computation that assigns each atom a cluster, fragment, or
aggregate ID.  Only atoms in the compute group are clustered and
assigned cluster IDs. Atoms not in the compute group are assigned an
ID = 0.

A cluster is defined as a set of atoms, each of which is within the
cutoff distance from one or more other atoms in the cluster.  If an
atom has no neighbors within the cutoff distance, then it is a 1-atom
cluster.

A fragment is similarly defined as a set of atoms, each of which has a
bond to another atom in the fragment.  Bonds can be defined initially
via the :doc:`data file <read_data>` or :doc:`create_bonds
<create_bonds>` commands, or dynamically by fixes which create or
break bonds like :doc:`fix bond/react <fix_bond_react>`, :doc:`fix
bond/create <fix_bond_create>`, :doc:`fix bond/swap <fix_bond_swap>`,
or :doc:`fix bond/break <fix_bond_break>`.  The cluster ID or fragment
ID of every atom in the cluster will be set to the smallest atom ID of
any atom in the cluster or fragment, respectively.

For the *fragment/atom* style, the *single* keyword determines whether
single atoms (not bonded to another atom) are treated as one-atom
fragments or not, based on the *yes* or *no* setting.  If the setting
is *no* (the default), their fragment IDs are set to 0.

An aggregate is defined by combining the rules for clusters and
fragments, i.e. a set of atoms, where each of it is within the cutoff
distance from one or more atoms within a fragment that is part of
the same cluster. This measure can be used to track molecular assemblies
like micelles.

For computes *cluster/atom* and *aggregate/atom* a neighbor list
needed to compute cluster IDs is constructed each time the compute is
invoked.  Thus it can be inefficient to compute/dump this quantity too
frequently or to have multiple *cluster/atom* or *aggregate/atom*
style computes.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included when computing the clusters. This
   does not apply when using long-range coulomb (\ *coul/long*, *coul/msm*,
   *coul/wolf* or similar.  One way to get around this would be to set
   special_bond scaling factors to very tiny numbers that are not exactly
   zero (e.g. 1.0e-50). Another workaround is to write a dump file, and
   use the :doc:`rerun <rerun>` command to compute the clusters for
   snapshots in the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

.. note::

   For the compute fragment/atom style, each fragment is identified
   using the current bond topology.  This will not account for bonds
   broken by the :doc:`bond_style quartic <bond_quartic>` command
   because it does not perform a full update of the bond topology data
   structures within LAMMPS.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values will be an ID > 0, as explained above.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`

Default
"""""""

The default for fragment/atom is single no.

