.. index:: special\_bonds

special\_bonds command
======================

Syntax
""""""


.. parsed-literal::

   special_bonds keyword values ...

* one or more keyword/value pairs may be appended
* keyword = *amber* or *charmm* or *dreiding* or *fene* or *lj/coul* or *lj* or *coul* or *angle* or *dihedral*
  
  .. parsed-literal::
  
       *amber* values = none
       *charmm* values = none
       *dreiding* values = none
       *fene* values = none
       *lj/coul* values = w1,w2,w3
         w1,w2,w3 = weights (0.0 to 1.0) on pairwise Lennard-Jones and Coulombic interactions
       *lj* values = w1,w2,w3
         w1,w2,w3 = weights (0.0 to 1.0) on pairwise Lennard-Jones interactions
       *coul* values = w1,w2,w3
         w1,w2,w3 = weights (0.0 to 1.0) on pairwise Coulombic interactions
       *angle* value = *yes* or *no*
       *dihedral* value = *yes* or *no*



Examples:


.. parsed-literal::

   special_bonds amber
   special_bonds charmm
   special_bonds fene dihedral no
   special_bonds lj/coul 0.0 0.0 0.5 angle yes dihedral yes
   special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.0 dihedral yes

Description
"""""""""""

Set weighting coefficients for pairwise energy and force contributions
between pairs of atoms that are also permanently bonded to each other,
either directly or via one or two intermediate bonds.  These weighting
factors are used by nearly all :doc:`pair styles <pair_style>` in LAMMPS
that compute simple pairwise interactions.  Permanent bonds between
atoms are specified by defining the bond topology in the data file
read by the :doc:`read_data <read_data>` command.  Typically a
:doc:`bond_style <bond_style>` command is also used to define a bond
potential.  The rationale for using these weighting factors is that
the interaction between a pair of bonded atoms is all (or mostly)
specified by the bond, angle, dihedral potentials, and thus the
non-bonded Lennard-Jones or Coulombic interaction between the pair of
atoms should be excluded (or reduced by a weighting factor).

.. note::

   These weighting factors are NOT used by :doc:`pair styles <pair_style>` that compute many-body interactions, since the
   "bonds" that result from such interactions are not permanent, but are
   created and broken dynamically as atom conformations change.  Examples
   of pair styles in this category are EAM, MEAM, Stillinger-Weber,
   Tersoff, COMB, AIREBO, and ReaxFF.  In fact, it generally makes no
   sense to define permanent bonds between atoms that interact via these
   potentials, though such bonds may exist elsewhere in your system,
   e.g. when using the :doc:`pair_style hybrid <pair_hybrid>` command.
   Thus LAMMPS ignores special\_bonds settings when many-body potentials
   are calculated.  Please note, that the existence of explicit bonds
   for atoms that are described by a many-body potential will alter the
   neighbor list and thus can render the computation of those interactions
   invalid, since those pairs are not only used to determine direct
   pairwise interactions but also neighbors of neighbors and more.
   The recommended course of action is to remove such bonds, or - if
   that is not possible - use a special bonds setting of 1.0 1.0 1.0.

.. note::

   Unlike some commands in LAMMPS, you cannot use this command
   multiple times in an incremental fashion: e.g. to first set the LJ
   settings and then the Coulombic ones.  Each time you use this command
   it sets all the coefficients to default values and only overrides the
   one you specify, so you should set all the options you need each time
   you use it.  See more details at the bottom of this page.

The Coulomb factors are applied to any Coulomb (charge interaction)
term that the potential calculates.  The LJ factors are applied to the
remaining terms that the potential calculates, whether they represent
LJ interactions or not.  The weighting factors are a scaling
pre-factor on the energy and force between the pair of atoms.  A value
of 1.0 means include the full interaction; a value of 0.0 means
exclude it completely.

The 1st of the 3 coefficients (LJ or Coulombic) is the weighting
factor on 1-2 atom pairs, which are pairs of atoms directly bonded to
each other.  The 2nd coefficient is the weighting factor on 1-3 atom
pairs which are those separated by 2 bonds (e.g. the two H atoms in a
water molecule).  The 3rd coefficient is the weighting factor on 1-4
atom pairs which are those separated by 3 bonds (e.g. the 1st and 4th
atoms in a dihedral interaction).  Thus if the 1-2 coefficient is set
to 0.0, then the pairwise interaction is effectively turned off for
all pairs of atoms bonded to each other.  If it is set to 1.0, then
that interaction will be at full strength.

.. note::

   For purposes of computing weighted pairwise interactions, 1-3
   and 1-4 interactions are not defined from the list of angles or
   dihedrals used by the simulation.  Rather, they are inferred
   topologically from the set of bonds specified when the simulation is
   defined from a data or restart file (see :doc:`read_data <read_data>` or
   :doc:`read_restart <read_restart>` commands).  Thus the set of
   1-2,1-3,1-4 interactions that the weights apply to is the same whether
   angle and dihedral potentials are computed or not, and remains the
   same even if bonds are constrained, or turned off, or removed during a
   simulation.

The two exceptions to this rule are (a) if the *angle* or *dihedral*
keywords are set to *yes* (see below), or (b) if the
:doc:`delete_bonds <delete_bonds>` command is used with the *special*
option that re-computes the 1-2,1-3,1-4 topologies after bonds are
deleted; see the :doc:`delete_bonds <delete_bonds>` command for more
details.

The *amber* keyword sets the 3 coefficients to 0.0, 0.0, 0.5 for LJ
interactions and to 0.0, 0.0, 0.8333 for Coulombic interactions, which
is the default for a commonly used version of the AMBER force field,
where the last value is really 5/6.  See :ref:`(Cornell) <Cornell>` for a
description of the AMBER force field.

The *charmm* keyword sets the 3 coefficients to 0.0, 0.0, 0.0 for both
LJ and Coulombic interactions, which is the default for a commonly
used version of the CHARMM force field.  Note that in pair styles
*lj/charmm/coul/charmm* and *lj/charmm/coul/long* the 1-4 coefficients
are defined explicitly, and these pairwise contributions are computed
as part of the charmm dihedral style - see the
:doc:`pair_coeff <pair_coeff>` and :doc:`dihedral_style <dihedral_style>`
commands for more information.  See :ref:`(MacKerell) <MacKerell>` for a
description of the CHARMM force field.

The *dreiding* keyword sets the 3 coefficients to 0.0, 0.0, 1.0 for both
LJ and Coulombic interactions, which is the default for the Dreiding
force field, as discussed in :ref:`(Mayo) <Mayo>`.

The *fene* keyword sets the 3 coefficients to 0.0, 1.0, 1.0 for both
LJ and Coulombic interactions, which is consistent with a
coarse-grained polymer model with :doc:`FENE bonds <bond_fene>`.  See
:ref:`(Kremer) <Kremer>` for a description of FENE bonds.

The *lj/coul*\ , *lj*\ , and *coul* keywords allow the 3 coefficients to
be set explicitly.  The *lj/coul* keyword sets both the LJ and
Coulombic coefficients to the same 3 values.  The *lj* and *coul*
keywords only set either the LJ or Coulombic coefficients.  Use both
of them if you wish to set the LJ coefficients to different values
than the Coulombic coefficients.

The *angle* keyword allows the 1-3 weighting factor to be ignored for
individual atom pairs if they are not listed as the first and last
atoms in any angle defined in the simulation or as 1,3 or 2,4 atoms in
any dihedral defined in the simulation.  For example, imagine the 1-3
weighting factor is set to 0.5 and you have a linear molecule with 4
atoms and bonds as follows: 1-2-3-4.  If your data file defines 1-2-3
as an angle, but does not define 2-3-4 as an angle or 1-2-3-4 as a
dihedral, then the pairwise interaction between atoms 1 and 3 will
always be weighted by 0.5, but different force fields use different
rules for weighting the pairwise interaction between atoms 2 and 4.
If the *angle* keyword is specified as *yes*\ , then the pairwise
interaction between atoms 2 and 4 will be unaffected (full weighting
of 1.0).  If the *angle* keyword is specified as *no* which is the
default, then the 2,4 interaction will also be weighted by 0.5.

The *dihedral* keyword allows the 1-4 weighting factor to be ignored
for individual atom pairs if they are not listed as the first and last
atoms in any dihedral defined in the simulation.  For example, imagine
the 1-4 weighting factor is set to 0.5 and you have a linear molecule
with 5 atoms and bonds as follows: 1-2-3-4-5.  If your data file
defines 1-2-3-4 as a dihedral, but does not define 2-3-4-5 as a
dihedral, then the pairwise interaction between atoms 1 and 4 will
always be weighted by 0.5, but different force fields use different
rules for weighting the pairwise interaction between atoms 2 and 5.
If the *dihedral* keyword is specified as *yes*\ , then the pairwise
interaction between atoms 2 and 5 will be unaffected (full weighting
of 1.0).  If the *dihedral* keyword is specified as *no* which is the
default, then the 2,5 interaction will also be weighted by 0.5.


----------


.. note::

   LAMMPS stores and maintains a data structure with a list of the
   1st, 2nd, and 3rd neighbors of each atom (within the bond topology of
   the system).  If new bonds are created (or molecules added containing
   atoms with more special neighbors), the size of this list needs to
   grow.  Note that adding a single bond always adds a new 1st neighbor
   but may also induce \*many\* new 2nd and 3rd neighbors, depending on the
   molecular topology of your system.  Using the *extra/special/per/atom*
   keyword to either :doc:`read_data <read_data>` or :doc:`create_box <create_box>`
   reserves empty space in the list for this N additional 1st, 2nd, or 3rd
   neighbors to be added.  If you do not do this, you may get an error
   when bonds (or molecules) are added.


----------


.. note::

   If you reuse this command in an input script, you should set all
   the options you need each time.  This command cannot be used a 2nd
   time incrementally.  E.g. these two commands:


.. parsed-literal::

   special_bonds lj 0.0 1.0 1.0
   special_bonds coul 0.0 0.0 1.0

are not the same as


.. parsed-literal::

   special_bonds lj 0.0 1.0 1.0 coul 0.0 0.0 1.0

In the first case you end up with (after the 2nd command):


.. parsed-literal::

   LJ: 0.0 0.0 0.0
   Coul: 0.0 0.0 1.0

while only in the second case, you get the desired settings of:


.. parsed-literal::

   LJ: 0.0 1.0 1.0
   Coul: 0.0 0.0 1.0

This happens because the LJ (and Coul) settings are reset to
their default values before modifying them, each time the
*special\_bonds* command is issued.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`delete_bonds <delete_bonds>`, :doc:`fix bond/create <fix_bond_create>`

Default
"""""""

All 3 Lennard-Jones and 3 Coulombic weighting coefficients = 0.0,
angle = no, dihedral = no.


----------


.. _Cornell:



**(Cornell)** Cornell, Cieplak, Bayly, Gould, Merz, Ferguson,
Spellmeyer, Fox, Caldwell, Kollman, JACS 117, 5179-5197 (1995).

.. _Kremer:



**(Kremer)** Kremer, Grest, J Chem Phys, 92, 5057 (1990).

.. _MacKerell:



**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).

.. _Mayo:



**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).


