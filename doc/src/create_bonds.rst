.. index:: create_bonds

create_bonds command
====================

Syntax
""""""

.. code-block:: LAMMPS

   create_bonds style args ... keyword value ...

* style = *many* or *single/bond* or *single/angle* or *single/dihedral* or *single/improper*

.. parsed-literal::

     *many* args = group-ID group2-ID btype rmin rmax
       group-ID = ID of first group
       group2-ID = ID of second group, bonds will be between atoms in the 2 groups
       btype = bond type of created bonds
       rmin = minimum distance between pair of atoms to bond together
       rmax = maximum distance between pair of atoms to bond together
     *single/bond* args = btype batom1 batom2
       btype = bond type of new bond
       batom1,batom2 = atom IDs for two atoms in bond
     *single/angle* args = atype aatom1 aatom2 aatom3
       atype = angle type of new angle
       aatom1,aatom2,aatom3 = atom IDs for three atoms in angle
     *single/dihedral* args = dtype datom1 datom2 datom3 datom4
       dtype = dihedral type of new dihedral
       datom1,datom2,datom3,datom4 = atom IDs for four atoms in dihedral
     *single/improper* args = itype iatom1 iatom2 iatom3 iatom4
       itype = improper type of new improper
       iatom1,iatom2,iatom3,iatom4 = atom IDs for four atoms in improper

* zero or more keyword/value pairs may be appended
* keyword = *special*

.. parsed-literal::

     *special* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   create_bonds many all all 1 1.0 1.2
   create_bonds many surf solvent 3 2.0 2.4
   create_bonds single/bond 1 1 2
   create_bonds single/angle 5 52 98 107 special no
   create_bonds single/dihedral 2 4 19 27 101
   create_bonds single/improper 3 23 26 31 57

Description
"""""""""""

Create bonds between pairs of atoms that meet a specified distance
criteria.  Or create a single bond, angle, dihedral or improper between 2, 3,
or 4 specified atoms.

The new bond (angle, dihedral, improper) interactions will then be computed
during a simulation by the bond (angle, dihedral, improper) potential defined by
the :doc:`bond_style <bond_style>`, :doc:`bond_coeff <bond_coeff>`,
:doc:`angle_style <angle_style>`, :doc:`angle_coeff <angle_coeff>`,
:doc:`dihedral_style <dihedral_style>`,
:doc:`dihedral_coeff <dihedral_coeff>`, :doc:`improper_style <improper_style>`,
:doc:`improper_coeff <improper_coeff>` commands.

The *many* style is useful for adding bonds to a system (e.g., between
nearest neighbors in a lattice of atoms) without having to enumerate
all the bonds in the data file read by the :doc:`read_data <read_data>`
command.

The *single* styles are useful for adding bonds, angles, dihedrals, and
impropers to a system incrementally, then continuing a simulation.

Note that this command does not auto-create any angle, dihedral, or improper
interactions when a bond is added, nor does it auto-create any bonds
when an angle, dihedral, or improper is added.  It also will not auto-create
any angles when a dihedral or improper is added.  Thus, the flexibility of this
command is limited.  It can be used several times to create different types of
bond at different distances, but it cannot typically auto-create all the
bonds or angles or dihedrals or impropers that would normally be defined in a
data file for a complex system of molecules.

.. note::

   If the system has no bonds (angles, dihedrals, impropers) to begin
   with, or if more bonds per atom are being added than currently exist,
   then you must ensure that the number of bond types and the maximum
   number of bonds per atom are set to large enough values, and
   similarly for angles, dihedrals, impropers, and special neighbors,
   otherwise an error may occur when too many bonds (angles, dihedrals,
   impropers) are added to an atom.  If the :doc:`read_data <read_data>`
   command is used to define the system, these parameters can be set via
   its optional *extra/bond/types*, *extra/bond/per/atom*, and similar
   keywords to the command.  If the :doc:`create_box <create_box>`
   command is used to define the system, these two parameters can be set
   via its optional *bond/types* and *extra/bond/per/atom* arguments,
   and similarly for angles, dihedrals, and impropers.  See the
   corresponding documentation pages for these two commands for details.

----------

The *many* style will create bonds between pairs of atoms :math:`I,J`,
where :math:`I` is in one of the two specified groups and :math:`J` is in the
other.  The two groups can be the same (e.g., group "all").  The created bonds
will be of bond type *btype*, where *btype* must be a value between 1 and the
number of bond types defined.

For a bond to be created, an :math:`I,J` pair of atoms must be a distance
:math:`D` apart such that :math:`r_\text{min} \le D \le r_\text{max}`.

The following settings must have been made in an input script before
the *many* style is used:

* special_bonds weight for 1--2 interactions must be 0.0
* a :doc:`pair_style <pair_style>` must be defined
* no :doc:`kspace_style <kspace_style>` defined
* minimum :doc:`pair_style <pair_style>` cutoff + :doc:`neighbor <neighbor>`
  skin :math:`\ge r_\text{max}`

These settings are required so that a neighbor list can be created to
search for nearby atoms.  Pairs of atoms that are already bonded cannot
appear in the neighbor list, to avoid creation of duplicate bonds.  The
neighbor list for all atom type pairs must also extend to a distance
that encompasses the *rmax* for new bonds to create.  When using
periodic boundary conditions, the box length in each periodic dimension
must be larger than *rmax*, so that no bonds are created between the
system and its own periodic image.

.. note::

   If you want to create bonds between pairs of 1--3 or 1--4 atoms in
   the current bond topology, then you need to use :doc:`special_bonds
   lj 0 1 1 <special_bonds>` to ensure those pairs appear in the
   neighbor list.  They will not appear with the default special_bonds
   settings, which are zero for 1--2, 1--3, and 1--4 atoms.  1--3 or 1--4
   atoms are those which are two hops or three hops apart in the bond
   topology.

An additional requirement for this style is that your system must be
ready to perform a simulation.  This means, for example, that all
:doc:`pair_style <pair_style>` coefficients be set via the
:doc:`pair_coeff <pair_coeff>` command.  A :doc:`bond_style <bond_style>`
command and all bond coefficients must also be set, even if no bonds
exist before this command is invoked.  This is because the building of
neighbor list requires initialization and setup of a simulation,
similar to what a :doc:`run <run>` command would require.

Note that you can change any of these settings after this command
executes (e.g., if you wish to use long-range Coulombic interactions)
via the :doc:`kspace_style <kspace_style>` command for your subsequent
simulation.

----------

The *single/bond* style creates a single bond of type *btype* between
two atoms with IDs *batom1* and *batom2*\ .  *Btype* must be a value
between 1 and the number of bond types defined.

The *single/angle* style creates a single angle of type *atype*
between three atoms with IDs *aatom1*, *aatom2*, and *aatom3*\ .  The
ordering of the atoms is the same as in the *Angles* section of a data
file read by the :doc:`read_data <read_data>` command (i.e., the three atoms
are ordered linearly within the angle; the central atom is *aatom2*).
*Atype* must be a value between 1 and the number of angle types
defined.

The *single/dihedral* style creates a single dihedral of type *dtype*
between four atoms with IDs *datom1*, *datom2*, *datom3*, and *datom4*\ .  The
ordering of the atoms is the same as in the *Dihedrals* section of a data file
read by the :doc:`read_data <read_data>` command.  I.e. the 4 atoms are ordered
linearly within the dihedral.  *dtype* must be a value between 1 and
the number of dihedral types defined.

The *single/improper* style creates a single improper of type *itype*
between four atoms with IDs *iatom1*, *iatom2*, *iatom3*, and *iatom4*\ .  The
ordering of the atoms is the same as in the *Impropers* section of a data file
read by the :doc:`read_data <read_data>` command.  I.e. the 4 atoms are ordered
linearly within the improper.  *itype* must be a value between 1 and
the number of improper types defined.

----------

The keyword *special* controls whether an internal list of special
bonds is created after one or more bonds, or a single angle, dihedral, or
improper is added to the system.

The default value is *yes*\ .  A value of *no* cannot be used
with the *many* style.

This is an expensive operation since the bond topology for the system
must be walked to find all 1--2, 1--3, and 1--4 interactions to store in an
internal list, which is used when pairwise interactions are weighted;
see the :doc:`special_bonds <special_bonds>` command for details.

Thus if you are adding a few bonds or a large list of angles all at
the same time, by using this command repeatedly, it is more efficient
to only trigger the internal list to be created once, after the last
bond (or angle, or dihedral, or improper) is added:

.. code-block:: LAMMPS

   create_bonds single/bond 5 52 98 special no
   create_bonds single/bond 5 73 74 special no
   ...
   create_bonds single/bond 5 17 386 special no
   create_bonds single/bond 4 112 183 special yes

Note that you **must** ensure the internal list is rebuilt after the last
bond (angle, dihedral, improper) is added, *before* performing a simulation.
Otherwise, pairwise interactions will not be properly excluded or
weighted.  LAMMPS does **not** check that you have done this correctly.

----------

Restrictions
""""""""""""

This command cannot be used with molecular systems defined using
molecule template files via the :doc:`molecule <molecule>` and
:doc:`atom_style template <atom_style>` commands.

For style *many*, no :doc:`kspace style <kspace_style>` must be
defined. Also, the *rmax* value must be smaller than any periodic box
length and the neighbor list cutoff (largest pair cutoff plus neighbor
skin).

Related commands
""""""""""""""""

:doc:`create_atoms <create_atoms>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

The keyword default is special = yes.
