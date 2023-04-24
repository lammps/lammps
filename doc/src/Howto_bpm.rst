Bonded particle models
======================

The BPM package implements bonded particle models which can be used to
simulate mesoscale solids.  Solids are constructed as a collection of
particles, which each represent a coarse-grained region of space much
larger than the atomistic scale.  Particles within a solid region are
then connected by a network of bonds to provide solid elasticity.

Unlike traditional bonds in molecular dynamics, the equilibrium bond
length can vary between bonds. Bonds store the reference state.  This
includes setting the equilibrium length equal to the initial distance
between the two particles, but can also include data on the bond
orientation for rotational models. This produces a stress-free initial
state. Furthermore, bonds are allowed to break under large strains,
producing fracture. The examples/bpm directory has sample input scripts
for simulations of the fragmentation of an impacted plate and the
pouring of extended, elastic bodies.

----------

Bonds can be created using a :doc:`read data <read_data>` or
:doc:`create bonds <create_bonds>` command. Alternatively, a
:doc:`molecule <molecule>` template with bonds can be used with
:doc:`fix deposit <fix_deposit>` or :doc:`fix pour <fix_pour>` to create
solid grains.

In this implementation, bonds store their reference state when they are
first computed in the setup of the first simulation run. Data is then
preserved across run commands and is written to :doc:`binary restart
files <restart>` such that restarting the system will not reset the
reference state of a bond. Bonds that are created midway into a run,
such as those created by pouring grains using :doc:`fix pour
<fix_pour>`, are initialized on that timestep.

----------

Currently, there are two types of bonds included in the BPM package. The
first bond style, :doc:`bond bpm/spring <bond_bpm_spring>`, only applies
pairwise, central body forces. Point particles must have :doc:`bond atom
style <atom_style>` and may be thought of as nodes in a spring
network. Alternatively, the second bond style, :doc:`bond bpm/rotational
<bond_bpm_rotational>`, resolves tangential forces and torques arising
with the shearing, bending, and twisting of the bond due to rotation or
displacement of particles.  Particles are similar to those used in the
:doc:`granular package <Howto_granular>`, :doc:`atom style sphere
<atom_style>`. However, they must also track the current orientation of
particles and store bonds, and therefore use a :doc:`bpm/sphere atom
style <atom_style>`.  This also requires a unique integrator :doc:`fix
nve/bpm/sphere <fix_nve_bpm_sphere>` which numerically integrates
orientation similar to :doc:`fix nve/asphere <fix_nve_asphere>`.

In addition to bond styles, a new pair style :doc:`pair bpm/spring
<pair_bpm_spring>` was added to accompany the bpm/spring bond
style. This pair style is simply a hookean repulsion with similar
velocity damping as its sister bond style.

----------

Bond data can be output using a combination of standard LAMMPS commands.
A list of IDs for bonded atoms can be generated using the
:doc:`compute property/local <compute_property_local>` command.
Various properties of bonds can be computed using the
:doc:`compute bond/local <compute_bond_local>` command. This
command allows one to access data saved to the bond's history,
such as the reference length of the bond. More information on
bond history data can be found on the documentation pages for the specific
BPM bond styles. Finally, this data can be output using a :doc:`dump local <dump>`
command. As one may output many columns from the same compute, the
:doc:`dump modify <dump_modify>` *colname* option may be used to provide
more helpful column names. An example of this procedure is found in
/examples/bpm/pour/. External software, such as OVITO, can read these dump
files to render bond data.

----------

As bonds can be broken between neighbor list builds, the
:doc:`special_bonds <special_bonds>` command works differently for BPM
bond styles. There are two possible settings which determine how pair
interactions work between bonded particles.  First, one can overlay
pair forces with bond forces such that all bonded particles also
feel pair interactions. This can be accomplished by using the *overlay/pair*
keyword present in all bpm bond styles and by using the following special
bond settings

   .. code-block:: LAMMPS

      special_bonds lj/coul 1 1 1

Alternatively, one can turn off all pair interactions between bonded
particles. Unlike :doc:`bond quartic <bond_quartic>`, this is not done
by subtracting pair forces during the bond computation, but rather by
dynamically updating the special bond list. This is the default behavior
of BPM bond styles and is done by updating the 1-2 special bond list as
bonds break.  To do this, LAMMPS requires :doc:`newton <newton>` bond off
such that all processors containing an atom know when a bond breaks.
Additionally, one must use the following special bond settings

   .. code-block:: LAMMPS

      special_bonds lj 0 1 1 coul 1 1 1

These settings accomplish two goals. First, they turn off 1-3 and 1-4
special bond lists, which are not currently supported for BPMs. As
BPMs often have dense bond networks, generating 1-3 and 1-4 special
bond lists is expensive.  By setting the lj weight for 1-2 bonds to
zero, this turns off pairwise interactions.  Even though there are no
charges in BPM models, setting a nonzero coul weight for 1-2 bonds
ensures all bonded neighbors are still included in the neighbor list
in case bonds break between neighbor list builds.

To monitor the fracture of bonds in the system, all BPM bond styles
have the ability to record instances of bond breakage to output using
the :doc:`dump local <dump>` command. Since one may frequently output
a list of broken bonds and the time they broke, the
:doc:`dump modify <dump_modify>` option *header no* may be useful to
avoid repeatedly printing the header of the dump file. An example of
this procedure is found in /examples/bpm/impact/. Additionally,
one can use :doc:`compute nbond/atom <compute_nbond_atom>` to tally the
current number of bonds per atom.

See the :doc:`Howto <Howto_broken_bonds>` page on broken bonds for
more information.

----------

While LAMMPS has many utilities to create and delete bonds, *only*
the following are currently compatible with BPM bond styles:

* :doc:`create_bonds <create_bonds>`
* :doc:`delete_bonds <delete_bonds>`
* :doc:`fix bond/create <fix_bond_create>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`fix bond/swap <fix_bond_swap>`

.. note::

   The :doc:`create_bonds <create_bonds>` command requires certain
   :doc:`special_bonds <special_bonds>` settings.  To subtract pair
   interactions, one will need to switch between different *special_bonds*
   settings in the input script. An example is found in
   ``examples/bpm/impact``.
