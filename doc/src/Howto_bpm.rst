Bonded particle models
===============

Bonded particle models are used to simulate mesoscale solids.
Solids are constructed as a collection of particles which each
represent a coarse-grained region of space much larger than the
atomistic scale. Particles within a solid region are then connected
by a network of bonds to provide solid elasticity.

Unlike traditional bonds in molecular dynamics, the equilibrium
bond length can vary between bonds. Bonds store the reference state.
This includes setting the equilibrium length equal to the initial
distance between the two particles but can also include data on the
bond orientation for rotational models. This produces a stress free
initial state. Furthermore, bonds are allowed to break under large
strains producing fracture.

Bonds can be created using a :doc:`read data <read_data>`
or :doc:`create bond <create_bond>` command. Alternatively, a
:doc:`molecule <molecule>` template with bonds can be used with
:doc:`fix deposit <fix_deposit>` or :doc:`fix pour <fix_pour>` to
create solid grains.
In this implementation, bonds store their reference state when they
are first computed in the setup of a simulation run. Data is then
preserved across run commands and is written to :doc:`binary restart files <restart>`
such that restarting the system will not reset the reference state of a bond.

As bonds can be broken between neighbor list builds, :doc:`special_bonds <special_bonds>`
work differently for BPM bond styles. There are two possible special
bond settings which determine how pair interactions work between bonded
particles. First, one can simply overlay pair interactions such that all
bonded particles also feel pair interactions. This can be accomplished by
simply turning off all special bonds by setting

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

Alternatively, one can censor all pair interactions between bonded particles.
Unlike :doc:`bond quartic <bond_quartic>`, this is not done by subtracting
pair forces during the bond computation but rather by dynamically updating
the special bond list. To do this, one must both define an instance of
:doc:`fix update/special/bonds <fix_update_special_bonds>` and have the special bond
settings

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

This fix ensures the 1-2 special bond list remains updated as bonds break. The fix
also requires :doc:`newton <newton>` bond off such that whena bond breaks between
atoms across multiple processors, all processors are aware of the event.
The special bond settings then accomplish two tasks. First, they turns off 1-3 and
1-4 special bond lists, which are not currently supported for BPMs. As BPMs often
have dense bond networks, generating 1-3 and 1-4 special bond lists is expensive.
By setting the lj weight for 1-2 bonds to zero, this censors pairwise interactions.
However, setting a nonzero coul weight for 1-2 bonds ensures all bonded
neighbors are included in the neighbor list. All bonded neighbors must be included
in neighbor lists as they could become unbonded at any timestep.

Currently there are two types of bonds included in this package. The first
bond style, :doc:`bond bpm/spring <bond_bpm_spring>`, only applies pairwise,
central body forces. Point particles must have :doc:`bond atom style <atom_style>`
and may be thought of as nodes in a spring network. Alternatively,
the second bond style, :doc:`bond bpm/rotational <bond_bpm_rotational>`,
resolves tangential forces and torques arising with the shearing, bending,
and twisting of the bond due to rotation or displacement of particles.
Particles are similar to those used in the :doc:`granular package <Howto_granular>`,
:doc:`atom style sphere <atom_style>`. However, they must also track the
current orientation of particles and therefore use a derived :doc:`atom style sphere/bpm <atom_style>`.
This also requires a unique integrator :doc:`fix nve/sphere/bpm <fix_nve_sphere_bpm>`
which numerically integrates orientation similar to :doc:`fix nve/asphere <fix_nve_asphere>`.

To monitor the fracture of bonds in the system, all BPM bond styles
can be associated with an instance of :doc:`fix store/local <fix_store_local>`
to record all instances of bond breakage for output. Additionally, one can use
:doc:`compute nbond/atom <compute_nbond_atom>` to tally the current number of bonds per atom.

In addition to bond styles, a new pair style :doc:`pair bpm/spring <pair_bpm_spring>` was added
to accompany the bpm/spring bond style. This pair style is simply a hookean repulsion with
similar velocity damping as its sister bond style.
