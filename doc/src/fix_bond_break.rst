.. index:: fix bond/break

fix bond/break command
======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID bond/break Nevery bondtype Rmax keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* bond/break = style name of this fix command
* Nevery = attempt bond breaking every this many steps
* bondtype = type of bonds to break
* Rmax = bond longer than Rmax can break (distance units)
* zero or more keyword/value pairs may be appended to args
* keyword = *prob*
  
  .. parsed-literal::
  
       *prob* values = fraction seed
         fraction = break a bond with this probability if otherwise eligible
         seed = random number seed (positive integer)



Examples
""""""""


.. parsed-literal::

   fix 5 all bond/break 10 2 1.2
   fix 5 polymer bond/break 1 1 2.0 prob 0.5 49829

Description
"""""""""""

Break bonds between pairs of atoms as a simulation runs according to
specified criteria.  This can be used to model the dissolution of a
polymer network due to stretching of the simulation box or other
deformations.  In this context, a bond means an interaction between a
pair of atoms computed by the :doc:`bond_style <bond_style>` command.
Once the bond is broken it will be permanently deleted, as will all
angle, dihedral, and improper interactions that bond is part of.

This is different than a :doc:`pairwise <pair_style>` bond-order
potential such as Tersoff or AIREBO which infers bonds and many-body
interactions based on the current geometry of a small cluster of atoms
and effectively creates and destroys bonds and higher-order many-body
interactions from timestep to timestep as atoms move.

A check for possible bond breakage is performed every *Nevery*
timesteps.  If two bonded atoms I,J are further than a distance *Rmax*
of each other, if the bond is of type *bondtype*\ , and if both I and J
are in the specified fix group, then I,J is labeled as a "possible"
bond to break.

If several bonds involving an atom are stretched, it may have multiple
possible bonds to break.  Every atom checks its list of possible bonds
to break and labels the longest such bond as its "sole" bond to break.
After this is done, if atom I is bonded to atom J in its sole bond,
and atom J is bonded to atom I in its sole bond, then the I,J bond is
"eligible" to be broken.

Note that these rules mean an atom will only be part of at most one
broken bond on a given timestep.  It also means that if atom I chooses
atom J as its sole partner, but atom J chooses atom K is its sole
partner (due to Rjk > Rij), then this means atom I will not be part of
a broken bond on this timestep, even if it has other possible bond
partners.

The *prob* keyword can effect whether an eligible bond is actually
broken.  The *fraction* setting must be a value between 0.0 and 1.0.
A uniform random number between 0.0 and 1.0 is generated and the
eligible bond is only broken if the random number < fraction.

When a bond is broken, data structures within LAMMPS that store bond
topology are updated to reflect the breakage.  Likewise, if the bond
is part of a 3-body (angle) or 4-body (dihedral, improper)
interaction, that interaction is removed as well.  These changes
typically affect pairwise interactions between atoms that used to be
part of bonds, angles, etc.

.. note::

   One data structure that is not updated when a bond breaks are
   the molecule IDs stored by each atom.  Even though one molecule
   becomes two molecules due to the broken bond, all atoms in both new
   molecules retain their original molecule IDs.

Computationally, each timestep this fix operates, it loops over all
the bonds in the system and computes distances between pairs of bonded
atoms.  It also communicates between neighboring processors to
coordinate which bonds are broken.  Moreover, if any bonds are broken,
neighbor lists must be immediately updated on the same timestep.  This
is to insure that any pairwise interactions that should be turned "on"
due to a bond breaking, because they are no longer excluded by the
presence of the bond and the settings of the
:doc:`special_bonds <special_bonds>` command, will be immediately
recognized.  All of these operations increase the cost of a timestep.
Thus you should be cautious about invoking this fix too frequently.

You can dump out snapshots of the current bond topology via the :doc:`dump local <dump>` command.

.. note::

   Breaking a bond typically alters the energy of a system.  You
   should be careful not to choose bond breaking criteria that induce a
   dramatic change in energy.  For example, if you define a very stiff
   harmonic bond and break it when 2 atoms are separated by a distance
   far from the equilibrium bond length, then the 2 atoms will be
   dramatically released when the bond is broken.  More generally, you
   may need to thermostat your system to compensate for energy changes
   resulting from broken bonds (and angles, dihedrals, impropers).


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes two statistics which it stores in a global vector of
length 2, which can be accessed by various :doc:`output commands <Howto_output>`.  The vector values calculated by this fix
are "intensive".

These are the 2 quantities:

* (1) # of bonds broken on the most recent breakage timestep
* (2) cumulative # of bonds broken

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix bond/create <fix_bond_create>`, :doc:`fix bond/react <fix_bond_react>`, :doc:`fix bond/swap <fix_bond_swap>`,
:doc:`dump local <dump>`, :doc:`special_bonds <special_bonds>`

Default
"""""""

The option defaults are prob = 1.0.


