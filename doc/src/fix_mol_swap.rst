.. index:: fix mol/swap

fix mol/swap command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID mol/swap N X itype jtype seed T keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* atom/swap = style name of this fix command
* N = invoke this fix every N steps
* X = number of swaps to attempt every N steps
* itype,jtype = two atom types (1-Ntypes or type label) to swap with each other
* seed = random # seed (positive integer)
* T = scaling temperature of the MC swaps (temperature units)
* zero or more keyword/value pairs may be appended to args
* keyword = *ke*

  .. parsed-literal::

       *ke* value = *no* or *yes*
         *no* = no conservation of kinetic energy after atom swaps
         *yes* = kinetic energy is conserved after atom swaps

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all mol/swap 100 1 2 3 29494 300.0 ke no
   fix mySwap fluid mol/swap 500 10 1 2 482798 1.0

   labelmap atom 1 A 2 B
   fix mySwap fluid mol/swap 500 10 A B 482798 1.0

Description
"""""""""""

This fix performs Monte Carlo swaps of two specified atom types within
a randomly selected molecule.  Two possible use cases are as follows.

First, consider a mixture of some molecules with atoms of itype and
other molecules with atoms of jtype.  The fix will select a random
molecule and attempt to swap all the itype atoms to jtype for the
first kind of molecule, or all the jtype atoms to itype for the second
kind.  Because the swap will only take place if it is energetically
favorable, the fix can be used to determine the miscibility of 2
different kinds of molecules much more quickly than just dynamics
would do it.

Second, consider diblock co-polymers with two types of monomers itype
and jtype.  The fix will select a random molecule and attempt to do a
itype <--> jtype swap of all those monomers within the molecule.  Thus
the fix can be used to find the energetically favorable fractions of
two flavors of diblock co-polymers.

Intra-molecular swaps of atom types are attempted every N timesteps.  On
that timestep, X swaps are attempted.  For each attempt a single
molecule ID is randomly selected.  The range of possible molecule IDs
from loID to hiID is pre-computed before each run begins.  The
loID/hiID is set for the molecule with the smallest/largest ID which
has any itype or jtype atoms in it.  Note that if you define a system
with many molecule IDs between loID and hiID which have no itype or
jtype atoms, then the fix will be inefficient at performing swaps.
Also note that if atoms with molecule ID = 0 exist, they are not
considered molecules by this fix; they are assumed to be solvent atoms
or molecules.

Candidate atoms for swapping must also be in the fix group.  Atoms
within the selected molecule which are not itype or jtype are ignored.

When an atom is swapped from itype to jtype (or vice versa), if
charges are defined, the charge values for itype versus jtype atoms
are also swapped.  This requires that all itype atoms in the system
have the same charge value.  Likewise all jtype atoms in the system
must have the same charge value.  If this is not the case, LAMMPS
issues a warning that it cannot swap charge values.

If the *ke* keyword is set to yes, which is the default, and the
masses of itype and jtype atoms are different, then when a swap
occurs, the velocity of the swapped atom is rescaled by the sqrt of
the mass ratio, so as to conserve the kinetic energy of the atom.

----------

The potential energy of the entire system is computed before and after
each swap is performed within a single molecule.  The specified
temperature T is used in the Metropolis criterion to accept or reject
the attempted swap.  If the swap is rejected all swapped values are
reversed.

The potential energy calculations can include systems and models with
the following features:

* manybody pair styles, including EAM
* hybrid pair styles
* long-range electrostatics (kspace)
* triclinic systems
* potential energy contributions from other fixes

For the last bullet point, fixes can have an associated potential
energy. Examples of such fixes include: :doc:`efield <fix_efield>`,
:doc:`gravity <fix_gravity>`, :doc:`addforce <fix_addforce>`,
:doc:`langevin <fix_langevin>`, :doc:`restrain <fix_restrain>`,
:doc:`temp/berendsen <fix_temp_berendsen>`, :doc:`temp/rescale
<fix_temp_rescale>`, and :doc:`wall fixes <fix_wall>`.  For that
energy to be included in the total potential energy of the system (the
quantity used for the swap accept/reject decision), you MUST enable
the :doc:`fix_modify <fix_modify>` *energy* option for that fix.  The
doc pages for individual :doc:`fix <fix>` commands specify if this
should be done.

.. note::

  One comment on computational efficiency.  If the cutoff lengths
  defined for the pair style are different for itype versus jtype
  atoms (for any of their interactions with any other atom type), then
  a new neighbor list needs to be generated for every attempted swap.
  This is potentially expensive if N is small or X is large.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the random number
generator seed, the next timestep for MC exchanges, the number of
exchange attempts and successes etc.  See the :doc:`read_restart
<read_restart>` command for info on how to re-specify a fix in an
input script that reads a restart file, so that the operation of the
fix continues in an uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep <reset_timestep>`.
   The fix will try to detect it and stop with an error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a global vector of length 2, which can be accessed
by various :doc:`output commands <Howto_output>`.  The vector values are
the following global cumulative quantities:

* 1 = swap attempts
* 2 = swap accepts

The vector values calculated by this fix are "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix atom/swap <fix_atom_swap>`, :doc:`fix gcmc <fix_gcmc>`

Default
"""""""

The option default is ke = yes.
