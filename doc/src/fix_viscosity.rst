.. index:: fix viscosity

fix viscosity command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID viscosity N vdim pdim Nbin keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* viscosity = style name of this fix command
* N = perform momentum exchange every N steps
* vdim = *x* or *y* or *z* = which momentum component to exchange
* pdim = *x* or *y* or *z* = direction of momentum transfer
* Nbin = # of layers in pdim direction (must be even number)
* zero or more keyword/value pairs may be appended
* keyword = *swap* or *target*
  
  .. parsed-literal::
  
       *swap* value = Nswap = number of swaps to perform every N steps
       *vtarget* value = V or INF = target velocity of swap partners (velocity units)



Examples
""""""""


.. parsed-literal::

   fix 1 all viscosity 100 x z 20
   fix 1 all viscosity 50 x z 20 swap 2 vtarget 1.5

Description
"""""""""""

Use the Muller-Plathe algorithm described in :ref:`this paper <Muller-Plathe2>` to exchange momenta between two particles in
different regions of the simulation box every N steps.  This induces a
shear velocity profile in the system.  As described below this enables
a viscosity of the fluid to be calculated.  This algorithm is
sometimes called a reverse non-equilibrium MD (reverse NEMD) approach
to computing viscosity.  This is because the usual NEMD approach is to
impose a shear velocity profile on the system and measure the response
via an off-diagonal component of the stress tensor, which is
proportional to the momentum flux.  In the Muller-Plathe method, the
momentum flux is imposed, and the shear velocity profile is the
system's response.

The simulation box is divided into *Nbin* layers in the *pdim*
direction, where the layer 1 is at the low end of that dimension and
the layer *Nbin* is at the high end.  Every N steps, Nswap pairs of
atoms are chosen in the following manner.  Only atoms in the fix group
are considered.  Nswap atoms in layer 1 with positive velocity
components in the *vdim* direction closest to the target value *V* are
selected.  Similarly, Nswap atoms in the "middle" layer (see below) with
negative velocity components in the *vdim* direction closest to the
negative of the target value *V* are selected.  The two sets of Nswap
atoms are paired up and their *vdim* momenta components are swapped
within each pair.  This resets their velocities, typically in opposite
directions.  Over time, this induces a shear velocity profile in the
system which can be measured using commands such as the following,
which writes the profile to the file tmp.profile:


.. parsed-literal::

   compute layers all chunk/atom bin/1d z lower 0.05 units reduced
   fix f1 all ave/chunk 100 10 1000 layers vx file tmp.profile

Note that by default, Nswap = 1 and vtarget = INF, though this can be
changed by the optional *swap* and *vtarget* keywords.  When vtarget =
INF, one or more atoms with the most positive and negative velocity
components are selected.  Setting these parameters appropriately, in
conjunction with the swap rate N, allows the momentum flux rate to be
adjusted across a wide range of values, and the momenta to be
exchanged in large chunks or more smoothly.

The "middle" layer for momenta swapping is defined as the *Nbin*\ /2 + 1
layer.  Thus if *Nbin* = 20, the two swapping layers are 1 and 11.
This should lead to a symmetric velocity profile since the two layers
are separated by the same distance in both directions in a periodic
sense.  This is why *Nbin* is restricted to being an even number.

As described below, the total momentum transferred by these velocity
swaps is computed by the fix and can be output.  Dividing this
quantity by time and the cross-sectional area of the simulation box
yields a momentum flux.  The ratio of momentum flux to the slope of
the shear velocity profile is proportional to the viscosity of the
fluid, in appropriate units.  See the :ref:`Muller-Plathe paper <Muller-Plathe2>` for details.

.. note::

   If your system is periodic in the direction of the momentum
   flux, then the flux is going in 2 directions.  This means the
   effective momentum flux in one direction is reduced by a factor of 2.
   You will see this in the equations for viscosity in the Muller-Plathe
   paper.  LAMMPS is simply tallying momentum which does not account for
   whether or not your system is periodic; you must use the value
   appropriately to yield a viscosity for your system.

.. note::

   After equilibration, if the velocity profile you observe is not
   linear, then you are likely swapping momentum too frequently and are
   not in a regime of linear response.  In this case you cannot
   accurately infer a viscosity and should try increasing the Nevery
   parameter.

An alternative method for calculating a viscosity is to run a NEMD
simulation, as described on the :doc:`Howto nemd <Howto_nemd>` doc page.
NEMD simulations deform the simulation box via the :doc:`fix deform <fix_deform>` command.  Thus they cannot be run on a charged
system using a :doc:`PPPM solver <kspace_style>` since PPPM does not
currently support non-orthogonal boxes.  Using fix viscosity keeps the
box orthogonal; thus it does not suffer from this limitation.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
momentum transferred between the bottom and middle of the simulation
box (in the *pdim* direction) is stored as a scalar quantity by this
fix.  This quantity is zeroed when the fix is defined and accumulates
thereafter, once every N steps.  The units of the quantity are
momentum = mass\*velocity.  The scalar value calculated by this fix is
"intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Swaps conserve both momentum and kinetic energy, even if the masses of
the swapped atoms are not equal.  Thus you should not need to
thermostat the system.  If you do use a thermostat, you may want to
apply it only to the non-swapped dimensions (other than *vdim*\ ).

LAMMPS does not check, but you should not use this fix to swap
velocities of atoms that are in constrained molecules, e.g. via :doc:`fix shake <fix_shake>` or :doc:`fix rigid <fix_rigid>`.  This is because
application of the constraints will alter the amount of transferred
momentum.  You should, however, be able to use flexible molecules.
See the :ref:`Maginn paper <Maginn>` for an example of using this algorithm
in a computation of alcohol molecule properties.

When running a simulation with large, massive particles or molecules
in a background solvent, you may want to only exchange momenta between
solvent particles.

Related commands
""""""""""""""""

:doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix thermal/conductivity <fix_thermal_conductivity>`

Default
"""""""

The option defaults are swap = 1 and vtarget = INF.


----------


.. _Muller-Plathe2:



**(Muller-Plathe)** Muller-Plathe, Phys Rev E, 59, 4894-4898 (1999).

.. _Maginn:



**(Maginn)** Kelkar, Rafferty, Maginn, Siepmann, Fluid Phase Equilibria,
260, 218-231 (2007).


