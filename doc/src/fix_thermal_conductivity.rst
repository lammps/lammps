.. index:: fix thermal/conductivity

fix thermal/conductivity command
================================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID thermal/conductivity N edim Nbin keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* thermal/conductivity = style name of this fix command
* N = perform kinetic energy exchange every N steps
* edim = *x* or *y* or *z* = direction of kinetic energy transfer
* Nbin = # of layers in edim direction (must be even number)
* zero or more keyword/value pairs may be appended
* keyword = *swap*
  
  .. parsed-literal::
  
       *swap* value = Nswap = number of swaps to perform every N steps



Examples
""""""""


.. parsed-literal::

   fix 1 all thermal/conductivity 100 z 20
   fix 1 all thermal/conductivity 50 z 20 swap 2

Description
"""""""""""

Use the Muller-Plathe algorithm described in :ref:`this paper <Muller-Plathe1>` to exchange kinetic energy between two particles
in different regions of the simulation box every N steps.  This
induces a temperature gradient in the system.  As described below this
enables the thermal conductivity of a material to be calculated.  This
algorithm is sometimes called a reverse non-equilibrium MD (reverse
NEMD) approach to computing thermal conductivity.  This is because the
usual NEMD approach is to impose a temperature gradient on the system
and measure the response as the resulting heat flux.  In the
Muller-Plathe method, the heat flux is imposed, and the temperature
gradient is the system's response.

See the :doc:`compute heat/flux <compute_heat_flux>` command for details
on how to compute thermal conductivity in an alternate way, via the
Green-Kubo formalism.

The simulation box is divided into *Nbin* layers in the *edim*
direction, where the layer 1 is at the low end of that dimension and
the layer *Nbin* is at the high end.  Every N steps, Nswap pairs of
atoms are chosen in the following manner.  Only atoms in the fix group
are considered.  The hottest Nswap atoms in layer 1 are selected.
Similarly, the coldest Nswap atoms in the "middle" layer (see below)
are selected.  The two sets of Nswap atoms are paired up and their
velocities are exchanged.  This effectively swaps their kinetic
energies, assuming their masses are the same.  If the masses are
different, an exchange of velocities relative to center of mass motion
of the 2 atoms is performed, to conserve kinetic energy.  Over time,
this induces a temperature gradient in the system which can be
measured using commands such as the following, which writes the
temperature profile (assuming z = edim) to the file tmp.profile:


.. parsed-literal::

   compute   ke all ke/atom
   variable  temp atom c_ke/1.5
   compute   layers all chunk/atom bin/1d z lower 0.05 units reduced
   fix       3 all ave/chunk 10 100 1000 layers v_temp file tmp.profile

Note that by default, Nswap = 1, though this can be changed by the
optional *swap* keyword.  Setting this parameter appropriately, in
conjunction with the swap rate N, allows the heat flux to be adjusted
across a wide range of values, and the kinetic energy to be exchanged
in large chunks or more smoothly.

The "middle" layer for velocity swapping is defined as the *Nbin*\ /2 +
1 layer.  Thus if *Nbin* = 20, the two swapping layers are 1 and 11.
This should lead to a symmetric temperature profile since the two
layers are separated by the same distance in both directions in a
periodic sense.  This is why *Nbin* is restricted to being an even
number.

As described below, the total kinetic energy transferred by these
swaps is computed by the fix and can be output.  Dividing this
quantity by time and the cross-sectional area of the simulation box
yields a heat flux.  The ratio of heat flux to the slope of the
temperature profile is proportional to the thermal conductivity of the
fluid, in appropriate units.  See the :ref:`Muller-Plathe paper <Muller-Plathe1>` for details.

.. note::

   If your system is periodic in the direction of the heat flux,
   then the flux is going in 2 directions.  This means the effective heat
   flux in one direction is reduced by a factor of 2.  You will see this
   in the equations for thermal conductivity (kappa) in the Muller-Plathe
   paper.  LAMMPS is simply tallying kinetic energy which does not
   account for whether or not your system is periodic; you must use the
   value appropriately to yield a kappa for your system.

.. note::

   After equilibration, if the temperature gradient you observe is
   not linear, then you are likely swapping energy too frequently and are
   not in a regime of linear response.  In this case you cannot
   accurately infer a thermal conductivity and should try increasing the
   Nevery parameter.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
kinetic energy transferred between the bottom and middle of the
simulation box (in the *edim* direction) is stored as a scalar
quantity by this fix.  This quantity is zeroed when the fix is defined
and accumulates thereafter, once every N steps.  The units of the
quantity are energy; see the :doc:`units <units>` command for details.
The scalar value calculated by this fix is "intensive".

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

LAMMPS does not check, but you should not use this fix to swap the
kinetic energy of atoms that are in constrained molecules, e.g. via
:doc:`fix shake <fix_shake>` or :doc:`fix rigid <fix_rigid>`.  This is
because application of the constraints will alter the amount of
transferred momentum.  You should, however, be able to use flexible
molecules.  See the :ref:`Zhang paper <Zhang2>` for a discussion and results
of this idea.

When running a simulation with large, massive particles or molecules
in a background solvent, you may want to only exchange kinetic energy
between solvent particles.

Related commands
""""""""""""""""

:doc:`fix ehex <fix_ehex>`, :doc:`fix heat <fix_heat>`, :doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix viscosity <fix_viscosity>`,
:doc:`compute heat/flux <compute_heat_flux>`

Default
"""""""

The option defaults are swap = 1.


----------


.. _Muller-Plathe1:



**(Muller-Plathe)** Muller-Plathe, J Chem Phys, 106, 6082 (1997).

.. _Zhang2:



**(Zhang)** Zhang, Lussetti, de Souza, Muller-Plathe, J Phys Chem B,
109, 15060-15067 (2005).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
