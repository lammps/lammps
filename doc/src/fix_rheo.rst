.. index:: fix rheo

fix rheo command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo cut kstyle zmin keyword values...

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo = style name of this fix command
* cut = cutoff for the kernel (distance)
* kstyle = *quintic* or *RK0* or *RK1* or *RK2*
* zmin = minimal number of neighbors for reproducing kernels
* zero or more keyword/value pairs may be appended to args
* keyword = *thermal* or *interface/reconstruct* or *surface/detection* or
            *shift* or *rho/sum* or *density* or *speed/sound*

  .. parsed-literal::

       *thermal* values = none, turns on thermal evolution
       *interface/reconstruct* values = none, reconstructs interfaces with solid particles
       *surface/detection* values = *sdstyle* *limit* *limit/splash*
         *sdstyle* = *coordination* or *divergence*
         *limit* = threshold for surface particles
         *limit/splash* = threshold for splash particles
       *shift* values = none, turns on velocity shifting
       *rho/sum* values = none, uses the kernel to compute the density of particles
       *density* values = *rho01*, ... *rho0N* (density)
       *speed/sound* values = *cs0*, ... *csN* (velocity)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo 1.0 quintic thermal density 0.1 speed/sound 10.0
   fix 1 all rheo 1.0 RK1 shift surface/detection coordination 40

Description
"""""""""""

.. versionadded:: TBD

Perform time integration for RHEO particles, updating positions, velocities,
and densities. For a detailed breakdown of the integration timestep and
numerical details, see :ref:`(Palermo) <howto_rheo_palermo>`. For an
overview of other features available in the RHEO package, see
:doc:`the RHEO howto <Howto_rheo>`.

The type of kernel is specified using *kstyle* and the cutoff is *cut*. Four
kernels are currently available. The *quintic* kernel is a standard quintic
spline function commonly used in SPH. The other options, *RK0*, *RK1*, and
*RK2*, are zeroth, first, and second order reproducing. To generate a reproducing kernel, a particle must have sufficient neighbors inside the
kernel cutoff distance (a coordination number) to accurately calculate
moments. This threshold is set by *zmin*. If reproducing kernels are
requested but a particle has fewer neighbors, then it will revert to a
non-reproducing quintic kernel until it gains more neighbors.

To model temperature evolution, one must specify the *thermal* keyword,
define a separate instance of :doc:`fix rheo/thermal <fix_rheo_thermal>`,
and use atom style rheo/thermal.

By default, the density of solid RHEO particles does not evolve and forces
with fluid particles are calculated using the current velocity of the solid
particle. If the *interface/reconstruct* keyword is used, then the density
and velocity of solid particles are alternatively reconstructed for every
fluid-solid interaction to ensure no-slip and pressure-balanced boundaries.
This is done by estimating the location of the fluid-solid interface and
extrapolating fluid particle properties across the interface to calculate a
temporary apparent density and velocity for a solid particle. The numerical
details are the same as those described in
:ref:`(Palermo) <howto_rheo_palermo>` except there is an additional
restriction that the reconstructed solid density cannot be less than the
equilibrium density. This prevents fluid particles from sticking to solid
surfaces.

A modified form of Fickian particle shifting can be enabled with the
*shift* keyword. This effectively shifts particle positions to generate a
more uniform spatial distribution. In systems with free surfaces, the
*surface/detection* keyword can be used to classify the location of
particles as being within the bulk fluid, on a free surface, or isolated
from other particles in a splash or droplet. Shifting is then disabled in
the direction away from the free surface to prevent it from diffusing
particles away from the bulk fluid. Surface detection can also be used
to control surface-nucleated effects like oxidation when used in combination
with :doc:`fix rheo/oxidation <fix_rheo_oxidation>`.

The *surface/detection* keyword takes three arguments: *sdstyle*, *limit*,
and *limi/splash*. The first, *sdstyle*, specifies whether surface particles
are identified using a coordination number (*coordination*) or the divergence
of the local particle positions (*divergence*). The threshold value for a
surface particle for either of these criteria is set by the numerical value
of *limit*. Additionally, if a particle's coordination number is too low,
i.e. if it has separated off from the bulk in a droplet, it is not possible
to define surfaces and a particle is classified as a splash. The coordination
threshold for this classification is set by the numerical value of
*limit/splash*.

By default, RHEO integrates particles' densities using a mass diffusion
equation. Alternatively, one can update densities every timestep by performing
a kernel summation of the masses of neighboring particles by specifying the *rho/sum* keyword.

The *density* is used to specify the equilbrium density of each of the N
particle types. It must be followed by N numerical values specifying each
type's equilibrium density *rho0*.

The *density* is used to specify the speed of sound of each of the N particle
types. It must be followed by N numerical values specifying each type's speed
of sound *cs*.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with atom style rheo or rheo/thermal.
This fix must be used in conjuction with
:doc:`fix rheo/pressure <fix_rheo_pressure>`. and
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`, If the *thermal*
setting is used, there must also be an instance of
:doc:`fix rheo/thermal <fix_rheo_thermal>`. The fix group must be
set to all. Only one instance of fix rheo may be defined and it
must be defined prior to all other RHEO fixes.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

*rho0* and *cs* are set to 1.0 for all atom types.

----------

.. _howto-howto_rheo_palermo:

**(Palermo)** Palermo, Clemmer, Wolf, O'Connor, in preparation.
