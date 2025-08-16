.. index:: fix rheo

fix rheo command
================

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
* keyword = *thermal* or *interface/reconstruct* or *surface/detection* or *shift* or *rho/sum* or *density* or *speed/sound*

  .. parsed-literal::

       *thermal* turns on thermal evolution
         values = none
       *interface/reconstruct* reconstructs interfaces with solid particles
         values = none
       *surface/detection* detects free-surfaces with an absence of particles
         values = *sdstyle* *limit* *limit/splash*
           *sdstyle* = *coordination* or *divergence*
           *limit* = threshold for surface particles
           *limit/splash* = threshold for splash particles (unitless)
       *shift* turns on velocity shifting
         values = none
         optional args = *exclude/type* or *scale/cross/type*
           *exclude/type* values = *types*
             *types* = list of types
           *scale/cross/type* values = *shiftscale* *cmin* *wmin*
             *shiftscale* = fraction of shifting in normal direction to preserve (unitless)
             *cmin* = minimum color function value required for scaling (unitless)
             *wmin* = minimum local same-type support required for any shifting (unitless)
       *rho/sum* density evolution performed by a kernel summation
         values = none
         optional args = *self/mass*
           *self/mass* values = none, a particle uses its own mass in summation
       *density* specify equilibrium densities for each atom type
         values = *rho01*, ... *rho0N* (density)
       *speed/sound* specify speeds of sound for each atom type
         values = *cs0*, ... *csN* (velocity)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo 3.0 quintic 0 thermal density 0.1 0.1 speed/sound 10.0 1.0
   fix 1 all rheo 3.0 RK1 10 shift surface/detection coordination 40
   fix 1 all rheo 3.0 RK1 10 shift exclude/type 2*4 scale/cross/type 0.05 0.02 0.5
   fix 1 all rheo 3.0 RK1 10 rhosum self/mass

Description
"""""""""""

.. versionadded:: 29Aug2024

Perform time integration for RHEO particles, updating positions, velocities,
and densities. For a detailed breakdown of the integration timestep and
numerical details, see :ref:`(Palermo) <fix_rheo_palermo>`. For an overview
and list of other features available in the RHEO package, see
:doc:`the RHEO howto <Howto_rheo>`.

The type of kernel is specified using *kstyle* and the cutoff is *cut*. Four
kernels are currently available. The *quintic* kernel is a standard quintic
spline function commonly used in SPH. The other options, *RK0*, *RK1*, and
*RK2*, are zeroth, first, and second order reproducing. To generate a
reproducing kernel, a particle must have sufficient neighbors inside the
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
:ref:`(Palermo) <fix_rheo_palermo>` except there is an additional
restriction that the reconstructed solid density cannot be less than the
equilibrium density. This prevents fluid particles from sticking to solid
surfaces.

A modified form of Fickian particle shifting can be enabled with the
*shift* keyword. This effectively shifts particle positions to generate a
more uniform spatial distribution. By default, shifting does not consider the
type of a particle and therefore may be inappropriate in systems consisting
of multiple atom types representing multiple fluid phases. However, two
optional sub-arguments can follow the *shift* keyword, *exclude/type* and
*scale/cross/type* to adjust shifting at fluid interfaces.

The *exclude/type* option lets the user specify a list of atom types which
are not shifted, *types*. A wild-card asterisk can be used in place
of or in conjunction with the *types* argument to toggle shifting for
multiple atom types.  This takes the form "\*" or "\*n" or "m\*"
or "m\*n".  If :math:`N` is the number of atom types, then an asterisk with
no numeric values means all types from 1 to :math:`N`.  A leading asterisk
means all types from 1 to n (inclusive).  A trailing asterisk means all types
from m to :math:`N` (inclusive).  A middle asterisk means all types from m to n
(inclusive).

The *scale/cross/type* option is designed to handle interfaces between fluids
made up of different atom types. Similar to the method by
:ref:`(Yang) <fix_rheo_yang>`, a color function is calculated and used to
estimate a local interfacial normal vector. Shifting along this normal direction
is rescaled by a factor of *scaleshift*, such that a value of *scaleshift* of
zero implies there is no shifting in the normal direction and a value of
*scaleshift* of one implies no change in behavior. This scaling is only applied
to atoms with a color function value greater than *cmin*. To handle scenarios
of a small inclusion of one fluid type (e.g. a single atom) inside another,
the degree of same-type support is calculated

.. math::
   W_{i,\mathrm{same}} = \sum_{j} W_{ij} \delta_{ij}

where :math:`\delta_{ij}` is zero if atoms :math:`i` and :math:`j` have different
types but unity otherwise. If :math:`W_{i,\mathrm{same}}` is ever less than the
specified value of *wmin*, shifting is turned off for particle :math:`i`

In systems with free surfaces (atom-vacuum), the *surface/detection* keyword
can classify the location of particles as being within the bulk fluid, on a
free surface, or isolated from other particles in a splash or droplet.
Shifting is then disabled in the normal direction away from the free surface
to prevent particles from diffusing away. Surface detection can also be used
to control surface-nucleated effects like oxidation when used in combination
with :doc:`fix rheo/oxidation <fix_rheo_oxidation>`. Surface detection is not
performed on solid bodies.

The *surface/detection* keyword takes three arguments: *sdstyle*, *limit*,
and *limit/splash*. The first, *sdstyle*, specifies whether surface particles
are identified using a coordination number (*coordination*) or the divergence
of the local particle positions (*divergence*). The threshold value for a
surface particle for either of these criteria is set by the numerical value
of *limit*. Additionally, if a particle's coordination number is too low,
i.e. if it has separated off from the bulk in a droplet, it is not possible
to define surfaces and the particle is classified as a splash. The coordination
threshold for this classification is set by the numerical value of
*limit/splash*.

By default, RHEO integrates particles' densities using a mass diffusion
equation. Alternatively, one can update densities every timestep by performing
a kernel summation of the masses of neighboring particles by specifying the *rho/sum*
keyword. Following this keyword, one may include the optional *self/mass* sub-argument
which modifies the behavior of the density summation. Typically, the density
:math:`\rho` of a particle is calculated as the sum over neighbors

.. math::
   \rho_i = \sum_{j} W_{ij} M_j

where :math:`W_{ij}` is the kernel, and :math:`M_j` is the mass of particle :math:`j`.
The *self/mass* keyword augments this expression by replacing :math:`M_j` with
:math:`M_i`. This may be useful in simulations of multiple fluid phases with large
differences in density, :ref:`(Hu) <fix_rheo_hu>`.

The *density* keyword is used to specify the equilibrium density of each of the N
particle types. It must be followed by N numerical values specifying each type's
equilibrium density *rho0*.

The *speed/sound* keyword is used to specify the speed of sound of each of the
N particle types. It must be followed by N numerical values specifying each type's
speed of sound *cs*. These values may be ignored if the pressure equation of
state has a non-constant speed of sound, as discussed further in
:doc:`fix rheo/pressure <fix_rheo_pressure>`.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with atom style rheo or rheo/thermal. This fix must
be used in conjunction with :doc:`fix rheo/pressure <fix_rheo_pressure>`.
and :doc:`fix rheo/viscosity <fix_rheo_viscosity>`. If the *thermal* setting
is used, there must also be an instance of
:doc:`fix rheo/thermal <fix_rheo_thermal>`. The fix group must be set to all.
Only one instance of fix rheo may be defined and it  must be defined prior
to all other RHEO fixes in the input script.

This fix is part of the RHEO package.  It is only enabled if LAMMPS was built
with that package. See the :doc:`Build package <Build_package>` page for more info.

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

.. _fix_rheo_palermo:

**(Palermo)** Palermo, Wolf, Clemmer, O'Connor, Phys. Fluids, 36, 113337 (2024).

.. _fix_rheo_yang:

**(Yang)** Yang, Rakhsha, Hu, Negrut, J. Comp. Physics, 458, 111079 (2022).

.. _fix_rheo_hu:

**(Hu)** Hu, and Adams, J. Comp. Physics, 213, 844-861 (2006).
