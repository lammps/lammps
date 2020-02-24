.. index:: pair\_style lubricateU

pair\_style lubricateU command
==============================

pair\_style lubricateU/poly command
===================================

Syntax
""""""


.. parsed-literal::

   pair_style style mu flaglog cutinner cutoff gdot flagHI flagVF

* style = *lubricateU* or *lubricateU/poly*
* mu = dynamic viscosity (dynamic viscosity units)
* flaglog = 0/1 to exclude/include log terms in the lubrication approximation
* cutinner = inner cut off distance (distance units)
* cutoff = outer cutoff for interactions (distance units)
* gdot = shear rate (1/time units)
* flagHI (optional) = 0/1 to exclude/include 1/r hydrodynamic interactions
* flagVF (optional) = 0/1 to exclude/include volume fraction corrections in the long-range isotropic terms

**Examples:** (all assume radius = 1)


.. parsed-literal::

   pair_style lubricateU 1.5 1 2.01 2.5 0.01 1 1
   pair_coeff 1 1 2.05 2.8
   pair_coeff \* \*

Description
"""""""""""

Styles *lubricateU* and *lubricateU/poly* compute velocities and
angular velocities for finite-size spherical particles such that the
hydrodynamic interaction balances the force and torque due to all
other types of interactions.

The interactions have 2 components.  The first is
Ball-Melrose lubrication terms via the formulas in :ref:`(Ball and Melrose) <Ball2>`

.. math::

   W & =  - a_{sq} | (v_1 - v_2) \bullet \mathbf{nn} |^2 - 
   a_{sh} | (\omega_1 + \omega_2) \bullet 
   (\mathbf{I} - \mathbf{nn}) - 2 \Omega_N |^2 - \\
   &  a_{pu} | (\omega_1 - \omega_2) \bullet (\mathbf{I} - \mathbf{nn}) |^2 -
   a_{tw} | (\omega_1 - \omega_2) \bullet \mathbf{nn} |^2  \qquad r < r_c \\
   & \\
   \Omega_N & = \mathbf{n} \times (v_1 - v_2) / r


which represents the dissipation W between two nearby particles due to
their relative velocities in the presence of a background solvent with
viscosity *mu*\ .  Note that this is dynamic viscosity which has units of
mass/distance/time, not kinematic viscosity.

The Asq (squeeze) term is the strongest and is included as long as
*flagHI* is set to 1 (default). It scales as 1/gap where gap is the
separation between the surfaces of the 2 particles. The Ash (shear)
and Apu (pump) terms are only included if *flaglog* is set to 1. They
are the next strongest interactions, and the only other singular
interaction, and scale as log(gap). Note that *flaglog* = 1 and
*flagHI* = 0 is invalid, and will result in a warning message, after
which *flagHI* will be set to 1. The Atw (twist) term is currently not
included. It is typically a very small contribution to the lubrication
forces.

The *flagHI* and *flagVF* settings are optional.  Neither should be
used, or both must be defined.

*Cutinner* sets the minimum center-to-center separation that will be
used in calculations irrespective of the actual separation.  *Cutoff*
is the maximum center-to-center separation at which an interaction is
computed.  Using a *cutoff* less than 3 radii is recommended if
*flaglog* is set to 1.

The other component is due to the Fast Lubrication Dynamics (FLD)
approximation, described in :ref:`(Kumar) <Kumar2>`.  The equation being
solved to balance the forces and torques is

.. math::

   -R_{FU}(U-U^{\infty}) = -R_{FE}E^{\infty} - F^{rest}


where U represents the velocities and angular velocities of the
particles, :math:`U^{\infty}` represents the velocities and the angular
velocities of the undisturbed fluid, and :math:`E^{\infty}` represents
the rate of strain tensor of the undisturbed fluid flow with viscosity
*mu*\ . Again, note that this is dynamic viscosity which has units of
mass/distance/time, not kinematic viscosity.  Volume fraction
corrections to R\_FU are included if *flagVF* is set to 1 (default).

F\ *rest* represents the forces and torques due to all other types of
interactions, e.g. Brownian, electrostatic etc.  Note that this
algorithm neglects the inertial terms, thereby removing the
restriction of resolving the small interial time scale, which may not
be of interest for colloidal particles.  This pair style solves for
the velocity such that the hydrodynamic force balances all other types
of forces, thereby resulting in a net zero force (zero inertia limit).
When defining this pair style, it must be defined last so that when
this style is invoked all other types of forces have already been
computed.  For the same reason, it won't work if additional non-pair
styles are defined (such as bond or Kspace forces) as they are
calculated in LAMMPS after the pairwise interactions have been
computed.

.. note::

   When using these styles, the these pair styles are designed to
   be used with implicit time integration and a correspondingly larger
   timestep.  Thus either :doc:`fix nve/noforce <fix_nve_noforce>` should
   be used for spherical particles defined via :doc:`atom_style sphere <atom_style>` or :doc:`fix nve/asphere/noforce <fix_nve_asphere_noforce>` should be used for
   spherical particles defined via :doc:`atom_style ellipsoid <atom_style>`.  This is because the velocity and angular
   momentum of each particle is set by the pair style, and should not be
   reset by the time integration fix.

Style *lubricateU* requires monodisperse spherical particles; style
*lubricateU/poly* allows for polydisperse spherical particles.

If the suspension is sheared via the :doc:`fix deform <fix_deform>`
command then the pair style uses the shear rate to adjust the
hydrodynamic interactions accordingly. Volume changes due to fix
deform are accounted for when computing the volume fraction
corrections to R\_FU.

When computing the volume fraction corrections to R\_FU, the presence
of walls (whether moving or stationary) will affect the volume
fraction available to colloidal particles. This is currently accounted
for with the following types of walls: :doc:`wall/lj93 <fix_wall>`,
:doc:`wall/lj126 <fix_wall>`, :doc:`wall/colloid <fix_wall>`, and
:doc:`wall/harmonic <fix_wall>`.  For these wall styles, the correct
volume fraction will be used when walls do not coincide with the box
boundary, as well as when walls move and thereby cause a change in the
volume fraction. To use these wall styles with pair\_style *lubricateU*
or *lubricateU/poly*\ , the *fld yes* option must be specified in the
fix wall command.

Since lubrication forces are dissipative, it is usually desirable to
thermostat the system at a constant temperature. If Brownian motion
(at a constant temperature) is desired, it can be set using the
:doc:`pair_style brownian <pair_brownian>` command. These pair styles
and the brownian style should use consistent parameters for *mu*\ ,
*flaglog*\ , *flagfld*\ , *cutinner*\ , *cutoff*\ , *flagHI* and *flagVF*\ .


----------


The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutinner (distance units)
* cutoff (distance units)

The two coefficients are optional.  If neither is specified, the two
cutoffs specified in the pair\_style command are used.  Otherwise both
must be specified.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the two cutoff distances for this
pair style can be mixed.  The default mix value is *geometric*\ .  See
the "pair\_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


These styles are part of the COLLOID package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Currently, these pair styles assume that all other types of
forces/torques on the particles have been already been computed when
it is invoked.  This requires this style to be defined as the last of
the pair styles, and that no fixes apply additional constraint forces.
One exception is the :doc:`fix wall/colloid <fix_wall>` commands, which
has an "fld" option to apply their wall forces correctly.

Only spherical monodisperse particles are allowed for pair\_style
lubricateU.

Only spherical particles are allowed for pair\_style lubricateU/poly.

For sheared suspensions, it is assumed that the shearing is done in
the xy plane, with x being the velocity direction and y being the
velocity-gradient direction. In this case, one must use :doc:`fix deform <fix_deform>` with the same rate of shear (erate).

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style lubricate <pair_lubricate>`

Default
"""""""

The default settings for the optional args are flagHI = 1 and flagVF =
1.


----------


.. _Ball2:



**(Ball)** Ball and Melrose, Physica A, 247, 444-472 (1997).

.. _Kumar2:



**(Kumar)** Kumar and Higdon, Phys Rev E, 82, 051401 (2010).
