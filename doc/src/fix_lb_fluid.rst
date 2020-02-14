.. index:: fix lb/fluid

fix lb/fluid command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID lb/fluid nevery LBtype viscosity density keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/fluid = style name of this fix command
* nevery = update the lattice-Boltzmann fluid every this many timesteps
* LBtype = 1 to use the standard finite difference LB integrator,
  2 to use the LB integrator of :ref:`Ollila et al. <Ollila>`
* viscosity = the fluid viscosity (units of mass/(time\*length)).
* density = the fluid density.
* zero or more keyword/value pairs may be appended
* keyword = *setArea* or *setGamma* or *scaleGamma* or *dx* or *dm* or *a0* or *noise* or *calcforce* or *trilinear* or *D3Q19* or *read\_restart* or *write\_restart* or *zwall\_velocity* or *bodyforce* or *printfluid*
  
  .. parsed-literal::
  
       *setArea* values = type node_area
           type = atom type (1-N)
           node_area = portion of the surface area of the composite object associated with the particular atom type (used when the force coupling constant is set by default).
       *setGamma* values = gamma
           gamma = user set value for the force coupling constant.
       *scaleGamma* values = type gammaFactor
           type = atom type (1-N)
           gammaFactor = factor to scale the *setGamma* gamma value by, for the specified atom type.
       *dx* values = dx_LB = the lattice spacing.
       *dm* values = dm_LB = the lattice-Boltzmann mass unit.
       *a0* values = a_0_real = the square of the speed of sound in the fluid.
       *noise* values = Temperature seed
           Temperature = fluid temperature.
           seed = random number generator seed (positive integer)
       *calcforce* values = N forcegroup-ID
           N = output the force and torque every N timesteps
           forcegroup-ID = ID of the particle group to calculate the force and torque of
       *trilinear* values = none (used to switch from the default Peskin interpolation stencil to the trilinear stencil).
       *D3Q19* values = none (used to switch from the default D3Q15, 15 velocity lattice, to the D3Q19, 19 velocity lattice).
       *read_restart* values = restart file = name of the restart file to use to restart a fluid run.
       *write_restart* values = N = write a restart file every N MD timesteps.
       *zwall_velocity* values = velocity_bottom velocity_top = velocities along the y-direction of the bottom and top walls (located at z=zmin and z=zmax).
       *bodyforce* values = bodyforcex bodyforcey bodyforcez = the x,y and z components of a constant body force added to the fluid.
       *printfluid* values = N = print the fluid density and velocity at each grid point every N timesteps.



Examples
""""""""


.. parsed-literal::

   fix 1 all lb/fluid 1 2 1.0 1.0 setGamma 13.0 dx 4.0 dm 10.0 calcforce sphere1
   fix 1 all lb/fluid 1 1 1.0 0.0009982071 setArea 1 1.144592082 dx 2.0 dm 0.3 trilinear noise 300.0 8979873

Description
"""""""""""

Implement a lattice-Boltzmann fluid on a uniform mesh covering the LAMMPS
simulation domain.  The MD particles described by *group-ID* apply a velocity
dependent force to the fluid.

The lattice-Boltzmann algorithm solves for the fluid motion governed by
the Navier Stokes equations,

.. math::
   
   \partial_t \rho + \partial_{\beta}\left(\rho u_{\beta}\right)= & 0 \\
   \partial_t\left(\rho u_{\alpha}\right) + \partial_{\beta}\left(\rho u_{\alpha} u_{\beta}\right) = & \partial_{\beta}\sigma_{\alpha \beta} + F_{\alpha} + \partial_{\beta}\left(\eta_{\alpha \beta \gamma \nu}\partial_{\gamma} u_{\nu}\right)


with,

.. math::

   \eta_{\alpha \beta \gamma \nu} = \eta\left[\delta_{\alpha \gamma}\delta_{\beta \nu} + \delta_{\alpha \nu}\delta_{\beta \gamma} - \frac{2}{3}\delta_{\alpha \beta}\delta_{\gamma \nu}\right] + \Lambda \delta_{\alpha \beta}\delta_{\gamma \nu}


where :math:`\rho` is the fluid density, *u* is the local
fluid velocity, :math:`\sigma` is the stress tensor, *F* is a local external
force, and :math:`\eta` and :math:`\Lambda` are the shear and bulk viscosities
respectively.  Here, we have implemented

.. math::

   \sigma_{\alpha \beta} = -P_{\alpha \beta} = -\rho a_0 \delta_{\alpha \beta}


with :math:`a_0` set to :math:`\frac{1}{3} \frac{dx}{dt}^2` by default.

The algorithm involves tracking the time evolution of a set of partial
distribution functions which evolve according to a velocity
discretized version of the Boltzmann equation,

.. math::

   \left(\partial_t + e_{i\alpha}\partial_{\alpha}\right)f_i = -\frac{1}{\tau}\left(f_i - f_i^{eq}\right) + W_i


where the first term on the right hand side represents a single time
relaxation towards the equilibrium distribution function, and :math:`\tau` is a
parameter physically related to the viscosity.  On a technical note,
we have implemented a 15 velocity model (D3Q15) as default; however,
the user can switch to a 19 velocity model (D3Q19) through the use of
the *D3Q19* keyword.  This fix provides the user with the choice of
two algorithms to solve this equation, through the specification of
the keyword *LBtype*\ .  If *LBtype* is set equal to 1, the standard
finite difference LB integrator is used.  If *LBtype* is set equal to
2, the algorithm of :ref:`Ollila et al. <Ollila>` is used.

Physical variables are then defined in terms of moments of the distribution
functions,

.. math::

   \rho = & \displaystyle\sum\limits_{i} f_i \\
   \rho u_{\alpha} = \displaystyle\sum\limits_{i} f_i e_{i\alpha}

Full details of the lattice-Boltzmann algorithm used can be found in
:ref:`Mackay et al. <fluid-Mackay>`.

The fluid is coupled to the MD particles described by *group-ID*
through a velocity dependent force.  The contribution to the fluid
force on a given lattice mesh site j due to MD particle alpha is
calculated as:

.. math::

   {\bf F}_{j \alpha} = \gamma \left({\bf v}_n - {\bf u}_f \right) \zeta_{j\alpha}


where :math:`\mathbf{v}_n` is the velocity of the MD particle,
:math:`\mathbf{u}_f` is the fluid
velocity interpolated to the particle location, and gamma is the force
coupling constant.  :math:`\zeta` is a weight assigned to the grid point,
obtained by distributing the particle to the nearest lattice sites.
For this, the user has the choice between a trilinear stencil, which
provides a support of 8 lattice sites, or the immersed boundary method
Peskin stencil, which provides a support of 64 lattice sites.  While
the Peskin stencil is seen to provide more stable results, the
trilinear stencil may be better suited for simulation of objects close
to walls, due to its smaller support.  Therefore, by default, the
Peskin stencil is used; however the user may switch to the trilinear
stencil by specifying the keyword, *trilinear*\ .

By default, the force coupling constant, :math:`\gamma`, is calculated
according to

.. math::

   \gamma = \frac{2m_um_v}{m_u+m_v}\left(\frac{1}{\Delta t_{collision}}\right)


Here, :math:`m_v` is the mass of the MD particle, :math:`m_u` is a
representative fluid mass at the particle location, and :math:`\Delta
t_{collision}` is a collision time, chosen such that
:math:`\frac{\tau}{\Delta t_{collision}} = 1` (see :ref:`Mackay and
Denniston <Mackay2>` for full details).  In order to calculate :math:`m_u`,
the fluid density is interpolated to the MD particle location, and
multiplied by a volume, node\_area\*dx\_lb, where node\_area
represents the portion of the surface area of the composite object
associated with a given MD particle.  By default, node\_area is set
equal to dx\_lb\*dx\_lb; however specific values for given atom types
can be set using the *setArea* keyword.

The user also has the option of specifying their own value for the
force coupling constant, for all the MD particles associated with the
fix, through the use of the *setGamma* keyword.  This may be useful
when modelling porous particles.  See :ref:`Mackay et al. <fluid-Mackay>` for a
detailed description of the method by which the user can choose an
appropriate gamma value.

.. note::

   while this fix applies the force of the particles on the fluid,
   it does not apply the force of the fluid to the particles.  When the
   force coupling constant is set using the default method, there is only
   one option to include this hydrodynamic force on the particles, and
   that is through the use of the :doc:`lb/viscous <fix_lb_viscous>` fix.
   This fix adds the hydrodynamic force to the total force acting on the
   particles, after which any of the built-in LAMMPS integrators can be
   used to integrate the particle motion.  However, if the user specifies
   their own value for the force coupling constant, as mentioned in
   :ref:`Mackay et al. <fluid-Mackay>`, the built-in LAMMPS integrators may prove to
   be unstable.  Therefore, we have included our own integrators :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>`, and :doc:`fix lb/pc <fix_lb_pc>`, to solve for the particle motion in these
   cases.  These integrators should not be used with the
   :doc:`lb/viscous <fix_lb_viscous>` fix, as they add hydrodynamic forces
   to the particles directly.  In addition, they can not be used if the
   force coupling constant has been set the default way.

.. note::

   if the force coupling constant is set using the default method,
   and the :doc:`lb/viscous <fix_lb_viscous>` fix is NOT used to add the
   hydrodynamic force to the total force acting on the particles, this
   physically corresponds to a situation in which an infinitely massive
   particle is moving through the fluid (since collisions between the
   particle and the fluid do not act to change the particle's velocity).
   Therefore, the user should set the mass of the particle to be
   significantly larger than the mass of the fluid at the particle
   location, in order to approximate an infinitely massive particle (see
   the dragforce test run for an example).


----------


Inside the fix, parameters are scaled by the lattice-Boltzmann
timestep, dt, grid spacing, dx, and mass unit, dm.  dt is set equal to
(nevery\*dt\_MD), where dt\_MD is the MD timestep.  By default, dm is set
equal to 1.0, and dx is chosen so that tau/(dt) =
(3\*eta\*dt)/(rho\*dx\^2) is approximately equal to 1.  However, the user
has the option of specifying their own values for dm, and dx, by using
the optional keywords *dm*\ , and *dx* respectively.

.. note::

   Care must be taken when choosing both a value for dx, and a
   simulation domain size.  This fix uses the same subdivision of the
   simulation domain among processors as the main LAMMPS program.  In
   order to uniformly cover the simulation domain with lattice sites, the
   lengths of the individual LAMMPS sub-domains must all be evenly
   divisible by dx.  If the simulation domain size is cubic, with equal
   lengths in all dimensions, and the default value for dx is used, this
   will automatically be satisfied.

Physical parameters describing the fluid are specified through
*viscosity*\ , *density*\ , and *a0*\ . If the force coupling constant is
set the default way, the surface area associated with the MD particles
is specified using the *setArea* keyword.  If the user chooses to
specify a value for the force coupling constant, this is set using the
*setGamma* keyword.  These parameters should all be given in terms of
the mass, distance, and time units chosen for the main LAMMPS run, as
they are scaled by the LB timestep, lattice spacing, and mass unit,
inside the fix.


----------


The *setArea* keyword allows the user to associate a surface area with
a given atom type.  For example if a spherical composite object of
radius R is represented as a spherical shell of N evenly distributed
MD particles, all of the same type, the surface area per particle
associated with that atom type should be set equal to 4\*pi\*R\^2/N.
This keyword should only be used if the force coupling constant,
gamma, is set the default way.

The *setGamma* keyword allows the user to specify their own value for
the force coupling constant, gamma, instead of using the default
value.

The *scaleGamma* keyword should be used in conjunction with the
*setGamma* keyword, when the user wishes to specify different gamma
values for different atom types.  This keyword allows the user to
scale the *setGamma* gamma value by a factor, gammaFactor, for a given
atom type.

The *dx* keyword allows the user to specify a value for the LB grid
spacing.

The *dm* keyword allows the user to specify the LB mass unit.

If the *a0* keyword is used, the value specified is used for the
square of the speed of sound in the fluid.  If this keyword is not
present, the speed of sound squared is set equal to (1/3)\*(dx/dt)\^2.
Setting a0 > (dx/dt)\^2 is not allowed, as this may lead to
instabilities.

If the *noise* keyword is used, followed by a positive temperature
value, and a positive integer random number seed, a thermal
lattice-Boltzmann algorithm is used.  If *LBtype* is set equal to 1
(i.e. the standard LB integrator is chosen), the thermal LB algorithm
of :ref:`Adhikari et al. <Adhikari>` is used; however if *LBtype* is set
equal to 2 both the LB integrator, and thermal LB algorithm described
in :ref:`Ollila et al. <Ollila>` are used.

If the *calcforce* keyword is used, both the fluid force and torque
acting on the specified particle group are printed to the screen every
N timesteps.

If the keyword *trilinear* is used, the trilinear stencil is used to
interpolate the particle nodes onto the fluid mesh.  By default, the
immersed boundary method, Peskin stencil is used.  Both of these
interpolation methods are described in :ref:`Mackay et al. <fluid-Mackay>`.

If the keyword *D3Q19* is used, the 19 velocity (D3Q19) lattice is
used by the lattice-Boltzmann algorithm.  By default, the 15 velocity
(D3Q15) lattice is used.

If the keyword *write\_restart* is used, followed by a positive
integer, N, a binary restart file is printed every N LB timesteps.
This restart file only contains information about the fluid.
Therefore, a LAMMPS restart file should also be written in order to
print out full details of the simulation.

.. note::

   When a large number of lattice grid points are used, the restart
   files may become quite large.

In order to restart the fluid portion of the simulation, the keyword
*read\_restart* is specified, followed by the name of the binary
lb\_fluid restart file to be used.

If the *zwall\_velocity* keyword is used y-velocities are assigned to
the lower and upper walls.  This keyword requires the presence of
walls in the z-direction.  This is set by assigning fixed boundary
conditions in the z-direction.  If fixed boundary conditions are
present in the z-direction, and this keyword is not used, the walls
are assumed to be stationary.

If the *bodyforce* keyword is used, a constant body force is added to
the fluid, defined by it's x, y and z components.

If the *printfluid* keyword is used, followed by a positive integer, N,
the fluid densities and velocities at each lattice site are printed to the
screen every N timesteps.


----------


For further details, as well as descriptions and results of several
test runs, see :ref:`Mackay et al. <fluid-Mackay>`.  Please include a citation to
this paper if the lb\_fluid fix is used in work contributing to
published research.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

Due to the large size of the fluid data, this fix writes it's own
binary restart files, if requested, independent of the main LAMMPS
:doc:`binary restart files <restart>`; no information about *lb\_fluid*
is written to the main LAMMPS :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the USER-LB package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix can only be used with an orthogonal simulation domain.

Walls have only been implemented in the z-direction.  Therefore, the
boundary conditions, as specified via the main LAMMPS boundary command
must be periodic for x and y, and either fixed or periodic for z.
Shrink-wrapped boundary conditions are not permitted with this fix.

This fix must be used before any of :doc:`fix lb/viscous <fix_lb_viscous>`, :doc:`fix lb/momentum <fix_lb_momentum>`, :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>`, and/ or :doc:`fix lb/pc <fix_lb_pc>` , as the fluid needs to be initialized before
any of these routines try to access its properties.  In addition, in
order for the hydrodynamic forces to be added to the particles, this
fix must be used in conjunction with the
:doc:`lb/viscous <fix_lb_viscous>` fix if the force coupling constant is
set by default, or either the :doc:`lb/viscous <fix_lb_viscous>` fix or
one of the :doc:`lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>` or
:doc:`lb/pc <fix_lb_pc>` integrators, if the user chooses to specify
their own value for the force coupling constant.

Related commands
""""""""""""""""

:doc:`fix lb/viscous <fix_lb_viscous>`, :doc:`fix lb/momentum <fix_lb_momentum>`, :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>`, :doc:`fix lb/pc <fix_lb_pc>`

Default
"""""""

By default, the force coupling constant is set according to

.. math::

   \gamma = \frac{2m_um_v}{m_u+m_v}\left(\frac{1}{\Delta t_{collision}}\right)


and an area of dx\_lb\^2 per node, used to calculate the fluid mass at
the particle node location, is assumed.

dx is chosen such that tau/(delta t\_LB) =
(3 eta dt\_LB)/(rho dx\_lb\^2) is approximately equal to 1.
dm is set equal to 1.0.
a0 is set equal to (1/3)\*(dx\_lb/dt\_lb)\^2.
The Peskin stencil is used as the default interpolation method.
The D3Q15 lattice is used for the lattice-Boltzmann algorithm.
If walls are present, they are assumed to be stationary.


----------


.. _Ollila:



**(Ollila et al.)** Ollila, S.T.T., Denniston, C., Karttunen, M., and Ala-Nissila, T., Fluctuating lattice-Boltzmann model for complex fluids, J. Chem. Phys. 134 (2011) 064902.

.. _fluid-Mackay:



**(Mackay et al.)** Mackay, F. E., Ollila, S.T.T., and Denniston, C., Hydrodynamic Forces Implemented into LAMMPS through a lattice-Boltzmann fluid, Computer Physics Communications 184 (2013) 2021-2031.

.. _Mackay2:



**(Mackay and Denniston)** Mackay, F. E., and Denniston, C., Coupling MD particles to a lattice-Boltzmann fluid through the use of conservative forces, J. Comput. Phys. 237 (2013) 289-298.

.. _Adhikari:



**(Adhikari et al.)** Adhikari, R., Stratford, K.,  Cates, M. E., and Wagner, A. J., Fluctuating lattice Boltzmann, Europhys. Lett. 71 (2005) 473-479.
