.. index:: fix lb/fluid

fix lb/fluid command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID lb/fluid nevery viscosity density keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/fluid = style name of this fix command
* nevery = update the lattice-Boltzmann fluid every this many timesteps (should normally be 1)
* viscosity = the fluid viscosity (units of mass/(time\*length)).
* density = the fluid density.
* zero or more keyword/value pairs may be appended
* keyword = *dx* or *dm* or *noise* or *stencil* or *read_restart* or *write_restart* or *zwall_velocity* or *pressurebcx* or *bodyforce* or *D3Q19* or *dumpxdmf* or *linearInit* or *dof* or *scaleGamma* or *a0* or *npits* or *wp* or *sw*

  .. parsed-literal::

       *dx* values = dx_LB = the lattice spacing.
       *dm* values = dm_LB = the lattice-Boltzmann mass unit.
       *noise* values = Temperature seed
           Temperature = fluid temperature.
           seed = random number generator seed (positive integer)
       *stencil* values = 2 (trilinear stencil, the default), 3 (3-point immersed boundary stencil), or 4 (4-point Keys' interpolation stencil)
       *read_restart* values = restart file = name of the restart file to use to restart a fluid run.
       *write_restart* values = N = write a restart file every N MD timesteps.
       *zwall_velocity* values = velocity_bottom velocity_top = velocities along the y-direction of the bottom and top walls (located at z=zmin and z=zmax).
       *pressurebcx* values = pgradav = imposes a pressure jump at the (periodic) x-boundary of pgradav*Lx*1000.
       *bodyforce* values = bodyforcex bodyforcey bodyforcez = the x,y and z components of a constant body force added to the fluid.
       *D3Q19* values = none (used to switch from the default D3Q15, 15 velocity lattice, to the D3Q19, 19 velocity lattice).
       *dumpxdmf* values = N file timeI
           N = output the force and torque every N timesteps
           file = output file name
           timeI = 1 (use simulation time to index xdmf file), 0 (use output frame number to index xdmf file)
       *linearInit* values = none = initialize density and velocity using linear interpolation (default is uniform density, no velocities)
       *dof* values = dof = specify the number of degrees of freedom for temperature calculation
       *scaleGamma* values = type gammaFactor
           type = atom type (1-N)
           gammaFactor = factor to scale the *setGamma* gamma value by, for the specified atom type.
       *a0* values = a_0_real = the square of the speed of sound in the fluid.
       *npits* values = npits h_p l_p l_pp l_e
           npits = number of pit regions
           h_p = z-height of pit regions (floor to bottom of slit)
           l_p = x-length of pit regions
           l_pp = x-length of slit regions between consecutive pits
           l_e = x-length of slit regions at ends
       *wp* values = w_p = y-width of slit regions (defaults to full width if not present or if sw active)
       *sw* values = none  (turns on y-sidewalls (in xz plane) if npits option active)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all lb/fluid 1 1.0 0.0009982071 dx 1.2 dm 0.001
   fix 1 all lb/fluid 1 1.0 0.0009982071 dx 1.2 dm 0.001 noise 300.0 2761
   fix 1 all lb/fluid 1 1.0 1.0 dx 4.0 dm 10.0 dumpxdmf 500 fflow 0 pressurebcx 0.01 npits 2 20 40 5 0 wp 30


Description
"""""""""""

.. versionchanged:: 24Mar2022

Implement a lattice-Boltzmann fluid on a uniform mesh covering the
LAMMPS simulation domain.  Note that this fix was updated in 2022 and is
not backward compatible with the previous version.  If you need the
previous version, please download an older version of LAMMPS.  The MD
particles described by *group-ID* apply a velocity dependent force to
the fluid.

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
You should not normally need to change this default.

The algorithm involves tracking the time evolution of a set of partial
distribution functions which evolve according to a velocity discretized
version of the Boltzmann equation,

.. math::

   \left(\partial_t + e_{i\alpha}\partial_{\alpha}\right)f_i = -\frac{1}{\tau}\left(f_i - f_i^{eq}\right) + W_i

where the first term on the right hand side represents a single time
relaxation towards the equilibrium distribution function, and
:math:`\tau` is a parameter physically related to the viscosity.  On a
technical note, we have implemented a 15 velocity model (D3Q15) as
default; however, the user can switch to a 19 velocity model (D3Q19)
through the use of the *D3Q19* keyword.  Physical variables are then
defined in terms of moments of the distribution functions,

.. math::

   \rho = & \displaystyle\sum\limits_{i} f_i \\
   \rho u_{\alpha} = & \displaystyle\sum\limits_{i} f_i e_{i\alpha}

Full details of the lattice-Boltzmann algorithm used can be found in
:ref:`Denniston et al. <fluid-Denniston>`.

The fluid is coupled to the MD particles described by *group-ID* through
a velocity dependent force.  The contribution to the fluid force on a
given lattice mesh site j due to MD particle :math:`\alpha` is
calculated as:

.. math::

   {\bf F}_{j \alpha} = \gamma \left({\bf v}_n - {\bf u}_f \right) \zeta_{j\alpha}

where :math:`\mathbf{v}_n` is the velocity of the MD particle,
:math:`\mathbf{u}_f` is the fluid velocity interpolated to the particle
location, and :math:`\gamma` is the force coupling constant.  This
force, as with most forces in LAMMPS, and hence the velocities, are
calculated at the half-time step. :math:`\zeta` is a weight assigned to
the grid point, obtained by distributing the particle to the nearest
lattice sites.

The force coupling constant, :math:`\gamma`, is calculated
according to

.. math::

   \gamma = \frac{2m_um_v}{m_u+m_v}\left(\frac{1}{\Delta t}\right)

Here, :math:`m_v` is the mass of the MD particle, :math:`m_u` is a
representative fluid mass at the particle location, and :math:`\Delta t`
is the time step.  The fluid mass :math:`m_u` that the MD particle
interacts with is calculated internally.  This coupling is chosen to
constrain the particle and associated fluid velocity to match at the end
of the time step.  As with other constraints, such as :doc:`shake
<fix_shake>`, this constraint can remove degrees of freedom from the
simulation which are accounted for internally in the algorithm.

.. note::

   While this fix applies the force of the particles on the fluid, it
   does not apply the force of the fluid to the particles.  There is
   only one option to include this hydrodynamic force on the particles,
   and that is through the use of the :doc:`lb/viscous <fix_lb_viscous>`
   fix.  This fix adds the hydrodynamic force to the total force acting
   on the particles, after which any of the built-in LAMMPS integrators
   can be used to integrate the particle motion.  If the
   :doc:`lb/viscous <fix_lb_viscous>` fix is NOT used to add the
   hydrodynamic force to the total force acting on the particles, this
   physically corresponds to a situation in which an infinitely massive
   particle is moving through the fluid (since collisions between the
   particle and the fluid do not act to change the particle's velocity).
   In this case, setting *scaleGamma* to -1 for the corresponding
   particle type will explicitly take this limit (of infinite particle
   mass) in computing the force coupling for the fluid force.

----------

Physical parameters describing the fluid are specified through
*viscosity* and *density*.  These parameters should all be given in
terms of the mass, distance, and time units chosen for the main LAMMPS
run, as they are scaled by the LB timestep, lattice spacing, and mass
unit, inside the fix.

The *dx* keyword allows the user to specify a value for the LB grid
spacing and the *dm* keyword allows the user to specify the LB mass
unit.  Inside the fix, parameters are scaled by the lattice-Boltzmann
timestep, :math:`dt_{LB}`, grid spacing, :math:`dx_{LB}`, and mass unit,
:math:`dm_{LB}`.  :math:`dt_{LB}` is set equal to
:math:`\mathrm{nevery}\cdot dt_{MD}`, where :math:`dt_{MD}` is the MD
timestep.  By default, :math:`dm_{LB}` is set equal to 1.0, and
:math:`dx_{LB}` is chosen so that :math:`\frac{\tau}{dt} = \frac{3\eta
dt}{\rho dx^2}` is approximately equal to 1.

 .. note::

   Care must be taken when choosing both a value for :math:`dx_{LB}`,
   and a simulation domain size.  This fix uses the same subdivision of
   the simulation domain among processors as the main LAMMPS program.  In
   order to uniformly cover the simulation domain with lattice sites, the
   lengths of the individual LAMMPS subdomains must all be evenly
   divisible by :math:`dx_{LB}`.  If the simulation domain size is cubic,
   with equal lengths in all dimensions, and the default value for
   :math:`dx_{LB}` is used, this will automatically be satisfied.

If the *noise* keyword is used, followed by a positive temperature
value, and a positive integer random number seed, the thermal LB algorithm
of :ref:`Adhikari et al. <Adhikari>` is used.

If the keyword *stencil* is used, the value sets the number of
interpolation points used in each direction.  For this, the user has the
choice between a trilinear stencil (*stencil* 2), which provides a
support of 8 lattice sites, or the 3-point immersed boundary method
stencil (*stencil* 3), which provides a support of 27 lattice sites, or
the 4-point Keys' interpolation stencil (stencil 4), which provides a
support of 64 lattice sites.  The trilinear stencil is the default as it
is better suited for simulation of objects close to walls or other
objects, due to its smaller support.  The 3-point stencil provides
smoother motion of the lattice and is suitable for particles not likely
to be to close to walls or other objects.

If the keyword *write_restart* is used, followed by a positive integer,
N, a binary restart file is printed every N LB timesteps.  This restart
file only contains information about the fluid.  Therefore, a LAMMPS
restart file should also be written in order to print out full details
of the simulation.

.. note::

   When a large number of lattice grid points are used, the restart
   files may become quite large.

In order to restart the fluid portion of the simulation, the keyword
*read_restart* is specified, followed by the name of the binary
lb_fluid restart file to be used.

If the *zwall_velocity* keyword is used y-velocities are assigned to
the lower and upper walls.  This keyword requires the presence of
walls in the z-direction.  This is set by assigning fixed boundary
conditions in the z-direction.  If fixed boundary conditions are
present in the z-direction, and this keyword is not used, the walls
are assumed to be stationary.

If the *pressurebcx* keyword is used, a pressure jump (implemented by a
step jump in density) is imposed at the (periodic) x-boundary.  The
value set specifies what would be the resulting equilibrium average
pressure gradient in the x-direction if the system had a constant
cross-section (i.e. resistance to flow).  It is converted to a pressure
jump by multiplication by the system size in the x-direction.  As this
value should normally be quite small, it is also assumed to be scaled
by 1000.

If the *bodyforce* keyword is used, a constant body force is added to
the fluid, defined by it's x, y and z components.

If the keyword *D3Q19* is used, the 19 velocity (D3Q19) lattice is
used by the lattice-Boltzmann algorithm.  By default, the 15 velocity
(D3Q15) lattice is used.

If the *dumpxdmf* keyword is used, followed by a positive integer, N,
and a file name, the fluid densities and velocities at each lattice site
are output to an xdmf file every N timesteps.  This is a binary file
format that can be read by visualization packages such as `Paraview
<https://www.paraview.org/>`_ .  The xdmf file format contains a time
index for each frame dump and the value timeI = 1 uses simulation time
while 0 uses the output frame number to index xdmf file.  The later can
be useful if the :doc:`dump vtk <dump_vtk>` command is used to output
the particle positions at the same timesteps and you want to visualize
both the fluid and particle data together in `Paraview
<https://www.paraview.org/>`_ .

The *scaleGamma* keyword allows the user to scale the :math:`\gamma`
value by a factor, gammaFactor, for a given atom type.  Setting
*scaleGamma* to -1 for the corresponding particle type will explicitly
take the limit of infinite particle mass in computing the force coupling
for the fluid force (see note above).

If the *a0* keyword is used, the value specified is used for the square
of the speed of sound in the fluid.  If this keyword is not present, the
speed of sound squared is set equal to
:math:`\frac{1}{3}\left(\frac{dx_{LB}}{dt_{LB}}\right)^2`.  Setting
:math:`a0 > (\frac{dx_{LB}}{dt_{LB}})^2` is not allowed, as this may
lead to instabilities.  As the speed of sound should usually be much
larger than any fluid velocity of interest, its value does not normally
have a significant impact on the results.  As such, it is usually best
to use the default for this option.

The *npits* keyword (followed by integer arguments: npits, h_p, l_p,
l_pp, l_e) sets the fluid domain to the pits geometry.  These arguments
should only be used if you actually want something more complex than a
rectangular/cubic geometry.  The npits value sets the number of pits
regions (arranged along x).  The remaining arguments are sizes measured
in multiples of dx_lb: h_p is the z-height of the pit regions, l_p is
the x-length of the pit regions, l_pp is the length of the region
between consecutive pits (referred to as a "slit" region), and l_e is
the x-length of the slit regions at each end of the channel.  The pit
geometry must fill the system in the x-direction but can be longer, in
which case it is truncated (which enables asymmetric entrance/exit end
sections).  The additional *wp* keyword allows the width (in
y-direction) of the pit to be specified (the default is full width) and
the *sw* keyword indicates that there should be sidewalls in the
y-direction (default is periodic in y-direction).  These parameters are
illustrated below::

    Sideview (in xz plane) of pit geometry:
    ______________________________________________________________________
      slit                          slit                          slit     ^
                                                                           |
    <---le---><---------lp-------><---lpp---><-------lp--------><---le---> hs = (Nbz-1) - hp
                                                                           |
    __________                    __________                    __________ v
              |                  |          |                  |           ^       z
              |                  |          |                  |           |       |
              |       pit        |          |       pit        |           hp      +-x
              |                  |          |                  |           |
              |__________________|          |__________________|           v

    Endview (in yz plane) of pit geometry (no sw so wp is active):
    _____________________
                          ^
                          |
                          hs
                          |
    _____________________ v
        |          |      ^
        |          |      |          z
        |<---wp--->|      hp         |
        |          |      |          +-y
        |__________|      v


----------

For further details, as well as descriptions and results of several test
runs, see :ref:`Denniston et al. <fluid-Denniston>`.  Please include a
citation to this paper if the lb_fluid fix is used in work contributing
to published research.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Due to the large size of the fluid data, this fix writes it's own
binary restart files, if requested, independent of the main LAMMPS
:doc:`binary restart files <restart>`; no information about *lb_fluid*
is written to the main LAMMPS :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

The fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the current
temperature of the group of particles described by *group-ID* along with
the fluid constrained to move with them. The temperature is computed via
the kinetic energy of the group and fluid constrained to move with them
and the total number of degrees of freedom (calculated internally).  If
the particles are not integrated independently (such as via :doc:`fix
NVE <fix_nve>`) but have additional constraints imposed on them (such as
via integration using :doc:`fix rigid <fix_rigid>`) the degrees of
freedom removed from these additional constraints will not be properly
accounted for.  In this case, the user can specify the total degrees of
freedom independently using the *dof* keyword.

The fix also computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.  There are 5 entries in
the array.  The first entry is the temperature of the fluid, the second
entry is the total mass of the fluid plus particles, the third through
fifth entries give the x, y, and z total momentum of the fluid plus
particles.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the LATBOLTZ package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This fix can only be used with an orthogonal simulation domain.

The boundary conditions for the fluid are specified independently to the
particles.  However, these should normally be specified consistently via
the main LAMMPS :doc:`boundary <boundary>` command (p p p, p p f, and p
f f are the only consistent possibilities).  Shrink-wrapped boundary
conditions are not permitted with this fix.

This fix must be used before any of :doc:`fix lb/viscous
<fix_lb_viscous>` and :doc:`fix lb/momentum <fix_lb_momentum>` as the
fluid needs to be initialized before any of these routines try to access
its properties.  In addition, in order for the hydrodynamic forces to be
added to the particles, this fix must be used in conjunction with the
:doc:`lb/viscous <fix_lb_viscous>` fix.

This fix needs to be used in conjunction with a standard LAMMPS
integrator such as :doc:`fix NVE <fix_nve>` or :doc:`fix rigid
<fix_rigid>`.

Related commands
""""""""""""""""

:doc:`fix lb/viscous <fix_lb_viscous>`, :doc:`fix lb/momentum <fix_lb_momentum>`

Default
"""""""

*dx* is chosen such that :math:`\frac{\tau}{dt_{LB}} = \frac{3\eta dt_{LB}}{\rho dx_{LB}^2}` is approximately equal to 1.
*dm* is set equal to 1.0.
*a0* is set equal to :math:`\frac{1}{3}\left(\frac{dx_{LB}}{dt_{LB}}\right)^2`.
The trilinear stencil is used as the default interpolation method.
The D3Q15 lattice is used for the lattice-Boltzmann algorithm.

----------

.. _fluid-Denniston:

**(Denniston et al.)** Denniston, C., Afrasiabian, N., Cole-Andre, M.G., Mackay, F. E., Ollila, S.T.T., and Whitehead, T., LAMMPS lb/fluid fix version 2: Improved Hydrodynamic Forces Implemented into LAMMPS through a lattice-Boltzmann fluid, Computer Physics Communications 275 (2022) `108318 <https://doi.org/10.1016/j.cpc.2022.108318>`_ .

.. _Mackay2:

**(Mackay and Denniston)** Mackay, F. E., and Denniston, C., Coupling MD particles to a lattice-Boltzmann fluid through the use of conservative forces, J. Comput. Phys. 237 (2013) 289-298.

.. _Adhikari:

**(Adhikari et al.)** Adhikari, R., Stratford, K.,  Cates, M. E., and Wagner, A. J., Fluctuating lattice Boltzmann, Europhys. Lett. 71 (2005) 473-479.

