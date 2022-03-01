.. index:: kspace_tild 

:doc:`kspace_style tild <kspace_tild>` command
====================================================

:doc:`kspace_modify tild <kspace_tild>` command
=====================================================

Syntax
""""""

.. code-block:: LAMMPS

   kspace_modify keyword value ...

* one or more keyword/value pairs may be listed
* keyword = *collective* or *compute* or *fftbench* or *force* or *gridres* or *kmax/ewald* or *mesh* or *minorder* or *mix/disp* or *order/disp* or *order* or *overlap* or *slab* or *splittol*

  .. parsed-literal::

       *collective* value = *yes* or *no*
       *compute* value = *yes* or *no*
       *fftbench* value = *yes* or *no*
       *mesh* value = x y z
         x,y,z = grid size in each dimension for long-range Coulombics
       *minorder* value = M
         M = min allowed extent of Gaussian when auto-adjusting to minimize grid communication
       *order* value = N
         N = extent of Gaussian for PM mapping of density to grid
       *overlap* = *yes* or *no* = whether the grid stencil is allowed to overlap into more than the nearest-neighbor processor
       *pressure/scalar* value = *yes* or *no*
       *shape* values = itype itype SHAPE parameters
         SHAPE = *gaussian* or *erfc* or *none*
         parameters = parameters for each function (see discussion below)
       *prefactor* values = itype jtype prefactor
         prefactor = magnitude of the function 
       *cross-interaction* values = itype jtype interaction_shape parameters
         interaction_shape = *gaussian* or *gaussian-erfc* or *erfc*
         parameters = parameters for each function (see discussion below)
       *set_rho0* value = rho0
         rho0 = total density of TILD particles
       *subtract_rho0* value = *yes* or *no*
       *normalize_by_rho0* value = *yes* or *no*
       *write_grid_data* values = freq filename
         freq = frequency of output. 0 or less disables output
         filename = name of file to output time averages to
       *ave/grid* values = Nevery Nrepeat Nfreq filename 
         Nevery = use input values every this many timesteps
         Nrepeat = # of times to use input values for calculating averages
         Nfreq = calculate averages every this many timesteps
         filename = name of file to output time averages to

       

Examples
""""""""

.. code-block:: LAMMPS

   kspace_modify mesh 24 24 30 order 6
   kspace_modify slab 3.0
   kspace_modify scafacos tolerance energy

Description
"""""""""""

The *tild* style is an implementation of 'theoretically informed Langevin Dynamics' method (previously known as `Dynamical Mean Field Theory`) :ref:`Chao <Chao>`: :ref:`Fredrickson<Fredrickson>`: :ref:`Grzetic<Grzetic>`. This interaction potential uses a particle-mesh scheme to calculate non-bonded pairwise forces indirectly through a gridded density representation. *tild* assigns a potentials for each simulation particle type, defined with :doc:`kspace_modify tild <kspace_modify>`. This potential does NOT calculate any Coulombic interactions.  

.. note::

   Unlike other KSpace solvers in LAMMPS, the kspace TILD accounts for
   non-bonded interactions, both short-range and long-range interactions through
   a "short-ranged" potentital. Therefore, there is no accompanying short range
   pair-style required. To fully implement the TILD methodology, use :doc:`fix
   langevin<fix_langevin>` with *tild*. (There is no warning produced if TILD is used without `fix langevin`. 


Set parameters used by the kspace solvers defined by the
:doc:`kspace_style <kspace_style>` command.  Not all parameters are
relevant to all kspace styles.

----------

The *collective* keyword applies only to TILD.  It is set to *no* by
default, except on IBM BlueGene machines.  If this option is set to
*yes*, LAMMPS will use MPI collective operations to remap data for
3d-FFT operations instead of the default point-to-point communication.
This is faster on IBM BlueGene machines, and may also be faster on
other machines if they have an efficient implementation of MPI
collective operations and adequate hardware.

----------

The *compute* keyword allows Kspace computations to be turned off,
even though a :doc:`kspace_style <kspace_style>` is defined.  This is
not useful for running a real simulation, but can be useful for
debugging purposes or for computing only partial forces that do not
include the Kspace contribution.  You can also do this by simply not
defining a :doc:`kspace_style <kspace_style>`, but a Kspace-compatible
:doc:`pair_style <pair_style>` requires a kspace style to be defined.
This keyword gives you that option.

----------

The *fftbench* keyword applies only to PPPM. It is off by default. If
this option is turned on, LAMMPS will perform a short FFT benchmark
computation and report its timings, and will thus finish some seconds
later than it would if this option were off.

----------

The *gridres* keyword overrides the current grid resolution parameter set by
the :doc:`kspace_style tild <kspace_tild>` command with an absolute 
accuracy.  The resolution determines the mesh grid for the long-range solver.

----------

The *mesh* keyword sets the grid size for kspace style *tild*\ .
In the case of TILD, this is the FFT mesh, and each dimension
must be factorizable into powers of 2, 3, and 5.  When this option is 
not set, the TILD solver chooses its own grid size, consistent with the
user-specified accuracy and pairwise cutoff.  Values for x,y,z of
0,0,0 unset the option.

----------

The *minorder* keyword allows LAMMPS to reduce the *order* setting if
necessary to keep the communication of ghost grid point limited to
exchanges between nearest-neighbor processors.  See the discussion of
the *overlap* keyword for details.  If the *overlap* keyword is set to
*yes*, which is the default, this is never needed.  If it set to *no*
and overlap occurs, then LAMMPS will reduce the order setting, one
step at a time, until the ghost grid overlap only extends to nearest
neighbor processors.  The *minorder* keyword limits how small the
*order* setting can become.  The minimum allowed value for PPPM is 2,
which is the default.  If *minorder* is set to the same value as
*order* then no reduction is allowed, and LAMMPS will generate an
error if the grid communication is non-nearest-neighbor and *overlap*
is set to *no*\ . The *minorder* keyword is not currently supported in
MSM.

----------

The *order* keyword determines how many grid spacings an atom's charge
extends when it is mapped to the grid in kspace style *tild*\ .
The default for this parameter is 5 for TILD, which
means each charge spans 5 grid cells in each dimension,
respectively.  For TILD, the minimum allowed
setting is 2 and the maximum allowed setting is 7. Note that there is an
inherent trade-off involved: a small grid will lower the cost of FFTs, but a larger order parameter will increase the cost
of interpolating charge/fields to/from the grid.

----------

The *overlap* keyword can be used in conjunction with the *minorder*
keyword with the TILD styles to adjust the amount of communication
that occurs when values on the FFT grid are exchanged between
processors.  This communication is distinct from the communication
inherent in the parallel FFTs themselves, and is required because
processors interpolate charge and field values using grid point values
owned by neighboring processors (i.e. ghost point communication).  If
the *overlap* keyword is set to *yes* then this communication is
allowed to extend beyond nearest-neighbor processors, e.g. when using
lots of processors on a small problem.  If it is set to *no* then the
communication will be limited to nearest-neighbor processors and the
*order* setting will be reduced if necessary, as explained by the
*minorder* keyword discussion. The *overlap* keyword is always set to
*yes* in MSM.

----------

The *tild shape* keywords specifies the shape potential of a given molecule
type. This is used to automatically generate interaction potentials between
particles of different types. There are two currently supported types:
`gaussian` and `erfc`. A `none` type is supported particles that do not have a
corresponding shape function. For interactions between two Gaussian particles,
we analytically convolve the two shape potentials together; for all other
interactions, we do a numerical convolution to get the proper convolved
interactions. 

The current shpae function styles used in *tild shape* are
.. math::

   U_{g} = & \frac{A}{\rho_0 (2\pi \sigma^2)^{3/2}} \exp(-r^2/2\sigma^2) \\
         = & \frac{A}{\rho_0} u_G (r) \\
   U_{erfc} = & - \frac{A}{\rho_0} \text{erfc} \left(\frac{\vert r \vert - R_p}{\xi}\right) \\ 
   U_{g-erfc} = & \frac{A}{\rho_0} u_G (r) * \text{erfc}
   \left(\frac{\vert r \vert - R_p}{\xi}\right)

where :math:`A` is the value set by `tild prefactor`\, :math:`\rho_0` is the total density of the TILD particles, :math:`\sigma`\ is the gaussian width, :math:`R_p` is the erfc particle radius and :math:`xi` is the erfc width.

The first required keyword for the *tild shape* option is the model. 
Currently supported options for shape function models
and their required arguments are:

1. *gaussian* : :math:`\sigma` (distance units)
2. *erfc* : :math:`R_p`, :math:`\xi` (both in distance units)

----------

The *tild prefactor* keyword sets the prefactor in front of a given shape. For
typical polymer represented by Gaussian monomers, the prefactors represents the
Flory-Higgins prefactor :math:`\chi` \ . See the :math:`A` prefactors in the
*tild shape* potentials.

----------

The *tild set_rho0* keyword is used when particles with a `tild shape` of `erfc`
exist within the simulation box and are used to ensure that the overall TILD
density of the box is the same as the user's input. Please note if the box
contains only `gaussian` shapes, this has no effect on the simulation. 

----------

The *tild normalize_by_rho0* keyword will divide the interactions by the
calcualted TILD :math:`\rho_0`\, the total density of the TILD particles. Please note this division will divide the
prefactors specified in `tild prefactor`\ .

----------

The *tild cross-interaction* keyword is used to override any specified interaction
from `tild shape`. At this time, we currently only support three non-zero
interaction styles (`gaussian`, `erfc`, `gaussian-erfc`), which model the
interactions between two gaussian potentials, two erfc potentials, or the
interaction between a gaussian particle and an erfc particle. There is also a
`none` style to force no-interactions between certain particle types and also a
`delete` command to remove any previously entered `tild cross-interaction`\ .

The current interaction styles used in *tild cross-interaction* are

.. math::

   U_{g} = & \frac{A}{\rho_0 (2\pi \sigma^2)^{3/2}} \exp(-r^2/2\sigma^2) \\
         = & \frac{A}{\rho_0} u_G (r) \\
   U_{erfc} = & - \frac{A}{\rho_0} \text{erfc} \left(\frac{\vert r \vert - R_p}{\xi}\right) \\ 
   U_{g-erfc} = & \frac{A}{\rho_0} u_G (r) * \text{erfc}
   \left(\frac{\vert r \vert - R_p}{\xi}\right)

where :math:`A` is the value set by `tild prefactor`\ , :math:`\rho_0` is the total density of the TILD particles, :math:`\sigma` is the gaussian width, :math:`R_p` is the erfc particle radius and :math:`\xi` is the erfc width.

The first required keyword for the *tild cross-interaction* option is the interaction model. 
Currently supported options for interaction models
and their required arguments are:

1. *gaussian* : :math:`\sigma` (distance units)
2. *gaussian-erfc* : :math:`\sigma`\ , :math:`R_p`, :math:`\xi` (all in distance units)
3. *erfc* : :math:`R_p`\ , :math:`\xi` (both in distance units)

----------

The *write_grid_data* writes the instantaneous gridded density to *filename*. Every $freq$ timesteps, the density is overwritten.

----------

The *ave/grid* keywords determines how freuently the density grids are averaged and 
output. The *Nevery*, *Nrepeat*, and *Nfreq* arguments specify on what
timesteps the input values will be used in order to contribute to the average.
The final averaged quantities are generated on timesteps that are a multiple of
*Nfreq*. The average is over *Nrepeat* quantities, computed in the preceding
portion of the simulation every *Nevery* timesteps. *Nfreq* must be a multiple
of *Nevery* and *Nevery* must be non-zero even if *Nrepeat* is 1. Also, the
timesteps contributing to the average value cannot overlap, i.e. Nrepeat*Nevery
can not exceed Nfreq.


----------

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`kspace_style <kspace_style>`

Default
"""""""

The option defaults are mesh = 0 0 0, order = 5 (TILD), minorder = 2, overlap = yes, gewald = gewald/disp = 0.0, slab = 1.0, fftbench = no (PPPM), mix/disp = pair, force/disp/real = -1.0, force/disp/kspace = -1.0, split = 0, tol = 1.0e-6, tild mix = convolution, tild subtract_rho0 = yes, tild normalize_by_rho0 = yes and disp/auto = no. 

----------

.. _Chao:

**(Chao)** Chao, H., Koski, J. & Riggleman, R. (2017)
"Solvent vapor annealing in block copolymer nanocomposite films: 
a dynamic mean field approach" Soft Matter, 13(1) 239-249.

.. _Fredrickson:

**(Fredrickson)** Fredrickson, G. H. and Orland, H.  (2017)
"Dynamics of polymers: A mean-field theory" The Journal of Chemical Physics 
140, 084902 (2014) https://doi.org/10.1063/1.4865911

.. _Grzetic:

**(Grzetic)** Grzetic, D. J., Wickman, R. A., and Shi, A.-C., "Statistical
dynamics of classical systems: A self-consistent field approach", The Journal of
Chemical Physics 140, 244907 (2014) https://doi.org/10.1063/1.4884825
