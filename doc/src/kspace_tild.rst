.. index:: kspace_tild 

:doc:`kspace_style tild <kspace_tild>` command
====================================================

:doc:`kspace_modify tild <kspace_tild>` command
=====================================================

Syntax
""""""

.. code-block:: LAMMPS

   kspace_style tild gridsize
   kspace_modify keyword value ...

* one or more keyword/value pairs may be listed
* keyword = *gridsize* or *mesh* or *minorder* or *order* or *overlap* or *shape* or *prefactor* or *cross-interaction* or *set_rho0* or *subtract_rho0* or *normalize_by_rho0* or *write_grid_data* or *ave/grid*

  .. parsed-literal::

       *gridsize* value = gridsize
       *mesh* value = x y z
         x,y,z = grid size in each dimension for TILD interaction
       *minorder* value = M
         M = min allowed extent of mapping when auto-adjusting to minimize grid communication
       *order* value = N
         N = extent of mapping when auto-adjusting to minimize grid communication
       *overlap* = *yes* or *no* = whether the grid stencil is allowed to overlap into more than the nearest-neighbor processor
       *shape* values = itype itype shape parameters
         shape  = *gaussian* or *erfc* or *none*
         parameters = parameters for each function (see discussion below)
       *prefactor* values = itype jtype prefactor
         prefactor = magnitude of the function 
       *cross-interaction* values = itype jtype interaction_shape parameters
         interaction_shape = *gaussian* or *gaussian_erfc* or *erfc*
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

   kspace_style  tild 0.5
   kspace_modify tild set_rho0 4
   kspace_modify tild shape 2 2 gaussian 3 shape 3 3 erfc 3.5 0.25 3.0  shape 4 4 none
   kspace_modify tild cross-interaction 2 2 gaussian 18
   kspace_modify tild cross-interaction 3 1 gaussian_erfc 0.25 3.5 0.25
   kspace_modify tild cross-interaction 3 3 erfc 3.5 0.25

Description
"""""""""""

The *tild* style is an implementation of "theoretically informed Langevin Dynamics" method (previously known as "Dynamical Mean Field Theory") :ref:`Chao<Chao>`: :ref:`Fredrickson<Fredrickson>`: :ref:`Grzetic<Grzetic>`. This interaction potential uses a particle-mesh scheme to calculate non-bonded pairwise forces indirectly through a gridded density representation. *tild* assigns a potential for each simulation particle type, defined with :doc:`kspace_modify tild <kspace_modify>`. This potential does NOT calculate any Coulombic or magnetic interactions and currently is incompatible with other `kspace_style` in LAMMPS.

.. note::

   Unlike other KSpace solvers in LAMMPS, the kspace TILD accounts for
   non-bonded interactions, both short-range and long-range interactions through
   a "short-ranged" potential. Therefore, there is no accompanying short-range
   pair style required. To fully implement the TILD methodology, use 
   :doc:`fix langevin<fix_langevin>` with *tild*\ . 
   (There is no warning produced if TILD is used without `fix langevin`.) 


Set parameters used by the kspace solvers defined by the
:doc:`kspace_style <kspace_style>` command.  Not all parameters are
relevant to all kspace styles.

----------

The *gridsize* keyword overrides the current grid resolution parameter set by
the `kspace_style tild` command with a new size in distance units. 
The grid size determines the mesh grid for the long-range solver.

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
*order* setting can become.  The minimum allowed value for TILD is 2,
which is the default.  If *minorder* is set to the same value as
*order* then no reduction is allowed, and LAMMPS will generate an
error if the grid communication is non-nearest-neighbor and *overlap*
is set to *no*\ . 

----------

The *order* keyword determines how many grid spacings an atom's charge
extends when it is mapped to the grid in kspace style *tild*\ .
The default for this parameter is 5 for TILD, which
means each charge spans 5 grid cells in each dimension,
respectively.  For TILD, the minimum allowed
setting is 2 and the maximum allowed setting is 7. Note that there is an
inherent trade-off involved: a small grid will lower the cost of FFTs, 
but a larger order parameter will increase the cost
of interpolating particles/fields to/from the grid.

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
*minorder* keyword discussion. 

----------

The *tild shape* keywords specifies the shape potential of a given molecule
type. This is used to automatically generate interaction potentials between
particles of different types. There are two currently supported types:
`gaussian` and `erfc`. A `none` type is supported particles that do not have a
corresponding shape function. For interactions between two Gaussian particles,
we analytically convolve the two shape potentials together; for all other
interactions, we do a numerical convolution to get the proper convolved
interactions. Therefore, it does not make sense to have a shape defined between
two particles using this keyword; one should instead use the `cross-interaction`
keyword if one wishes to specify the cross-interaction. The code will *NOT*
error out if you use two different types. The input keeps two types as its input
to maintain a familiar interface as the `pair_style` keyword. 

The current shape function styles used in *tild shape* are

.. math::

   U_{g} = & \frac{A}{\rho_0 (2\pi \sigma^2)^{3/2}} \exp(-r^2/2\sigma^2) \\
         = & \frac{A}{\rho_0} u_G (r) \\
   U_{erfc} = & \frac{A}{2\rho_0}  \rho_{NP} \text{erfc} \left(\frac{\vert r \vert - R_p}{\xi}\right) \\ 

where :math:`A` is the value set by `tild prefactor`\, :math:`\rho_0` is the total density of the TILD particles, :math:`\rho_{NP}` is the density of the TILD erfc nanoparticle, :math:`\sigma`\ is the gaussian width, :math:`R_p` is the erfc particle radius and :math:`\xi` is the erfc width.

The first required keyword for the *tild shape* option is the model. 
Currently supported options for shape function models
and their required arguments are:

1. *gaussian* : :math:`\sigma` (distance units)
2. *erfc* : :math:`R_p`, :math:`\xi` and :math:`\rho_{NP}` (:math:`R_p` and :math:`\xi` in distance units; :math:`\rho_{NP}` in inverse volume units)

----------

The *tild prefactor* keyword sets the prefactor in front of a given shape. For
typical polymer represented by Gaussian monomers, the prefactors represents the
Flory-Higgins prefactor :math:`\chi`. See the :math:`A` prefactors in the
*tild shape* potentials.

----------

The *tild set_rho0* keyword is used to set the TILD density which is calculated
separately from any other density in LAMMPS. Each defined `gaussian` shape
The *tild set_rho0* keyword is used to set the TILD number density which is calculated
 separately from any other density in LAMMPS. Each defined `gaussian` shape
 particle contributes 1, while each defined `erfc` shape has a density of 
:math:`4/3 \pi r^3 \rho_{NP}`\ . Particles without any defined shape functions do not contribute to the
overall density, even if they are included in a `cross-interaction`. 
Defining a *rho0* for a system without any shape functions (but with `cross-interaction` functions) will
accept the value as is (provided it is non-negative) and use that for
normalization purposes. Similarly, a function consisting of whose only defined
shapes are purely `gaussian` will also accept the user specified *rho0* as is.
For simulations with shape defined `erfc` particles, the *rho0* of all the
nanoparticles will be adjusted so that the overall density of the system matches
the user specified density. 

----------

The *tild normalize_by_rho0* keyword will divide the interactions by the
calculated TILD :math:`\rho_0`\, the total box density of the TILD particles. 
Please note this division will divide the prefactors specified in `tild prefactor`\ .

----------

The *tild cross-interaction* keyword is used to override any specified interaction
from `tild shape`. At this time, we currently only support three non-zero
interaction styles (`gaussian`, `erfc` and `gaussian_erfc`), which model the
interactions between two gaussian potentials, two erfc potentials, or the
interaction between a gaussian particle and an erfc particle. There is also a
`none` style to force no-interactions between certain particle types and also a
`delete` command to remove any previously entered `tild cross-interaction`\ .

The current interaction styles used in *tild cross-interaction* are

.. math::

   U_{g} = & \frac{A\exp(-r^2/2\sigma^2)}{\rho_0 (2\pi \sigma^2)^{3/2}}  \\
         = & \frac{A u_G (r)}{\rho_0} \\
   U_{erfc} = & \frac{A}{\rho_0} \text{erfc} \left(\frac{\vert r \vert - R_p}{\xi}\right) \\ 
   U_{g-erfc} = & \frac{A}{\rho_0} u_G (r) * \text{erfc}
   \left(\frac{\vert r \vert - R_p}{\xi}\right)

where :math:`A` is the value set by `tild prefactor`\ , :math:`\rho_0` is the TILD density of the simulation box, :math:`\sigma` is the gaussian width, :math:`R_p` is the erfc particle radius and :math:`\xi` is the erfc width, which controls how quickly the particle density drops from :math:`\rho_0`` to zero. :math:`U_{g-erfc}` involves convoluting the :math:`U_{g}` and :math:`U_{erfc}` functions.

The first required keyword for the *tild cross-interaction* option is the interaction model. 
Currently supported options for interaction models
and their required arguments are:

1. *gaussian* : :math:`\sigma` (distance units)
2. *gaussian_erfc* : :math:`\sigma`\ , :math:`R_p`, :math:`\xi` (all in distance units)
1. *gaussian* : :math:`\sigma^2` (distance:math:`^2` units)
 2. *gaussian_erfc* : :math:`\sigma^2`\ , :math:`R_p`, :math:`\xi` (distance units except for :math:`\sigma^2` )
 3. *erfc* : :math:`R_p`\ , :math:`\xi` (both in distance units)

.. note::

   ``cross-interaction`` and `shape` definitions have slightly different input parameters and so mapping is explicitly laid out.
   For the ``gaussian`` `shape`, the input parameter is :math:`\sigma_{i}`\ ; the code will square this automatically. 
   For interactions between two ``gaussian`` defined `shape` particles, the code analytically and behind the scenes performs the convolution so that the interaction potential uses :math:`\sigma^2_{12} = \sigma^2_{1} + \sigma^2_{2}`. For the convolution between a ``gaussian`` `shape` and a ``erfc`` `shape`, the code convolves the ``gaussian`` and ``erfc`` `shape` potentials computationally; this is also true for interactions between two ``erfc`` `shape` particles. 

   However, for ``cross-interaction``, the code treats the user input for ``gaussian`` as :math:`\sigma^2` so the user should manually calculate their own :math:`\sigma_{12}^2` before the run. For ``gaussian_erfc``, the code takes in :math:`\sigma^2` instead of :math:`\sigma`. Additionally, the ``gaussian_erfc`` and ``erfc`` commands do not take into account the :math:`\rho_{NP}` since ``cross-interactions`` assume to know nothing about this. Thus, if you have a value of :math:`\rho_{NP}` that is not 1, you should multiply it by :math:`\rho_{NP}` or :math:`(\rho_{NP})^2` for ``gaussian_erfc`` and ``erfc``, respectively.
   
   Identical simulations defined both ways can be found in examples/tild.

----------

The *write_grid_data* writes the instantaneous gridded density to *filename*\ . Every *freq* timesteps, the density is overwritten.

----------

The *ave/grid* keywords determines how frequently the density grids are averaged and 
output. The *Nevery*, *Nrepeat*, and *Nfreq* arguments specify on what
timesteps the input values will be used in order to contribute to the average.
The final averaged quantities are generated on timesteps that are a multiple of
*Nfreq*\ . The average is over *Nrepeat* quantities, computed in the preceding
portion of the simulation every *Nevery* timesteps. *Nfreq* must be a multiple
of *Nevery* and *Nevery* must be non-zero even if *Nrepeat* is 1. Also, the
timesteps contributing to the average value cannot overlap, i.e. *Nrepeat* * *Nevery*
can not exceed *Nfreq*\ .

----------

Examples using both input types for potentials can be found in examples/tild. 

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`kspace_style <kspace_style>`
:doc:`kspace_modify <kspace_modify>`

Default
"""""""

The option defaults are mesh = 0 0 0, order = 5 (TILD), minorder = 2, overlap = yes, tild subtract_rho0 = yes, and tild normalize_by_rho0 = yes.

----------

.. _Chao:

**(Chao)** Chao, H., Koski, J. & Riggleman, R. 
"Solvent vapor annealing in block copolymer nanocomposite films: 
a dynamic mean field approach" Soft Matter, 13(1), 239-249 (2017) 
https://doi.org/10.1039/c6sm00770h

.. _Fredrickson:

**(Fredrickson)** Fredrickson, G. H. and Orland, H. 
"Dynamics of polymers: A mean-field theory" J Chem Phys,
140, 084902 (2014) https://doi.org/10.1063/1.4865911

.. _Grzetic:

**(Grzetic)** Grzetic, D. J., Wickman, R. A., and Shi, A.-C., "Statistical
dynamics of classical systems: A self-consistent field approach", J Chem Phys, 
140, 244907 (2014) https://doi.org/10.1063/1.4884825
