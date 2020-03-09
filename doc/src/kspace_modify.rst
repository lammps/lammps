.. index:: kspace_modify

kspace_modify command
=====================

Syntax
""""""


.. code-block:: LAMMPS

   kspace_modify keyword value ...

* one or more keyword/value pairs may be listed
* keyword = *collective* or *compute* or *cutoff/adjust* or *diff* or *disp/auto* or *fftbench* or *force/disp/kspace* or *force/disp/real* or *force* or *gewald/disp* or *gewald* or *kmax/ewald* or *mesh* or *minorder* or *mix/disp* or *order/disp* or *order* or *overlap* or *scafacos* or *slab* or *splittol*
  
  .. parsed-literal::
  
       *collective* value = *yes* or *no*
       *compute* value = *yes* or *no*
       *cutoff/adjust* value = *yes* or *no*
       *diff* value = *ad* or *ik* = 2 or 4 FFTs for PPPM in smoothed or non-smoothed mode
       *disp/auto* value = yes or no
       *fftbench* value = *yes* or *no*
       *force/disp/real* value = accuracy (force units)
       *force/disp/kspace* value = accuracy (force units)
       *force* value = accuracy (force units)
       *gewald* value = rinv (1/distance units)
         rinv = G-ewald parameter for Coulombics
       *gewald/disp* value = rinv (1/distance units)
         rinv = G-ewald parameter for dispersion
       *kmax/ewald* value = kx ky kz
         kx,ky,kz = number of Ewald sum kspace vectors in each dimension
       *mesh* value = x y z
         x,y,z = grid size in each dimension for long-range Coulombics
       *mesh/disp* value = x y z
         x,y,z = grid size in each dimension for 1/r\^6 dispersion
       *minorder* value = M
         M = min allowed extent of Gaussian when auto-adjusting to minimize grid communication
       *mix/disp* value = *pair* or *geom* or *none*
       *order* value = N
         N = extent of Gaussian for PPPM or MSM mapping of charge to grid
       *order/disp* value = N
         N = extent of Gaussian for PPPM mapping of dispersion term to grid
       *overlap* = *yes* or *no* = whether the grid stencil for PPPM is allowed to overlap into more than the nearest-neighbor processor
       *pressure/scalar* value = *yes* or *no*
       *scafacos* values = option value1 value2 ...
         option = *tolerance*
           value = *energy* or *energy_rel* or *field* or *field_rel* or *potential* or *potential_rel*
         option = *fmm_tuning*
           value = *0* or *1*
       *slab* value = volfactor or *nozforce*
         volfactor = ratio of the total extended volume used in the
           2d approximation compared with the volume of the simulation domain
         *nozforce* turns off kspace forces in the z direction
       *splittol* value = tol
         tol = relative size of two eigenvalues (see discussion below)



Examples
""""""""


.. code-block:: LAMMPS

   kspace_modify mesh 24 24 30 order 6
   kspace_modify slab 3.0
   kspace_modify scafacos tolerance energy

Description
"""""""""""

Set parameters used by the kspace solvers defined by the
:doc:`kspace_style <kspace_style>` command.  Not all parameters are
relevant to all kspace styles.


----------


The *collective* keyword applies only to PPPM.  It is set to *no* by
default, except on IBM BlueGene machines.  If this option is set to
*yes*\ , LAMMPS will use MPI collective operations to remap data for
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


The *cutoff/adjust* keyword applies only to MSM. If this option is
turned on, the Coulombic cutoff will be automatically adjusted at the
beginning of the run to give the desired estimated error. Other
cutoffs such as LJ will not be affected. If the grid is not set using
the *mesh* command, this command will also attempt to use the optimal
grid that minimizes cost using an estimate given by
:ref:`(Hardy) <Hardy1>`. Note that this cost estimate is not exact, somewhat
experimental, and still may not yield the optimal parameters.


----------


The *diff* keyword specifies the differentiation scheme used by the
PPPM method to compute forces on particles given electrostatic
potentials on the PPPM mesh.  The *ik* approach is the default for
PPPM and is the original formulation used in :ref:`(Hockney) <Hockney1>`.  It
performs differentiation in Kspace, and uses 3 FFTs to transfer each
component of the computed fields back to real space for total of 4
FFTs per timestep.

The analytic differentiation *ad* approach uses only 1 FFT to transfer
information back to real space for a total of 2 FFTs per timestep.  It
then performs analytic differentiation on the single quantity to
generate the 3 components of the electric field at each grid point.
This is sometimes referred to as "smoothed" PPPM.  This approach
requires a somewhat larger PPPM mesh to achieve the same accuracy as
the *ik* method. Currently, only the *ik* method (default) can be
used for a triclinic simulation cell with PPPM. The *ad* method is
always used for MSM.

.. note::

   Currently, not all PPPM styles support the *ad* option.  Support
   for those PPPM variants will be added later.


----------


The *disp/auto* option controls whether the pppm/disp is allowed to
generate PPPM parameters automatically. If set to *no*\ , parameters have
to be specified using the *gewald/disp*\ , *mesh/disp*\ ,
*force/disp/real* or *force/disp/kspace* keywords, or
the code will stop with an error message. When this option is set to
*yes*\ , the error message will not appear and the simulation will start.
For a typical application, using the automatic parameter generation
will provide simulations that are either inaccurate or slow. Using this
option is thus not recommended. For guidelines on how to obtain good
parameters, see the :doc:`How-To <Howto_dispersion>` discussion.


----------


The *fftbench* keyword applies only to PPPM. It is off by default. If
this option is turned on, LAMMPS will perform a short FFT benchmark
computation and report its timings, and will thus finish a some seconds
later than it would if this option were off.


----------


The *force/disp/real* and *force/disp/kspace* keywords set the force
accuracy for the real and space computations for the dispersion part
of pppm/disp. As shown in :ref:`(Isele-Holder) <Isele-Holder1>`, optimal
performance and accuracy in the results is obtained when these values
are different.


----------


The *force* keyword overrides the relative accuracy parameter set by
the :doc:`kspace_style <kspace_style>` command with an absolute
accuracy.  The accuracy determines the RMS error in per-atom forces
calculated by the long-range solver and is thus specified in force
units.  A negative value for the accuracy setting means to use the
relative accuracy parameter.  The accuracy setting is used in
conjunction with the pairwise cutoff to determine the number of
K-space vectors for style *ewald*\ , the FFT grid size for style
*pppm*\ , or the real space grid size for style *msm*\ .


----------


The *gewald* keyword sets the value of the Ewald or PPPM G-ewald
parameter for charge as *rinv* in reciprocal distance units.  Without
this setting, LAMMPS chooses the parameter automatically as a function
of cutoff, precision, grid spacing, etc.  This means it can vary from
one simulation to the next which may not be desirable for matching a
KSpace solver to a pre-tabulated pairwise potential.  This setting can
also be useful if Ewald or PPPM fails to choose a good grid spacing
and G-ewald parameter automatically.  If the value is set to 0.0,
LAMMPS will choose the G-ewald parameter automatically.  MSM does not
use the *gewald* parameter.


----------


The *gewald/disp* keyword sets the value of the Ewald or PPPM G-ewald
parameter for dispersion as *rinv* in reciprocal distance units.  It
has the same meaning as the *gewald* setting for Coulombics.


----------


The *kmax/ewald* keyword sets the number of kspace vectors in each
dimension for kspace style *ewald*\ .  The three values must be positive
integers, or else (0,0,0), which unsets the option.  When this option
is not set, the Ewald sum scheme chooses its own kspace vectors,
consistent with the user-specified accuracy and pairwise cutoff. In
any case, if kspace style *ewald* is invoked, the values used are
printed to the screen and the log file at the start of the run.


----------


The *mesh* keyword sets the grid size for kspace style *pppm* or
*msm*\ .  In the case of PPPM, this is the FFT mesh, and each dimension
must be factorizable into powers of 2, 3, and 5.  In the case of MSM,
this is the finest scale real-space mesh, and each dimension must be
factorizable into powers of 2.  When this option is not set, the PPPM
or MSM solver chooses its own grid size, consistent with the
user-specified accuracy and pairwise cutoff.  Values for x,y,z of
0,0,0 unset the option.


----------


The *mesh/disp* keyword sets the grid size for kspace style
*pppm/disp*\ .  This is the FFT mesh for long-range dispersion and ach
dimension must be factorizable into powers of 2, 3, and 5.  When this
option is not set, the PPPM solver chooses its own grid size,
consistent with the user-specified accuracy and pairwise cutoff.
Values for x,y,z of 0,0,0 unset the option.


----------


The *minorder* keyword allows LAMMPS to reduce the *order* setting if
necessary to keep the communication of ghost grid point limited to
exchanges between nearest-neighbor processors.  See the discussion of
the *overlap* keyword for details.  If the *overlap* keyword is set to
*yes*\ , which is the default, this is never needed.  If it set to *no*
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


The *mix/disp* keyword selects the mixing rule for the dispersion
coefficients.  With *pair*\ , the dispersion coefficients of unlike
types are computed as indicated with :doc:`pair_modify <pair_modify>`.
With *geom*\ , geometric mixing is enforced on the dispersion
coefficients in the kspace coefficients. When using the arithmetic
mixing rule, this will speed-up the simulations but introduces some
error in the force computations, as shown in :ref:`(Wennberg) <Wennberg>`.
With *none*\ , it is assumed that no mixing rule is
applicable. Splitting of the dispersion coefficients will be performed
as described in :ref:`(Isele-Holder) <Isele-Holder1>`.

This splitting can be influenced with the *splittol* keywords.  Only
the eigenvalues that are larger than tol compared to the largest
eigenvalues are included. Using this keywords the original matrix of
dispersion coefficients is approximated. This leads to faster
computations, but the accuracy in the reciprocal space computations of
the dispersion part is decreased.


----------


The *order* keyword determines how many grid spacings an atom's charge
extends when it is mapped to the grid in kspace style *pppm* or *msm*\ .
The default for this parameter is 5 for PPPM and 8 for MSM, which
means each charge spans 5 or 8 grid cells in each dimension,
respectively.  For the LAMMPS implementation of MSM, the order can
range from 4 to 10 and must be even. For PPPM, the minimum allowed
setting is 2 and the maximum allowed setting is 7.  The larger the
value of this parameter, the smaller that LAMMPS will set the grid
size, to achieve the requested accuracy.  Conversely, the smaller the
order value, the larger the grid size will be.  Note that there is an
inherent trade-off involved: a small grid will lower the cost of FFTs
or MSM direct sum, but a larger order parameter will increase the cost
of interpolating charge/fields to/from the grid.

The PPPM order parameter may be reset by LAMMPS when it sets up the
FFT grid if the implied grid stencil extends beyond the grid cells
owned by neighboring processors.  Typically this will only occur when
small problems are run on large numbers of processors.  A warning will
be generated indicating the order parameter is being reduced to allow
LAMMPS to run the problem. Automatic adjustment of the order parameter
is not supported in MSM.


----------


The *order/disp* keyword determines how many grid spacings an atom's
dispersion term extends when it is mapped to the grid in kspace style
*pppm/disp*\ .  It has the same meaning as the *order* setting for
Coulombics.


----------


The *overlap* keyword can be used in conjunction with the *minorder*
keyword with the PPPM styles to adjust the amount of communication
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


The *pressure/scalar* keyword applies only to MSM. If this option is
turned on, only the scalar pressure (i.e. (Pxx + Pyy + Pzz)/3.0) will
be computed, which can be used, for example, to run an isotropic barostat.
Computing the full pressure tensor with MSM is expensive, and this option
provides a faster alternative. The scalar pressure is computed using a
relationship between the Coulombic energy and pressure :ref:`(Hummer) <Hummer>`
instead of using the virial equation. This option cannot be used to access
individual components of the pressure tensor, to compute per-atom virial,
or with suffix kspace/pair styles of MSM, like OMP or GPU.


----------


The *scafacos* keyword is used for settings that are passed to the
ScaFaCoS library when using :doc:`kspace_style scafacos <kspace_style>`.

The *tolerance* option affects how the *accuracy* specified with the
:doc:`kspace_style <kspace_style>` command is interpreted by ScaFaCoS.
The following values may be used:

* energy = absolute accuracy in total Coulombic energy
* energy\_rel = relative accuracy in total Coulombic energy
* potential = absolute accuracy in total Coulombic potential
* potential\_rel = relative accuracy in total Coulombic potential
* field = absolute accuracy in electric field
* field\_rel = relative accuracy in electric field

The values with suffix \_rel indicate the tolerance is a relative
tolerance; the other values impose an absolute tolerance on the given
quantity. Absolute tolerance in this case means, that for a given
quantity q and a given absolute tolerance of t\_a the result should
be between q-t\_a and q+t\_a. For a relative tolerance t\_r the relative
error should not be greater than t\_r, i.e. abs(1 - (result/q)) < t\_r.
As a consequence of this, the tolerance type should be checked, when
performing computations with a high absolute field / energy. E.g.
if the total energy in the system is 1000000.0 an absolute tolerance
of 1e-3 would mean that the result has to be between 999999.999 and
1000000.001, which would be equivalent to a relative tolerance of
1e-9.

The energy and energy\_rel values, set a tolerance based on the total
Coulombic energy of the system.  The potential and potential\_rel set a
tolerance based on the per-atom Coulombic energy.  The field and
field\_rel tolerance types set a tolerance based on the electric field
values computed by ScaFaCoS.  Since per-atom forces are derived from
the per-atom electric field, this effectively sets a tolerance on the
forces, similar to other LAMMPS KSpace styles, as explained on the
:doc:`kspace_style <kspace_style>` doc page.

Note that not all ScaFaCoS solvers support all tolerance types.
These are the allowed values for each method:

* fmm = energy and energy\_rel
* p2nfft = field (1d-,2d-,3d-periodic systems) or potential (0d-periodic)
* p3m = field
* ewald = field
* direct = has no tolerance tuning

If the tolerance type is not changed, the default values for the
tolerance type are the first values in the above list, e.g. energy
is the default tolerance type for the fmm solver.

The *fmm\_tuning* option is only relevant when using the FMM method.
It activates (value=1) or deactivates (value=0) an internal tuning
mechanism for the FMM solver.  The tuning operation runs sequentially
and can be very time-consuming.  Usually it is not needed for systems
with a homogeneous charge distribution. The default for this option is
therefore *0*\ . The FMM internal tuning is performed once, when the
solver is set up.


----------


The *slab* keyword allows an Ewald or PPPM solver to be used for a
systems that are periodic in x,y but non-periodic in z - a
:doc:`boundary <boundary>` setting of "boundary p p f".  This is done by
treating the system as if it were periodic in z, but inserting empty
volume between atom slabs and removing dipole inter-slab interactions
so that slab-slab interactions are effectively turned off.  The
volfactor value sets the ratio of the extended dimension in z divided
by the actual dimension in z.  The recommended value is 3.0.  A larger
value is inefficient; a smaller value introduces unwanted slab-slab
interactions.  The use of fixed boundaries in z means that the user
must prevent particle migration beyond the initial z-bounds, typically
by providing a wall-style fix.  The methodology behind the *slab*
option is explained in the paper by :ref:`(Yeh) <Yeh>`.  The *slab* option
is also extended to non-neutral systems :ref:`(Ballenegger) <Ballenegger>`.

An alternative slab option can be invoked with the *nozforce* keyword
in lieu of the volfactor.  This turns off all kspace forces in the z
direction.  The *nozforce* option is not supported by MSM. For MSM,
any combination of periodic, non-periodic, or shrink-wrapped
boundaries can be set using :doc:`boundary <boundary>` (the slab
approximation in not needed).  The *slab* keyword is not currently
supported by Ewald or PPPM when using a triclinic simulation cell. The
slab correction has also been extended to point dipole interactions
:ref:`(Klapp) <Klapp>` in :doc:`kspace_style <kspace_style>` *ewald/disp*\ ,
*ewald/dipole*\ , and *pppm/dipole*\ .

.. note::

   If you wish to apply an electric field in the Z-direction, in
   conjunction with the *slab* keyword, you should do it by adding
   explicit charged particles to the +/- Z surfaces.  If you do it via
   the :doc:`fix efield <fix_efield>` command, it will not give the correct
   dielectric constant due to the Yeh/Berkowitz :ref:`(Yeh) <Yeh>` correction
   not being compatible with how :doc:`fix efield <fix_efield>` works.


----------


The *force/disp/real* and *force/disp/kspace* keywords set the force
accuracy for the real and space computations for the dispersion part
of pppm/disp. As shown in :ref:`(Isele-Holder) <Isele-Holder1>`, optimal
performance and accuracy in the results is obtained when these values
are different.

The *disp/auto* option controls whether the pppm/disp is allowed to
generate PPPM parameters automatically. If set to *no*\ , parameters
have to be specified using the *gewald/disp*\ , *mesh/disp*\ ,
*force/disp/real* or *force/disp/kspace* keywords, or the code will
stop with an error message. When this option is set to *yes*\ , the
error message will not appear and the simulation will start.  For a
typical application, using the automatic parameter generation will
provide simulations that are either inaccurate or slow. Using this
option is thus not recommended.  For guidelines on how to obtain good
parameters, see the :doc:`Howto dispersion <Howto_dispersion>` doc page.


----------


Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`kspace_style <kspace_style>`, :doc:`boundary <boundary>`

Default
"""""""

The option defaults are mesh = mesh/disp = 0 0 0, order = order/disp =
5 (PPPM), order = 10 (MSM), minorder = 2, overlap = yes, force = -1.0,
gewald = gewald/disp = 0.0, slab = 1.0, compute = yes, cutoff/adjust =
yes (MSM), pressure/scalar = yes (MSM), fftbench = no (PPPM), diff =
ik (PPPM), mix/disp = pair, force/disp/real = -1.0, force/disp/kspace
= -1.0, split = 0, tol = 1.0e-6, and disp/auto = no. For pppm/intel,
order = order/disp = 7.  For scafacos settings, the scafacos tolerance
option depends on the method chosen, as documented above.  The
scafacos fmm\_tuning default = 0.


----------


.. _Hockney1:



**(Hockney)** Hockney and Eastwood, Computer Simulation Using Particles,
Adam Hilger, NY (1989).

.. _Yeh:



**(Yeh)** Yeh and Berkowitz, J Chem Phys, 111, 3155 (1999).

.. _Ballenegger:



**(Ballenegger)** Ballenegger, Arnold, Cerda, J Chem Phys, 131, 094107
(2009).

.. _Klapp:



**(Klapp)** Klapp, Schoen, J Chem Phys, 117, 8050 (2002).

.. _Hardy1:



**(Hardy)** David Hardy thesis: Multilevel Summation for the Fast
Evaluation of Forces for the Simulation of Biomolecules, University of
Illinois at Urbana-Champaign, (2006).

.. _Hummer:



**(Hummer)** Hummer, Gronbech-Jensen, Neumann, J Chem Phys, 109, 2791 (1998)

.. _Isele-Holder1:



**(Isele-Holder)** Isele-Holder, Mitchell, Hammond, Kohlmeyer, Ismail, J
Chem Theory Comput, 9, 5412 (2013).

.. _Wennberg:



**(Wennberg)** Wennberg, Murtola, Hess, Lindahl, J Chem Theory Comput,
9, 3527 (2013).
