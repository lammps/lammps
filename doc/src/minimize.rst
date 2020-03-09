.. index:: minimize

minimize command
================

minimize/kk command
===================

Syntax
""""""


.. parsed-literal::

   minimize etol ftol maxiter maxeval

* etol = stopping tolerance for energy (unitless)
* ftol = stopping tolerance for force (force units)
* maxiter = max iterations of minimizer
* maxeval = max number of force/energy evaluations

Examples
""""""""


.. parsed-literal::

   minimize 1.0e-4 1.0e-6 100 1000
   minimize 0.0 1.0e-8 1000 100000

Description
"""""""""""

Perform an energy minimization of the system, by iteratively adjusting
atom coordinates.  Iterations are terminated when one of the stopping
criteria is satisfied.  At that point the configuration will hopefully
be in local potential energy minimum.  More precisely, the
configuration should approximate a critical point for the objective
function (see below), which may or may not be a local minimum.

The minimization algorithm used is set by the
:doc:`min_style <min_style>` command.  Other options are set by the
:doc:`min_modify <min_modify>` command.  Minimize commands can be
interspersed with :doc:`run <run>` commands to alternate between
relaxation and dynamics.  The minimizers bound the distance atoms move
in one iteration, so that you can relax systems with highly overlapped
atoms (large energies and forces) by pushing the atoms off of each
other.

Alternate means of relaxing a system are to run dynamics with a small
or :doc:`limited timestep <fix_nve_limit>`.  Or dynamics can be run
using :doc:`fix viscous <fix_viscous>` to impose a damping force that
slowly drains all kinetic energy from the system.  The :doc:`pair_style soft <pair_soft>` potential can be used to un-overlap atoms while
running dynamics.
un-overlap atoms while running dynamics.

Note that you can minimize some atoms in the system while holding the
coordinates of other atoms fixed by applying :doc:`fix setforce
<fix_setforce>` to the other atoms.  See a fuller discussion of using
fixes while minimizing below.

The :doc:`minimization styles <min_style>` *cg*\ , *sd*\ , and *hftn*
involves an outer iteration loop which sets the search direction along
which atom coordinates are changed.  An inner iteration is then
performed using a line search algorithm.  The line search typically
evaluates forces and energies several times to set new coordinates.
Currently, a backtracking algorithm is used which may not be optimal
in terms of the number of force evaluations performed, but appears to
be more robust than previous line searches we've tried.  The
backtracking method is described in Nocedal and Wright's Numerical
Optimization (Procedure 3.1 on p 41).

The :doc:`minimization styles <min_style>` *quickmin*\ , *fire* and
*fire/old* perform damped dynamics using an Euler integration step.  Thus
they require a :doc:`timestep <timestep>` be defined.

.. note::

   The damped dynamic minimizers use whatever timestep you have
   defined via the :doc:`timestep <timestep>` command.  Often they
   will converge more quickly if you use a timestep about 10x larger
   than you would normally use for dynamics simulations.


----------


In all cases, the objective function being minimized is the total
potential energy of the system as a function of the N atom
coordinates:

.. math::
 E(r_1,r_2, \ldots ,r_N)  = & \sum_{i,j} E_{\it pair}(r_i,r_j) +
                              \sum_{ij} E_{\it bond}(r_i,r_j) +
                              \sum_{ijk} E_{\it angle}(r_i,r_j,r_k) + \\
                            & \sum_{ijkl} E_{\it dihedral}(r_i,r_j,r_k,r_l) +
                              \sum_{ijkl} E_{\it improper}(r_i,r_j,r_k,r_l) +
                              \sum_i E_{\it fix}(r_i)


where the first term is the sum of all non-bonded :doc:`pairwise
interactions <pair_style>` including :doc:`long-range Coulombic
interactions <kspace_style>`, the 2nd through 5th terms are :doc:`bond
<bond_style>`, :doc:`angle <angle_style>`, :doc:`dihedral
<dihedral_style>`, and :doc:`improper <improper_style>` interactions
respectively, and the last term is energy due to :doc:`fixes <fix>`
which can act as constraints or apply force to atoms, such as through
interaction with a wall.  See the discussion below about how fix
commands affect minimization.

The starting point for the minimization is the current configuration
of the atoms.


----------


The minimization procedure stops if any of several criteria are met:

* the change in energy between outer iterations is less than *etol*
* the 2-norm (length) of the global force vector is less than the *ftol*
* the line search fails because the step distance backtracks to 0.0
* the number of outer iterations or timesteps exceeds *maxiter*
* the number of total force evaluations exceeds *maxeval*

.. note::

   the :doc:`minimization style <min_style>` *spin*\ ,
   *spin/cg*\ , and *spin/lbfgs* replace
   the force tolerance *ftol* by a torque tolerance.
   The minimization procedure stops if the 2-norm (length) of the torque vector on atom
   (defined as the cross product between the
   atomic spin and its precession vectors omega) is less than *ftol*\ ,
   or if any of the other criteria are met. Torque have the same units as the energy.

.. note::

   You can also use the :doc:`fix halt <fix_halt>` command to specify
   a general criterion for exiting a minimization, that is a
   calculation performed on the state of the current system, as
   defined by an :doc:`equal-style variable <variable>`.

For the first criterion, the specified energy tolerance *etol* is
unitless; it is met when the energy change between successive
iterations divided by the energy magnitude is less than or equal to
the tolerance.  For example, a setting of 1.0e-4 for *etol* means an
energy tolerance of one part in 10\^4.  For the damped dynamics
minimizers this check is not performed for a few steps after
velocities are reset to 0, otherwise the minimizer would prematurely
converge.

For the second criterion, the specified force tolerance *ftol* is in
force units, since it is the length of the global force vector for all
atoms, e.g. a vector of size 3N for N atoms.  Since many of the
components will be near zero after minimization, you can think of
*ftol* as an upper bound on the final force on any component of any
atom.  For example, a setting of 1.0e-4 for *ftol* means no x, y, or z
component of force on any atom will be larger than 1.0e-4 (in force
units) after minimization.

Either or both of the *etol* and *ftol* values can be set to 0.0, in
which case some other criterion will terminate the minimization.

During a minimization, the outer iteration count is treated as a
timestep.  Output is triggered by this timestep, e.g. thermodynamic
output or dump and restart files.

Using the :doc:`thermo_style custom <thermo_style>` command with the
*fmax* or *fnorm* keywords can be useful for monitoring the progress
of the minimization.  Note that these outputs will be calculated only
from forces on the atoms, and will not include any extra degrees of
freedom, such as from the :doc:`fix box/relax <fix_box_relax>` command.

Following minimization, a statistical summary is printed that lists
which convergence criterion caused the minimizer to stop, as well as
information about the energy, force, final line search, and iteration
counts.  An example is as follows:


.. parsed-literal::

   Minimization stats:
     Stopping criterion = max iterations
     Energy initial, next-to-last, final =
          -0.626828169302     -2.82642039062     -2.82643549739
     Force two-norm initial, final = 2052.1 91.9642
     Force max component initial, final = 346.048 9.78056
     Final line search alpha, max atom move = 2.23899e-06 2.18986e-05
     Iterations, force evaluations = 2000 12724

The 3 energy values are for before and after the minimization and on
the next-to-last iteration.  This is what the *etol* parameter checks.

The two-norm force values are the length of the global force vector
before and after minimization.  This is what the *ftol* parameter
checks.

The max-component force values are the absolute value of the largest
component (x,y,z) in the global force vector, i.e. the infinity-norm
of the force vector.

The alpha parameter for the line-search, when multiplied by the max
force component (on the last iteration), gives the max distance any
atom moved during the last iteration.  Alpha will be 0.0 if the line
search could not reduce the energy.  Even if alpha is non-zero, if the
"max atom move" distance is tiny compared to typical atom coordinates,
then it is possible the last iteration effectively caused no atom
movement and thus the evaluated energy did not change and the
minimizer terminated.  Said another way, even with non-zero forces,
it's possible the effect of those forces is to move atoms a distance
less than machine precision, so that the energy cannot be further
reduced.

The iterations and force evaluation values are what is checked by the
*maxiter* and *maxeval* parameters.


----------


.. note::

   There are several force fields in LAMMPS which have
   discontinuities or other approximations which may prevent you from
   performing an energy minimization to high tolerances.  For example,
   you should use a :doc:`pair style <pair_style>` that goes to 0.0 at the
   cutoff distance when performing minimization (even if you later change
   it when running dynamics).  If you do not do this, the total energy of
   the system will have discontinuities when the relative distance
   between any pair of atoms changes from cutoff+epsilon to
   cutoff-epsilon and the minimizer may behave poorly.  Some of the
   many-body potentials use splines and other internal cutoffs that
   inherently have this problem.  The :doc:`long-range Coulombic styles <kspace_style>` (PPPM, Ewald) are approximate to within the
   user-specified tolerance, which means their energy and forces may not
   agree to a higher precision than the Kspace-specified tolerance.  In
   all these cases, the minimizer may give up and stop before finding a
   minimum to the specified energy or force tolerance.

Note that a cutoff Lennard-Jones potential (and others) can be shifted
so that its energy is 0.0 at the cutoff via the
:doc:`pair_modify <pair_modify>` command.  See the doc pages for
individual :doc:`pair styles <pair_style>` for details.  Note that
Coulombic potentials always have a cutoff, unless versions with a
long-range component are used (e.g. :doc:`pair_style lj/cut/coul/long <pair_lj>`).  The CHARMM potentials go to 0.0 at
the cutoff (e.g. :doc:`pair_style lj/charmm/coul/charmm <pair_charmm>`),
as do the GROMACS potentials (e.g. :doc:`pair_style lj/gromacs <pair_gromacs>`).

If a soft potential (:doc:`pair_style soft <pair_soft>`) is used the
Astop value is used for the prefactor (no time dependence).

The :doc:`fix box/relax <fix_box_relax>` command can be used to apply an
external pressure to the simulation box and allow it to shrink/expand
during the minimization.

Only a few other fixes (typically those that add forces) are invoked
during minimization.  See the doc pages for individual :doc:`fix <fix>`
commands to see which ones are relevant.  Current examples of fixes
that can be used include:

* :doc:`fix addforce <fix_addforce>`
* :doc:`fix addtorque <fix_addtorque>`
* :doc:`fix efield <fix_efield>`
* :doc:`fix enforce2d <fix_enforce2d>`
* :doc:`fix indent <fix_indent>`
* :doc:`fix lineforce <fix_lineforce>`
* :doc:`fix planeforce <fix_planeforce>`
* :doc:`fix setforce <fix_setforce>`
* :doc:`fix spring <fix_spring>`
* :doc:`fix spring/self <fix_spring_self>`
* :doc:`fix viscous <fix_viscous>`
* :doc:`fix wall <fix_wall>`
* :doc:`fix wall/region <fix_wall_region>`

.. note::

   Some fixes which are invoked during minimization have an
   associated potential energy.  For that energy to be included in the
   total potential energy of the system (the quantity being minimized),
   you MUST enable the :doc:`fix_modify <fix_modify>` *energy* option for
   that fix.  The doc pages for individual :doc:`fix <fix>` commands
   specify if this should be done.

.. note::

   The minimizers in LAMMPS do not allow for bonds (or angles, etc)
   to be held fixed while atom coordinates are being relaxed, e.g. via
   :doc:`fix shake <fix_shake>` or :doc:`fix rigid <fix_rigid>`.  See more
   info in the Restrictions section below.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


Features that are not yet implemented are listed here, in case someone
knows how they could be coded:

It is an error to use :doc:`fix shake <fix_shake>` with minimization
because it turns off bonds that should be included in the potential
energy of the system.  The effect of a fix shake can be approximated
during a minimization by using stiff spring constants for the bonds
and/or angles that would normally be constrained by the SHAKE
algorithm.

:doc:`Fix rigid <fix_rigid>` is also not supported by minimization.  It
is not an error to have it defined, but the energy minimization will
not keep the defined body(s) rigid during the minimization.  Note that
if bonds, angles, etc internal to a rigid body have been turned off
(e.g. via :doc:`neigh_modify exclude <neigh_modify>`), they will not
contribute to the potential energy which is probably not what is
desired.

Pair potentials that produce torque on a particle (e.g. :doc:`granular potentials <pair_gran>` or the :doc:`GayBerne potential <pair_gayberne>` for ellipsoidal particles) are not
relaxed by a minimization.  More specifically, radial relaxations are
induced, but no rotations are induced by a minimization, so such a
system will not fully relax.

Related commands
""""""""""""""""

:doc:`min_modify <min_modify>`, :doc:`min_style <min_style>`,
:doc:`run_style <run_style>`

**Default:** none
