.. index:: fix qeq/point
.. index:: fix qeq/point/omp
.. index:: fix qeq/point/xlmd
.. index:: fix qeq/point/xlmd/omp
.. index:: fix qeq/shielded
.. index:: fix qeq/shielded/omp
.. index:: fix qeq/shielded/xlmd
.. index:: fix qeq/shielded/xlmd/omp
.. index:: fix qeq/slater
.. index:: fix qeq/slater/omp
.. index:: fix qeq/slater/xlmd
.. index:: fix qeq/slater/xlmd/omp
.. index:: fix qeq/ctip
.. index:: fix qeq/ctip/omp
.. index:: fix qeq/dynamic
.. index:: fix qeq/dynamic/omp
.. index:: fix qeq/fire
.. index:: fix qeq/fire/omp

fix qeq/point command
=====================

Accelerator Variants: *qeq/point/omp*

fix qeq/point/xlmd command
==========================

Accelerator Variants: *qeq/point/xlmd/omp*

fix qeq/shielded command
========================

Accelerator Variants: *qeq/shielded/omp*

fix qeq/shielded/xlmd command
=============================

Accelerator Variants: *qeq/shielded/xlmd/omp*

fix qeq/slater command
======================

Accelerator Variants: *qeq/slater/omp*

fix qeq/slater/xlmd command
===========================

Accelerator Variants: *qeq/slater/xlmd/omp*

fix qeq/ctip command
====================

Accelerator Variants: *qeq/ctip/omp*

fix qeq/dynamic command
=======================

Accelerator Variants: *qeq/dynamic/omp*

fix qeq/fire command
====================

Accelerator Variants: *qeq/fire/omp*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style Nevery cutoff tolerance maxiter qfile keyword ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *qeq/point* or *qeq/point/xlmd* or *qeq/shielded* or *qeq/shielded/xlmd* or *qeq/slater* or *qeq/slater/xlmd* or *qeq/ctip* or *qeq/dynamic* or *qeq/fire*
* Nevery = perform charge equilibration every this many steps
* cutoff = global cutoff for charge-charge interactions (distance unit)
* tolerance = precision to which charges will be equilibrated
* maxiter = maximum iterations to perform charge equilibration
* qfile = a filename with QEq parameters or *coul/streitz* or *coul/ctip* or *reaxff*
* zero or more keyword/value pairs may be appended
* keyword = *alpha* or *cdamp* or *maxrepeat* or *qdamp* or *qstep* or *warn* or *xlcg* or *xldamp* or *xlkappa*

  .. parsed-literal::

       *alpha* value = Slater type orbital exponent (qeq/slater only). Can be followed by optional arguments:
         *wolf* value = width of taper to terminate Coulomb integrals for the Wolf summation (default value is zero)
         *dsf* value = width of taper to terminate Coulomb integrals for the Fennell-Gezelter summation (default value is zero)
       *cdamp* value = damping parameter for Coulomb interactions (qeq/ctip only)
       *maxrepeat* value = number of equilibration cycles allowed to ensure no atoms cross charge bounds (qeq/ctip only)
       *qdamp* value = damping factor for damped dynamics charge solver (qeq/dynamic and qeq/fire only)
       *qstep* value = time step size for damped dynamics charge solver (qeq/dynamic and qeq/fire only)
       *warn* value = do (=yes) or do not (=no) print a warning when the maximum number of iterations is reached
       *xlcg* value = number of conjugate gradient iterations per step (\ *xlmd* styles only)
       *xldamp* value = order of the dissipation term for the auxiliary variables, 0 or 3 to 9 (\ *xlmd* styles only)
       *xlkappa* value = coupling constant :math:`\kappa = \omega^2 \delta t^2` of the undamped propagation (\ *xlmd* styles only)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all qeq/point 1 10 1.0e-6 200 param.qeq1
   fix 1 qeq qeq/shielded 1 8 1.0e-6 100 param.qeq2
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2 wolf
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2 wolf 2.0
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2 dsf
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2 dsf 2.0
   fix 1 all qeq/ctip 1 12 1.0e-8 100 coul/ctip cdamp 0.30 maxrepeat 10
   fix 1 qeq qeq/dynamic 1 12 1.0e-3 100 my_qeq
   fix 1 all qeq/fire 1 10 1.0e-3 100 my_qeq qdamp 0.2 qstep 0.1
   fix 1 all qeq/shielded/xlmd 1 10 1.0e-6 400 reaxff
   fix 1 all qeq/shielded/xlmd 1 10 1.0e-6 400 reaxff xlcg 2
   fix 1 all qeq/point/xlmd 1 10 1.0e-6 200 param.qeq1 xldamp 0 xlkappa 2.0

Description
"""""""""""

Perform the charge equilibration (QEq) method as described in
:ref:`(Rappe) <Rappe1>` and formulated in :ref:`(Nakano)
<Nakano1>` (also known as the matrix inversion method) and in
:ref:`(Rick) <Rick1>` (also known as the extended Lagrangian
method) based on the electronegativity equilization principle.

These fixes can be used with any :doc:`pair style <pair_style>` in
LAMMPS, so long as per-atom charges are defined.  The most typical
use-case is in conjunction with a :doc:`pair style <pair_style>` that
performs charge equilibration periodically (e.g. every timestep), such
as the ReaxFF or Streitz-Mintmire potential.  But these fixes can also
be used with potentials that normally assume per-atom charges are fixed,
e.g. a :doc:`Buckingham <pair_buck>` or :doc:`LJ/Coulombic <pair_lj>`
potential.

Because the charge equilibration calculation is effectively independent
of the pair style, these fixes can also be used to perform a one-time
assignment of charges to atoms.  For example, you could define the QEq
fix, perform a zero-timestep run via the :doc:`run <run>` command
without any pair style defined which would set per-atom charges (based
on the current atom configuration), then remove the fix via the
:doc:`unfix <unfix>` command before performing further dynamics.

.. note::

   Computing and using charge values different from published
   values defined for a fixed-charge potential like Buckingham or CHARMM
   or AMBER, can have a strong effect on energies and forces, and
   produces a different model than the published versions.

.. note::

   The :doc:`fix qeq/comb <fix_qeq_comb>` command must still be used to
   perform charge equilibration with the :doc:`COMB potential
   <pair_comb>`.  The :doc:`fix qeq/reaxff <fix_qeq_reaxff>` command can be
   used to perform charge equilibration with the :doc:`ReaxFF force
   field <pair_reaxff>`, although fix qeq/shielded yields the same
   results as fix qeq/reaxff if *Nevery*, *cutoff*, and *tolerance*
   are the same.  Eventually the fix qeq/reaxff command will be
   deprecated.

The QEq method minimizes the electrostatic energy of the system (or
equalizes the derivative of energy with respect to charge of all the
atoms) by adjusting the partial charge on individual atoms based on
interactions with their neighbors within *cutoff*\ .  It requires a few
parameters in the appropriate units for each atom type which are read
from a file specified by *qfile*\ .  The file has the following format:

.. parsed-literal::

   1 chi eta gamma zeta qcore
   2 chi eta gamma zeta qcore
   ...
   Ntype chi eta gamma zeta qcore

except for fix style *qeq/ctip* where the format is:

.. parsed-literal::

   1 chi eta gamma zeta qcore qmin qmax omega
   2 chi eta gamma zeta qcore qmin qmax omega
   ...
   Ntype chi eta gamma zeta qcore qmin qmax omega

There have to be parameters given for every atom type. Wildcard entries
are possible using the same type range syntax as for "coeff" commands
(i.e., n\*m, n\*, \*m, \*). Later entries will overwrite previous ones.
Empty lines or any text following the pound sign (#) are ignored.  Each
line starts with the atom type followed by eight parameters.  Only a
subset of the parameters is used by each QEq style as described below,
thus the others can be set to 0.0 if desired, but all eight entries per
line are required.

* *chi* = electronegativity in energy units
* *eta* = self-Coulomb potential in energy units
* *gamma* = shielded Coulomb constant defined by :ref:`ReaxFF force field <vanDuin>` in distance units
* *zeta* = Slater type orbital exponent defined by the :ref:`Streitz <Streitz1>` potential in reverse distance units
* *qcore* = charge of the nucleus defined by the :ref:`Streitz-Mintmire potential <Streitz1>` potential in charge units
* *qmin* = lower bound on the allowed charge defined by the :ref:`CTIP <CTIP1>` potential in charge units
* *qmax* = upper bound on the allowed charge defined by the :ref:`CTIP <CTIP1>` potential in charge units
* *omega* = penalty parameter used to enforce charge bounds defined by the :ref:`CTIP <CTIP1>` potential in energy units

The fix qeq styles will print a warning if the charges are not
equilibrated within *tolerance* by *maxiter* steps, unless the
*warn* keyword is used with "no" as argument.  This latter option
may be useful for testing and benchmarking purposes, as it allows
to use a fixed number of QEq iterations when *tolerance* is set
to a small enough value to always reach the *maxiter* limit.  Turning
off warnings will avoid the excessive output in that case.

The *qeq/point* style describes partial charges on atoms as point
charges.  Interaction between a pair of charged particles is 1/r,
which is the simplest description of the interaction between charges.
Only the *chi* and *eta* parameters from the *qfile* file are used.
Note that Coulomb catastrophe can occur if repulsion between the pair
of charged particles is too weak.  This style solves partial charges
on atoms via the matrix inversion method.  A tolerance of 1.0e-6 is
usually a good number.

The *qeq/shielded* style describes partial charges on atoms also as
point charges, but uses a shielded Coulomb potential to describe the
interaction between a pair of charged particles.  Interaction through
the shielded Coulomb is given by equation (13) of the :ref:`ReaxFF force
field <vanDuin>` paper.  The shielding accounts for charge overlap
between charged particles at small separation.  This style is the same
as :doc:`fix qeq/reaxff <fix_qeq_reaxff>`, and can be used with
:doc:`pair_style reaxff <pair_reaxff>`.  Only the *chi*, *eta*, and
*gamma* parameters from the *qfile* file are used. When using the string
*reaxff* as filename, these parameters are extracted directly from an
active *reaxff* pair style.  This style solves partial charges on atoms
via the matrix inversion method.  A tolerance of 1.0e-6 is usually a
good number.

The *qeq/slater* style describes partial charges on atoms as spherical
charge densities centered around atoms via the Slater 1\ *s* orbital, so
that the interaction between a pair of charged particles is the product
of two Slater 1\ *s* orbitals.  The expression for the Slater 1\ *s*
orbital is given under equation (6) of the :ref:`Streitz
<Streitz1>` paper.  Only the *chi*, *eta*, *zeta*, and *qcore*
parameters from the *qfile* file are used. When using the string
*coul/streitz* as filename, these parameters are extracted directly from
an active *coul/streitz* pair style.  This style solves partial charges
on atoms via the matrix inversion method.  A tolerance of 1.0e-6 is
usually a good number.  Keyword *alpha* can be used to change the Slater
type orbital exponent.

.. versionadded:: 19Nov2024

The *qeq/ctip* style describes partial charges on atoms in the same way
as style *qeq/shielded* but also enables the definition of charge
bounds.  Only the *chi*, *eta*, *gamma*, *qmin*, *qmax*, and *omega*
parameters from the *qfile* file are used.  When using the string
*coul/ctip* as filename, these parameters are extracted directly from an
active *coul/ctip* pair style.  This style solves partial charges on
atoms via the matrix inversion method.  Keyword *cdamp* can be used to
change the damping parameter used to calculate Coulomb interactions.
Keyword *maxrepeat* can be used to adjust the number of equilibration
cycles allowed to ensure no atoms have crossed the charge bounds.  A
value of 10 is usually a good choice.  A tolerance between 1.0e-6 and
1.0e-8 is usually a good choice but should be checked in conjunction
with the timestep for adequate energy conservation during dynamic runs.

The *qeq/dynamic* style describes partial charges on atoms as point
charges that interact through 1/r, but the extended Lagrangian method is
used to solve partial charges on atoms.  Only the *chi* and *eta*
parameters from the *qfile* file are used.  Note that Coulomb
catastrophe can occur if repulsion between the pair of charged particles
is too weak.  A tolerance of 1.0e-3 is usually a good number.  Keyword
*qdamp* can be used to change the damping factor, while keyword *qstep*
can be used to change the time step size.

The :ref:`\ *qeq/fire*\ <Shan>` style describes the same charge model
and charge solver as the *qeq/dynamic* style, but employs a FIRE
minimization algorithm to solve for equilibrium charges.  Keyword
*qdamp* can be used to change the damping factor, while keyword *qstep*
can be used to change the time step size.

Note that *qeq/point*, *qeq/shielded*, *qeq/slater*, and *qeq/ctip* describe
different charge models, whereas the matrix inversion method and the
extended Lagrangian method (\ *qeq/dynamic* and *qeq/fire*\ ) are
different solvers.

Note that *qeq/point*, *qeq/dynamic* and *qeq/fire* styles all
describe charges as point charges that interact through 1/r
relationship, but solve partial charges on atoms using different
solvers.  These three styles should yield comparable results if the QEq
parameters and *Nevery*, *cutoff*, and *tolerance* are the same.
Style *qeq/point* is typically faster, *qeq/dynamic* scales better on
larger sizes, and *qeq/fire* is faster than *qeq/dynamic*\ .

.. note::

   In order to solve the self-consistent equations for electronegativity
   equalization, LAMMPS imposes the additional constraint that all the
   charges in the fix group must add up to zero.  The initial charge
   assignments should also satisfy this constraint.  LAMMPS will print a
   warning if that is not the case.

.. note::

   Developing QEq parameters (chi, eta, gamma, zeta, and qcore) is
   non-trivial.  Charges on atoms are not guaranteed to equilibrate with
   arbitrary choices of these parameters.  We do not develop these QEq
   parameters.  See the examples/qeq directory for some examples.

.. versionadded:: 11Feb2026

In previous versions of LAMMPS, the real-space summations of Coulomb
interactions were done by replacing *1/r* using a damped potential
*erfc(alpha*r)/r* with the parameter *alpha* controlling the rate of
decay. However, any finite value of *alpha* leads to a jump at the
cutoff, which interferes with equilibration if atoms move across the
cutoff. The charge-neutralized potential of :ref:`(Wolf et al.) <Wolf5>`
(*wolf*) and its extension by :ref:`(Fennell and Gezelter) <Fennell3>`
(*dsf*) solve this problem. An extension was implemented to specify the
width of taper (see :ref:`(Mei et al.) <Mei2>`) to smoothly terminate the
Coulomb integrals at the cutoff. This is done by specifying the optional
arguments *wolf* and *dsf* with the value representing the width of
taper that smoothly terminates the Coulomb integrals. For example, if
the cutoff is 8 A and the taper width is 2 A, the Coulomb integrals are
smoothly rescaled from their actual value at r=6 A to zero at r=8 A. For
backward compatibility, the default taper width is zero.

.. versionadded:: TBD

The *qeq/point/xlmd*, *qeq/shielded/xlmd*, and *qeq/slater/xlmd* styles
use the same charge models as their parent styles, but employ the
extended-Lagrangian scheme of :ref:`(Nomura) <Nomura2015>` instead of
iterating the charges to self-consistency on every step.  Auxiliary
per-atom variables, which follow the self-consistent solution through a
harmonic coupling, are propagated time-reversibly by a Verlet
integrator alongside the atoms.  On each timestep they provide the
initial guess for the iterative solver, which then applies only a small,
fixed number of conjugate gradient iterations (keyword *xlcg*, default
2) instead of iterating to *tolerance*\ .  This reduces the cost of the
charge equilibration several-fold (typically 5x to 8x) while conserving
the total energy far better than simply truncating or loosening the
converged solve, which can heat up the system or become unstable within
a few hundred steps.

By default, the auxiliary variables are propagated including a weak
dissipation term of order 5 following :ref:`(Niklasson) <Niklasson2009>`.
The order can be changed with the *xldamp* keyword: higher orders damp
more weakly (and thus perturb the time-reversible dynamics less), lower
orders damp more strongly (and thus bias the charges more).  A value of
0 selects the original, fully time-reversible scheme of :ref:`(Nomura)
<Nomura2015>` without dissipation; in this case the *xlkappa* keyword
may be used to change the coupling constant :math:`\kappa = \omega^2
\delta t^2` from its default value of 2.0.  The undamped propagation
slowly accumulates numerical noise in the auxiliary variables, which
makes long simulations (more than a few thousand steps) unstable, so
*xldamp* 0 is mainly useful for validation studies.  With the default
settings the energy conservation was close to that of fully converged
solves in our tests.  Reducing *xlcg* to 1 reproduces the setting of
:ref:`(Nomura) <Nomura2015>` and further reduces the cost, but with a
smaller stability margin: deviations may accumulate unnoticed over
many picoseconds before degrading the dynamics.

The extended-Lagrangian propagation requires a valid history of
previous solutions.  Whenever such a history is not available -- at the
beginning of every :doc:`run <run>`, after reading a :doc:`restart file
<read_restart>`, and during :doc:`energy minimization <minimize>` (which
always uses fully converged solves) -- the charges are equilibrated to
*tolerance* as with the parent fix styles, and the extended-Lagrangian
propagation switches on automatically after the first few steps.
Because the auxiliary dynamics is tied to the timestep, these styles
require *Nevery* = 1.  The charges obtained from the truncated solves
are not exactly at the energy minimum; the deviations remain small and
bounded, but it is recommended to verify energy conservation against a
fully converged reference run when applying these styles to a new
system.  The average number of solver iterations reported by the fix
(see below) shows the warm-up and the savings directly.

.. note::

   For charge equilibration with :doc:`pair_style reaxff <pair_reaxff>`
   the recommended combination is *checkqeq no* in the pair style and
   ``fix qeq/shielded/xlmd 1 10.0 1.0e-6 400 reaxff``, which extracts
   the QEq parameters directly from the pair style.  This is the
   extended-Lagrangian replacement for :doc:`fix qeq/reaxff
   <fix_qeq_reaxff>` requested in issue `#507
   <https://github.com/lammps/lammps/issues/507>`_.  Note that,
   different from fix qeq/reaxff, the QEQ package fixes ignore
   :doc:`fix efield <fix_efield>`.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about these fixes is written to :doc:`binary restart
files <restart>`.  These fixes compute a global scalar with the number
of iterations used by the charge solver on the current timestep, which
can be accessed by various :doc:`output commands <Howto_output>`.  No
parameter of these fixes can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

These fixes are invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

These fixes are part of the QEQ package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

These qeq fixes will ignore electric field contributions from
:doc:`fix efield <fix_efield>`.

Related commands
""""""""""""""""

:doc:`fix qeq/reaxff <fix_qeq_reaxff>`, :doc:`fix qeq/comb <fix_qeq_comb>`

Default
"""""""

warn yes; for the *xlmd* styles additionally: xlcg 2, xldamp 5, xlkappa 2.0

----------

.. _Rappe1:

**(Rappe)** A. K. Rappe and W. A. Goddard III, J Physical
Chemistry, 95, 3358-3363 (1991).

.. _Nakano1:

**(Nakano)** A. Nakano, Computer Physics Communications, 104, 59-69 (1997).

.. _Rick1:

**(Rick)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys 101,
6141 (1994).

.. _Streitz1:

**(Streitz)** F. H. Streitz, J. W. Mintmire, Physical Review B, 50,
16, 11996 (1994)

.. _CTIP1:

**(CTIP)** G. Plummer, J. P. Tavenner, M. I. Mendelev, Z. Wu, J. W. Lawson,
J Chemical Physics, 162, 054709 (2025)

.. _vanDuin:

**(ReaxFF)** A. C. T. van Duin, S. Dasgupta, F. Lorant, W. A. Goddard III, J
Physical Chemistry, 105, 9396-9049 (2001)

.. _Shan:

**(QEq/Fire)** T.-R. Shan, A. P. Thompson, S. J. Plimpton, in preparation

.. _Wolf5:

**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J. Chem. Phys. 110, 8254 (1999).

.. _Fennell3:

**(Fennell)** J. Fennell, J. D. Gezelter, J. Chem. Phys. 124, 234104 (2006).

.. _Mei2:

**(Mei)** J. Mei, J. W. Davenport, G. W. Fernando, Phys. Rev. B 43, 4653 (1991).

.. _Nomura2015:

**(Nomura)** K. Nomura, P. E. Small, R. K. Kalia, A. Nakano, P. Vashishta,
Computer Physics Communications, 192, 91-96 (2015).

.. _Niklasson2009:

**(Niklasson)** A. M. N. Niklasson, P. Steneteg, A. Odell, N. Bock,
M. Challacombe, C. J. Tymczak, E. Holmstrom, G. Zheng, V. Weber,
J Chemical Physics, 130, 214109 (2009).
