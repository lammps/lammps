.. index:: fix qeq/reaxff
.. index:: fix qeq/reaxff/kk
.. index:: fix qeq/reaxff/omp

fix qeq/reaxff command
======================

Accelerator Variants: *qeq/reaxff/kk*, *qeq/reaxff/omp*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID qeq/reaxff Nevery cutlo cuthi tolerance params args

* ID, group-ID are documented in :doc:`fix <fix>` command
* qeq/reaxff = style name of this fix command
* Nevery = perform QEq every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reaxff or a filename
* one or more keywords or keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *dual* or *maxiter* or *nowarn* or *matfree* or *xlmd* or *xlcg* or *xldamp* or *xlkappa* or *xlguard*
       *dual* = process S and T matrix in parallel (only for qeq/reaxff/omp)
       *maxiter* N = limit the number of iterations to *N*
       *nowarn* = do not print a warning message if the maximum number of iterations was reached
       *matfree* = use a matrix-free approach for applying the H matrix (only for qeq/reaxff/kk)
       *xlmd* value = *yes* or *no* = whether or not to use the extended-Lagrangian propagation of the charges
       *xlcg* value = number of conjugate gradient iterations per step (requires *xlmd yes*)
       *xldamp* value = order of the dissipation term for the auxiliary variables, 0 or 3 to 9 (requires *xlmd yes*)
       *xlkappa* value = coupling constant :math:`\kappa = \omega^2 \delta t^2` of the undamped propagation (requires *xlmd yes*)
       *xlguard* value = relative residual of the truncated solves above which the extended-Lagrangian propagation is restarted, <= 0 disables (requires *xlmd yes*)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff
   fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 param.qeq maxiter 500
   fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff xlmd yes

Description
"""""""""""

Perform the charge equilibration (QEq) method as described in
:ref:`(Rappe and Goddard) <Rappe2>` and formulated in :ref:`(Nakano)
<Nakano2>`.  It is typically used in conjunction with the ReaxFF force
field model as implemented in the :doc:`pair_style reaxff <pair_reaxff>`
command, but it can be used with any potential in LAMMPS, so long as it
defines and uses charges on each atom.  The :doc:`fix qeq/comb
<fix_qeq_comb>` command should be used to perform charge equilibration
with the :doc:`COMB potential <pair_comb>`.  For more technical details
about the charge equilibration performed by fix qeq/reaxff, see the
:ref:`(Aktulga) <qeq-Aktulga>` paper.

The QEq method minimizes the electrostatic energy of the system by
adjusting the partial charge on individual atoms based on interactions
with their neighbors.  It requires some parameters for each atom type.
If the *params* setting above is the word "reaxff", then these are
extracted from the :doc:`pair_style reaxff <pair_reaxff>` command and
the ReaxFF force field file it reads in.  If a file name is specified
for *params*, then the parameters are taken from the specified file and
the file must contain one line for each atom type.  The latter form must
be used when performing QEq with a non-ReaxFF potential.  Each line
should be formatted as follows:

.. parsed-literal::

   itype chi eta gamma

where *itype* is the atom type from 1 to Ntypes, *chi* denotes the
electronegativity in eV, *eta* denotes the self-Coulomb potential in eV,
and *gamma* denotes the valence orbital exponent.  Note that these 3
quantities are also in the ReaxFF potential file, except that eta is
defined here as twice the eta value in the ReaxFF file.  Note that
unlike the rest of LAMMPS, the units of this fix are hard-coded to be A,
eV, and electronic charge.

The optional *dual* keyword allows to perform the optimization of the S
and T matrices in parallel.  This is only supported for the
*qeq/reaxff/omp* style.  Otherwise they are processed separately.  The
*qeq/reaxff/kk* style always solves the S and T matrices in parallel.

The optional *maxiter* keyword allows changing the max number of
iterations in the linear solver.  The default value is 200.

The optional *nowarn* keyword silences the warning message printed when
the maximum number of iterations was reached.  This can be useful for
comparing serial and parallel results where having the same fixed number
of QEq iterations is desired, which can be achieved by using a very
small tolerance and setting *maxiter* to the desired number of
iterations.

.. versionadded:: TBD

The optional *xlmd* keyword with the value *yes* enables the
extended-Lagrangian propagation of the charges, which is explained in
more detail for the :doc:`fix qeq/shielded/xlmd <fix_qeq>` command.
Auxiliary per-atom variables are propagated time-reversibly alongside
the atoms and provide the initial guess for the iterative solver, which
then applies only a small, fixed number of conjugate gradient iterations
(keyword *xlcg*, default 2) instead of iterating to *tolerance*.  This
makes the charge equilibration part of a time step about 1.5 to 3 times
faster (see the estimates below), while the energy conservation usually
remains close to that of fully converged solves.  The *xldamp* and
*xlkappa* keywords set the dissipation and the coupling of the auxiliary
variables in the same way as for fix qeq/shielded/xlmd (default:
dissipation of order 5).  The *xlguard* keyword sets the threshold of a
safeguard against occasional instabilities of the propagation, which is
also explained there: when the relative residual of a truncated solve
exceeds this value (default 1.0e-3), the solve is completed to
*tolerance*, a warning is printed (unless the *nowarn* keyword is used),
and the propagation is restarted from the converged charges.  A value <=
0 disables the safeguard.  At the beginning of every run, after the
number of atoms changed, and during :doc:`energy minimization
<minimize>`, the charges are equilibrated to *tolerance* until a valid
history of solutions from consecutive steps is available.  The
extended-Lagrangian propagation requires *Nevery* = 1.  It is available
for all variants of this fix and for the :doc:`fix qtpie/reaxff
<fix_qtpie_reaxff>` and :doc:`fix qeq/rel/reaxff <fix_qeq_rel_reaxff>`
commands, but not for :doc:`fix acks2/reaxff <fix_acks2_reaxff>`.  It is
recommended to verify the energy conservation of a simulation against a
run with fully converged charges.

The optional *matfree* keyword replaces the sequence of explicitly
constructing the H matrix, then (repeatedly) applying it with a
matrix-free approach where the H matrix is effectively regenerated each
time it is applied.  This trades performance for reduced memory
requirements because it avoids the overheads of storing the matrix.
This is only supported for the *qeq/reaxff/kk* style, with both full and
half qeq neighbor lists supported.

.. versionadded:: TBD

The following estimates for using the *xlmd* and *matfree* keywords are
from the PETN crystal benchmark in the ``bench/POTENTIALS`` folder for a
system of about 1M atoms.  With *xlmd yes*, only 3 conjugate gradient
iterations per time step are used instead of an average of about 12.
This reduced the time per step by about 15% on 8 CPU cores, by about 20%
on one GPU, and by 35% on four GPUs.  The savings are larger on multiple
GPUs, since every iteration requires communication between them.  With
the *matfree* keyword, the GPU memory needed was about a third smaller,
so that up to 75% more atoms fit on one GPU, but the simulation took
about 50% longer with GPUs and more than twice as much time on the CPU.
Combined with *xlmd yes*, which applies the matrix only a few times per
step, the simulation was still up to 35% faster on GPUs and about as
fast on the CPU.  The actual savings depend on the individual system,
the *tolerance*, the hardware, and the number of MPI processes or GPUs.

.. versionchanged:: TBD

If the fix group does not contain all atoms, only the charges of the
atoms in the fix group are equilibrated.  The charges of all other atoms
are not changed, but they still contribute to the electrostatic
potential at the atoms in the fix group, i.e. they act as fixed external
charges.  All accelerated variants of this fix follow the same
convention.

.. note::

   In order to solve the self-consistent equations for electronegativity
   equalization, LAMMPS imposes the additional constraint that all the
   charges in the fix group must add up to zero.  The initial charge
   assignments should also satisfy this constraint.  LAMMPS will print a
   warning if that is not the case.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  This fix computes a global scalar (the number of
iterations) for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This fix is part of the REAXFF package.  It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This fix does not correctly handle interactions involving multiple
periodic images of the same atom.  Hence, it should not be used for
periodic cell dimensions smaller than the non-bonded cutoff radius,
which is typically :math:`10~\AA` for ReaxFF simulations.

This fix may be used in combination with :doc:`fix efield <fix_efield>`
and will apply the external electric field during charge equilibration,
but there may be only one fix efield instance used and the electric
field vector may only have components in non-periodic directions.
Equal-style variables can be used for electric field vector components
without any further settings.  Atom-style variables can be used for
spatially-varying electric field vector components, but the resulting
electric potential must be specified as an atom-style variable using the
*potential* keyword for `fix efield`.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix qeq/shielded <fix_qeq>`,
:doc:`fix acks2/reaxff <fix_acks2_reaxff>`, :doc:`fix qtpie/reaxff <fix_qtpie_reaxff>`,
:doc:`fix qeq/rel/reaxff <fix_qeq_rel_reaxff>`

Default
"""""""

maxiter 200, xlmd no, xlcg 2, xldamp 5, xlguard 1.0e-3. The last three
settings only apply for xlmd yes.

----------

.. _Rappe2:

**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).

.. _Nakano2:

**(Nakano)** Nakano, Computer Physics Communications, 104, 59-69 (1997).

.. _qeq-Aktulga:

**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).
