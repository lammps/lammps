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

     keyword = *dual* or *maxiter* or *nowarn*
       *dual* = process S and T matrix in parallel (only for qeq/reaxff/omp)
       *maxiter* N = limit the number of iterations to *N*
       *nowarn* = do not print a warning message if the maximum number of iterations was reached

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff
   fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 param.qeq maxiter 500

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
for *params*, then the parameters are taken from the specified file
and the file must contain one line for each atom type.  The latter
form must be used when performing QeQ with a non-ReaxFF potential.
Each line should be formatted as follows:

.. parsed-literal::

   itype chi eta gamma

where *itype* is the atom type from 1 to Ntypes, *chi* denotes the
electronegativity in eV, *eta* denotes the self-Coulomb
potential in eV, and *gamma* denotes the valence orbital
exponent.  Note that these 3 quantities are also in the ReaxFF
potential file, except that eta is defined here as twice the eta value
in the ReaxFF file. Note that unlike the rest of LAMMPS, the units
of this fix are hard-coded to be A, eV, and electronic charge.

The optional *dual* keyword allows to perform the optimization
of the S and T matrices in parallel. This is only supported for
the *qeq/reaxff/omp* style. Otherwise they are processed separately.
The *qeq/reaxff/kk* style always solves the S and T matrices in
parallel.

The optional *maxiter* keyword allows changing the max number
of iterations in the linear solver. The default value is 200.

The optional *nowarn* keyword silences the warning message printed
when the maximum number of iterations was reached.  This can be
useful for comparing serial and parallel results where having the
same fixed number of QEq iterations is desired, which can be achieved
by using a very small tolerance and setting *maxiter* to the desired
number of iterations.

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

This fix is part of the REAXFF package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This fix does not correctly handle interactions involving multiple
periodic images of the same atom.  Hence, it should not be used for
periodic cell dimensions less than 10 Angstroms.

This fix may be used in combination with :doc:`fix efield <fix_efield>`
and will apply the external electric field during charge equilibration,
but there may be only one fix efield instance used and the electric field
vector may only have components in non-periodic directions. Equal-style
variables can be used for electric field vector components without any further
settings. Atom-style variables can be used for spatially-varying electric field
vector components, but the resulting electric potential must be specified
as an atom-style variable using the *potential* keyword for `fix efield`.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix qeq/shielded <fix_qeq>`

Default
"""""""

maxiter 200

----------

.. _Rappe2:

**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).

.. _Nakano2:

**(Nakano)** Nakano, Computer Physics Communications, 104, 59-69 (1997).

.. _qeq-Aktulga:

**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).
