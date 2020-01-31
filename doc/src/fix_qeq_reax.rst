.. index:: fix qeq/reax

fix qeq/reax command
====================

fix qeq/reax/kk command
=======================

fix qeq/reax/omp command
========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID qeq/reax Nevery cutlo cuthi tolerance params args

* ID, group-ID are documented in :doc:`fix <fix>` command
* qeq/reax = style name of this fix command
* Nevery = perform QEq every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reax/c or a filename
* args   = *dual* (optional)

Examples
""""""""


.. parsed-literal::

   fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
   fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 param.qeq

Description
"""""""""""

Perform the charge equilibration (QEq) method as described in :ref:`(Rappe and Goddard) <Rappe2>` and formulated in :ref:`(Nakano) <Nakano2>`.  It is
typically used in conjunction with the ReaxFF force field model as
implemented in the :doc:`pair_style reax/c <pair_reaxc>` command, but
it can be used with any potential in LAMMPS, so long as it defines and
uses charges on each atom.  The :doc:`fix qeq/comb <fix_qeq_comb>`
command should be used to perform charge equilibration with the :doc:`COMB potential <pair_comb>`.  For more technical details about the
charge equilibration performed by fix qeq/reax, see the
:ref:`(Aktulga) <qeq-Aktulga>` paper.

The QEq method minimizes the electrostatic energy of the system by
adjusting the partial charge on individual atoms based on interactions
with their neighbors.  It requires some parameters for each atom type.
If the *params* setting above is the word "reax/c", then these are
extracted from the :doc:`pair_style reax/c <pair_reaxc>` command and
the ReaxFF force field file it reads in.  If a file name is specified
for *params*\ , then the parameters are taken from the specified file
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
the *qeq/reax/omp* style. Otherwise they are processed separately.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  No global scalar or vector or per-atom
quantities are stored by this fix for access by various :doc:`output commands <Howto_output>`.  No parameter of this fix can be used
with the *start/stop* keywords of the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.


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


This fix is part of the USER-REAXC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix does not correctly handle interactions
involving multiple periodic images of the same atom. Hence, it should not
be used for periodic cell dimensions less than 10 angstroms.

Related commands
""""""""""""""""

:doc:`pair_style reax/c <pair_reaxc>`

**Default:** none


----------


.. _Rappe2:



**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).

.. _Nakano2:



**(Nakano)** Nakano, Computer Physics Communications, 104, 59-69 (1997).

.. _qeq-Aktulga:



**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).


