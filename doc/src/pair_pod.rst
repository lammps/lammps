.. index:: pair_style pod
.. index:: pair_style pod/kk

pair_style pod command
========================

Accelerator Variants: *pod/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style pod

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pod
   pair_coeff * * Ta_param.pod Ta_coefficients.pod Ta

Description
"""""""""""

.. versionadded:: 22Dec2022

Pair style *pod* defines the proper orthogonal descriptor (POD)
potential :ref:`(Nguyen and Rohskopf) <Nguyen20222b>`,
:ref:`(Nguyen2023) <Nguyen20232b>`, :ref:`(Nguyen2024) <Nguyen20242b>`,
and :ref:`(Nguyen and Sema) <Nguyen20243b>`.  The :doc:`fitpod
<fitpod_command>` is used to fit the POD potential.

Only a single pair_coeff command is used with the *pod* style which
specifies a POD parameter file followed by a coefficient file, a
projection matrix file, and a centroid file.

The POD parameter file (``Ta_param.pod``) can contain blank and comment
lines (start with #) anywhere. Each non-blank non-comment line must
contain one keyword/value pair. See :doc:`fitpod <fitpod_command>` for
the description of all the keywords that can be assigned in the
parameter file.

The coefficient file (``Ta_coefficients.pod``) contains coefficients for
the POD potential. The top of the coefficient file can contain any
number of blank and comment lines (start with #), but follows a strict
format after that. The first non-blank non-comment line must contain:

* model_coefficients: *ncoeff* *nproj* *ncentroid*

This is followed by *ncoeff* coefficients, *nproj* projection matrix entries,
and *ncentroid* centroid coordinates, one per line. The coefficient
file is generated after training the POD potential using :doc:`fitpod
<fitpod_command>`.

As an example, if a LAMMPS indium phosphide simulation has 4 atoms
types, with the first two being indium and the third and fourth being
phophorous, the pair_coeff command would look like this:

.. code-block:: LAMMPS

   pair_coeff * * pod InP_param.pod InP_coefficients.pod In In P P

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The two filenames are for the parameter and coefficient files, respectively.
The two trailing 'In' arguments map LAMMPS atom types 1 and 2 to the
POD 'In' element. The two trailing 'P' arguments map LAMMPS atom types
3 and 4 to the POD 'P' element.

If a POD mapping value is specified as NULL, the mapping is not
performed.  This can be used when a *pod* potential is used as part of
the *hybrid* pair style.  The NULL values are placeholders for atom
types that will be used with other potentials.

Examples about training and using POD potentials are found in the
directory lammps/examples/PACKAGES/pod and the Github repo https://github.com/cesmix-mit/pod-examples.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential files.  Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This style is part of the ML-POD package.  It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fitpod <fitpod_command>`,
:doc:`compute pod/atom <compute_pod_atom>`,
:doc:`compute podd/atom <compute_pod_atom>`,
:doc:`compute pod/local <compute_pod_atom>`,
:doc:`compute pod/global <compute_pod_atom>`

Default
"""""""

none

----------

.. _Nguyen20222b:

**(Nguyen and Rohskopf)** Nguyen and Rohskopf,  Journal of Computational Physics, 480, 112030, (2023).

.. _Nguyen20232b:

**(Nguyen2023)** Nguyen, Physical Review B, 107(14), 144103, (2023).

.. _Nguyen20242b:

**(Nguyen2024)** Nguyen, Journal of Computational Physics, 113102, (2024).

.. _Nguyen20243b:

**(Nguyen and Sema)** Nguyen and Sema, https://arxiv.org/abs/2405.00306, (2024).


