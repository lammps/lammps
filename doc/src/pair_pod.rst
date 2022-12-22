.. index:: pair_style pod

pair_style pod command
========================

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
potential :ref:`(Nguyen) <Nguyen20221>`.  The mathematical definition of
the POD potential is described from :doc:`fitpod <fitpod_command>`, which is
used to fit the POD potential to *ab initio* energy and force data.

Only a single pair_coeff command is used with the *pod* style which
specifies a POD parameter file followed by a coefficient file.

The coefficient file (``Ta_coefficients.pod``) contains coefficients for the
POD potential. The top of the coefficient file can contain any number of
blank and comment lines (start with #), but follows a strict format
after that. The first non-blank non-comment line must contain:

* POD_coefficients: *ncoeff*

This is followed by *ncoeff* coefficients, one per line. The coefficient
file is generated after training the POD potential using :doc:`fitpod
<fitpod_command>`.

The POD parameter file (``Ta_param.pod``) can contain blank and comment lines
(start with #) anywhere. Each non-blank non-comment line must contain
one keyword/value pair. See :doc:`fitpod <fitpod_command>` for the description
of all the keywords that can be assigned in the parameter file.

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
directory lammps/examples/PACKAGES/pod.

----------

Restrictions
""""""""""""

This style is part of the ML-POD package.  It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This pair style does not compute per-atom energies and per-atom stresses.

Related commands
""""""""""""""""

:doc:`fitpod <fitpod_command>`,

Default
"""""""

none

----------

.. _Nguyen20221:

**(Nguyen)** Nguyen and Rohskopf, arXiv preprint arXiv:2209.02362 (2022).
