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
   pair_coeff * * pod.txt coefficient.txt

Description
"""""""""""

.. versionadded:: TBD

Pair style *pod* defines the proper orthogonal descriptor (POD)
potential :ref:`(Nguyen) <Nguyen20221>`.  The mathematical definition of
the POD potential is described from :doc:`fitpod <fitpod>`, which is
used to fit the POD potential to *ab initio* energy and force data.

Only a single pair_coeff command is used with the *pod* style which
specifies a POD parameter file followed by a coefficient file.

The coefficient file (``coefficient.txt``) contains coefficients for the
POD potential. The top of the coefficient file can contain any number of
blank and comment lines (start with #), but follows a strict format
after that. The first non-blank non-comment line must contain:

* POD_coefficients: *ncoeff*

This is followed by *ncoeff* coefficients, one per line. The coefficient
file is generated after training the POD potential using :doc:`fitpod
<fitpod>`.

The POD parameter file (``pod.txt``) can contain blank and comment lines
(start with #) anywhere. Each non-blank non-comment line must contain
one keyword/value pair. See :doc:`fitpod <fitpod>` for the description
of all the keywords that can be assigned in the parameter file.

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

:doc:`fitpod <fitpod>`,

Default
"""""""

none

----------

.. _Nguyen20221:

**(Nguyen)** Nguyen and Rohskopf, arXiv preprint arXiv:2209.02362 (2022).
