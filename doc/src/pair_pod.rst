.. index:: pair_style pod

pair_style pod command
======================

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

Pair style *pod* defines the proper orthogonal descriptor (POD) potential. The mathematical 
definition of the POD potential is described from :doc:`compute podfit <compute_podfit>`, which is 
used to fit the POD potential to *ab initio* energy and force data. More details on the POD potential
are given in :ref:`(Nguyen) <Nguyen2022_2>`

Only a single pair_coeff command is used with the *pod* style which
specifies a POD parameter file followed by a coefficient file.

The coefficient file contains coefficients for the POD potential. The top of the coefficient 
file can contain any number of blank and comment lines (start with #), but follows a 
strict format after that. The first non-blank non-comment line must contain:

* POD_coefficients: *ncoeff*

This is followed by *ncoeff* coefficients, one per line. The coefficient file
is generated after training the POD potential using :doc:`compute podfit <compute_podfit>`.  

The POD parameter file can contain blank and comment lines (start
with #) anywhere. Each non-blank non-comment line must contain one
keyword/value pair. See :doc:`compute podfit <compute_podfit>` for the description 
of all the keywords that can be assigned in the parameter file. 

----------

Restrictions
""""""""""""

This style is part of the ML-POD package.  It is only enabled if LAMMPS
was built with that package by setting -D PKG_ML-POD=on. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute podfit <compute_podfit>`,

Default
"""""""

none

----------

.. _Nguyen2022_2:

**(Nguyen)** Ngueyn, et. al., arxiv (2022).
