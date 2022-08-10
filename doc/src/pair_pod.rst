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

Pair style *pod* defines the proper orthogonal descriptor (POD) potential, 
:ref:`(Thompson) <Thompson20142>`.  The mathematical definition of the POD potential
is described from :doc:`compute podfit <compute_podfit>`, which is used to fit the POD
potential to *ab initio* energy and force data.  

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

.. _Thompson20142:

**(Thompson)** Thompson, Swiler, Trott, Foiles, Tucker, J Comp Phys, 285, 316 (2015).

.. _Bartok20102:

**(Bartok2010)** Bartok, Payne, Kondor, Csanyi, Phys Rev Lett, 104, 136403 (2010).

.. _Wood20182:

**(Wood)** Wood and Thompson, J Chem Phys, 148, 241721, (2018)

.. _Cusentino20202:

**(Cusentino)** Cusentino, Wood, and Thompson, J Phys Chem A, xxx, xxxxx, (2020)
