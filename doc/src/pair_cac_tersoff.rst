.. index:: pair_style cac/sw

pair_style cac/sw command
==========================

Syntax
""""""

.. parsed-literal::

   pair_style style args

* style = *cac/tersoff*
* args = list of arguments for a particular style

.. parsed-literal::

     *cac/tersoff* args = quadrature_flag
       quadrature_flag = left blank or set to the *one* keyword

Examples
""""""""

.. code-block:: LAMMPS

   pair_style cac/tersoff
   pair_coeff * * Si.tersoff Si
   pair_coeff * * Si.tersoff Si C Si

   pair_style cac/sw one
   pair_coeff * * Si.tersoff Si

Description
"""""""""""

The *cac/tersoff* style computes 3-body forces using the Tersoff
potential for a CAC model. Other than the allotment of the
*one* keyword, the rest of the input parameters and pair_coeff settings
are as found in :doc:`pair_style tersoff <pair_tersoff>`. For more information on
the quadrature scheme and the effect of the *one* keyword see :doc:`Howto CAC <Howto_cac>`
for the details.

.. note::

   CAC Pair Styles do not currently support being substyles of pair_style
   hybrid

.. note::

   this pair style does not require the MANYBODY package to be installed.

Restrictions
""""""""""""

CAC Pair styles require the definition of a CAC atom style.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`atom_style <atom_style>`

**Default:** none
