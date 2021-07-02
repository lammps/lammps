.. index:: pair_style cac/buck

pair_style cac/buck command
============================

Syntax
""""""

.. parsed-literal::

   pair_style style args

* style = *cac/buck*
* args = list of arguments for a particular style

.. parsed-literal::

     *cac/buck* args = cutoff quadrature_flag
       cutoff = global cutoff for Buckingham interactions (distance units)
       quadrature_flag = left blank or set to the *one* keyword

Examples
""""""""

.. code-block:: LAMMPS

   pair_style cac/buck 8.0
   pair_coeff * * 100.0 1.5 200.0
   pair_coeff * * 100.0 1.5 200.0 3.0

   pair_style cac/buck 8.0 one
   pair_coeff * * 100.0 1.5 200.0
   pair_coeff * * 100.0 1.5 200.0 3.0

Description
"""""""""""

The *cac/buck* style computes a computes pair forces using the Buckingham 
potential for a CAC model. Other than the allotment of the *one* keyword,
the rest of the input parameters and pair_coeff settings are as found in
:doc:`pair_style buck <pair_buck>`. For more information on the quadrature
scheme and the effect of the *one* keyword see :doc:`Howto CAC <Howto_cac>` 
for the details.

.. note::

   CAC Pair Styles do not currently support being substyles of pair_style
   hybrid

Restrictions
""""""""""""

CAC Pair styles require the definition of a CAC atom style.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`atom_style <atom_style>`

**Default:** none
