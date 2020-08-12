.. index:: pair_style cac/lj

pair_style cac/lj command
==========================

Syntax
""""""

.. parsed-literal::

   pair_style style args

* style = *cac/lj*
* args = list of arguments for a particular style

.. parsed-literal::

     *cac/lj* args = cutoff quadrature_flag
       cutoff = global cutoff for Lennard Jones interactions (distance units)
       quadrature_flag = left blank or set to the *one*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style cac/lj 8.0
   pair_coeff   1 1 0.583 2.27

   pair_style cac/lj 8.0 one
   pair_coeff   1 1 0.583 2.27

Description
"""""""""""

The *cac/lj* style computes pair forces using the standard 12/6
Lennard Jones potential for a CAC model. Other than the allotment of the 
*one* keyword, the rest of the input parameters and pair_coeff settings 
are as found in :doc:`pair_style lj/cut <pair_lj>`. For more information on 
the quadrature scheme and the effect of the *one* keyword see :doc:`Howto CAC <Howto_cac>` 
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
