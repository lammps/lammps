.. index:: pair_style cac/coul/wolf

pair_style cac/coul/wolf command
=================================

Syntax
""""""

.. parsed-literal::

   pair_style style args

* style = *cac/coul/wolf*
* args = list of arguments for a particular style

.. parsed-literal::

     *cac/coul/wolf* args = alpha cutoff quadrature_flag
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for Coulombic interactions
       quadrature_flag = left blank or set to the *one* keyword

Examples
""""""""

.. code-block:: LAMMPS

   pair_style cac/coul/wolf
   pair_coeff * *

   pair_style cac/coul/wolf one
   pair_coeff * *

Description
"""""""""""

The *cac/coul/wof* style computes coulomb pair interactions using the Wolf
summation method for a CAC model. Other than the allotment of the 
*one* keyword, the rest of the input parameters and pair_coeff settings 
are as found in :doc:`pair_style coul <pair_coul>`. For more information on 
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
