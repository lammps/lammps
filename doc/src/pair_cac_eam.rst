.. index:: pair_style cac/eam

pair_style cac/eam command
===========================

pair_style cac/eam/alloy command
=================================

Syntax
""""""

.. parsed-literal::

   pair_style style args

* style = *cac/eam* or *cac/eam/interp* or *cac/eam/alloy* or *cac/eam/alloy/interp*
* args = list of arguments for a particular style


.. parsed-literal::

     *cac/eam* args = quadrature_flag
       quadrature_flag = left blank or set to the *one* keyword

Examples
""""""""

.. code-block:: LAMMPS

   pair_style cac/eam
   pair_coeff     * * Cu_u3.eam

   pair_style cac/eam/alloy
   pair_coeff     * * CoAl.eam.alloy

   pair_style cac/eam one
   pair_coeff     * * Cu_u3.eam

   pair_style cac/eam/interp
   pair_coeff     * * Cu_u3.eam

Description
"""""""""""

The *cac/eam* style computes pairwise forces using
the Embedded-Atom-Method (EAM) potential for a CAC model. Other than the allotment
of the *one* keyword, the rest of the input parameters and pair_coeff settings
are as found in :doc:`pair_style eam <pair_eam>`. For more information on
the quadrature scheme and the effect of the *one* keyword see :doc:`Howto CAC <Howto_cac>`
for the details.

----------

Style *cac/eam/alloy* computes pairwise interactions using the same
formula as style *cac/eam* for a CAC model.  However the associated
:doc:`pair_coeff <pair_coeff>` command reads a DYNAMO *setfl* file
instead of a *funcfl* file.  *Setfl* files can be used to model a
single-element or alloy system. See :doc:`pair_style eam <pair_eam>` for
more information.

----------

The *cac/eam/interp* style computes pairwise forces using the Embedded-Atom-Method (EAM)
potential for a CAC model. The distinction between this style and the *cac/eam* style
is the use of finite element interpolation to compute electron densities.
As a result the potential calculation is substantially faster than that
performed in *cac/eam* at a cost to the accuracy.

----------

The distinction between *cac/eam/interp/alloy* and the *cac/eam/alloy* style
is the use of finite element interpolation to compute electron densities.
As a result the potential calculation is substantially faster than that
performed in *cac/eam/alloy* at a cost to the accuracy.

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
