.. index:: pair\_style smd/tlsph

pair\_style smd/tlsph command
=============================

Syntax
""""""


.. parsed-literal::

   pair_style smd/tlsph args

Examples
""""""""

pair\_style smd/tlsph

Description
"""""""""""

The *smd/tlsph* style computes particle interactions according to
continuum mechanics constitutive laws and a Total-Lagrangian
Smooth-Particle Hydrodynamics algorithm.

This pair style is invoked with the following command:


.. parsed-literal::

   pair_style smd/tlsph
   pair_coeff i j \*COMMON rho0 E nu Q1 Q2 hg Cp &
                  \*END

Here, *i* and *j* denote the *LAMMPS* particle types for which this
pair style is defined. Note that *i* and *j* must be equal, i.e., no
*tlsph* cross interactions between different particle types are
allowed.  In contrast to the usual *LAMMPS* *pair coeff* definitions,
which are given solely a number of floats and integers, the *tlsph*
*pair coeff* definition is organized using keywords. These keywords
mark the beginning of different sets of parameters for particle
properties, material constitutive models, and damage models. The *pair
coeff* line must be terminated with the *\*END* keyword. The use the
line continuation operator *&* is recommended. A typical invocation of
the *tlsph* for a solid body would consist of an equation of state for
computing the pressure (the diagonal components of the stress tensor),
and a material model to compute shear stresses (the off-diagonal
components of the stress tensor). Damage and failure models can also
be added.

Please see the `SMD user guide <PDF/SMD_LAMMPS_userguide.pdf>`_ for a
complete listing of the possible keywords and material models.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

No mixing is performed automatically.  Currently, no part of USER-SMD
supports restarting nor minimization.  rRESPA does not apply to this
pair style.


----------


Restrictions
""""""""""""


This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------
