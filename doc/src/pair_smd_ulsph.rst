.. index:: pair_style smd/ulsph

pair_style smd/ulsph command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style smd/ulsph args

* these keywords must be given

.. parsed-literal::

   keyword = *\*DENSITY_SUMMATION* or *\*DENSITY_CONTINUITY* and *\*VELOCITY_GRADIENT* or *\*NO_VELOCITY_GRADIENT* and *\*GRADIENT_CORRECTION* or *\*NO_GRADIENT_CORRECTION*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smd/ulsph *DENSITY_CONTINUITY *VELOCITY_GRADIENT *NO_GRADIENT_CORRECTION

Description
"""""""""""

The *smd/ulsph* style computes particle interactions according to
continuum mechanics constitutive laws and an updated Lagrangian
Smooth-Particle Hydrodynamics algorithm.

This pair style is invoked similar to the following command:

.. code-block:: LAMMPS

   pair_style smd/ulsph *DENSITY_CONTINUITY *VELOCITY_GRADIENT *NO_GRADIENT_CORRECTION
   pair_coeff i j *COMMON rho0 c0 Q1 Cp hg &
                  *END

Here, *i* and *j* denote the *LAMMPS* particle types for which this
pair style is defined. Note that *i* and *j* can be different, i.e.,
*ulsph* cross interactions between different particle types are
allowed. However, *i*\ --\ *i* respectively *j*\ --\ *j* pair_coeff lines have
to precede a cross interaction.  In contrast to the usual *LAMMPS*
*pair coeff* definitions, which are given solely a number of floats
and integers, the *ulsph* *pair coeff* definition is organized using
keywords. These keywords mark the beginning of different sets of
parameters for particle properties, material constitutive models, and
damage models. The *pair coeff* line must be terminated with the
*\*END* keyword. The use the line continuation operator *&* is
recommended. A typical invocation of the *ulsph* for a solid body
would consist of an equation of state for computing the pressure (the
diagonal components of the stress tensor), and a material model to
compute shear stresses (the off-diagonal components of the stress
tensor).

Note that the use of \*GRADIENT_CORRECTION can lead to severe numerical
instabilities. For a general fluid simulation, \*NO_GRADIENT_CORRECTION
is recommended.

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
