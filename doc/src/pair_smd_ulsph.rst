.. index:: pair_style smd/ulsph

pair_style smd/ulsph command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style smd/ulsph density gradient correction

* density = *\*DENSITY_SUMMATION* or *\*DENSITY_CONTINUITY*
* gradient = *\*VELOCITY_GRADIENT* or *\*NO_VELOCITY_GRADIENT*
* correction = *\*GRADIENT_CORRECTION* or *\*NO_GRADIENT_CORRECTION*

.. code-block:: LAMMPS

   pair_coeff I I *COMMON rho0 c0 Q1 Cp hg keyword args ... *END
   pair_coeff I J *CROSS

* I, J = atom types
* rho0 = reference mass density (mass/volume units)
* c0 = reference speed of sound (velocity units)
* Q1 = linear artificial viscosity coefficient (unitless)
* Cp = ignored, but must be present
* hg = hourglass control coefficient (unitless)
* one or more keyword/args groups must be appended, with exactly one equation of state and at most one strength or viscosity model
* keyword = *\*EOS_LINEAR* or *\*EOS_TAIT* or *\*EOS_PERFECT_GAS* or *\*STRENGTH_LINEAR* or *\*STRENGTH_LINEAR_PLASTIC* or *\*VISCOSITY_NEWTON* or *\*ARTIFICIAL_PRESSURE* or *\*ARTIFICIAL_STRESS*

  .. parsed-literal::

       *\*EOS_LINEAR* args = none
       *\*EOS_TAIT* args = n
         n = Tait exponent (unitless)
       *\*EOS_PERFECT_GAS* args = gamma
         gamma = heat capacity ratio (unitless)
       *\*STRENGTH_LINEAR* args = G
         G = shear modulus (pressure units)
       *\*STRENGTH_LINEAR_PLASTIC* args = G sigma_y H
         G = shear modulus (pressure units)
         sigma_y = initial yield stress (pressure units)
         H = hardening parameter (pressure units)
       *\*VISCOSITY_NEWTON* args = mu
         mu = dynamic viscosity (mass/(distance time) units)
       *\*ARTIFICIAL_PRESSURE* args = amplitude
         amplitude = strength of the artificial pressure (unitless)
       *\*ARTIFICIAL_STRESS* args = amplitude
         amplitude = strength of the artificial stress (unitless)

* *\*END* = keyword that terminates the list of parameters

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smd/ulsph *DENSITY_CONTINUITY *VELOCITY_GRADIENT *NO_GRADIENT_CORRECTION
   pair_coeff 1 1 *COMMON 1000.0 10.0 1.0 1.0 0.0 &
              *EOS_TAIT 7.0 &
              *END

   pair_style smd/ulsph *DENSITY_CONTINUITY *VELOCITY_GRADIENT *NO_GRADIENT_CORRECTION
   pair_coeff 1 1 *COMMON 1000.0 10.0 1.0 1.0 0.0 &
              *EOS_TAIT 7.0 &
              *END
   pair_coeff 2 2 *COMMON 1000.0 10.0 1.0 1.0 0.0 &
              *EOS_TAIT 7.0 &
              *VISCOSITY_NEWTON 0.1 &
              *END
   pair_coeff 1 2 *CROSS

Description
"""""""""""

The *smd/ulsph* style computes particle interactions according to
continuum mechanics constitutive laws and an updated Lagrangian smoothed
particle hydrodynamics (ULSPH) algorithm.  In this formulation, all SPH
expressions are evaluated in the current configuration, so that the
particles interacting with a given particle change during a simulation.
It is best suited for fluids and for solids undergoing arbitrarily large
deformations.  The method and the material models are explained on the
:doc:`Howto machdyn <Howto_machdyn>` page.

The three arguments of the pair style select how the mass density is
computed, whether the velocity gradient is computed, and whether the
kernel gradients are corrected to first order:

- *\*DENSITY_SUMMATION* recomputes the mass density every time step by
  kernel summation, while *\*DENSITY_CONTINUITY* integrates the mass
  density in time from the rate of deformation.
- *\*VELOCITY_GRADIENT* enables the computation of the velocity
  gradient.  This is required for *\*DENSITY_CONTINUITY* and for
  strength and viscosity models.  *\*NO_VELOCITY_GRADIENT* saves
  computation time and is only useful for inviscid fluids.
- *\*GRADIENT_CORRECTION* computes first-order corrected kernel
  gradients.  This makes the method more accurate, but can lead to
  severe numerical instabilities if the neighborhood of a particle is
  ill-defined.  For a general fluid simulation, *\*NO_GRADIENT_CORRECTION*
  is recommended.

In contrast to most other pair styles, the parameters of the
:doc:`pair_coeff <pair_coeff>` command are organized in groups that
start with a keyword.  The *\*COMMON* keyword with its 5 parameters is
mandatory and must come first.  The bulk modulus of the material is
derived from the reference density and speed of sound as :math:`K = c_0^2
\rho_0`.  The *\*COMMON* group must be followed by an equation of state,
which determines the pressure, and optionally by either a strength model
or a viscosity model, which determine the deviatoric (shear) stress.
The list of parameters must be terminated with the *\*END* keyword.
Using the line continuation operator "&" to put each keyword group on
its own line is recommended.  The available models and their equations
are described in the :ref:`material models section
<howto_smd_materials>` of the Howto machdyn page.

The *\*ARTIFICIAL_PRESSURE* and *\*ARTIFICIAL_STRESS* keywords add
Monaghan's artificial pressure or artificial stress terms, which
counteract the clumping of particles under tension (tensile
instability).

Interactions between particles of different types *I* and *J*, e.g.
between a highly viscous and a low viscosity fluid, must be requested
with the *\*CROSS* keyword, after the parameters for both types *I* and
*J* have been set.  The amplitudes of the artificial pressure and
artificial stress for the cross interaction are the averages of the
values for the two types, if both are larger than zero.

Recommended values for *Q1* are between 0.05 and 1.0.  The hourglass
control coefficient *hg* is usually set to 0 for fluids.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No mixing is performed automatically.  Currently, no part of MACHDYN
supports restarting nor minimization.  rRESPA does not apply to this
pair style.

----------

Restrictions
""""""""""""

This pair style is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires the :doc:`atom style smd <atom_style>`, and
the particles interacting with it must be integrated with :doc:`fix
smd/integrate_ulsph <fix_smd_integrate_ulsph>`.

A strength model and a viscosity model cannot be used for the same
particle type.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style smd/tlsph <pair_smd_tlsph>`,
:doc:`pair_style smd/hertz <pair_smd_hertz>`,
:doc:`fix smd/integrate_ulsph <fix_smd_integrate_ulsph>`,
:doc:`fix smd/adjust_dt <fix_smd_adjust_dt>`

Default
"""""""

none
