.. index:: pair_style smd/tlsph

pair_style smd/tlsph command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style smd/tlsph keyword value

* zero or one keyword/value pair may be appended
* keyword = *\*UPDATE_CONSTANT* or *\*UPDATE_PAIRWISE*

  .. parsed-literal::

       *\*UPDATE_CONSTANT* value = dmax
         dmax = update the reference configuration when the change of a pair distance exceeds dmax (distance units)
       *\*UPDATE_PAIRWISE* value = ratio
         ratio = update the reference configuration when the change of a pair distance divided by the kernel range exceeds ratio (unitless)

.. code-block:: LAMMPS

   pair_coeff I I *COMMON rho0 E nu Q1 Q2 hg Cp keyword args ... *END

* I = atom type (both types must be the same)
* rho0 = reference mass density (mass/volume units)
* E = Young's modulus (pressure units)
* nu = Poisson ratio (unitless)
* Q1 = linear artificial viscosity coefficient (unitless)
* Q2 = quadratic artificial viscosity coefficient (unitless)
* hg = hourglass control coefficient (unitless)
* Cp = ignored, but must be present
* zero or more keyword/args groups may be appended, with at most one equation of state and one strength model
* keyword = *\*EOS_NONE* or *\*EOS_LINEAR* or *\*EOS_SHOCK* or *\*EOS_POLYNOMIAL* or *\*STRENGTH_LINEAR* or *\*STRENGTH_LINEAR_PLASTIC* or *\*JOHNSON_COOK* or *\*FAILURE_MAX_PLASTIC_STRAIN* or *\*FAILURE_MAX_PAIRWISE_STRAIN* or *\*FAILURE_ENERGY_RELEASE_RATE*

  .. parsed-literal::

       *\*EOS_NONE* args = none
       *\*EOS_LINEAR* args = none
       *\*EOS_SHOCK* args = c0 S Gamma
         c0 = reference speed of sound (velocity units)
         S = slope of the linear relation between shock and particle velocity (unitless)
         Gamma = Grueneisen parameter (unitless)
       *\*EOS_POLYNOMIAL* args = C0 C1 C2 C3 C4 C5 C6
         C0, C1, C2, C3 = polynomial coefficients (pressure units)
         C4, C5, C6 = ignored, but must be present
       *\*STRENGTH_LINEAR* args = none
       *\*STRENGTH_LINEAR_PLASTIC* args = sigma_y H
         sigma_y = initial yield stress (pressure units)
         H = hardening parameter (pressure units)
       *\*JOHNSON_COOK* args = A B a C epdot0 T0 Tmelt m
         A = initial yield stress (pressure units)
         B = strain hardening factor (pressure units)
         a = strain hardening exponent (unitless)
         C = strain rate factor (unitless)
         epdot0 = reference plastic strain rate (1/time units)
         T0, Tmelt, m = ignored, but must be present
       *\*FAILURE_MAX_PLASTIC_STRAIN* args = epsmax
         epsmax = equivalent plastic strain at which a particle fails (unitless)
       *\*FAILURE_MAX_PAIRWISE_STRAIN* args = epsmax
         epsmax = pairwise strain at which the degradation of a pair interaction starts (unitless)
       *\*FAILURE_ENERGY_RELEASE_RATE* args = Gc
         Gc = critical energy release rate (energy/area units)

* *\*END* = keyword that terminates the list of parameters

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smd/tlsph
   pair_coeff 1 1 *COMMON 1.0 1.0 0.3 0.06 0.0 10.0 1.0 &
              *STRENGTH_LINEAR &
              *EOS_LINEAR &
              *END

   pair_style smd/tlsph
   pair_coeff 1 1 *COMMON 1.0 70.0 0.3 0.56 0.0 10.0 1.0 &
              *EOS_POLYNOMIAL 0.0 74.2 60.5 36.5 1.96 0.0 0.0 &
              *JOHNSON_COOK 0.3241 0.1138 0.42 0.0 1.0e-3 298.15 500.15 1.34 &
              *FAILURE_MAX_PLASTIC_STRAIN 1.2 &
              *END

   pair_style hybrid/overlay smd/tlsph smd/hertz 1.0
   pair_coeff 1 1 smd/tlsph *COMMON 1.0 1.0 0.3 0.06 0.0 10.0 1.0 &
              *STRENGTH_LINEAR_PLASTIC 0.051 0.0 &
              *EOS_LINEAR &
              *END

Description
"""""""""""

The *smd/tlsph* style computes particle interactions according to
continuum mechanics constitutive laws and a total Lagrangian smoothed
particle hydrodynamics (TLSPH) algorithm.  In this formulation, all SPH
expressions are evaluated in the undeformed reference configuration of
the material, which avoids the tensile instability of SPH.  It is best
suited for solids.  The method, including the hourglass control and the
material models, is explained on the :doc:`Howto machdyn
<Howto_machdyn>` page.

By default, the reference configuration is the initial configuration
and is never changed.  With the *\*UPDATE_CONSTANT* or
*\*UPDATE_PAIRWISE* keyword, the reference configuration is updated to
the current configuration whenever the change of the distance of any
pair of interacting particles relative to their reference distance
exceeds *dmax*, or exceeds the fraction *ratio* of the kernel range,
respectively.  This is needed when the neighbor relations of the
particles change, e.g. due to plastic flow.  Please see the
:doc:`Howto machdyn <Howto_machdyn>` page for what happens during an
update.  A warning is printed if *dmax* exceeds the communication
cutoff.

In contrast to most other pair styles, the parameters of the
:doc:`pair_coeff <pair_coeff>` command are organized in groups that
start with a keyword.  The *\*COMMON* keyword with its 7 parameters is
mandatory and must come first.  From Young's modulus and the Poisson
ratio, the shear modulus, the bulk modulus, and the speed of sound of
the material are derived, which are also used by the equations of state
and strength models.  The *\*COMMON* group can be followed by keywords
for an equation of state, which determines the pressure, a strength
model, which determines the deviatoric (shear) stress, and failure
models.  The list of parameters must be terminated with the *\*END*
keyword.  Using the line continuation operator "&" to put each keyword
group on its own line is recommended.  If no equation of state or no
strength model is given, the pressure or the deviatoric stress is zero,
respectively.  The available models and their equations are described
in the :ref:`material models section <howto_smd_materials>` of the
Howto machdyn page.

Recommended values for *Q1* are between 0.05 and 1.0 and for *hg*
between about 1 and 20.

Only interactions between particles of the same type can be defined,
i.e. *I* and *J* must be the same.  Particles of different types or
different bodies can interact through contact forces computed with the
:doc:`smd/hertz <pair_smd_hertz>` pair style, which can be combined with
this pair style through :doc:`pair_style hybrid/overlay <pair_hybrid>`.

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

This pair style requires the :doc:`atom style smd <atom_style>`, an
atom map (see the :doc:`atom_modify <atom_modify>` command), and the
:doc:`newton <newton>` setting *off*.  All particles interacting
with this pair style must be members of a :doc:`group <group>` named
*tlsph*, and they must be integrated with :doc:`fix smd/integrate_tlsph
<fix_smd_integrate_tlsph>`.

The *\*FAILURE_MAX_PAIRWISE_STRAIN* keyword cannot be combined with the
*\*UPDATE_CONSTANT* or *\*UPDATE_PAIRWISE* keywords.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style smd/ulsph <pair_smd_ulsph>`,
:doc:`pair_style smd/hertz <pair_smd_hertz>`,
:doc:`fix smd/integrate_tlsph <fix_smd_integrate_tlsph>`,
:doc:`fix smd/adjust_dt <fix_smd_adjust_dt>`

Default
"""""""

The reference configuration is not updated.
