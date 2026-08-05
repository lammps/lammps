.. index:: bond_style bpm/peri

bond_style bpm/peri command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/peri keyword value attribute1 attribute2 ...

* zero or more keyword/value pairs may be appended
* keyword = *store/local*

  .. parsed-literal::

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of two atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the two atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the two atoms (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/peri
   bond_coeff 1 pmb 1.6863e22 0.0015001 0.0005 0.25
   bond_coeff 1 lps 14.9e9 14.9e9 0.0015001 0.0005 0.25
   bond_coeff 1 ves 14.9e9 14.9e9 0.5 0.001 0.0015001 0.0005 0.25
   bond_coeff 1 eps 14.9e9 14.9e9 118.43 0.0015001 0.0005 0.25

Description
"""""""""""

.. versionadded:: TBD

The *bpm/peri* bond style implements the four peridynamic constitutive
models of the :doc:`PERI package <pair_peri>` --- bond-based prototype
microelastic brittle (*pmb*), state-based linear peridynamic solid
(*lps*), state-based viscoelastic solid (*ves*), and state-based
elastic-plastic solid (*eps*) --- recast as a bond style within the
:doc:`BPM package <Howto_bpm>`.  Each bond stores its reference
length and breaks individually on a stretch criterion.  The model is
selected by the first :doc:`bond_coeff <bond_coeff>` argument.

Like the other :doc:`BPM bond styles <Howto_bpm>`, the reference state is
stored by each bond when it is first computed in the setup of a run, is
preserved across run commands, and is written to :doc:`binary restart
files <restart>`.

A short-range repulsive contact force between non-bonded near neighbors
is provided by the companion :doc:`pair_style bpm/peri <pair_bpm_peri>`,
which must also be defined (see the Restrictions section).

----------

Per-atom nodal volume
"""""""""""""""""""""""

Every peridynamic node carries a nodal volume :math:`V` (the discretized
material volume it represents).  Unlike the internal bookkeeping fields
the bond style creates automatically, the nodal volume is *input* data
and so it must be declared by the user with a :doc:`fix property/atom
<fix_property_atom>` named *vfrac*, using group *all* and with ghost
communication enabled, **before** the bond style is defined:

.. code-block:: LAMMPS

   fix vol all property/atom d_vfrac ghost yes
   set group all d_vfrac 1.25e-10            # uniform nodal volume

For a body with non-uniform nodal volume the values may instead be read
from a data file (see :doc:`read_data <read_data>` and :doc:`fix
property/atom <fix_property_atom>`).  The nodal mass is supplied
independently, either as a uniform per-type :doc:`mass <mass>` (equal to
density times nodal volume) or, for variable mass, by using
:doc:`atom_style hybrid sphere bond <atom_style>` and setting the
per-atom mass.  See the :doc:`peridynamics Howto <Howto_peri>` for a
complete worked example.

----------

Constitutive models and coefficients
"""""""""""""""""""""""""""""""""""""

The model keyword is the first :doc:`bond_coeff <bond_coeff>` argument,
followed by that model's coefficients:

For *pmb* (bond-based prototype microelastic brittle):

* *c* (energy/distance/volume\^2 units)
* horizon :math:`\delta` (distance units)
* s00 (unitless)
* :math:`\alpha` (unitless)

The bond force is :math:`F = c\, s\, V_j` with stretch :math:`s = (r -
r_0)/r_0`, where *c* is the micromodulus :math:`c = 18 K / (\pi
\delta^4)` for bulk modulus *K*.

For *lps* (state-based linear peridynamic solid):

* *K* (force/area units), bulk modulus
* *G* (force/area units), shear modulus
* horizon :math:`\delta` (distance units)
* s00 (unitless)
* :math:`\alpha` (unitless)

For *ves* (state-based viscoelastic solid), as *lps* plus the
viscoelastic parameters inserted after *G*:

* *K*, *G* (force/area units)
* :math:`\lambda` (unitless), relaxation parameter in [0,1]
* :math:`\tau` (time units), relaxation time
* horizon, s00, :math:`\alpha`

For small :math:`\lambda` the response approaches the linear-elastic
*lps* model.

For *eps* (state-based elastic-plastic solid), as *lps* with the yield
stress inserted after *G*:

* *K*, *G* (force/area units)
* yield stress (force/area units)
* horizon, s00, :math:`\alpha`

The *lps*, *ves*, and *eps* models are state-based: the force on a bond
depends on the dilatation :math:`\theta` and weighted volume of *both*
endpoints, which the bond style accumulates and communicates internally
each step.  The canonical references are :ref:`(Silling 2000)
<Silling2000-bpm>`, :ref:`(Silling 2007) <Silling2007-bpm>`, and
:ref:`(Parks) <Parks-bpm>`.  The *ves* and *eps* formulations are from
:ref:`(Mitchell) <Mitchell2011-bpm>` and :ref:`(Mitchell2)
<Mitchell2011a-bpm>`; the underlying state-based viscoplasticity theory is
:ref:`(Foster 2010) <Foster2010-bpm>`, and the original PDLAMMPS
implementations of these two models by Rahman and Foster are described in
:ref:`(Rahman) <Rahman-bpm>`.

.. note::

   The *eps* plastic deviatoric extension is updated by a radial return
   onto the yield surface that is consistent with the returned force
   state.  This differs from the update in the legacy :doc:`pair_style
   peri/eps <pair_peri>`, which is only conditionally stable under
   sustained plastic flow; the yield-surface norm here also uses the same
   deviatoric force state the model applies.

----------

Bond breaking criterion
"""""""""""""""""""""""

A peridynamic bond between particles *i* and *j* breaks irreversibly once
its stretch :math:`s = (r - r_0)/r_0` exceeds a critical stretch.
Following :ref:`(Parks) <Parks-bpm>` (eq. 9), the per-bond critical
stretch is :math:`s_{00} - \alpha\, \max(s_{min,i}, s_{min,j})`, where
:math:`s_{min}` is the minimum (most compressive) stretch over a
particle's bonds on the previous step.  The criterion is evaluated *per
bond* using that bond's own s00 and :math:`\alpha`, so type-dependent
coefficients behave as intended (for example a weaker s00 at a material
interface initiates a crack there).  Broken bonds are removed from the
bond topology by setting the bond type to 0.

The :doc:`compute property/atom <compute_property_atom>` *s0* and *smin*
per-atom properties report the per-particle critical stretch and minimum
stretch; for *eps*, *lambda* reports the accumulated plastic multiplier.

For the state-based models (*lps*, *ves*, *eps*) the per-step dilatation
:math:`\theta` can optionally be exposed for visualization: declare a
per-atom *theta* property before the bond style,

.. code-block:: LAMMPS

   fix dil all property/atom d_theta ghost yes

and the bond style writes the dilatation of each owned node into it every
step, readable with :doc:`compute property/atom <compute_property_atom>`
*theta*.  (It must be declared up front, like *vfrac*, so a compute that
references it can be defined from the start of the input.)

Similarly, a volume-weighted bond damage of each peridynamics node will
be calculated by declaring a per-atom *damage* property before the bond
style,

.. code-block:: LAMMPS

   fix dmg all property/atom d_damage

The damage of node *i* is

.. math::

   d_i = 1 - \frac{\sum_{j\,\mathrm{(intact)}} V_j}{\sum_{j\,\mathrm{(initial)}} V_j}

where the numerator sums the nodal volumes :math:`V_j` of the bonds that
are still intact and the denominator is the reference interaction volume
(the same sum over all bonds present in the initial configuration).
Damage is 0 for a fully intact node and approaches 1 as a node loses
bonds, so it is a convenient measure for visualizing cracks and fracture
surfaces.  It is the BPM-framework equivalent of the legacy :doc:`compute
damage/atom <compute_damage_atom>`.

----------

Restart and other info
"""""""""""""""""""""""

This bond style writes the reference state of each bond and its
per-type coefficients to :doc:`binary restart files <restart>`.  Loading
a restart file restores bonds and their reference state.  The reference
state is NOT written to data files.

If the *store/local* option is used, an internal fix records data for
each breaking bond into a local vector or array accessible through a
:doc:`dump local <dump>` command, as for the other :doc:`BPM bond styles
<bond_bpm_spring>`.

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This bond style requires the companion :doc:`pair_style bpm/peri
<pair_bpm_peri>` (or :doc:`pair_style zero <pair_zero>` for a deliberate
contact-free test), a per-atom *vfrac* property as described above, and a
cubic lattice with equal spacing in *x*, *y*, and *z* (used by the
partial-volume horizon correction).

As with the other BPM bond styles, :doc:`newton <newton>` must be set to
*bond off* and the special bond weights must be

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

See the :doc:`peridynamics Howto <Howto_peri>` for a complete input
script and a guide to migrating PERI inputs to the BPM framework.

Related commands
""""""""""""""""

:doc:`pair_style bpm/peri <pair_bpm_peri>`, :doc:`bond_coeff
<bond_coeff>`, :doc:`pair_style peri/pmb <pair_peri>`, :doc:`Howto BPM
<Howto_bpm>`, :doc:`Howto peridynamics <Howto_peri>`

Default
"""""""

none

----------

.. _Parks-bpm:

**(Parks)** Parks, Lehoucq, Plimpton, Silling, Comp Phys Comm, 179(11),
777-783 (2008).

.. _Silling2000-bpm:

**(Silling 2000)** Silling, J Mech Phys Solids, 48, 175-209 (2000).

.. _Silling2007-bpm:

**(Silling 2007)** Silling, Epton, Weckner, Xu, Askari, J Elasticity,
88, 151-184 (2007).

.. _Mitchell2011-bpm:

**(Mitchell)** Mitchell. A non-local, ordinary-state-based
viscoelasticity model for peridynamics. Sandia National Lab Report,
8064:1-28 (2011).

.. _Mitchell2011a-bpm:

**(Mitchell2)** Mitchell. A Nonlocal, Ordinary, State-Based
Plasticity Model for Peridynamics. Sandia National Lab Report,
3166:1-34 (2011).

.. _Foster2010-bpm:

**(Foster 2010)** Foster, Silling, Chen, Viscoplasticity using
peridynamics, Int J Numer Methods Eng, 81(10), 1242-1258 (2010).

.. _Rahman-bpm:

**(Rahman)** Rahman, Foster, Implementation of linear viscoelasticity model
in PDLAMMPS (2013) and Implementation of elastic-plastic model in PDLAMMPS,
technical reports, University of Texas at San Antonio (included in the LAMMPS
distribution as `PDLammps_VES.pdf <PDF/PDLammps_VES.pdf>`_ and
`PDLammps_EPS.pdf <PDF/PDLammps_EPS.pdf>`_; see also the PDLAMMPS overview
`PDLammps_overview.pdf <PDF/PDLammps_overview.pdf>`_).
