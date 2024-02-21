.. index:: compute pressure

compute pressure command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID pressure temp-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure = style name of this compute command
* temp-ID = ID of compute that calculates temperature, can be NULL if not needed
* zero or more keywords may be appended
* keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or *improper* or *kspace* or *fix* or *virial* or *pair/hybrid*

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pressure thermo_temp
   compute 1 all pressure NULL pair bond
   compute 1 all pressure NULL pair/hybrid lj/cut

Description
"""""""""""

Define a computation that calculates the pressure of the entire system
of atoms.  The specified group must be "all".  See the
:doc:`compute stress/atom <compute_stress_atom>` command if you want per-atom
pressure (stress).  These per-atom values could be summed for a group
of atoms via the :doc:`compute reduce <compute_reduce>` command.

The pressure is computed by the formula

.. math::

   P = \frac{N k_B T}{V} + \frac{1}{V d}\sum_{i=1}^{N'} \vec r_i \cdot \vec f_i

where *N* is the number of atoms in the system (see discussion of DOF
below), :math:`k_B` is the Boltzmann constant, :math:`T` is the
temperature, *d* is the dimensionality of the system (2 for 2d, 3 for
3d), and *V* is the system volume (or area in 2d).  The second term is
the virial, equal to :math:`-dU/dV`, computed for all pairwise as well
as 2-body, 3-body, 4-body, many-body, and long-range interactions, where
:math:`\vec r_i` and :math:`\vec f_i` are the position and force vector
of atom *i*, and the dot indicates the dot product (scalar product).
This is computed in parallel for each subdomain and then summed over
all parallel processes. Thus :math:`N'` necessarily includes atoms from
neighboring subdomains (so-called ghost atoms) and the position and
force vectors of ghost atoms are thus included in the summation.  Only
when running in serial and without periodic boundary conditions is
:math:`N' = N` the number of atoms in the system.  :doc:`Fixes <fix>`
that impose constraints (e.g., the :doc:`fix shake <fix_shake>` command)
may also contribute to the virial term.

A symmetric pressure tensor, stored as a 6-element vector, is also
calculated by this compute.  The six components of the vector are
ordered :math:`xx,` :math:`yy,` :math:`zz,` :math:`xy,` :math:`xz,` :math:`yz.`
The equation for the :math:`(I,J)` components (where :math:`I` and :math:`J`
are :math:`x`, :math:`y`, or :math:`z`) is similar to the above formula,
except that the first term uses components of the kinetic energy tensor and the
second term uses components of the virial tensor:

.. math::

   P_{IJ} = \frac{1}{V}\sum_{k=1}^{N} m_k v_{k_I} v_{k_J} +
   \frac{1}{V}\sum_{k=1}^{N'} r_{k_I} f_{k_J}.

If no extra keywords are listed, the entire equations above are
calculated.  This includes a kinetic energy (temperature) term and the
virial as the sum of pair, bond, angle, dihedral, improper, kspace
(long-range), and fix contributions to the force on each atom.  If any
extra keywords are listed, then only those components are summed to
compute temperature or ke and/or the virial.  The *virial* keyword
means include all terms except the kinetic energy *ke*\ .

The *pair/hybrid* keyword means to only include contribution
from a sub-style in a *hybrid* or *hybrid/overlay* pair style.

Details of how LAMMPS computes the virial efficiently for the entire
system, including for many-body potentials and accounting for the
effects of periodic boundary conditions are discussed in
:ref:`(Thompson) <Thompson1>`.

The temperature and kinetic energy tensor is not calculated by this
compute, but rather by the temperature compute specified with the
command.  If the kinetic energy is not included in the pressure, than
the temperature compute is not used and can be specified as NULL.
Normally the temperature compute used by compute pressure should
calculate the temperature of all atoms for consistency with the virial
term, but any compute style that calculates temperature can be used
(e.g., one that excludes frozen atoms or other degrees of freedom).

Note that if desired the specified temperature compute can be one that
subtracts off a bias to calculate a temperature using only the thermal
velocity of the atoms (e.g., by subtracting a background streaming
velocity).
See the doc pages for individual :doc:`compute commands <compute>` to determine
which ones include a bias.

Also note that the :math:`N` in the first formula above is really
degrees-of-freedom divided by :math:`d` = dimensionality, where the DOF value
is calculated by the temperature compute.
See the various :doc:`compute temperature <compute>` styles for details.

A compute of this style with the ID of thermo_press is created when
LAMMPS starts up, as if this command were in the input script:

.. code-block:: LAMMPS

   compute thermo_press all pressure thermo_temp

where thermo_temp is the ID of a similarly defined compute of style
"temp".  See the :doc:`thermo_style <thermo_style>` command for more details.

----------

.. include:: accel_styles.rst

----------

Output info
"""""""""""

This compute calculates a global scalar (the pressure) and a global
vector of length 6 (pressure tensor), which can be accessed by indices
1--6.  These values can be used by any command that uses global scalar
or vector values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

The ordering of values in the symmetric pressure tensor is as follows:
:math:`p_{xx},` :math:`p_{yy},` :math:`p_{zz},` :math:`p_{xy},`
:math:`p_{xz},` :math:`p_{yz}.`

The scalar and vector values calculated by this compute are
"intensive".  The scalar and vector values will be in pressure
:doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute stress/atom <compute_stress_atom>`,
:doc:`thermo_style <thermo_style>`, :doc:`fix numdiff/virial <fix_numdiff_virial>`,

Default
"""""""

By default the compute includes contributions from the keywords:
``ke pair bond angle dihedral improper kspace fix``

----------

.. _Thompson1:

**(Thompson)** Thompson, Plimpton, Mattson, J Chem Phys, 131, 154107 (2009).
