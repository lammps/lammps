.. index:: pair_style rebomos
.. index:: pair_style rebomos/omp

pair_style rebomos command
==========================

Accelerator Variants: *rebomos/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rebomos

* rebomos = name of this pair style

Examples
""""""""

.. code-block:: LAMMPS

   pair_style rebomos
   pair_coeff * * ../potentials/MoS.rebomos Mo S

Example input scripts available: examples/threebody/

Description
"""""""""""

.. versionadded:: TBD

The *rebomos* pair style computes the interactions between molybdenum
and sulfur atoms :ref:`(Stewart) <Stewart>` utilizing an adaptive
interatomic reactive empirical bond order potential that is similar in
form to the AIREBO potential :ref:`(Stuart) <Stuart2>`.  The potential
is based on an earlier parameterizations for :math:`\text{MoS}_2`
developed by :ref:`(Liang) <Liang>`.

The REBOMoS potential consists of two terms:

.. math::

   E & = \frac{1}{2} \sum_i \sum_{j \neq i}
   \left[ E^{\text{REBO}}_{ij} + E^{\text{LJ}}_{ij}  \right] \\

The :math:`E^{\text{REBO}}` term describes the covalently bonded
interactions between Mo and S atoms while the :math:`E^{\text{LJ}}` term
describes longer range dispersion forces between layers.  A cubic spline
function is applied to smoothly switch between covalent bonding at short
distances to dispersion interactions at longer distances. This allows
the model to capture bond formation and breaking events which may occur
between adjacent MoS2 layers, edges, defects, and more.

----------

Only a single pair_coeff command is used with the *rebomos* pair style
which specifies an REBOMoS potential file with parameters for Mo and S.
These are mapped to LAMMPS atom types by specifying N additional
arguments after the filename in the pair_coeff command, where N is the
number of LAMMPS atom types:

* filename
* :math:`N` element names = mapping of REBOMoS elements to atom types

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the potential file.

As an example, if your LAMMPS simulation has three atom types and you want
the first two to be Mo, and the third to be S, you would use the following
pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * MoS.rebomos Mo Mo S

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first two Mo arguments map LAMMPS atom types 1 and 2 to the Mo
element in the REBOMoS file.  The final S argument maps LAMMPS atom type
3 to the S element in the REBOMoS file.  If a mapping value is specified
as NULL, the mapping is not performed.  This can be used when a
*rebomos* potential is used as part of the *hybrid* pair style.  The
NULL values are placeholders for atom types that will be used with other
potentials.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write their information to :doc:`binary restart
files <restart>`, since it is stored in potential files.  Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

These pair potentials require the :doc:`newton <newton>` setting to be
"on" for pair interactions.

The MoS.rebomos potential file provided with LAMMPS (see the potentials
directory) is parameterized for metal :doc:`units <units>`.  You can use
the *rebomos* pair style with any LAMMPS units setting, but you would
need to create your own REBOMoS potential file with coefficients listed
in the appropriate units.

The pair style provided here **only** supports potential files parameterized
for the elements molybdenum and sulfur (designated with "Mo" and "S" in the
*pair_coeff* command.  Using potential files for other elements will trigger
an error.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair style rebo <pair_airebo>`

Default
"""""""

none

----------

.. _Stewart:

**(Steward)**  Stewart, Spearot, Modelling Simul. Mater. Sci. Eng. 21, 045003, (2013).

.. _Stuart2:

**(Stuart)** Stuart, Tutein, Harrison, J Chem Phys, 112, 6472-6486, (2000).

.. _Liang:

**(Liang)**  Liang, Phillpot, Sinnott Phys. Rev. B79 245110, (2009), Erratum: Phys. Rev. B85 199903(E), (2012)
