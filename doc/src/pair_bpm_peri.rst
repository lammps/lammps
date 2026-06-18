.. index:: pair_style bpm/peri

pair_style bpm/peri command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style bpm/peri

Examples
""""""""

.. code-block:: LAMMPS

   pair_style bpm/peri
   pair_coeff * * 1.6863e22 0.0015001

Description
"""""""""""

.. versionadded:: TBD

Style *bpm/peri* supplies the short-range repulsive contact force
between non-bonded peridynamic nodes for the :doc:`bond_style bpm/peri
<bond_bpm_peri>` model.  It is the BPM-framework equivalent of the
contact term built into the legacy :doc:`peridynamic pair styles
<pair_peri>`.  Bonded pairs are handled entirely by the bond style and
are excluded here through the 1-2 :doc:`special_bonds <special_bonds>`
weight; this pair style therefore acts only between nodes that are *not*
bonded (for example new contacts that form after fracture or large
deformation).

For a pair of non-bonded nodes within the contact distance
:math:`d = 1.35\, l_c` (where :math:`l_c` is the lattice spacing) the
force has magnitude

.. math::

   F = 15\, k\, \bar{V}\, \frac{r - d}{\delta}

where *k* is the contact stiffness, :math:`\delta` the horizon,
:math:`r` the current distance, and :math:`\bar{V} = (V_i + V_j)/2` the
mean nodal volume of the two nodes (read from the per-atom *vfrac*
property declared for the bond style).  The factor 15 follows the EMU
Theory Manual convention used by the legacy peridynamic styles.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the example above,
or in the restart files read by the :doc:`read_restart <read_restart>`
command:

* *k* (energy/distance/volume\^2 units), contact stiffness
* horizon :math:`\delta` (distance units)

Typically *k* and the horizon are set equal to the corresponding
:doc:`bond_coeff <bond_bpm_peri>` values.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Coefficients for all I,J pairs
must be specified explicitly.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, or tail options.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires a per-atom *vfrac* property (a :doc:`fix
property/atom <fix_property_atom>` *d_vfrac* with ghost communication)
and a cubic lattice with equal spacing in *x*, *y*, and *z*.  It is
intended to be used together with :doc:`bond_style bpm/peri
<bond_bpm_peri>`; see the :doc:`peridynamics Howto <Howto_peri>`.

Related commands
""""""""""""""""

:doc:`bond_style bpm/peri <bond_bpm_peri>`, :doc:`pair_coeff
<pair_coeff>`, :doc:`Howto peridynamics <Howto_peri>`

Default
"""""""

none
