.. index:: compute temp/uvt

compute temp/uvt command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/uvt fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/uvt = style name of this compute command
* fix-ID = ID of the :doc:`fix uvt <fix_uvt>` supplying the electronic state

Examples
""""""""

.. code-block:: LAMMPS

   variable dEdN equal 0.0
   fix cp all uvt temp 300 300 100 mu 0 0 100 ne 1 ne_velocity 0 dedn v_dEdN
   compute combined all temp/uvt cp
   thermo_style custom step f_cp[1] temp f_cp[2] c_combined

Description
"""""""""""

Calculate a combined temperature for the atomic velocities and the
single electron-number velocity controlled by :doc:`fix uvt <fix_uvt>`:

.. math::

   T = \frac{2 K_\mathrm{nuc} + W_\mathrm{e}\dot{N}_\mathrm{e}^{\,2}}
            {(g+1) k_B}.

Here :math:`g` is the nuclear number of degrees of freedom, determined
as for :doc:`compute temp <compute_temp>`, and :math:`K_\mathrm{nuc}` is
the nuclear kinetic energy.  The electronic mass :math:`W_\mathrm{e}` and
velocity :math:`\dot{N}_\mathrm{e}` are read from the referenced fix.
The reported number of degrees of freedom is :math:`g+1`.

The nuclear contribution uses the same per-type or per-atom masses,
group selection, unit conversion, and MPI reduction as compute temp.
The replicated electronic contribution is added once after the nuclear
reduction, independent of the number of MPI ranks.

The six-element vector is inherited from compute temp and contains only
the nuclear kinetic-energy tensor, in the order xx, yy, zz, xy, xz, yz.
The electron-number velocity has no Cartesian direction and is not
included in that tensor.  Thus its trace is not the numerator of the
combined scalar temperature.

The *extra/dof* and *dynamic/dof* options of
:doc:`compute_modify <compute_modify>` apply to the nuclear DOF calculation;
the one electronic DOF is added afterward.  Use *dynamic/dof yes* when
the atom count changes during a run.  Nuclear DOF must not be negative
for a nonempty group, and combined DOF must be positive.

Output info
"""""""""""

This compute calculates an intensive global scalar in temperature units
and an extensive six-element global vector in energy units.  The vector
uses the compute temp convention without the factor of one half.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  Using it also requires
the EXTRA-FIX package for the referenced fix uvt, which must use the same
group.  It must exist when the run is initialized;
the electronic mass must be initialized before the combined temperature
is evaluated.  A normal run setup initializes that mass in fix uvt.

Fix uvt creates a compute of this style automatically, with ID
*fix-ID_temp*, and uses it for the combined thermostat.  An alternative
temp/uvt compute can be selected with :doc:`fix_modify temp <fix_modify>`
if it references the same UVT fix and uses the same group.

Fix uvt also exposes this temperature directly as its second vector entry
(``f_cp[2]`` for fix ID ``cp``), so reporting it does not require an explicit
compute command. The scalar can alternatively report the combined temperature.  Other
atomic thermostats do not scale the electron-number velocity and must
not use this combined temperature to control the coupled system.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`fix uvt <fix_uvt>`,
:doc:`compute_modify <compute_modify>`

Default
"""""""

none
