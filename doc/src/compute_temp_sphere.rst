.. index:: compute temp/sphere

compute temp/sphere command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/sphere keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/sphere = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *bias* or *dof*

  .. parsed-literal::

       *bias* value = bias-ID
         bias-ID = ID of a temperature compute that removes a velocity bias
       *dof* value = *all* or *rotate*
         all = compute temperature of translational and rotational degrees of freedom
         rotate = compute temperature of just rotational degrees of freedom

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all temp/sphere
   compute myTemp mobile temp/sphere bias tempCOM
   compute myTemp mobile temp/sphere dof rotate

Description
"""""""""""

Define a computation that calculates the temperature of a group of
spherical particles, including a contribution from both their
translational and rotational kinetic energy.  This differs from the
usual :doc:`compute temp <compute_temp>` command, which assumes point
particles with only translational kinetic energy.

Both point and finite-size particles can be included in the group.
Point particles do not rotate, so they have only three translational
degrees of freedom.  For 3d spherical particles, each has six degrees of
freedom (three translational, three rotational).  For 2d spherical particles,
each has three degrees of freedom (two translational, one rotational).

.. note::

   This choice for degrees of freedom (DOF) assumes that all finite-size
   spherical particles in your model will freely rotate, sampling all
   their rotational DOF.  It is possible to use a combination of
   interaction potentials and fixes that induce no torque or otherwise
   constrain some of all of your particles so that this is not the case.
   Then there are less DOF and you should use the :doc:`compute_modify
   extra/dof <compute_modify>` command to adjust the DOF accordingly.

The translational kinetic energy is computed the same as is described
by the :doc:`compute temp <compute_temp>` command.  The rotational
kinetic energy is computed as :math:`\frac12 I \omega^2`, where :math:`I` is
the moment of inertia for a sphere and :math:`\omega` is the particle's angular
velocity.

.. note::

   For :doc:`2d models <dimension>`, particles are treated as
   spheres, not disks, meaning their moment of inertia will be the same
   as in 3d.

A kinetic energy tensor, stored as a six-element vector, is also
calculated by this compute.  The formula for the components of the
tensor is the same as the above formulas, except that :math:`v^2` and
:math:`\omega^2` are replaced by :math:`v_x v_y` and :math:`\omega_x
\omega_y` for the :math:`xy` component.  The six components of the
vector are ordered :math:`xx`, :math:`yy`, :math:`zz`, :math:`xy`,
:math:`xz`, :math:`yz`.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

This compute subtracts out translational degrees-of-freedom due to fixes
that constrain molecular motion, such as :doc:`fix shake <fix_shake>`
and :doc:`fix rigid <fix_rigid>`.  This means the temperature of groups
of atoms that include these constraints will be computed correctly.  If
needed, the subtracted degrees of freedom can be altered using the
*extra/dof* option of the :doc:`compute_modify <compute_modify>`
command.

See the :doc:`Howto thermostat <Howto_thermostat>` page for a
discussion of different ways to compute temperature and perform
thermostatting.

----------

The keyword/value option pairs are used in the following ways.

For the *bias* keyword, *bias-ID* refers to the ID of a temperature
compute that removes a "bias" velocity from each atom.  This allows
compute temp/sphere to compute its thermal temperature after the
translational kinetic energy components have been altered in a
prescribed way (e.g., to remove a flow velocity profile).  Thermostats
that use this compute will work with this bias term.  See the doc
pages for individual computes that calculate a temperature and the doc
pages for fixes that perform thermostatting for more details.

For the *dof* keyword, a setting of *all* calculates a temperature
that includes both translational and rotational degrees of freedom.
A setting of *rotate* calculates a temperature that includes only
rotational degrees of freedom.

----------

Output info
"""""""""""

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1--6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS
output options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This fix requires that atoms store torque and angular velocity (omega)
and a radius as defined by the :doc:`atom_style sphere <atom_style>`
command.

All particles in the group must be finite-size spheres, or point
particles with radius = 0.0.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute temp/asphere <compute_temp>`

Default
"""""""

The option defaults are no bias and dof = all.
