.. index:: compute temp/body

compute temp/body command
=========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp/body keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/body = style name of this compute command
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

.. parsed-literal::

   compute 1 all temp/body
   compute myTemp mobile temp/body bias tempCOM
   compute myTemp mobile temp/body dof rotate

Description
"""""""""""

Define a computation that calculates the temperature of a group of
body particles, including a contribution from both their
translational and rotational kinetic energy.  This differs from the
usual :doc:`compute temp <compute_temp>` command, which assumes point
particles with only translational kinetic energy.

Only body particles can be included in the group.  For 3d particles,
each has 6 degrees of freedom (3 translational, 3 rotational).  For 2d
body particles, each has 3 degrees of freedom (2 translational, 1
rotational).

.. note::

   This choice for degrees of freedom (dof) assumes that all body
   particles in your model will freely rotate, sampling all their
   rotational dof.  It is possible to use a combination of interaction
   potentials and fixes that induce no torque or otherwise constrain some
   of all of your particles so that this is not the case.  Then there are
   less dof and you should use the :doc:`compute_modify extra <compute_modify>` command to adjust the dof accordingly.

The translational kinetic energy is computed the same as is described
by the :doc:`compute temp <compute_temp>` command.  The rotational
kinetic energy is computed as 1/2 I w\^2, where I is the inertia tensor
for the aspherical particle and w is its angular velocity, which is
computed from its angular momentum.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute.  The formula for the components of the
tensor is the same as the above formula, except that v\^2 and w\^2 are
replaced by vx\*vy and wx\*wy for the xy component, and the appropriate
elements of the inertia tensor are used.  The 6 components of the
vector are ordered xx, yy, zz, xy, xz, yz.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

This compute subtracts out translational degrees-of-freedom due to
fixes that constrain molecular motion, such as :doc:`fix shake <fix_shake>` and :doc:`fix rigid <fix_rigid>`.  This means the
temperature of groups of atoms that include these constraints will be
computed correctly.  If needed, the subtracted degrees-of-freedom can
be altered using the *extra* option of the
:doc:`compute_modify <compute_modify>` command.

See the :doc:`Howto thermostat <Howto_thermostat>` doc page for a
discussion of different ways to compute temperature and perform
thermostatting.

----------

The keyword/value option pairs are used in the following ways.

For the *bias* keyword, *bias-ID* refers to the ID of a temperature
compute that removes a "bias" velocity from each atom.  This allows
compute temp/sphere to compute its thermal temperature after the
translational kinetic energy components have been altered in a
prescribed way, e.g. to remove a flow velocity profile.  Thermostats
that use this compute will work with this bias term.  See the doc
pages for individual computes that calculate a temperature and the doc
pages for fixes that perform thermostatting for more details.

For the *dof* keyword, a setting of *all* calculates a temperature
that includes both translational and rotational degrees of freedom.  A
setting of *rotate* calculates a temperature that includes only
rotational degrees of freedom.

----------

**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1-6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the BODY package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This compute requires that atoms store angular momentum and a
quaternion as defined by the :doc:`atom_style body <atom_style>`
command.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`

**Default:** none
