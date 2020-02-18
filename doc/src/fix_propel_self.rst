.. index:: fix propel/self

fix propel/self command
=======================

Syntax
""""""

fix ID group-ID propel/self mode magnitude keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* propel/self = style name of this fix command
* mode = velocity or quat
* magnitude = magnitude of the active force
* one or more keyword/value pairs may be appended to args
* keyword = *types*

  *types* values = one or more atom types



Examples
""""""""


.. parsed-literal::

   fix active_group all propel/self velocity 1.0
   fix constant_velocity all viscous 1.0

   fix active_group all propel/self quat 1.0

   fix active all propel/self quat 1.0 types 1 2 4

Description
"""""""""""

Adds a force of a constant magnitude to each atom in the group. The nature in
which the force is added depends on the mode.

For *mode* = *velocity*, the active force acts along the velocity vector of
each atom. This can be interpreted as a velocity-dependent friction,
such as proposed by :ref:`(Erdmann) <Erdmann>`.

For *mode* = *quat* the force is applied along the axis obtained
by rotating the x-axis along the atom's quaternion. In other words, the
force is along the x-axis in the atom's body frame. This mode requires
all atoms in the group to have a quaternion, so atom\_style should
either be ellipsoid or body.  In combination with Langevin thermostat
for translation and rotation in the overdamped regime, the quaternion
mode corresponds to the active Brownian particle model introduced by
:ref:`(Henkes) <Henkes>`, :ref:`(Bialke) <Bialke>` and :ref:`(Fily)
<Fily>`.

By default, this fix is applied to all atoms in the group. You can
override this behavior by specifying the atom types the fix should work
on through the *types* keyword.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix is not imposed  during minimization.

Restrictions
""""""""""""


In quat mode, this fix makes use of per-atom quaternions to take
into account the fact that the orientation can rotate and hence the
direction of the active force can change. The quat mode
of this fix only works with atom\_style ellipsoid.

Related commands
""""""""""""""""

:doc:`fix setforce <fix_setforce>`, :doc:`fix addforce <fix_addforce>`

.. _Erdmann:



**(Erdmann)** U. Erdmann , W. Ebeling, L. Schimansky-Geier, and F. Schweitzer,
Eur. Phys. J. B 15, 105-113, 2000.

.. _Henkes:



**(Henkes)** Henkes, S, Fily, Y., and Marchetti, M. C. Phys. Rev. E, 84, 040301(R), 2011.

.. _Bialke:



**(Bialke)** J. Bialke, T. Speck, and H Loewen, Phys. Rev. Lett. 108, 168301, 2012.

.. _Fily:



**(Fily)** Y. Fily and M.C. Marchetti, Phys. Rev. Lett. 108, 235702, 2012.

**Default:** types
