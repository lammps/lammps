.. index:: fix selfpropel

fix selfpropel command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID selfpropel selfpropulsionforce

* ID, group-ID are documented in :doc:`fix <fix>` command
* selfpropel = style name of this fix command



Examples
""""""""

.. code-block:: LAMMPS

   fix selfpropel all 40.0

Description
"""""""""""

Add a force to each atom in the group due to an self-propulsion. The
force is given by

.. math::

   F_i = f_P e_i

where *i* is the particle the force is being applied to, :math:`f_P`
is the self-propulsion force (selfpropulsionforce), and :math:`e_i`
is the dipole orientation of particle *i*.

.. note::

   If another command changes the magnitude of the dipole, this force will
   change accordingly (since :math:`|e_i|` will change, which is physically
   equivalent to re-scaling :math:`f_P` while keeping :math:`|e_i|` constant),
   and no warning will be provided by LAMMPS. This is almost never what you
   want, so ensure you aren't changing dipole magnitudes with another LAMMPS
   fix or pair style. Furthermore, self-propulsion forces (almost) always
   set :math:`e_i`  to be a unit vector for all times, so it's best to set
   all the dipole magnitudes to 1.0 unless you have a good reason not to
   (see the :doc:`set <set>` command on how to do this).

Along with adding a force contribution, this fix defaults to adding
a contribution to the pressure :math:`f_P \sum_i <e_i . r_i>/(d V)`,
where :math:`r_i` is the *unwrapped* coordinate of particle i in
the case of periodic boundary conditions.
See :ref:`(Winkler) <Winkler1>` for a discussion of this active
pressure contribution.

.. note::

   In contrast to equilibrium systems, pressure of active systems
   in general depends on the geometry of the container.
   The active pressure contribution as calculated in this fix
   is only valid for certain boundary conditions (spherical
   walls, rectangular walls, or periodic boundary conditions).
   For other geometries, the pressure must be measured via
   explicit calculation of the force per unit area on a wall,
   and so one must not calculate it using this fix.
   (Use :doc:`fix_modify <fix_modify>` as described below
   to turn off the virial contribution of this fix). Again,
   see :ref:`(Winkler) <Winkler1>` for discussion of why this
   is the case.
   
   Furthermore, when dealing with active systems, the temperature
   is no longer well defined. Therefore, one should ensure that
   the *virial* flag is used if computing pressure with the
   :doc:`compute pressure` command (turning off temperature
   contributions).

This command is extremely similar to :doc:`fix propel/self <fix_propel_self>`
(when the latter has the *quat* keyword invoked). There are two differences
between this command and :doc:`fix propel/self <fix_propel_self>`. Firstly,
this command uses dipoles to indicate the direction of self-propulsion (vs
ellipsoid quaternions). There is no real advantage/disadvantage of using
dipole or quaternions that the author can see. The second difference is
that this command allows for the calculation of the active pressure. This
is the reason why this new command exists, as it is useful to be able
to compute the pressure of the system.
   
----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to the
system's virial as part of :doc:`thermodynamic output <thermo_style>`.
The default is *virial yes*


No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.


Restrictions
""""""""""""

This fix only works when the DIPOLE package is enabled.
See the :doc:`Build package <Build_package>` doc page for more info.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.


Related commands
""""""""""""""""

:doc:`fix propel/self <fix_propel_self>`, :doc:`fix efield <fix_efield>` ,
:doc:`compute pressure`     

Default
"""""""

none

----------

.. _Winkler1:

**(Winkler)** Winkler, Wysocki, and Gompper, Soft Matter, 11, 6680 (2015).
