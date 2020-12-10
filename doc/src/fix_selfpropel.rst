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

Related commands
""""""""""""""""

:doc:`fix efield <fix_efield>`

Default
"""""""

none
