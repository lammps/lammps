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
is the dipole orientation (unit vector) of particle *i*.

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
