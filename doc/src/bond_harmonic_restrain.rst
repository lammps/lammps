.. index:: bond_style harmonic/restrain

bond_style harmonic/restrain command
====================================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style harmonic/restrain

Examples
""""""""

.. code-block:: LAMMPS

   bond_style harmonic
   bond_coeff 5 80.0

Description
"""""""""""

.. versionadded:: TBD

The *harmonic/restrain* bond style uses the potential

.. math::

   E = K (r - r_{t=0})^2

where :math:`r_{t=0}` is the distance between the bonded atoms at the
beginning of the first :doc:`run <run>` or :doc:`minimize <minimize>`
command after the bond style has been defined (*t=0*).  Note that the
usual 1/2 factor is included in :math:`K`.  This will effectively
restrain bonds to their initial length, whatever that is.  This is where
this bond style differs from :doc:`bond style harmonic <bond_harmonic>`
where the bond length is set through the per bond type coefficients.

The following coefficient must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands

* :math:`K` (energy/distance\^2)

This bond style differs from other options to add harmonic restraints
like :doc:`fix restrain <fix_restrain>` or :doc:`pair style list
<pair_list>` or :doc:`fix colvars <fix_colvars>` in that it requires a
bond topology, and thus the defined bonds will trigger exclusion of
special neighbors from the neighbor list according to the
:doc:`special_bonds <special_bonds>` settings.

Restart info
""""""""""""

This bond style supports the :doc:`write_restart <write_restart>` and
:doc:`read_restart <read_restart>` commands. The state of the initial
bond lengths is stored with restart files and read back.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>`
page for more info.

This bond style does NOT write its coefficients to a data file since the
:math:`r_{t=0}` values are for individual bonds and not bond types.

Since it uses :doc:`fix property/atom <fix_property_atom>` internally,
you should use the "nofix" argument to the :doc:`write_data command
<write_data>` to avoid writing a section to the data file that cannot be
read back by the bond style.

This bond style cannot be used with :doc:`fix shake or fix rattle
<fix_shake>`, with :doc:`fix filter/corotate <fix_filter_corotate>`, or
any :doc:`tip4p pair style <pair_lj_cut_tip4p>` since there is no specific
equilibrium distance for a given bond type.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`,  :doc:`bond_harmonic <bond_harmonic>`,
:doc:`fix restrain <fix_restrain>`, :doc:`pair style list <pair_list>`

Default
"""""""

none
