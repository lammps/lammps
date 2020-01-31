.. index:: bond_style special

bond_style special command
=================================

Syntax
""""""


.. code-block:: LAMMPS

   bond_style special

Examples
""""""""


.. code-block:: LAMMPS

   bond_style special
   bond_coeff 0.5 0.5

Description
"""""""""""

The *special* bond style can be used to impose weighted Lennard Jones and/or
Coulombic interactions on any two particles in the system. It can be used for
cases that cannot be handled in :doc:`special_bonds <special_bonds>`, such as
1-5 interactions. It is a potential of the form:

.. math::

   E =  w_{LJ} E_{LJ} + w_{Coul}E_{Coul}

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`w_{LJ}` weight (0.0 to 1.0) on pairwise Lennard-Jones interactions

* :math:`w_{Coul}` weight (0.0 to 1.0) on pairwise Coulombic interactions

----------

This style has strict requirements on the :doc:`special_bonds <special_bonds>`
setting. 1-2 interactions must have weights of zero. 1-3 interactions must
either have weights of one or the *angle* setting must be turned on. 1-4
interactions must have weights of one or the *dihedral* setting must be turned
on. These requirements ensure that the new bonds created by this style do not
create spurious 1-2, 1-3 or 1-4 interactions.

This style should be used in conjunction with a regular bond style via
:doc:`bond_style hybrid <bond_hybrid>`. Since it can be used to create
bonded interactions between particles that are further away than usual
(e.g. 1-5 or 1-6 interactions), this style may require an increase in the
communication cutoff via the :doc:`neigh_modify <neigh_modify>` command.


Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`special_bonds <special_bonds>`

**Default:** none
