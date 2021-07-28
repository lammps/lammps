.. index:: bond_style special

bond_style special command
==========================

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

The *special* bond style can be used to create conceptual bonds which
effectively impose weightings on the pairwise Lennard Jones and/or
Coulombic interactions between selected pairs of particles in the
system.  The form of the pairwise interaction will be whatever is
computed by the :doc:`pair_style <pair_style>` command defined for the
system; this command defines the weightings for its two terms.

This command can thus be useful to apply weightings that cannot be
handled by the :doc:`special_bonds <special_bonds>` command, such as
on 1-5 or 1-6 interactions.  Or it can be used to add pairwise forces
between one or more pairs of atoms that otherwise would not be include
in the :doc:`pair_style <pair_style>` computation.

The potential for this bond style has the form

.. math::

   E =  w_{LJ} E_{LJ} + w_{Coul} E_{Coul}

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`w_{LJ}` weight (0.0 to 1.0) on pairwise Lennard-Jones interactions

* :math:`w_{Coul}` weight (0.0 to 1.0) on pairwise Coulombic interactions

----------

Normally this bond style should be used in conjunction with one (or
more) other bond styles which compute forces between atoms directly
bonded to each other in a molecule.  This means the :doc:`bond_style
hybrid <bond_hybrid>` command should be used with bond_style special
as one of its sub-styles.

Note that the same as for any other bond style, pairs of bonded atoms
must be enumerated in the data file read by the :doc:`read_data
<read_data>` command.  Thus if this command is used to weight all 1-5
interactions in the system, all the 1-5 pairs of atoms must be listed
in the "Bonds" section of the data file.

This bond style imposes strict requirements on settings made with the
:doc:`special_bonds <special_bonds>` command.  These requirements
ensure that the new bonds created by this style do not create spurious
1-2, 1-3, or 1-4 interactions within the molecular topology.

Specifically 1-2 interactions must have weights of zero, 1-3
interactions must either have weights of unity or :doc:`special_bonds
angle yes <special_bonds>` must be used, and 1-4 interactions must
have weights of unity or :doc:`special_bonds dihedral yes <special_bonds>`
must be used.

If this command is used to create bonded interactions between
particles that are further apart than usual (e.g. 1-5 or 1-6
interactions), this style may require an increase in the communication
cutoff via the :doc:`comm_modify cutoff <comm_modify>` command.  If
LAMMPS cannot find a partner atom in a bond, an error will be issued.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

This bond style requires the use of a :doc:`pair_style <pair_style>` which
computes a pairwise additive interaction and provides the ability to
compute interactions for individual pairs of atoms.  Manybody potentials
are not compatible in general, but also some other pair styles are missing
the required functionality and thus will cause an error.

This command is not compatible with long-range Coulombic interactions. If a
`kspace_style <kspace_style>` is declared, an error will be issued.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`special_bonds <special_bonds>`

Default
"""""""

none
