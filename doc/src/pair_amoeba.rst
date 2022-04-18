.. index:: pair_style amoeba
.. index:: pair_style hippo

pair_style amoeba command
=========================

pair_style hippo command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *amoeba* or *hippo*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style amoeba
   pair_coeff * * protein.prm.amoeba protein.key.amoeba

.. code-block:: LAMMPS

   pair_style hippo
   pair_coeff * * water.prm.hippo water.key.hippo


Additional info
"""""""""""""""

* doc:`Howto amoeba <howto_ameoba>`
* examples/amoeba
* tools/amoeba
* potentials/*.amoeba
* potentials/*.hippo

Description
"""""""""""

The *amoeba* style computes the AMOEBA polarizeable field formulated
by Jay Ponder's group at the U Washington at St Louis :ref:`(Ren)
<amoeba-Ren>`, :ref:`(Shi) <amoeba-Shi>`.  The *hippo* style computes
the HIPPO polarizeable force field, an extension to AMOEBA, formulated
by Josh Rackers and collaborators in the Ponder group :ref:`(Rackers)
<amoeba-Rackers>`.

These force fields can be used when polarization effects are desired
in simulations of water, organic molecules, and biomolecules including
proteins, provided that parameterizations (Tinker PRM force field
files) are available for the systems you are interested in.  Files in
the LAMMPS potentials directory with a "amoeba" or "hippo" suffix can
be used.  The Tinker distribution and website have additional force
field files as well.

As discussed on the doc:`Howto amoeba <howto_ameoba>` doc page, the
intermolecular (non-bonded) portion of the AMOEBA force field contains
these terms:

.. math::

   U_{amoeba} = U_{hal} + U_{multipole} + U_{polar}

while the HIPPO force field contains these terms:

.. math::

   U_{hippo} = U_{repulsion} + U_{dispersion} + U_{multipole} + U_{polar} + U_{qxfer}

Conceptually, these terms compute the following interactions:
 
* :math:`U_{hal}` = buffered 14-7 van der Waals with offsets applied to hydrogen atoms
* :math:`U_{repulsion}` = Pauli repulsion due to rearrangement of electron density
* :math:`U_{dispersion}` = dispersion between correlated, instantaneous induced dipole moments
* :math:`U_{multipole}` = electrostatics between permanent point charges, dipoles, and quadrupoles
* :math:`U_{polar}` = electronic polarization bewteen induced point dipoles
* :math:`U_{qxfer}` = charge transfer effects

Note that the AMOEBA versus HIPPO force fields typically compute the
same term differently using their own formulas.  The references on
this doc page give full details for both force fields.

.. note::

  The AMOEBA and HIPPO force fields compute long-range charge, dipole,
  and quadrupole interactions as well as long-range dispersion
  effects.  However, unlike other models with long-range interactions
  in LAMMPS, this does not require use of a KSpace style via the
  :doc:`kspace_style <kspace_style>` command.  That is because for
  AMOEBA and HIPPO the long-range computations are intertwined with
  the pairwise computations.  So these pair style include both short-
  and long-range computations.  This means the energy and virial
  computed by the pair style as well as the "Pair" timing reported by
  LAMMPS will include the long-range calculations.

The implementation of the AMOEBA and HIPPO force fields in LAMMPS was
done using F90 code provided by the Ponder group from their `Tinker MD
code <https://dasher.wustl.edu/tinker/>`_.

NOTE: what version of AMOEBA and HIPPO does LAMMPS implement?

----------

Only a single pair_coeff command is used with either the *amoeba* and
*hippo* styles which specifies two Tinker files, a PRM and KEY file.

.. code-block:: LAMMPS

   pair_coeff * * ../potentials/protein.prm.amoeba ../potentials/protein.key.amoeba
   pair_coeff * * ../potentials/water.prm.hippo ../potentials/water.key.hippo

Examples of the PRM files are in the potentials directory with an
*.amoeba or *.hippo suffix.  The examples/amoeba directory has
examples of both PRM and KEY files.

A Tinker PRM file composed of sections, each of which has multiple
lines.  A Tinker KEY file is composed of lines, each of which has a keyword,
which can be followed by zero or more parameters.  

The list of PRM sections and KEY keywords which LAMMPS recognizes are
listed on the doc:`Howto amoeba <howto_ameoba>` doc page.  If not
recognized, the section or keyword is skipped.

Note that if the KEY file is specified as NULL, then no file is
required; default values for various AMOEBA/HIPPO settings are used.
The doc:`Howto amoeba <howto_ameoba>` doc page also gives the default
settings.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These pair styles do not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

These pair styles do not write their information to :doc:`binary
restart files <restart>`, since it is stored in potential files.
Thus, you need to re-specify the pair_style and pair_coeff commands in
an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

These pair styles are part of the AMOEBA package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build
package <Build_package>` doc page for more info.

The AMOEBA and HIPPO potential (PRM) and KEY files provided with
LAMMPS in the potentials and examples/amoeva directories are Tinker
files parameterized for Tinker units.  Their numeric parameters are
converted by LAMMPS to its real units :doc:`units <units>`.  Thus uou
can only use these pair styles with real units.

These potentials do not yet calculate per-atom energy or virial
contributions.

As explained on the :doc:`Howto amoeba <howto_ameoba>` doc page, use
of these pair styles to run a simulation with the AMOEBA or HIPPO
force fields requires several things.

The first is a data file generated by the tools/tinker/tinker2lmp.py
conversion script which uses Tinker file force field file input to
create a data file compatible with LAMMPS.

The second is use of these commands:

* :doc:`atom_style amoeba <atom_style>`
* :doc:`fix property/atom <fix_property_atom>`
* :doc:`special_bonds one/five <special_bonds>`

And third, depending on the model being simulated, these
commands for intramolecular interactions may also be required:

* :doc:`bond_style class2 <bond_class2>`
* :doc:`angle_style amoeba <angle_amoeba>`
* :doc:`dihedral_style fourier <dihedral_fourier>`
* :doc:`improper_style amoeba <improper_amoeba>`
* :doc:`fix amoeba/pitorsion <fix_pitorsion>`
* :doc:`fix amoeba/bitorsion <fix_bitorsion>`

----------

Related commands
""""""""""""""""

:doc:`atom_style amoeba <atom_style>`, `bond_style amoeba
:doc:<bond_amoeba>`, `angle_style amoeba <angle_amoeba>`,
:doc:`dihedral_style amoeba <dihedral_amoeba>`, `special_bonds
:doc:one/five <special_bonds>`

Default
"""""""

none

----------

.. _amoeba-Ren:

**(Ren)** Ren and Ponder, J Phys Chem B, 107, 5933 (2003).

.. _amoeba-Shi:

**(Shi)** Shi, Xiz, Znahg, Best, Wu, Ponder, Ren, J Chem Theory Comp,
 9, 4046, 2013.

.. _amoeba-Rackers:

**(Rackers)** Rackers and Ponder, J Chem Phys, 150, 084104 (2010).
