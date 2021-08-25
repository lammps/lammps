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

doc:`Howto amoeba <howto_ameoba>`
doc:`bond_style amoeba <bond_amoeba>`
doc:`angle_style amoeba <angle_amoeba>`
doc:`dihedral_style amoeba <dihedral_amoeba>`
examples/amoeba
tools/amoeba
potentials/\*.prm.ameoba
potentials/\*.key.ameoba
potentials/\*.prm.hippo
potentials/\*.key.hippo

Description
"""""""""""

NOTE: edit this paragraph however you wish including adding or changing
citations (see bottom of this file)

The *amoeba* style computes the AMOEBA polarizeable field formulated
by Jay Ponder's group at U Washington at St Louis :ref:`(Ren)
<amoeba-Ren>`, :ref:`(Shi) <amoeba-Shi>`.  The *hippo* style
computes the HIPPO polarizeable force field , an extension to AMOEBA,
formulated by Josh Rackers and collaborators in the Ponder group
:ref:`(Rackers) <amoeba-Rackers>`.  

These force fields can be used when polarization effects are desired
in simulations of water, organic molecules, and biomolecules including
proteins, provided that parameterizations (force field files) are
available for the systems you are interested in.  Files in the LAMMPS
potentials with a "amoeba" or "hippo" suffix can be used.  The Tinker
distribution and website may have other such files.

NOTE: replace these LaTeX formulas with a list of AMOEBA and HIPPO
terms in some simple notation.  If desired, a more detailed
mathematical expression for each term could also be included below
the initial 2 formulas.

The AMOEBA force field can be written as a collection of terms

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff, blah is the dipole, etc

The HIPPO force field is similar but alters some of the terms

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff, blah is the dipole, etc

NOTE: Add a sentence for each term to explain what physical effects
the FFs are encoding.  E.g. The polar term iteratively computes an
induced dipole for each atom, then calculates dipole-dipole
interactions ...

See the AMOEBA and HIPPO papers for further details.

.. note::

  The AMOEBA and HIPPO force fields compute long-range charge, dipole,
  and quadrupole interactions (NOTE: also long-range dispersion?).
  However, unlike other models with long-range interactions in LAMMPS,
  this does not require use of a KSpace style via the
  :doc:`kspace_style <kspace_style>` command.  That is because for
  AMOEBA and HIPPO the long-range computations are intertwined with
  the pairwise compuations.  So these pair style include both short-
  and long-range computations.  This means the energy and virial
  computed by the pair style as well as the "Pair" timing reported by
  LAMMPS will include the long-range components.

.. note::

  As explained on the :doc:`Howto amoeba <howto_ameoba>` doc page, use
  of these pair styles to run a simulation with the AMOEBA or HIPPO
  force fields requires your input script to use the :doc:`atom_style
  hybrid full amoeba <atom_style>` atom style, AMOEBA versions of
  bond/angle/dihedral styles, the :doc:`special_bonds one/five
  <special_bonds>` option, and the :doc:`fix property/atom one/five
  <fix_property>` command to define several additional per-atom
  properties.  The latter requires a "Tinker Types" section be
  included in the LAMMPS data file.  This can be auto-generated using
  the tools/amoeba/tinker2lmp.py Python script.  See the :doc:`Howto
  amoeba <howto_ameoba>` doc page and tools/amoeba/README file for
  details on using that tool.

The implementation of the AMOEBA and HIPPO force fields in LAMMPS was
done using code provided by the Ponder group from their `Tinker MD
code <https://dasher.wustl.edu/tinker/>`_ written in F90.

NOTE: what version of AMOEBA and HIPPO does LAMMPS implement?

----------

Only a single pair_coeff command is used with either the *amoeba* and
*hippo* styles which specifies two Tinker files, a PRM and KEY file.

.. code-block:: LAMMPS

   pair_coeff * * ../potentials/protein.prm.amoeba ../potentials/protein.key.amoeba
   pair_coeff * * ../potentials/water.prm.hippo ../potentials/water.key.hippo

See examples of these files in the potentials directory.

The format of a PRM file is a collection of sections, each with
multiple lines.  These are the sections which LAMMPS recognizes:

NOTE: need to list these, possibly there are some LAMMPS skips
      or which need to be removed for LAMMPS to use the PRM file?

The format of a KEY file is a series of lines, with one parameter and
its value per line.  These are the parameters which LAMMPS recognizes:

NOTE: need to list these, possibly there are some keywords LAMMPS skips
      or which need to be removed for LAMMPS to use the KEY file?

NOTE: Any other info about PRM and KEY files we should explain
for LAMMPS users here?

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
LAMMPS in the potentials directory are Tinker files parameterized for
Tinker units.  Their numeric parameters are converted by LAMMPS to its
real units :doc:`units <units>`.  You can only use these pair styles
with real units.

These potentials do not yet calculate per-atom energy or virial
contributions.

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
