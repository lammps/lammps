CHARMM, AMBER, COMPASS, and DREIDING force fields
=================================================

A force field has 2 parts: the formulas that define it and the
coefficients used for a particular system.  Here we only discuss
formulas implemented in LAMMPS that correspond to formulas commonly
used in the CHARMM, AMBER, COMPASS, and DREIDING force fields.  Setting
coefficients is done either from special sections in an input data file
via the :doc:`read\_data <read_data>` command or in the input script with
commands like :doc:`pair\_coeff <pair_coeff>` or
:doc:`bond\_coeff <bond_coeff>` and so on.  See the :doc:`Tools <Tools>` doc
page for additional tools that can use CHARMM, AMBER, or Materials
Studio generated files to assign force field coefficients and convert
their output into LAMMPS input.

See :ref:`(MacKerell) <howto-MacKerell>` for a description of the CHARMM force
field.  See :ref:`(Cornell) <howto-Cornell>` for a description of the AMBER
force field.  See :ref:`(Sun) <howto-Sun>` for a description of the COMPASS
force field.

.. _charmm: http://www.scripps.edu/brooks



.. _amber: http://amber.scripps.edu



The interaction styles listed below compute force field formulas that
are consistent with common options in CHARMM or AMBER.  See each
command's documentation for the formula it computes.

* :doc:`bond\_style <bond_harmonic>` harmonic
* :doc:`angle\_style <angle_charmm>` charmm
* :doc:`dihedral\_style <dihedral_charmm>` charmmfsh
* :doc:`dihedral\_style <dihedral_charmm>` charmm
* :doc:`pair\_style <pair_charmm>` lj/charmmfsw/coul/charmmfsh
* :doc:`pair\_style <pair_charmm>` lj/charmmfsw/coul/long
* :doc:`pair\_style <pair_charmm>` lj/charmm/coul/charmm
* :doc:`pair\_style <pair_charmm>` lj/charmm/coul/charmm/implicit
* :doc:`pair\_style <pair_charmm>` lj/charmm/coul/long

* :doc:`special\_bonds <special_bonds>` charmm
* :doc:`special\_bonds <special_bonds>` amber

.. note::

   For CHARMM, newer *charmmfsw* or *charmmfsh* styles were released
   in March 2017.  We recommend they be used instead of the older *charmm*
   styles.  See discussion of the differences on the :doc:`pair charmm <pair_charmm>` and :doc:`dihedral charmm <dihedral_charmm>` doc
   pages.

COMPASS is a general force field for atomistic simulation of common
organic molecules, inorganic small molecules, and polymers which was
developed using ab initio and empirical parameterization techniques.
See the :doc:`Tools <Tools>` doc page for the msi2lmp tool for creating
LAMMPS template input and data files from BIOVIA's Materials Studio
files.  Please note that the msi2lmp tool is very old and largely
unmaintained, so it does not support all features of Materials Studio
provided force field files, especially additions during the last decade.
You should watch the output carefully and compare results, where
possible.  See :ref:`(Sun) <howto-Sun>` for a description of the COMPASS force
field.

These interaction styles listed below compute force field formulas that
are consistent with the COMPASS force field.  See each command's
documentation for the formula it computes.

* :doc:`bond\_style <bond_class2>` class2
* :doc:`angle\_style <angle_class2>` class2
* :doc:`dihedral\_style <dihedral_class2>` class2
* :doc:`improper\_style <improper_class2>` class2

* :doc:`pair\_style <pair_class2>` lj/class2
* :doc:`pair\_style <pair_class2>` lj/class2/coul/cut
* :doc:`pair\_style <pair_class2>` lj/class2/coul/long

* :doc:`special\_bonds <special_bonds>` lj/coul 0 0 1

DREIDING is a generic force field developed by the `Goddard group <http://www.wag.caltech.edu>`_ at Caltech and is useful for
predicting structures and dynamics of organic, biological and main-group
inorganic molecules. The philosophy in DREIDING is to use general force
constants and geometry parameters based on simple hybridization
considerations, rather than individual force constants and geometric
parameters that depend on the particular combinations of atoms involved
in the bond, angle, or torsion terms. DREIDING has an :doc:`explicit hydrogen bond term <pair_hbond_dreiding>` to describe interactions involving a
hydrogen atom on very electronegative atoms (N, O, F).

See :ref:`(Mayo) <howto-Mayo>` for a description of the DREIDING force field

The interaction styles listed below compute force field formulas that
are consistent with the DREIDING force field.  See each command's
documentation for the formula it computes.

* :doc:`bond\_style <bond_harmonic>` harmonic
* :doc:`bond\_style <bond_morse>` morse

* :doc:`angle\_style <angle_harmonic>` harmonic
* :doc:`angle\_style <angle_cosine>` cosine
* :doc:`angle\_style <angle_cosine_periodic>` cosine/periodic

* :doc:`dihedral\_style <dihedral_charmm>` charmm
* :doc:`improper\_style <improper_umbrella>` umbrella

* :doc:`pair\_style <pair_buck>` buck
* :doc:`pair\_style <pair_buck>` buck/coul/cut
* :doc:`pair\_style <pair_buck>` buck/coul/long
* :doc:`pair\_style <pair_lj>` lj/cut
* :doc:`pair\_style <pair_lj>` lj/cut/coul/cut
* :doc:`pair\_style <pair_lj>` lj/cut/coul/long

* :doc:`pair\_style <pair_hbond_dreiding>` hbond/dreiding/lj
* :doc:`pair\_style <pair_hbond_dreiding>` hbond/dreiding/morse

* :doc:`special\_bonds <special_bonds>` dreiding


----------


.. _howto-MacKerell:



**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).

.. _howto-Cornell:



**(Cornell)** Cornell, Cieplak, Bayly, Gould, Merz, Ferguson,
Spellmeyer, Fox, Caldwell, Kollman, JACS 117, 5179-5197 (1995).

.. _howto-Sun:



**(Sun)** Sun, J. Phys. Chem. B, 102, 7338-7364 (1998).

.. _howto-Mayo:



**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
