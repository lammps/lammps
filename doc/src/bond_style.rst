.. index:: bond_style

bond_style command
==================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style style args

* style = *none* or *zero* or *hybrid* or *bpm/rotational* or *bpm/spring* or *class2* or *fene* or *fene/expand* or *fene/nm* or *gaussian* or *gromos* or *harmonic* or *harmonic/shift* or *harmonic/shift/cut* or *morse* or *nonlinear* or *oxdna/fene* or *oxdena2/fene* or *oxrna2/fene* or *quartic* or *special* or *table*

* args = none for any style except *hybrid*

  * *hybrid* args = list of one or more styles

Examples
""""""""

.. code-block:: LAMMPS

   bond_style harmonic
   bond_style fene
   bond_style hybrid harmonic fene

Description
"""""""""""

Set the formula(s) LAMMPS uses to compute bond interactions between
pairs of atoms.  In LAMMPS, a bond differs from a pairwise
interaction, which are set via the :doc:`pair_style <pair_style>`
command.  Bonds are defined between specified pairs of atoms and
remain in force for the duration of the simulation (unless the bond
breaks which is possible in some bond potentials).  The list of bonded
atoms is read in by a :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` command from a data or restart file.
By contrast, pair potentials are typically defined between all pairs
of atoms within a cutoff distance and the set of active interactions
changes over time.

Hybrid models where bonds are computed using different bond potentials
can be setup using the *hybrid* bond style.

The coefficients associated with a bond style can be specified in a
data or restart file or via the :doc:`bond_coeff <bond_coeff>` command.

All bond potentials store their coefficient data in binary restart
files which means bond_style and :doc:`bond_coeff <bond_coeff>` commands
do not need to be re-specified in an input script that restarts a
simulation.  See the :doc:`read_restart <read_restart>` command for
details on how to do this.  The one exception is that bond_style
*hybrid* only stores the list of sub-styles in the restart file; bond
coefficients need to be re-specified.

.. note::

   When both a bond and pair style is defined, the
   :doc:`special_bonds <special_bonds>` command often needs to be used to
   turn off (or weight) the pairwise interaction that would otherwise
   exist between 2 bonded atoms.

In the formulas listed for each bond style, *r* is the distance
between the 2 atoms in the bond.

----------

Here is an alphabetic list of bond styles defined in LAMMPS.  Click on
the style to display the formula it computes and coefficients
specified by the associated :doc:`bond_coeff <bond_coeff>` command.

Click on the style to display the formula it computes, any additional
arguments specified in the bond_style command, and coefficients
specified by the associated :doc:`bond_coeff <bond_coeff>` command.

There are also additional accelerated pair styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :doc:`Commands bond <Commands_bond>`
doc page are followed by one or more of (g,i,k,o,t) to indicate which
accelerated styles exist.

* :doc:`none <bond_none>` - turn off bonded interactions
* :doc:`zero <bond_zero>` - topology but no interactions
* :doc:`hybrid <bond_hybrid>` - define multiple styles of bond interactions

* :doc:`bpm/rotational <bond_bpm_rotational>` - breakable bond with forces and torques based on deviation from reference state
* :doc:`bpm/spring <bond_bpm_spring>` - breakable bond with forces based on deviation from reference length
* :doc:`class2 <bond_class2>` - COMPASS (class 2) bond
* :doc:`fene <bond_fene>` - FENE (finite-extensible non-linear elastic) bond
* :doc:`fene/expand <bond_fene_expand>` - FENE bonds with variable size particles
* :doc:`fene/nm <bond_fene>` - FENE bonds with a generalized Lennard-Jones potential
* :doc:`gaussian <bond_gaussian>` - multicentered Gaussian-based bond potential
* :doc:`gromos <bond_gromos>` - GROMOS force field bond
* :doc:`harmonic <bond_harmonic>` - harmonic bond
* :doc:`harmonic/shift <bond_harmonic_shift>` - shifted harmonic bond
* :doc:`harmonic/shift/cut <bond_harmonic_shift_cut>` - shifted harmonic bond with a cutoff
* :doc:`mesocnt <bond_mesocnt>` - Harmonic bond wrapper with parameterization presets for nanotubes
* :doc:`mm3 <bond_mm3>` - MM3 anharmonic bond
* :doc:`morse <bond_morse>` - Morse bond
* :doc:`nonlinear <bond_nonlinear>` - nonlinear bond
* :doc:`oxdna/fene <bond_oxdna>` - modified FENE bond suitable for DNA modeling
* :doc:`oxdna2/fene <bond_oxdna>` - same as oxdna but used with different pair styles
* :doc:`oxrna2/fene <bond_oxdna>` - modified FENE bond suitable for RNA modeling
* :doc:`quartic <bond_quartic>` - breakable quartic bond
* :doc:`special <bond_special>` - enable special bond exclusions for 1-5 pairs and beyond
* :doc:`table <bond_table>` - tabulated by bond length

----------

Restrictions
""""""""""""

Bond styles can only be set for atom styles that allow bonds to be
defined.

Most bond styles are part of the MOLECULE package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.  The doc pages for
individual bond potentials tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

.. code-block:: LAMMPS

   bond_style none
