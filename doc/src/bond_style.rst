.. index:: bond\_style

bond\_style command
===================

Syntax
""""""


.. parsed-literal::

   bond_style style args

* style = *none* or *hybrid* or *class2* or *fene* or *fene/expand* or         *harmonic* or *morse* or *nonlinear* or *quartic*


.. parsed-literal::

     args = none for any style except *hybrid*
     *hybrid* args = list of one or more styles

Examples
""""""""


.. parsed-literal::

   bond_style harmonic
   bond_style fene
   bond_style hybrid harmonic fene

Description
"""""""""""

Set the formula(s) LAMMPS uses to compute bond interactions between
pairs of atoms.  In LAMMPS, a bond differs from a pairwise
interaction, which are set via the :doc:`pair\_style <pair_style>`
command.  Bonds are defined between specified pairs of atoms and
remain in force for the duration of the simulation (unless the bond
breaks which is possible in some bond potentials).  The list of bonded
atoms is read in by a :doc:`read\_data <read_data>` or
:doc:`read\_restart <read_restart>` command from a data or restart file.
By contrast, pair potentials are typically defined between all pairs
of atoms within a cutoff distance and the set of active interactions
changes over time.

Hybrid models where bonds are computed using different bond potentials
can be setup using the *hybrid* bond style.

The coefficients associated with a bond style can be specified in a
data or restart file or via the :doc:`bond\_coeff <bond_coeff>` command.

All bond potentials store their coefficient data in binary restart
files which means bond\_style and :doc:`bond\_coeff <bond_coeff>` commands
do not need to be re-specified in an input script that restarts a
simulation.  See the :doc:`read\_restart <read_restart>` command for
details on how to do this.  The one exception is that bond\_style
*hybrid* only stores the list of sub-styles in the restart file; bond
coefficients need to be re-specified.

.. note::

   When both a bond and pair style is defined, the
   :doc:`special\_bonds <special_bonds>` command often needs to be used to
   turn off (or weight) the pairwise interaction that would otherwise
   exist between 2 bonded atoms.

In the formulas listed for each bond style, *r* is the distance
between the 2 atoms in the bond.


----------


Here is an alphabetic list of bond styles defined in LAMMPS.  Click on
the style to display the formula it computes and coefficients
specified by the associated :doc:`bond\_coeff <bond_coeff>` command.

Click on the style to display the formula it computes, any additional
arguments specified in the bond\_style command, and coefficients
specified by the associated :doc:`bond\_coeff <bond_coeff>` command.

There are also additional accelerated pair styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :doc:`Commands bond <Commands_bond>`
doc page are followed by one or more of (g,i,k,o,t) to indicate which
accelerated styles exist.

* :doc:`none <bond_none>` - turn off bonded interactions
* :doc:`zero <bond_zero>` - topology but no interactions
* :doc:`hybrid <bond_hybrid>` - define multiple styles of bond interactions

* :doc:`class2 <bond_class2>` - COMPASS (class 2) bond
* :doc:`fene <bond_fene>` - FENE (finite-extensible non-linear elastic) bond
* :doc:`fene/expand <bond_fene_expand>` - FENE bonds with variable size particles
* :doc:`gromos <bond_gromos>` - GROMOS force field bond
* :doc:`harmonic <bond_harmonic>` - harmonic bond
* :doc:`harmonic/shift <bond_harmonic_shift>` - shifted harmonic bond
* :doc:`harmonic/shift/cut <bond_harmonic_shift_cut>` - shifted harmonic bond with a cutoff
* :doc:`mm3 <bond_mm3>` - MM3 anharmonic bond
* :doc:`morse <bond_morse>` - Morse bond
* :doc:`nonlinear <bond_nonlinear>` - nonlinear bond
* :doc:`oxdna/fene <bond_oxdna>` - modified FENE bond suitable for DNA modeling
* :doc:`oxdna2/fene <bond_oxdna>` - same as oxdna but used with different pair styles
* :doc:`quartic <bond_quartic>` - breakable quartic bond
* :doc:`table <bond_table>` - tabulated by bond length


----------


Restrictions
""""""""""""


Bond styles can only be set for atom styles that allow bonds to be
defined.

Most bond styles are part of the MOLECULE package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.  The doc pages for
individual bond potentials tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`

Default
"""""""

bond\_style none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
