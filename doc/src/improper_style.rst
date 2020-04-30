.. index:: improper_style

improper_style command
======================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style style

* style = *none* or *hybrid* or *class2* or *cvff* or *harmonic*

Examples
""""""""

.. code-block:: LAMMPS

   improper_style harmonic
   improper_style cvff
   improper_style hybrid cvff harmonic

Description
"""""""""""

Set the formula(s) LAMMPS uses to compute improper interactions
between quadruplets of atoms, which remain in force for the duration
of the simulation.  The list of improper quadruplets is read in by a
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` command
from a data or restart file.  Note that the ordering of the 4 atoms in
an improper quadruplet determines the definition of the improper
angle used in the formula for each style.  See the doc pages of
individual styles for details.

Hybrid models where impropers are computed using different improper
potentials can be setup using the *hybrid* improper style.

The coefficients associated with an improper style can be specified in
a data or restart file or via the :doc:`improper_coeff <improper_coeff>`
command.

All improper potentials store their coefficient data in binary restart
files which means improper_style and
:doc:`improper_coeff <improper_coeff>` commands do not need to be
re-specified in an input script that restarts a simulation.  See the
:doc:`read_restart <read_restart>` command for details on how to do
this.  The one exception is that improper_style *hybrid* only stores
the list of sub-styles in the restart file; improper coefficients need
to be re-specified.

.. note::

   When both an improper and pair style is defined, the
   :doc:`special_bonds <special_bonds>` command often needs to be used to
   turn off (or weight) the pairwise interaction that would otherwise
   exist between a group of 4 bonded atoms.

----------

Here is an alphabetic list of improper styles defined in LAMMPS.
Click on the style to display the formula it computes and coefficients
specified by the associated :doc:`improper_coeff <improper_coeff>`
command.

Click on the style to display the formula it computes, any additional
arguments specified in the improper_style command, and coefficients
specified by the associated :doc:`improper_coeff <improper_coeff>`
command.

There are also additional accelerated pair styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :ref:`Commands improper <improper>` doc page are followed by one or
more of (g,i,k,o,t) to indicate which accelerated styles exist.

* :doc:`none <improper_none>` - turn off improper interactions
* :doc:`zero <improper_zero>` - topology but no interactions
* :doc:`hybrid <improper_hybrid>` - define multiple styles of improper interactions

* :doc:`class2 <improper_class2>` - COMPASS (class 2) improper
* :doc:`cossq <improper_cossq>` - improper with a cosine squared term
* :doc:`cvff <improper_cvff>` - CVFF improper
* :doc:`distance <improper_distance>` - improper based on distance between atom planes
* :doc:`distharm <improper_distharm>` - improper that is harmonic in the out-of-plane distance
* :doc:`fourier <improper_fourier>` - improper with multiple cosine terms
* :doc:`harmonic <improper_harmonic>` - harmonic improper
* :doc:`inversion/harmonic <improper_inversion_harmonic>` - harmonic improper with Wilson-Decius out-of-plane definition
* :doc:`ring <improper_ring>` - improper which prevents planar conformations
* :doc:`umbrella <improper_umbrella>` - DREIDING improper

:doc:`sqdistharm <improper_sqdistharm>` - improper that is harmonic in the square of the out-of-plane distance

----------

Restrictions
""""""""""""

Improper styles can only be set for atom_style choices that allow
impropers to be defined.

Most improper styles are part of the MOLECULE package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.  The doc pages for
individual improper potentials tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

.. code-block:: LAMMPS

   improper_style none
