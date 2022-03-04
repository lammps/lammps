.. index:: angle_style

angle_style command
===================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style style

* style = *none* or *hybrid* or *charmm* or *class2* or *cosine* or         *cosine/squared* or *harmonic*

Examples
""""""""

.. code-block:: LAMMPS

   angle_style harmonic
   angle_style charmm
   angle_style hybrid harmonic cosine

Description
"""""""""""

Set the formula(s) LAMMPS uses to compute angle interactions between
triplets of atoms, which remain in force for the duration of the
simulation.  The list of angle triplets is read in by a
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` command
from a data or restart file.

Hybrid models where angles are computed using different angle
potentials can be setup using the *hybrid* angle style.

The coefficients associated with a angle style can be specified in a
data or restart file or via the :doc:`angle_coeff <angle_coeff>` command.

All angle potentials store their coefficient data in binary restart
files which means angle_style and :doc:`angle_coeff <angle_coeff>`
commands do not need to be re-specified in an input script that
restarts a simulation.  See the :doc:`read_restart <read_restart>`
command for details on how to do this.  The one exception is that
angle_style *hybrid* only stores the list of sub-styles in the restart
file; angle coefficients need to be re-specified.

.. note::

   When both an angle and pair style is defined, the
   :doc:`special_bonds <special_bonds>` command often needs to be used to
   turn off (or weight) the pairwise interaction that would otherwise
   exist between 3 bonded atoms.

In the formulas listed for each angle style, *theta* is the angle
between the 3 atoms in the angle.

----------

Here is an alphabetic list of angle styles defined in LAMMPS.  Click on
the style to display the formula it computes and coefficients
specified by the associated :doc:`angle_coeff <angle_coeff>` command.

Click on the style to display the formula it computes, any additional
arguments specified in the angle_style command, and coefficients
specified by the associated :doc:`angle_coeff <angle_coeff>` command.

There are also additional accelerated pair styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :ref:`Commands angle <angle>` page are followed by one or more
of (g,i,k,o,t) to indicate which accelerated styles exist.

* :doc:`none <angle_none>` - turn off angle interactions
* :doc:`zero <angle_zero>` - topology but no interactions
* :doc:`hybrid <angle_hybrid>` - define multiple styles of angle interactions

* :doc:`charmm <angle_charmm>` - CHARMM angle
* :doc:`class2 <angle_class2>` - COMPASS (class 2) angle
* :doc:`class2/p6 <angle_class2>` - COMPASS (class 2) angle expanded to 6th order
* :doc:`cosine <angle_cosine>` - angle with cosine term
* :doc:`cosine/buck6d <angle_cosine_buck6d>` - same as cosine with Buckingham term between 1-3 atoms
* :doc:`cosine/delta <angle_cosine_delta>` - angle with difference of cosines
* :doc:`cosine/periodic <angle_cosine_periodic>` - DREIDING angle
* :doc:`cosine/shift <angle_cosine_shift>` - angle cosine with a shift
* :doc:`cosine/shift/exp <angle_cosine_shift_exp>` - cosine with shift and exponential term in spring constant
* :doc:`cosine/squared <angle_cosine_squared>` - angle with cosine squared term
* :doc:`cross <angle_cross>` - cross term coupling angle and bond lengths
* :doc:`dipole <angle_dipole>` - angle that controls orientation of a point dipole
* :doc:`fourier <angle_fourier>` - angle with multiple cosine terms
* :doc:`fourier/simple <angle_fourier_simple>` - angle with a single cosine term
* :doc:`gaussian <angle_gaussian>` - multicentered Gaussian-based angle potential
* :doc:`harmonic <angle_harmonic>` - harmonic angle
* :doc:`mm3 <angle_mm3>` - anharmonic angle
* :doc:`quartic <angle_quartic>` - angle with cubic and quartic terms
* :doc:`sdk <angle_sdk>` - harmonic angle with repulsive SDK pair style between 1-3 atoms
* :doc:`table <angle_table>` - tabulated by angle

----------

Restrictions
""""""""""""

Angle styles can only be set for atom_styles that allow angles to be
defined.

Most angle styles are part of the MOLECULE package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.  The doc pages for
individual bond potentials tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

.. code-block:: LAMMPS

   angle_style none
