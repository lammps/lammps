.. index:: dynamical_matrix

dynamical_matrix command
========================

Syntax
""""""

.. code-block:: LAMMPS

   dynamical_matrix group-ID style gamma args keyword value ...

* group-ID = ID of group of atoms to displace
* style = *regular* or *eskm*
* gamma = finite different displacement length (distance units)
* one or more keyword/arg pairs may be appended

  .. parsed-literal::

       keyword = *file* or *binary*
         *file* name = name of output file for the dynamical matrix
         *binary* arg = *yes* or *no* or *gzip*

Examples
""""""""

.. code-block:: LAMMPS

   dynamical_matrix 1 regular 0.000001
   dynamical_matrix 1 eskm 0.000001
   dynamical_matrix 3 regular 0.00004 file dynmat.dat
   dynamical_matrix 5 eskm 0.00000001 file dynamical.dat binary yes

Description
"""""""""""

Calculate the dynamical matrix by finite difference of the selected group,

.. math::

   D = \frac{\Phi_{ij}^{\alpha\beta}}{\sqrt{M_i M_j}}

where D is the dynamical matrix and :math:`\Phi` is the force constant
matrix defined by

.. math::

   \Phi_{ij}^{\alpha\beta} = \frac{\partial^2 U}{\partial x_{i,\alpha} \partial x_{j,\beta}}

The output for the dynamical matrix is printed three elements at a time.
The three elements are the three :math:`\beta` elements for a respective
i/:math:`\alpha`/j combination.  Each line is printed in order of j
increasing first, :math:`\alpha` second, and i last.

If the style eskm is selected, the dynamical matrix will be in units of
inverse squared femtoseconds. These units will then conveniently leave
frequencies in THz.

Restrictions
""""""""""""

The command collects an array of nine times the number of atoms in a group
on every single MPI rank, so the memory requirements can be very significant
for large systems.

This command is part of the USER-PHONON package.  It is only enabled if
LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix phonon <fix_phonon>`

:doc:`compute hma <compute_hma>` uses an analytic formulation of the
Hessian provided by a pair_style's Pair::single_hessian() function,
if implemented.

Default
"""""""

The default settings are file = "dynmat.dyn", binary = no
