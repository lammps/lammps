.. index:: pair_style sw/mod

pair_style sw/mod command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style sw/mod

Examples
""""""""

.. code-block:: LAMMPS

   pair_style sw/mod
   pair_coeff * * tmd.sw.mod Mo S S

Description
"""""""""""

The *sw/mod* style computes the energy E of a system of atoms, whose
formula is the same as the Stillinger-Weber potential :doc:`pair_style sw <pair_sw>`.
The only modification is in the three-body term, where 
:math:`\delta = \cos \theta_{ijk} - \cos \theta_{0ijk}`
is modified with the following function:

.. math::

  g_C(\delta) & = \left\{ \begin{array} {r@{\quad:\quad}l}
    1 & \delta < \delta_1 \\
    \frac{1}{2} + \frac{1}{2} \cos \left( \pi \frac{\delta-\delta_1}{\delta_2 - \delta_1} \right) &
      \delta_1 < \delta < \delta_2 \\
    0 & \delta > \delta_2
    \end{array} \right. \\


This potential is designed for simulations of materials when
distinguishing three-body angles are necessary, such as borophene
and transition metal dichalcogenide, which cannot be described 
by the original code for the Stillinger-Weber potential. Validation, 
benchmark tests, and applications of the *modify* keyword can be found in
:ref:`(Jiang_1) <Jiang1>` and :ref:`(Jiang_2) <Jiang2>`.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above from values in the potential file.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Stillinger-Weber potential files provided with LAMMPS (see the
potentials directory) are parameterized for metal :doc:`units <units>`.
You can use the SW potential with any LAMMPS units, but you would need
to create your own SW potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" units.

Related commands
""""""""""""""""
:doc:`pair_coeff <pair_coeff>`

:doc:`pair_style sw <pair_sw>`

Default
"""""""

none

----------

.. _Jiang1:

**(Jiang_1)** J.-W. Jiang, Nanotechnology 26, 315706 (2015).

.. _Jiang2:

**(Jiang_2)** J.-W. Jiang, Acta Mech. Solida. Sin 32, 17 (2019).
