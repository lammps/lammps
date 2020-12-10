.. index:: pair_style tersoff/shift

Syntax
""""""

.. code-block:: LAMMPS

   pair_style tersoff/shift delta

* delta = the shift applied to the equilibrium bond length of the tersoff potential.

Examples
""""""""

.. code-block:: LAMMPS

   pair_style tersoff/shift delta
   pair_coeff * * Si.tersoff Si
   pair_coeff * * BNC.tersoff B N

Description
"""""""""""

The *tersoff/shift* style computes the energy E of a system of atoms, whose formula
is the same as a 3-body Tersoff potential :ref:`(Tersoff_1) <Tersoff_11>`. The only
modification is that the original equilibrium bond length (:math: `r_0`) of the 
system is shifted to :math:`r_0-\delta`.

.. note::

   The sign of :math:`\delta` determines whether the resulting equilibrium bond length will be elongated 
   or shrinked relative to the original value. Specifically, values of :math:`\delta < 0` will result in 
   elongation of the bond length, while values of :math:`\delta > 0` will result in shrinking of the bond length.

This style is designed for simulations of closely matched van der Waals heterostructures. For instance, let's 
consider the case of a system with few-layers graphene atop a thick hexagonal boron nitride (h-BN) substrate 
simulated using periodic boundary conditions. The experimental lattice mismatch of ~1.8% between graphene and h-BN
is not well captured by the equilibrium lattice constants of available potentials, thus a small in-plane strain
will be introduced in the system when building a periodic supercell.To minimize the effect of strain on
simulation results, the *tersoff/shift* style is proposed that allows adjusting the equilibrium bond length
of one of the two materials (e.g., h-BN). Validitation, benchmark tests and applications of the *tersoff/shift* style
can be found in :ref:`(Mandelli_1) <Mandelli1>` and :ref:`(Ouyang_1) <Ouyang5>`.

For the specific case discussed above, the force field can be defined as

.. code-block:: LAMMPS

   pair_style  hybrid/overlay rebo tersoff/shift -4.07e-3 ilp/graphene/hbn 16.0 coul/shield 16.0
   pair_coeff  * * rebo              CH.rebo     NULL NULL C
   pair_coeff  * * tersoff/shift     BNC.tersoff B    N    NULL
   pair_coeff  * * ilp/graphene/hbn  BNCH.ILP    B    N    C
   pair_coeff  1 1 coul/shield 0.70
   pair_coeff  1 2 coul/shield 0.695
   pair_coeff  2 2 coul/shield 0.69

Except for the shift, the usage of the *tersoff/shift* style is the same as that for the
*tersoff* style :doc:`tersoff <pair_tersoff>`.


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
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the USER-MISC package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Tersoff potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.  You can
use the Tersoff potential with any LAMMPS units, but you would need to
create your own Tersoff potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_style tersoff <pair_tersoff>`,
:doc:`pair_style ilp/graphene/hbn <pair_ilp_graphene_hbn>`.

Default
"""""""

none

----------

.. _Mandelli1:

**(Mandelli_1)** D. Mandelli, W. Ouyang, M. Urbakh, and O. Hod, ACS Nano 13(7), 7603-7609 (2019).


.. _Ouyang5:

**(Ouyang_1)** W. Ouyang et al., J. Chem. Theory Comput. 16(1), 666-676 (2020).
