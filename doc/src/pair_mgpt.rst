.. index:: pair_style mgpt

pair_style mgpt command
=======================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style mgpt

Examples
""""""""


.. code-block:: LAMMPS

   pair_style mgpt
   pair_coeff * * Ta6.8x.mgpt.parmin Ta6.8x.mgpt.potin Omega
   cp ~/lammps/potentials/Ta6.8x.mgpt.parmin parmin
   cp ~/lammps/potentials/Ta6.8x.mgpt.potin potin
   pair_coeff * * parmin potin Omega volpress yes nbody 1234 precision double
   pair_coeff * * parmin potin Omega volpress yes nbody 12

Description
"""""""""""

Within DFT quantum mechanics, generalized pseudopotential theory (GPT)
(:ref:`Moriarty1 <Moriarty1>`) provides a first-principles approach to
multi-ion interatomic potentials in d-band transition metals, with a
volume-dependent, real-space total-energy functional for the N-ion
elemental bulk material in the form

.. math::

   E_{\rm tot}({\bf R}_1 \ldots {\bf R}_N) = NE_{\rm vol}(\Omega )
   + \frac{1}{2} \sum _{i,j} \mbox{}^\prime \ v_2(ij;\Omega )
   + \frac{1}{6} \sum _{i,j,k} \mbox{}^\prime \ v_3(ijk;\Omega )
   + \frac{1}{24} \sum _{i,j,k,l} \mbox{}^\prime \ v_4(ijkl;\Omega )


where the prime on each summation sign indicates the exclusion of all
self-interaction terms from the summation.  The leading volume term
E\_vol as well as the two-ion central-force pair potential v\_2 and the
three- and four-ion angular-force potentials, v\_3 and v\_4, depend
explicitly on the atomic volume Omega, but are structure independent
and transferable to all bulk ion configurations, either ordered or
disordered, and with of without the presence of point and line
defects.  The simplified model GPT or MGPT (:ref:`Moriarty2 <Moriarty2>`,
:ref:`Moriarty3 <Moriarty3>`), which retains the form of E\_tot and permits
more efficient large-scale atomistic simulations, derives from the GPT
through a series of systematic approximations applied to E\_vol and the
potentials v\_n that are valid for mid-period transition metals with
nearly half-filled d bands.

Both analytic (:ref:`Moriarty2 <Moriarty2>`) and matrix
(:ref:`Moriarty3 <Moriarty3>`) representations of MGPT have been developed.
In the more general matrix representation, which can also be applied
to f-band actinide metals and permits both canonical and non-canonical
d/f bands, the multi-ion potentials are evaluated on the fly during a
simulation through d- or f-state matrix multiplication, and the forces
that move the ions are determined analytically.  Fast matrix-MGPT
algorithms have been developed independently by Glosli
(:ref:`Glosli <Glosli>`, :ref:`Moriarty3 <Moriarty3>`) and by Oppelstrup
(:ref:`Oppelstrup <Oppelstrup>`)

The *mgpt* pair style calculates forces, energies, and the total
energy per atom, E\_tot/N, using the Oppelstrup matrix-MGPT algorithm.
Input potential and control data are entered through the
:doc:`pair_coeff <pair_coeff>` command.  Each material treated requires
input parmin and potin potential files, as shown in the above
examples, as well as specification by the user of the initial atomic
volume Omega through pair\_coeff.  At the beginning of a time step in
any simulation, the total volume of the simulation cell V should
always be equal to Omega\*N, where N is the number of metal ions
present, taking into account the presence of any vacancies and/or
interstitials in the case of a solid.  In a constant-volume
simulation, which is the normal mode of operation for the *mgpt* pair
style, Omega, V and N all remain constant throughout the simulation
and thus are equal to their initial values.  In a constant-stress
simulation, the cell volume V will change (slowly) as the simulation
proceeds.  After each time step, the atomic volume should be updated
by the code as Omega = V/N.  In addition, the volume term E\_vol and
the potentials v\_2, v\_3 and v\_4 have to be removed at the end of the
time step, and then respecified at the new value of Omega.  In all
simulations, Omega must remain within the defined volume range for
E\_vol and the potentials for the given material.

The default option volpress yes in the :doc:`pair_coeff <pair_coeff>`
command includes all volume derivatives of E\_tot required to calculate
the stress tensor and pressure correctly.  The option volpress no
disregards the pressure contribution resulting from the volume term
E\_vol, and can be used for testing and analysis purposes.  The
additional optional variable nbody controls the specific terms in
E\_tot that are calculated.  The default option and the normal option
for mid-period transition and actinide metals is nbody 1234 for which
all four terms in E\_tot are retained.  The option nbody 12, for
example, retains only the volume term and the two-ion pair potential
term and can be used for GPT series-end transition metals that can be
well described without v\_3 and v\_4.  The nbody option can also be used
to test or analyze the contribution of any of the four terms in E\_tot
to a given calculated property.

The *mgpt* pair style makes extensive use of matrix algebra and
includes optimized kernels for the BlueGene/Q architecture and the
Intel/AMD (x86) architectures.  When compiled with the appropriate
compiler and compiler switches (-msse3 on x86, and using the IBM XL
compiler on BG/Q), these optimized routines are used automatically.
For BG/Q machines, building with the default Makefile for that
architecture (e.g., "make bgq") should enable the optimized algebra
routines.  For x-86 machines, there is a provided Makefile.mgptfast
which enables the fast algebra routines, i.e. build LAMMPS with "make
mgptfast".  The user will be informed in the output files of the
matrix kernels in use. To further improve speed, on x86 the option
precision single can be added to the :doc:`pair_coeff <pair_coeff>`
command line, which improves speed (up to a factor of two) at the cost
of doing matrix calculations with 7 digit precision instead of the
default 16. For consistency the default option can be specified
explicitly by the option precision double.

All remaining potential and control data are contained with the parmin
and potin files, including cutoffs, atomic mass, and other basic MGPT
variables.  Specific MGPT potential data for the transition metals
tantalum (Ta4 and Ta6.8x potentials), molybdenum (Mo5.2 potentials),
and vanadium (V6.1 potentials) are contained in the LAMMPS potentials
directory.  The stored files are, respectively, Ta4.mgpt.parmin,
Ta4.mgpt.potin, Ta6.8x.mgpt.parmin, Ta6.8x.mgpt.potin,
Mo5.2.mgpt.parmin, Mo5.2.mgpt.potin, V6.1.mgpt.parmin, and
V6.1.mgpt.potin .  Useful corresponding informational "README" files
on the Ta4, Ta6.8x, Mo5.2 and V6.1 potentials are also included in the
potentials directory.  These latter files indicate the volume mesh and
range for each potential and give appropriate references for the
potentials.  It is expected that MGPT potentials for additional
materials will be added over time.

Useful example MGPT scripts are given in the examples/USER/mgpt
directory.  These scripts show the necessary steps to perform
constant-volume calculations and simulations.  It is strongly
recommended that the user work through and understand these examples
before proceeding to more complex simulations.

.. note::

   For good performance, LAMMPS should be built with the compiler
   flags "-O3 -msse3 -funroll-loops" when including this pair style.  The
   src/MAKE/OPTIONS/Makefile.mgptfast is an example machine Makefile with
   these options included as part of a standard MPI build.  Note that it
   as provided, it will build with whatever low-level compiler (g++, icc,
   etc) is the default for your MPI installation.


----------


**Mixing, shift, table tail correction, restart**\ :

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
needs to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This pair style is part of the USER-MGPT package and is only enabled
if LAMMPS is built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

The MGPT potentials require the :doc:`newtion <newton>` setting to be
"on" for pair style interactions.

The stored parmin and potin potential files provided with LAMMPS in
the "potentials" directory are written in Rydberg atomic units, with
energies in Rydbergs and distances in Bohr radii. The *mgpt* pair
style converts Rydbergs to Hartrees to make the potential files
compatible with LAMMPS electron :doc:`units <units>`.

The form of E\_tot used in the *mgpt* pair style is only appropriate
for elemental bulk solids and liquids.  This includes solids with
point and extended defects such as vacancies, interstitials, grain
boundaries and dislocations.  Alloys and free surfaces, however,
require significant modifications, which are not included in the
*mgpt* pair style.  Likewise, the *hybrid* pair style is not allowed,
where MGPT would be used for some atoms but not for others.

Electron-thermal effects are not included in the standard MGPT
potentials provided in the "potentials" directory, where the
potentials have been constructed at zero electron temperature.
Physically, electron-thermal effects may be important in 3d (e.g., V)
and 4d (e.g., Mo) transition metals at high temperatures near melt and
above.  It is expected that temperature-dependent MGPT potentials for
such cases will be added over time.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

The options defaults for the :doc:`pair_coeff <pair_coeff>` command are
volpress yes, nbody 1234, and precision double.


----------


.. _Moriarty1:



**(Moriarty1)** Moriarty, Physical Review B, 38, 3199 (1988).

.. _Moriarty2:



**(Moriarty2)** Moriarty, Physical Review B, 42, 1609 (1990).
Moriarty, Physical Review B 49, 12431 (1994).

.. _Moriarty3:



**(Moriarty3)** Moriarty, Benedict, Glosli, Hood, Orlikowski, Patel, Soderlind, Streitz, Tang, and Yang,
Journal of Materials Research, 21, 563 (2006).

.. _Glosli:



**(Glosli)** Glosli, unpublished, 2005.
Streitz, Glosli, Patel, Chan, Yates, de Supinski, Sexton and Gunnels, Journal of Physics: Conference
Series, 46, 254 (2006).

.. _Oppelstrup:



**(Oppelstrup)** Oppelstrup, unpublished, 2015.
Oppelstrup and Moriarty, to be published.
