.. index:: pair_style ilp/tmd
.. index:: pair_style ilp/tmd/opt

pair_style ilp/tmd command
===================================

Accelerator Variant: *ilp/tmd/opt*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style [hybrid/overlay ...] ilp/tmd cutoff tap_flag

* cutoff = global cutoff (distance units)
* tap_flag = 0/1 to turn off/on the taper function

Examples
""""""""

.. code-block:: LAMMPS

   pair_style  hybrid/overlay ilp/tmd 16.0 1
   pair_coeff  * * ilp/tmd  TMD.ILP Mo S S

   pair_style  hybrid/overlay sw/mod sw/mod ilp/tmd 16.0
   pair_coeff  * * sw/mod 1  tmd.sw.mod Mo S S NULL NULL NULL
   pair_coeff  * * sw/mod 2  tmd.sw.mod NULL NULL NULL Mo S S
   pair_coeff  * * ilp/tmd   TMD.ILP    Mo S S Mo S S

Description
"""""""""""

.. versionadded:: 17Feb2022

The *ilp/tmd* style computes the registry-dependent interlayer
potential (ILP) potential for transition metal dichalcogenides (TMD)
as described in :ref:`(Ouyang7) <Ouyang7>`.

.. math::

   E  = & \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
   V_{ij}  = & {\rm Tap}(r_{ij})\left \{ e^{-\alpha (r_{ij}/\beta -1)}
                \left [ \epsilon + f(\rho_{ij}) + f(\rho_{ji})\right ] -
                 \frac{1}{1+e^{-d\left [ \left ( r_{ij}/\left (s_R \cdot r^{eff} \right ) \right )-1 \right ]}}
                 \cdot \frac{C_6}{r^6_{ij}} \right \}\\
   \rho_{ij}^2 = & r_{ij}^2 - ({\bf r}_{ij} \cdot {\bf n}_i)^2 \\
   \rho_{ji}^2  = & r_{ij}^2 - ({\bf r}_{ij} \cdot {\bf n}_j)^2 \\
   f(\rho)  = &  C e^{ -( \rho / \delta )^2 } \\
   {\rm Tap}(r_{ij})  = & 20\left ( \frac{r_{ij}}{R_{cut}} \right )^7 -
                           70\left ( \frac{r_{ij}}{R_{cut}} \right )^6 +
                           84\left ( \frac{r_{ij}}{R_{cut}} \right )^5 -
                           35\left ( \frac{r_{ij}}{R_{cut}} \right )^4 + 1

Where :math:`\mathrm{Tap}(r_{ij})` is the taper function which provides
a continuous cutoff (up to third derivative) for interatomic separations
larger than :math:`r_c` :doc:`pair_style ilp_graphene_hbn <pair_ilp_graphene_hbn>`.

It is important to include all the pairs to build the neighbor list for
calculating the normals.

.. note::

   Since each MX2 (M = Mo, W and X = S, Se Te) layer contains two
   sub-layers of X atoms and one sub-layer of M atoms, the definition of the
   normal vectors used for graphene and h-BN is no longer valid for TMDs.
   In :ref:`(Ouyang7) <Ouyang7>`, a new definition is proposed, where for
   each atom `i`, its six nearest neighboring atoms belonging to the same
   sub-layer are chosen to define the normal vector `{\bf n}_i`.

The parameter file (e.g. TMD.ILP), is intended for use with *metal*
:doc:`units <units>`, with energies in meV. Two additional parameters,
*S*, and *rcut* are included in the parameter file. *S* is designed to
facilitate scaling of energies. *rcut* is designed to build the neighbor
list for calculating the normals for each atom pair.

.. note::

   The parameters presented in the parameter file (e.g. TMD.ILP),
   are fitted with taper function by setting the cutoff equal to 16.0
   Angstrom.  Using different cutoff or taper function should be careful.
   These parameters provide a good description in both short- and long-range
   interaction regimes. This feature is essential for simulations in high pressure
   regime (i.e., the interlayer distance is smaller than the equilibrium
   distance). The benchmark tests and comparison of these parameters can
   be found in :ref:`(Ouyang7) <Ouyang7>`.

This potential must be used in combination with hybrid/overlay.
Other interactions can be set to zero using pair_style *none*\ .

This pair style tallies a breakdown of the total interlayer potential
energy into sub-categories, which can be accessed via the :doc:`compute pair <compute_pair>` command as a vector of values of length 2.
The 2 values correspond to the following sub-categories:

1. *E_vdW* = vdW (attractive) energy
2. *E_Rep* = Repulsive energy

To print these quantities to the log file (with descriptive column
headings) the following commands could be included in an input script:

.. code-block:: LAMMPS

   compute 0 all pair ilp/tmd
   variable Evdw  equal c_0[1]
   variable Erep  equal c_0[2]
   thermo_style custom step temp epair v_Erep v_Evdw

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the pair_modify mix, shift, table, and
tail options.

This pair style does not write their information to binary restart
files, since it is stored in potential files. Thus, you need to
re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

Restrictions
""""""""""""

This pair style is part of the INTERLAYER package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires the newton setting to be *on* for pair
interactions.

The TMD.ILP potential file provided with LAMMPS (see the potentials
directory) are parameterized for *metal* units.  You can use this
potential with any LAMMPS units, but you would need to create your own
custom TMD.ILP potential file with coefficients listed in the appropriate
units, if your simulation does not use *metal* units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_none <pair_none>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`,
:doc:`pair_style drip <pair_drip>`,
:doc:`pair_style saip_metal <pair_saip_metal>`,
:doc:`pair_style ilp_graphene_hbn <pair_ilp_graphene_hbn>`,
:doc:`pair_style pair_kolmogorov_crespi_z <pair_kolmogorov_crespi_z>`,
:doc:`pair_style pair_kolmogorov_crespi_full <pair_kolmogorov_crespi_full>`,
:doc:`pair_style pair_lebedeva_z <pair_lebedeva_z>`,
:doc:`pair_style pair_coul_shield <pair_coul_shield>`.

Default
"""""""

tap_flag = 1


----------

.. _Ouyang7:

**(Ouyang7)** W. Ouyang, et al., J. Chem. Theory Comput. 17, 7237 (2021).
