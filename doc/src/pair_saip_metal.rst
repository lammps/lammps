.. index:: pair_style saip/metal
.. index:: pair_style saip/metal/opt
.. index:: pair_style saip/metal/tmd
.. index:: pair_style saip/metal/tmd/opt

pair_style saip/metal command
===================================

Accelerator Variant: *saip/metal/opt*

pair_style saip/metal/tmd command
===================================

Accelerator Variant: *saip/metal/tmd/opt*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style [hybrid/overlay ...] style cutoff tap_flag

* style = *saip/metal* or *saip/metal/tmd*
* cutoff = global cutoff (distance units)
* tap_flag = 0/1 to turn off/on the taper function

Examples
""""""""

.. code-block:: LAMMPS

   pair_style  hybrid/overlay saip/metal 16.0 1
   pair_coeff  * * saip/metal CHAu.ILP Au C H

   pair_style  hybrid/overlay eam rebo saip/metal 16.0
   pair_coeff  1 1 eam  Au_u3.eam  Au NULL NULL
   pair_coeff  * * rebo CH.rebo    NULL  C H
   pair_coeff  * * saip/metal      CHAu.ILP  Au C H

   pair_style  hybrid/overlay eam sw/mod saip/metal/tmd 16.0
   pair_coeff  1 1 eam  Au_u3.eam    Au NULL NULL
   pair_coeff  * * sw/mod tmd.sw.mod NULL  S Mo S
   pair_coeff  * * saip/metal/tmd    TMDAu.SAIP  Au S Mo S

Description
"""""""""""

.. versionadded:: 17Feb2022

The *saip/metal* style computes the semi-anisotropic interfacial
potential (SAIP) potential for hetero-junctions formed with hexagonal
2D materials and metal surfaces, as described in :ref:`(Ouyang6) <Ouyang6>` and :ref:`(Yao1) <Yao1>`.

.. versionadded:: 10Dec2025

The *saip/metal/tmd* style computes the semi-anisotropic interfacial
potential (SAIP) potential for hetero-junctions formed with transition
metal dichalcogenides (TMDCs) and metal surfaces, as described in :ref:`(Yao2) <Yao2>`.

.. math::

   E  = & \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
   V_{ij}  = & \mathrm{Tap}(r_{ij})\left \{ e^{-\alpha (r_{ij}/\beta -1)}
                \left [ \epsilon + f(\rho_{ij}) + f(\rho_{ji})\right ] -
                 \frac{1}{1+e^{-d\left [ \left ( r_{ij}/\left (s_R \cdot r^{eff} \right ) \right )-1 \right ]}}
                 \cdot \frac{C_6}{r^6_{ij}} \right \}\\
   \rho_{ij}^2 = & r_{ij}^2 - (\mathbf{r}_{ij} \cdot \mathbf{n}_i)^2 \\
   \rho_{ji}^2  = & r_{ij}^2 - (\mathbf{r}_{ij} \cdot \mathbf{n}_j)^2 \\
   f(\rho)  = &  C e^{ -( \rho / \delta )^2 } \\
   \mathrm{Tap}(r_{ij})  = & 20\left ( \frac{r_{ij}}{R_{cut}} \right )^7 -
                           70\left ( \frac{r_{ij}}{R_{cut}} \right )^6 +
                           84\left ( \frac{r_{ij}}{R_{cut}} \right )^5 -
                           35\left ( \frac{r_{ij}}{R_{cut}} \right )^4 + 1

Where :math:`\mathrm{Tap}(r_{ij})` is the taper function which provides
a continuous cutoff (up to third derivative) for interatomic separations
larger than :math:`r_c` :doc:`pair_style ilp_graphene_hbn <pair_ilp_graphene_hbn>`.

It is important to include all the pairs to build the neighbor list for
calculating the normals.

.. note::

   To account for the isotropic nature of the isolated gold atom
   electron cloud, their corresponding normal vectors (`\mathbf{n}_i`) are
   assumed to lie along the interatomic vector `\mathbf{r}_ij`. Notably, this
   assumption is suitable for many bulk material surfaces, for
   example, for systems possessing s-type valence orbitals or
   metallic surfaces, whose valence electrons are mostly
   delocalized, such that their Pauli repulsion with the electrons
   of adjacent surfaces are isotropic. Caution should be used in
   the case of very small gold contacts, for example, nano-clusters,
   where edge effects may become relevant.

The parameter file (e.g. CHAu.ILP), is intended for use with *metal*
:doc:`units <units>`, with energies in meV. Two additional parameters,
*S*, and *rcut* are included in the parameter file. *S* is designed to
facilitate scaling of energies. *rcut* is designed to build the neighbor
list for calculating the normals for each atom pair.

.. note::

   The parameters presented in the parameter file (e.g. TMDAu.SAIP),
   are fitted with taper function by setting the cutoff equal to 16.0
   Angstrom.  Using different cutoff or taper function should be careful.

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

   compute 0 all pair saip/metal
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

The CHAu.ILP and TMDAu.SAIP potential file provided with LAMMPS (see the potentials
directory) are parameterized for *metal* units.  You can use this
potential with any LAMMPS units, but you would need to create your own
custom CHAu.ILP potential file with coefficients listed in the appropriate
units, if your simulation does not use *metal* units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_none <pair_none>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`,
:doc:`pair_style drip <pair_drip>`,
:doc:`pair_style ilp_tmd <pair_ilp_tmd>`,
:doc:`pair_style ilp_graphene_hbn <pair_ilp_graphene_hbn>`,
:doc:`pair_style pair_kolmogorov_crespi_z <pair_kolmogorov_crespi_z>`,
:doc:`pair_style pair_kolmogorov_crespi_full <pair_kolmogorov_crespi_full>`,
:doc:`pair_style pair_lebedeva_z <pair_lebedeva_z>`,
:doc:`pair_style pair_coul_shield <pair_coul_shield>`.

Default
"""""""

tap_flag = 1


----------

.. _Ouyang6:

**(Ouyang6)** W. Ouyang, O. Hod, and R. Guerra, J. Chem. Theory Comput. 17, 7215 (2021).

.. _Yao1:

**(Yao1)** Y. Yao, ..., W. Ouyang, J. Phys. Chem. C, 128, 6836 (2024).

.. _Yao2:

**(Yao2)** Y. Yao, ..., W. Ouyang, Adv. Sci. 12, 2415884 (2025).
