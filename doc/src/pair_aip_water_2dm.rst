.. index:: pair_style aip/water/2dm
.. index:: pair_style aip/water/2dm/opt

pair_style aip/water/2dm command
===================================

Accelerator Variant: *aip/water/2dm/opt*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style [hybrid/overlay ...] aip/water/2dm cutoff tap_flag

* cutoff = global cutoff (distance units)
* tap_flag = 0/1 to turn off/on the taper function

Examples
""""""""

.. code-block:: LAMMPS

   pair_style  hybrid/overlay aip/water/2dm 16.0 1
   pair_coeff  * * aip/water/2dm  COH.aip.water.2dm C Ow Hw

   pair_style  hybrid/overlay aip/water/2dm 16.0 lj/cut/tip4p/long 2 3 1 1 0.1546 10 8.5
   pair_coeff  2 2   lj/cut/tip4p/long    8.0313e-3  3.1589  # O-O
   pair_coeff  2 3   lj/cut/tip4p/long    0.0        0.0     # O-H
   pair_coeff  3 3   lj/cut/tip4p/long    0.0        0.0     # H-H
   pair_coeff  * *   aip/water/2dm        COH.aip.water.2dm    C Ow Hw

Description
"""""""""""

.. versionadded:: TBD

The *aip/water/2dm* style computes the anisotropic interfacial potential
(AIP) potential for interfaces of water with two-dimensional (2D)
materials as described in :ref:`(Feng) <Feng>`.

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
larger than :math:`r_c` :doc:`pair_style ilp_graphene_hbn
<pair_ilp_graphene_hbn>`.

.. note::

   This pair style uses the atomic normal vector definition from
   :ref:`(Feng) <Feng>`), where the atomic normal vectors of the
   hydrogen atoms are assumed to lie along the corresponding
   oxygen-hydrogen bonds and the normal vector of the central oxygen
   atom is defined as their average.

The provided parameter file, ``COH.aip.water.2dm``, is intended for use
with *metal* :doc:`units <units>`, with energies in meV.  Two additional
parameters, *S*, and *rcut* are included in the parameter file. *S* is
designed to facilitate scaling of energies; *rcut* is the cutoff for an
internal, short distance neighbor list that is generated for speeding up
the calculation of the normals for all atom pairs.

.. note::

   The parameters presented in the provided parameter file,
   ``COH.aip.water.2dm``, are fitted with the taper function enabled by
   setting the cutoff equal to 16.0 Angstrom.  Using a different cutoff
   or taper function setting should be carefully checked as they can
   lead to significant errors.  These parameters provide a good
   description in both short- and long-range interaction regimes. This
   is essential for simulations in high pressure regime (i.e., the
   interlayer distance is smaller than the equilibrium distance).

This potential must be used in combination with hybrid/overlay.  Other
interactions can be set to zero using :doc:`pair_coeff settings
<pair_coeff>` with the pair style set to *none*\ .

This pair style tallies a breakdown of the total interlayer potential
energy into sub-categories, which can be accessed via the :doc:`compute
pair <compute_pair>` command as a vector of values of length 2.  The 2
values correspond to the following sub-categories:

1. *E_vdW* = vdW (attractive) energy
2. *E_Rep* = Repulsive energy

To print these quantities to the log file (with descriptive column
headings) the following commands could be included in an input script:

.. code-block:: LAMMPS

   compute 0 all pair aip/water/2dm
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

The ``COH.aip.water.2dm`` potential file provided with LAMMPS is
parameterized for *metal* units.  You can use this pair style with any
LAMMPS units, but you would need to create your own potential file with
parameters in the appropriate units, if your simulation does not use
*metal* units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_none <pair_none>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`,
:doc:`pair_style drip <pair_drip>`,
:doc:`pair_style ilp_tmd <pair_ilp_tmd>`,
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

.. _Feng:

**(Feng)** Z. Feng and W. Ouyang et al., J. Phys. Chem. C. 127, 8704-8713 (2023).
