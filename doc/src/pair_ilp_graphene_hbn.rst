.. index:: pair_style ilp/graphene/hbn

pair_style ilp/graphene/hbn command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style [hybrid/overlay ...] ilp/graphene/hbn cutoff tap_flag

* cutoff = global cutoff (distance units)
* tap\_flag = 0/1 to turn off/on the taper function

Examples
""""""""

.. code-block:: LAMMPS

   pair_style  hybrid/overlay ilp/graphene/hbn 16.0 1
   pair_coeff  * * ilp/graphene/hbn  BNCH.ILP B N C

   pair_style  hybrid/overlay rebo tersoff ilp/graphene/hbn 16.0 coul/shield 16.0
   pair_coeff  * * rebo              CH.rebo     NULL NULL C
   pair_coeff  * * tersoff           BNC.tersoff B    N    NULL
   pair_coeff  * * ilp/graphene/hbn  BNCH.ILP    B    N    C
   pair_coeff  1 1 coul/shield 0.70
   pair_coeff  1 2 coul/shield 0.695
   pair_coeff  2 2 coul/shield 0.69

Description
"""""""""""

The *ilp/graphene/hbn* style computes the registry-dependent interlayer
potential (ILP) potential as described in :ref:`(Leven1) <Leven1>`,
:ref:`(Leven2) <Leven2>` and :ref:`(Maaravi) <Maaravi2>`.
The normals are calculated in the way as described
in :ref:`(Kolmogorov) <Kolmogorov2>`.

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
larger than :math:`r_c` :ref:`(Maaravi) <Maaravi2>`. The definitions of
each parameter in the above equation can be found in :ref:`(Leven1)
<Leven1>` and :ref:`(Maaravi) <Maaravi2>`.

It is important to include all the pairs to build the neighbor list for
calculating the normals.

.. note::

   This potential (ILP) is intended for interlayer interactions between two
   different layers of graphene, hexagonal boron nitride (h-BN) and their hetero-junction.
   To perform a realistic simulation, this potential must be used in combination with
   intralayer potential, such as :doc:`AIREBO <pair_airebo>` or :doc:`Tersoff <pair_tersoff>` potential.
   To keep the intralayer properties unaffected, the interlayer interaction
   within the same layers should be avoided. Hence, each atom has to have a layer
   identifier such that atoms residing on the same layer interact via the
   appropriate intralayer potential and atoms residing on different layers
   interact via the ILP. Here, the molecule id is chosen as the layer identifier,
   thus a data file with the "full" atom style is required to use this potential.

The parameter file (e.g. BNCH.ILP), is intended for use with *metal*
:doc:`units <units>`, with energies in meV. Two additional parameters,
*S*\ , and *rcut* are included in the parameter file. *S* is designed to
facilitate scaling of energies. *rcut* is designed to build the neighbor
list for calculating the normals for each atom pair.

.. note::

   The parameters presented in the parameter file (e.g. BNCH.ILP),
   are fitted with taper function by setting the cutoff equal to 16.0
   Angstrom.  Using different cutoff or taper function should be careful.
   The parameters for atoms pairs between Boron and Nitrogen are fitted with
   a screened Coulomb interaction :doc:`coul/shield <pair_coul_shield>`. Therefore,
   to simulated the properties of h-BN correctly, this potential must be used in
   combination with the pair style :doc:`coul/shield <pair_coul_shield>`.

.. note::

   Four new sets of parameters of ILP for 2D layered Materials with bilayer and
   bulk configurations are presented in :ref:`(Ouyang1) <Ouyang1>` and :ref:`(Ouyang2) <Ouyang2>`, respectively.
   These parameters provide a good description in both short- and long-range interaction regimes.
   While the old ILP parameters published in :ref:`(Leven2) <Leven2>` and
   :ref:`(Maaravi) <Maaravi2>` are only suitable for long-range interaction
   regime. This feature is essential for simulations in high pressure
   regime (i.e., the interlayer distance is smaller than the equilibrium
   distance). The benchmark tests and comparison of these parameters can
   be found in :ref:`(Ouyang1) <Ouyang1>` and :ref:`(Ouyang2) <Ouyang2>`.

This potential must be used in combination with hybrid/overlay.
Other interactions can be set to zero using pair\_style *none*\ .

This pair style tallies a breakdown of the total interlayer potential
energy into sub-categories, which can be accessed via the :doc:`compute pair <compute_pair>` command as a vector of values of length 2.
The 2 values correspond to the following sub-categories:

1. *E\_vdW* = vdW (attractive) energy
2. *E\_Rep* = Repulsive energy

To print these quantities to the log file (with descriptive column
headings) the following commands could be included in an input script:

.. code-block:: LAMMPS

   compute 0 all pair ilp/graphene/hbn
   variable Evdw  equal c_0[1]
   variable Erep  equal c_0[2]
   thermo_style custom step temp epair v_Erep v_Evdw

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the pair\_modify mix, shift, table, and
tail options.

This pair style does not write their information to binary restart
files, since it is stored in potential files. Thus, you need to
re-specify the pair\_style and pair\_coeff commands in an input script
that reads a restart file.

Restrictions
""""""""""""

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair potential requires the newton setting to be *on* for pair
interactions.

The BNCH.ILP potential file provided with LAMMPS (see the potentials
directory) are parameterized for *metal* units.  You can use this
potential with any LAMMPS units, but you would need to create your
BNCH.ILP potential file with coefficients listed in the appropriate
units, if your simulation does not use *metal* units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_none <pair_none>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`,
:doc:`pair_style drip <pair_drip>`,
:doc:`pair_style pair\_kolmogorov\_crespi\_z <pair_kolmogorov_crespi_z>`,
:doc:`pair_style pair\_kolmogorov\_crespi\_full <pair_kolmogorov_crespi_full>`,
:doc:`pair_style pair\_lebedeva\_z <pair_lebedeva_z>`,
:doc:`pair_style pair\_coul\_shield <pair_coul_shield>`.

**Default:** tap\_flag = 1

----------

.. _Leven1:

**(Leven1)** I. Leven, I. Azuri, L. Kronik and O. Hod, J. Chem. Phys. 140, 104106 (2014).

.. _Leven2:

**(Leven2)** I. Leven et al, J. Chem.Theory Comput. 12, 2896-905 (2016).

.. _Maaravi2:

**(Maaravi)** T. Maaravi et al, J. Phys. Chem. C 121, 22826-22835 (2017).

.. _Kolmogorov2:

**(Kolmogorov)** A. N. Kolmogorov, V. H. Crespi, Phys. Rev. B 71, 235415 (2005).

.. _Ouyang1:

**(Ouyang1)** W. Ouyang, D. Mandelli, M. Urbakh and O. Hod, Nano Lett. 18, 6009-6016 (2018).

.. _Ouyang2:

**(Ouyang2)** W. Ouyang et al., J. Chem. Theory Comput. 16(1), 666-676 (2020).
