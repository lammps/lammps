.. index:: pair\_style

pair\_style command
===================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = one of the styles from the list below
* args = arguments used by a particular style

Examples
""""""""


.. parsed-literal::

   pair_style lj/cut 2.5
   pair_style eam/alloy
   pair_style hybrid lj/charmm/coul/long 10.0 eam
   pair_style table linear 1000
   pair_style none

Description
"""""""""""

Set the formula(s) LAMMPS uses to compute pairwise interactions.  In
LAMMPS, pair potentials are defined between pairs of atoms that are
within a cutoff distance and the set of active interactions typically
changes over time.  See the :doc:`bond\_style <bond_style>` command to
define potentials between pairs of bonded atoms, which typically
remain in place for the duration of a simulation.

In LAMMPS, pairwise force fields encompass a variety of interactions,
some of which include many-body effects, e.g. EAM, Stillinger-Weber,
Tersoff, REBO potentials.  They are still classified as "pairwise"
potentials because the set of interacting atoms changes with time
(unlike molecular bonds) and thus a neighbor list is used to find
nearby interacting atoms.

Hybrid models where specified pairs of atom types interact via
different pair potentials can be setup using the *hybrid* pair style.

The coefficients associated with a pair style are typically set for
each pair of atom types, and are specified by the
:doc:`pair\_coeff <pair_coeff>` command or read from a file by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands.

The :doc:`pair\_modify <pair_modify>` command sets options for mixing of
type I-J interaction coefficients and adding energy offsets or tail
corrections to Lennard-Jones potentials.  Details on these options as
they pertain to individual potentials are described on the doc page
for the potential.  Likewise, info on whether the potential
information is stored in a :doc:`restart file <write_restart>` is listed
on the potential doc page.

In the formulas listed for each pair style, *E* is the energy of a
pairwise interaction between two atoms separated by a distance *r*\ .
The force between the atoms is the negative derivative of this
expression.

If the pair\_style command has a cutoff argument, it sets global
cutoffs for all pairs of atom types.  The distance(s) can be smaller
or larger than the dimensions of the simulation box.

Typically, the global cutoff value can be overridden for a specific
pair of atom types by the :doc:`pair\_coeff <pair_coeff>` command.  The
pair style settings (including global cutoffs) can be changed by a
subsequent pair\_style command using the same style.  This will reset
the cutoffs for all atom type pairs, including those previously set
explicitly by a :doc:`pair\_coeff <pair_coeff>` command.  The exceptions
to this are that pair\_style *table* and *hybrid* settings cannot be
reset.  A new pair\_style command for these styles will wipe out all
previously specified pair\_coeff values.


----------


Here is an alphabetic list of pair styles defined in LAMMPS.  They are
also listed in more compact form on the :doc:`Commands pair <Commands_pair>` doc page.

Click on the style to display the formula it computes, any additional
arguments specified in the pair\_style command, and coefficients
specified by the associated :doc:`pair\_coeff <pair_coeff>` command.

There are also additional accelerated pair styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :doc:`Commands pair <Commands_pair>`
doc page are followed by one or more of (g,i,k,o,t) to indicate which
accelerated styles exist.

* :doc:`none <pair_none>` - turn off pairwise interactions
* :doc:`hybrid <pair_hybrid>` - multiple styles of pairwise interactions
* :doc:`hybrid/overlay <pair_hybrid>` - multiple styles of superposed pairwise interactions
* :doc:`zero <pair_zero>` - neighbor list but no interactions

* :doc:`adp <pair_adp>` - angular dependent potential (ADP) of Mishin
* :doc:`agni <pair_agni>` - machine learned potential mapping atomic environment to forces
* :doc:`airebo <pair_airebo>` - AIREBO potential of Stuart
* :doc:`airebo/morse <pair_airebo>` - AIREBO with Morse instead of LJ
* :doc:`atm <pair_atm>` - Axilrod-Teller-Muto potential
* :doc:`awpmd/cut <pair_awpmd>` - Antisymmetrized Wave Packet MD potential for atoms and electrons
* :doc:`beck <pair_beck>` - Beck potential
* :doc:`body/nparticle <pair_body_nparticle>` - interactions between body particles
* :doc:`body/rounded/polygon <pair_body_rounded_polygon>` - granular-style 2d polygon potential
* :doc:`body/rounded/polyhedron <pair_body_rounded_polyhedron>` - granular-style 3d polyhedron potential
* :doc:`bop <pair_bop>` - BOP potential of Pettifor
* :doc:`born <pair_born>` - Born-Mayer-Huggins potential
* :doc:`born/coul/dsf <pair_born>` - Born with damped-shifted-force model
* :doc:`born/coul/dsf/cs <pair_cs>` - Born with damped-shifted-force and core/shell model
* :doc:`born/coul/long <pair_born>` - Born with long-range Coulombics
* :doc:`born/coul/long/cs <pair_cs>` - Born with long-range Coulombics and core/shell
* :doc:`born/coul/msm <pair_born>` - Born with long-range MSM Coulombics
* :doc:`born/coul/wolf <pair_born>` - Born with Wolf potential for Coulombics
* :doc:`born/coul/wolf/cs <pair_cs>` - Born with Wolf potential for Coulombics and core/shell model
* :doc:`brownian <pair_brownian>` - Brownian potential for Fast Lubrication Dynamics
* :doc:`brownian/poly <pair_brownian>` - Brownian potential for Fast Lubrication Dynamics with polydispersity
* :doc:`buck <pair_buck>` - Buckingham potential
* :doc:`buck/coul/cut <pair_buck>` - Buckingham with cutoff Coulomb
* :doc:`buck/coul/long <pair_buck>` - Buckingham with long-range Coulombics
* :doc:`buck/coul/long/cs <pair_cs>` - Buckingham with long-range Coulombics and core/shell
* :doc:`buck/coul/msm <pair_buck>` - Buckingham with long-range MSM Coulombics
* :doc:`buck/long/coul/long <pair_buck_long>` - long-range Buckingham with long-range Coulombics
* :doc:`buck/mdf <pair_mdf>` - Buckingham with a taper function
* :doc:`buck6d/coul/gauss/dsf <pair_buck6d_coul_gauss>` - dispersion-damped Buckingham with damped-shift-force model
* :doc:`buck6d/coul/gauss/long <pair_buck6d_coul_gauss>` - dispersion-damped Buckingham with long-range Coulombics
* :doc:`colloid <pair_colloid>` - integrated colloidal potential
* :doc:`comb <pair_comb>` - charge-optimized many-body (COMB) potential
* :doc:`comb3 <pair_comb>` - charge-optimized many-body (COMB3) potential
* :doc:`cosine/squared <pair_cosine_squared>` - Cooke-Kremer-Deserno membrane model potential
* :doc:`coul/cut <pair_coul>` - cutoff Coulombic potential
* :doc:`coul/cut/soft <pair_fep_soft>` - Coulombic potential with a soft core
* :doc:`coul/debye <pair_coul>` - cutoff Coulombic potential with Debye screening
* :doc:`coul/diel <pair_coul_diel>` - Coulomb potential with dielectric permittivity
* :doc:`coul/dsf <pair_coul>` - Coulombics with damped-shifted-force model
* :doc:`coul/long <pair_coul>` - long-range Coulombic potential
* :doc:`coul/long/cs <pair_cs>` - long-range Coulombic potential and core/shell
* :doc:`coul/long/soft <pair_fep_soft>` - long-range Coulombic potential with a soft core
* :doc:`coul/msm <pair_coul>` - long-range MSM Coulombics
* :doc:`coul/shield <pair_coul_shield>` - Coulombics for boron nitride for use with :doc:`ilp/graphene/hbn <pair_ilp_graphene_hbn>` potential
* :doc:`coul/streitz <pair_coul>` - Coulombics via Streitz/Mintmire Slater orbitals
* :doc:`coul/wolf <pair_coul>` - Coulombics via Wolf potential
* :doc:`coul/wolf/cs <pair_cs>` - ditto with core/shell adjustments
* :doc:`dpd <pair_dpd>` - dissipative particle dynamics (DPD)
* :doc:`dpd/fdt <pair_dpd_fdt>` - DPD for constant temperature and pressure
* :doc:`dpd/fdt/energy <pair_dpd_fdt>` - DPD for constant energy and enthalpy
* :doc:`dpd/tstat <pair_dpd>` - pair-wise DPD thermostatting
* :doc:`dsmc <pair_dsmc>` - Direct Simulation Monte Carlo (DSMC)
* :doc:`e3b <pair_e3b>` - Explicit-three body (E3B) water model
* :doc:`drip <pair_drip>` - Dihedral-angle-corrected registry-dependent interlayer potential (DRIP)
* :doc:`eam <pair_eam>` - embedded atom method (EAM)
* :doc:`eam/alloy <pair_eam>` - alloy EAM
* :doc:`eam/cd <pair_eam>` - concentration-dependent EAM
* :doc:`eam/cd/old <pair_eam>` - older two-site model for concentration-dependent EAM
* :doc:`eam/fs <pair_eam>` - Finnis-Sinclair EAM
* :doc:`edip <pair_edip>` - three-body EDIP potential
* :doc:`edip/multi <pair_edip>` - multi-element EDIP potential
* :doc:`edpd <pair_meso>` - eDPD particle interactions
* :doc:`eff/cut <pair_eff>` - electron force field with a cutoff
* :doc:`eim <pair_eim>` - embedded ion method (EIM)
* :doc:`exp6/rx <pair_exp6_rx>` - reactive DPD potential
* :doc:`extep <pair_extep>` - extended Tersoff potential
* :doc:`gauss <pair_gauss>` - Gaussian potential
* :doc:`gauss/cut <pair_gauss>` - generalized Gaussian potential
* :doc:`gayberne <pair_gayberne>` - Gay-Berne ellipsoidal potential
* :doc:`gran/hertz/history <pair_gran>` - granular potential with Hertzian interactions
* :doc:`gran/hooke <pair_gran>` - granular potential with history effects
* :doc:`gran/hooke/history <pair_gran>` - granular potential without history effects
* :doc:`gw <pair_gw>` - Gao-Weber potential
* :doc:`gw/zbl <pair_gw>` - Gao-Weber potential with a repulsive ZBL core
* :doc:`hbond/dreiding/lj <pair_hbond_dreiding>` - DREIDING hydrogen bonding LJ potential
* :doc:`hbond/dreiding/morse <pair_hbond_dreiding>` - DREIDING hydrogen bonding Morse potential
* :doc:`ilp/graphene/hbn <pair_ilp_graphene_hbn>` - registry-dependent interlayer potential (ILP)
* :doc:`kim <pair_kim>` - interface to potentials provided by KIM project
* :doc:`kolmogorov/crespi/full <pair_kolmogorov_crespi_full>` - Kolmogorov-Crespi (KC) potential with no simplifications
* :doc:`kolmogorov/crespi/z <pair_kolmogorov_crespi_z>` - Kolmogorov-Crespi (KC) potential with normals along z-axis
* :doc:`lcbop <pair_lcbop>` - long-range bond-order potential (LCBOP)
* :doc:`lebedeva/z <pair_lebedeva_z>` - Lebedeva interlayer potential for graphene with normals along z-axis
* :doc:`lennard/mdf <pair_mdf>` - LJ potential in A/B form with a taper function
* :doc:`line/lj <pair_line_lj>` - LJ potential between line segments
* :doc:`list <pair_list>` - potential between pairs of atoms explicitly listed in an input file
* :doc:`lj/charmm/coul/charmm <pair_charmm>` - CHARMM potential with cutoff Coulomb
* :doc:`lj/charmm/coul/charmm/implicit <pair_charmm>` - CHARMM for implicit solvent
* :doc:`lj/charmm/coul/long <pair_charmm>` - CHARMM with long-range Coulomb
* :doc:`lj/charmm/coul/long/soft <pair_fep_soft>` - CHARMM with long-range Coulomb and a soft core
* :doc:`lj/charmm/coul/msm <pair_charmm>` - CHARMM with long-range MSM Coulombics
* :doc:`lj/charmmfsw/coul/charmmfsh <pair_charmm>` - CHARMM with force switching and shifting
* :doc:`lj/charmmfsw/coul/long <pair_charmm>` - CHARMM with force switching and long-rnage Coulombics
* :doc:`lj/class2 <pair_class2>` - COMPASS (class 2) force field with no Coulomb
* :doc:`lj/class2/coul/cut <pair_class2>` - COMPASS with cutoff Coulomb
* :doc:`lj/class2/coul/cut/soft <pair_fep_soft>` - COMPASS with cutoff Coulomb with a soft core
* :doc:`lj/class2/coul/long <pair_class2>` - COMPASS with long-range Coulomb
* :doc:`lj/class2/coul/long/soft <pair_fep_soft>` - COMPASS with long-range Coulomb with a soft core
* :doc:`lj/class2/soft <pair_fep_soft>` - COMPASS (class 2) force field with no Coulomb with a soft core
* :doc:`lj/cubic <pair_lj_cubic>` - LJ with cubic after inflection point
* :doc:`lj/cut <pair_lj>` - cutoff Lennard-Jones potential with no Coulomb
* :doc:`lj/cut/coul/cut <pair_lj>` - LJ with cutoff Coulomb
* :doc:`lj/cut/coul/cut/soft <pair_fep_soft>` - LJ with cutoff Coulomb with a soft core
* :doc:`lj/cut/coul/debye <pair_lj>` - LJ with Debye screening added to Coulomb
* :doc:`lj/cut/coul/dsf <pair_lj>` - LJ with Coulombics via damped shifted forces
* :doc:`lj/cut/coul/long <pair_lj>` - LJ with long-range Coulombics
* :doc:`lj/cut/coul/long/cs <pair_cs>` - ditto with core/shell adjustments
* :doc:`lj/cut/coul/long/soft <pair_fep_soft>` - LJ with long-range Coulombics with a soft core
* :doc:`lj/cut/coul/msm <pair_lj>` - LJ with long-range MSM Coulombics
* :doc:`lj/cut/coul/wolf <pair_lj>` - LJ with Coulombics via Wolf potential
* :doc:`lj/cut/dipole/cut <pair_dipole>` - point dipoles with cutoff
* :doc:`lj/cut/dipole/long <pair_dipole>` - point dipoles with long-range Ewald
* :doc:`lj/cut/soft <pair_fep_soft>` - LJ with a soft core
* :doc:`lj/cut/thole/long <pair_thole>` - LJ with Coulombics with thole damping
* :doc:`lj/cut/tip4p/cut <pair_lj>` - LJ with cutoff Coulomb for TIP4P water
* :doc:`lj/cut/tip4p/long <pair_lj>` - LJ with long-range Coulomb for TIP4P water
* :doc:`lj/cut/tip4p/long/soft <pair_fep_soft>` - LJ with cutoff Coulomb for TIP4P water with a soft core
* :doc:`lj/expand <pair_lj_expand>` - Lennard-Jones for variable size particles
* :doc:`lj/expand/coul/long <pair_lj_expand>` - Lennard-Jones for variable size particles with long-range Coulombics
* :doc:`lj/gromacs <pair_gromacs>` - GROMACS-style Lennard-Jones potential
* :doc:`lj/gromacs/coul/gromacs <pair_gromacs>` - GROMACS-style LJ and Coulombic potential
* :doc:`lj/long/coul/long <pair_lj_long>` - long-range LJ and long-range Coulombics
* :doc:`lj/long/dipole/long <pair_dipole>` - long-range LJ and long-range point dipoles
* :doc:`lj/long/tip4p/long <pair_lj_long>` - long-range LJ and long-range Coulombics for TIP4P water
* :doc:`lj/mdf <pair_mdf>` - LJ potential with a taper function
* :doc:`lj/sdk <pair_sdk>` - LJ for SDK coarse-graining
* :doc:`lj/sdk/coul/long <pair_sdk>` - LJ for SDK coarse-graining with long-range Coulombics
* :doc:`lj/sdk/coul/msm <pair_sdk>` - LJ for SDK coarse-graining with long-range Coulombics via MSM
* :doc:`lj/sf/dipole/sf <pair_dipole>` - LJ with dipole interaction with shifted forces
* :doc:`lj/smooth <pair_lj_smooth>` - smoothed Lennard-Jones potential
* :doc:`lj/smooth/linear <pair_lj_smooth_linear>` - linear smoothed LJ potential
* `lj/switch3/coulgauss <pair_lj_switch3_coulgauss>`_ - smoothed LJ vdW potential with Gaussian electrostatics
* :doc:`lj96/cut <pair_lj96>` - Lennard-Jones 9/6 potential
* :doc:`local/density <pair_local_density>` - generalized basic local density potential
* :doc:`lubricate <pair_lubricate>` - hydrodynamic lubrication forces
* :doc:`lubricate/poly <pair_lubricate>` - hydrodynamic lubrication forces with polydispersity
* :doc:`lubricateU <pair_lubricateU>` - hydrodynamic lubrication forces for Fast Lubrication Dynamics
* :doc:`lubricateU/poly <pair_lubricateU>` - hydrodynamic lubrication forces for Fast Lubrication with polydispersity
* :doc:`mdpd <pair_meso>` - mDPD particle interactions
* :doc:`mdpd/rhosum <pair_meso>` - mDPD particle interactions for mass density
* :doc:`meam/c <pair_meamc>` - modified embedded atom method (MEAM) in C
* :doc:`meam/spline <pair_meam_spline>` - splined version of MEAM
* :doc:`meam/sw/spline <pair_meam_sw_spline>` - splined version of MEAM with a Stillinger-Weber term
* :doc:`mgpt <pair_mgpt>` - simplified model generalized pseudopotential theory (MGPT) potential
* :doc:`mie/cut <pair_mie>` - Mie potential
* :doc:`mm3/switch3/coulgauss/long <pair_mm3_switch3_coulgauss_long>` - smoothed MM3 vdW potential with Gaussian electrostatics
* :doc:`momb <pair_momb>` - Many-Body Metal-Organic (MOMB) force field
* :doc:`morse <pair_morse>` - Morse potential
* :doc:`morse/smooth/linear <pair_morse>` - linear smoothed Morse potential
* :doc:`morse/soft <pair_morse>` - Morse potential with a soft core
* :doc:`multi/lucy <pair_multi_lucy>` - DPD potential with density-dependent force
* :doc:`multi/lucy/rx <pair_multi_lucy_rx>` - reactive DPD potential with density-dependent force
* :doc:`nb3b/harmonic <pair_nb3b_harmonic>` - non-bonded 3-body harmonic potential
* :doc:`nm/cut <pair_nm>` - N-M potential
* :doc:`nm/cut/coul/cut <pair_nm>` - N-M potential with cutoff Coulomb
* :doc:`nm/cut/coul/long <pair_nm>` - N-M potential with long-range Coulombics
* :doc:`oxdna/coaxstk <pair_oxdna>` -
* :doc:`oxdna/excv <pair_oxdna>` -
* :doc:`oxdna/hbond <pair_oxdna>` -
* :doc:`oxdna/stk <pair_oxdna>` -
* :doc:`oxdna/xstk <pair_oxdna>` -
* :doc:`oxdna2/coaxstk <pair_oxdna2>` -
* :doc:`oxdna2/dh <pair_oxdna2>` -
* :doc:`oxdna2/excv <pair_oxdna2>` -
* :doc:`oxdna2/hbond <pair_oxdna2>` -
* :doc:`oxdna2/stk <pair_oxdna2>` -
* :doc:`oxdna2/xstk <pair_oxdna2>` -
* :doc:`peri/eps <pair_peri>` - peridynamic EPS potential
* :doc:`peri/lps <pair_peri>` - peridynamic LPS potential
* :doc:`peri/pmb <pair_peri>` - peridynamic PMB potential
* :doc:`peri/ves <pair_peri>` - peridynamic VES potential
* :doc:`polymorphic <pair_polymorphic>` - polymorphic 3-body potential
* :doc:`python <pair_python>` -
* :doc:`quip <pair_quip>` -
* :doc:`reax/c <pair_reaxc>` - ReaxFF potential in C
* :doc:`rebo <pair_airebo>` - 2nd generation REBO potential of Brenner
* :doc:`resquared <pair_resquared>` - Everaers RE-Squared ellipsoidal potential
* :doc:`sdpd/taitwater/isothermal <pair_sdpd_taitwater_isothermal>` - smoothed dissipative particle dynamics for water at isothermal conditions
* :doc:`smd/hertz <pair_smd_hertz>` -
* :doc:`smd/tlsph <pair_smd_tlsph>` -
* :doc:`smd/tri\_surface <pair_smd_triangulated_surface>` -
* :doc:`smd/ulsph <pair_smd_ulsph>` -
* :doc:`smtbq <pair_smtbq>` -
* :doc:`snap <pair_snap>` - SNAP quantum-accurate potential
* :doc:`soft <pair_soft>` - Soft (cosine) potential
* :doc:`sph/heatconduction <pair_sph_heatconduction>` -
* :doc:`sph/idealgas <pair_sph_idealgas>` -
* :doc:`sph/lj <pair_sph_lj>` -
* :doc:`sph/rhosum <pair_sph_rhosum>` -
* :doc:`sph/taitwater <pair_sph_taitwater>` -
* :doc:`sph/taitwater/morris <pair_sph_taitwater_morris>` -
* :doc:`spin/dipole/cut <pair_spin_dipole>` -
* :doc:`spin/dipole/long <pair_spin_dipole>` -
* :doc:`spin/dmi <pair_spin_dmi>` -
* :doc:`spin/exchange <pair_spin_exchange>` -
* :doc:`spin/magelec <pair_spin_magelec>` -
* :doc:`spin/neel <pair_spin_neel>` -
* :doc:`srp <pair_srp>` -
* :doc:`sw <pair_sw>` - Stillinger-Weber 3-body potential
* :doc:`table <pair_table>` - tabulated pair potential
* :doc:`table/rx <pair_table_rx>` -
* :doc:`tdpd <pair_meso>` - tDPD particle interactions
* :doc:`tersoff <pair_tersoff>` - Tersoff 3-body potential
* :doc:`tersoff/mod <pair_tersoff_mod>` - modified Tersoff 3-body potential
* :doc:`tersoff/mod/c <pair_tersoff_mod>` -
* :doc:`tersoff/table <pair_tersoff>` -
* :doc:`tersoff/zbl <pair_tersoff_zbl>` - Tersoff/ZBL 3-body potential
* :doc:`thole <pair_thole>` - Coulomb interactions with thole damping
* :doc:`tip4p/cut <pair_coul>` - Coulomb for TIP4P water w/out LJ
* :doc:`tip4p/long <pair_coul>` - long-range Coulombics for TIP4P water w/out LJ
* :doc:`tip4p/long/soft <pair_fep_soft>` -
* :doc:`tri/lj <pair_tri_lj>` - LJ potential between triangles
* :doc:`ufm <pair_ufm>` -
* :doc:`vashishta <pair_vashishta>` - Vashishta 2-body and 3-body potential
* :doc:`vashishta/table <pair_vashishta>` -
* :doc:`yukawa <pair_yukawa>` - Yukawa potential
* :doc:`yukawa/colloid <pair_yukawa_colloid>` - screened Yukawa potential for finite-size particles
* :doc:`zbl <pair_zbl>` - Ziegler-Biersack-Littmark potential


----------


Restrictions
""""""""""""


This command must be used before any coefficients are set by the
:doc:`pair\_coeff <pair_coeff>`, :doc:`read\_data <read_data>`, or
:doc:`read\_restart <read_restart>` commands.

Some pair styles are part of specific packages.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.  The doc pages for
individual pair potentials tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`read\_data <read_data>`,
:doc:`pair\_modify <pair_modify>`, :doc:`kspace\_style <kspace_style>`,
:doc:`dielectric <dielectric>`, :doc:`pair\_write <pair_write>`

Default
"""""""


.. parsed-literal::

   pair_style none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
