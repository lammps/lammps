Available Packages
==================

This is the list of packages included in LAMMPS.  The link for each
package name gives more details.

Packages are supported by either the LAMMPS developers or the
contributing authors and written in a syntax and style consistent with
the rest of LAMMPS.

The "Examples" column is a subdirectory in the examples directory of the
distribution which has one or more input scripts that use the package.
E.g. "peptide" refers to the examples/peptide directory; PACKAGES/atc refers
to the examples/PACKAGES/atc directory.  The "Lib" column indicates
whether an extra library is needed to build and use the package:

* no  = no library
* sys = system library: you likely have it on your machine
* int = internal library: provided with LAMMPS, but you may need to build it
* ext = external library: you will need to download and install it on your machine

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Package
     - Description
     - Doc page
     - Examples
     - Lib
   * - :ref:`ADIOS <PKG-ADIOS>`
     - dump output via ADIOS
     - :doc:`dump adios <dump_adios>`
     - PACKAGES/adios
     - ext
   * - :ref:`AMOEBA <PKG-AMOEBA>`
     - AMOEBA and HIPPO force fields
     - :doc:`AMOEBA and HIPPO howto <Howto_amoeba>`
     - amoeba
     - no
   * - :ref:`ASPHERE <PKG-ASPHERE>`
     - aspherical particle models
     - :doc:`Howto spherical <Howto_spherical>`
     - ellipse
     - no
   * - :ref:`ATC <PKG-ATC>`
     - Atom-to-Continuum coupling
     - :doc:`fix atc <fix_atc>`
     - PACKAGES/atc
     - int
   * - :ref:`AWPMD <PKG-AWPMD>`
     - wave packet MD
     - :doc:`pair_style awpmd/cut <pair_awpmd>`
     - PACKAGES/awpmd
     - int
   * - :ref:`BOCS <PKG-BOCS>`
     - BOCS bottom up coarse graining
     - :doc:`fix bocs <fix_bocs>`
     - PACKAGES/bocs
     - no
   * - :ref:`BODY <PKG-BODY>`
     - body-style particles
     - :doc:`Howto body <Howto_body>`
     - body
     - no
   * - :ref:`BPM <PKG-BPM>`
     - bonded particle models
     - :doc:`Howto bpm <Howto_bpm>`
     - bpm
     - no
   * - :ref:`BROWNIAN <PKG-BROWNIAN>`
     - Brownian dynamics, self-propelled particles
     - :doc:`fix brownian <fix_brownian>`, :doc:`fix propel/self <fix_propel_self>`
     - PACKAGES/brownian
     - no
   * - :ref:`CG-DNA <PKG-CG-DNA>`
     - coarse-grained DNA force fields
     - src/CG-DNA/README
     - PACKAGES/cgdna
     - no
   * - :ref:`CG-SPICA <PKG-CG-SPICA>`
     - SPICA (SDK) coarse-graining model
     - :doc:`pair_style lj/spica <pair_spica>`
     - PACKAGES/cgspica
     - no
   * - :ref:`CLASS2 <PKG-CLASS2>`
     - class 2 force fields
     - :doc:`pair_style lj/class2 <pair_class2>`
     - n/a
     - no
   * - :ref:`COLLOID <PKG-COLLOID>`
     - colloidal particles
     - :doc:`atom_style colloid <atom_style>`
     - colloid
     - no
   * - :ref:`COLVARS <PKG-COLVARS>`
     - `Colvars collective variables library <https://colvars.github.io/>`_
     - :doc:`fix colvars <fix_colvars>`
     - PACKAGES/colvars
     - int
   * - :ref:`COMPRESS <PKG-COMPRESS>`
     - I/O compression
     - :doc:`dump \*/gz <dump>`
     - n/a
     - sys
   * - :ref:`CORESHELL <PKG-CORESHELL>`
     - adiabatic core/shell model
     - :doc:`Howto coreshell <Howto_coreshell>`
     - coreshell
     - no
   * - :ref:`DIELECTRIC <PKG-DIELECTRIC>`
     - dielectric boundary solvers and force styles
     - :doc:`compute efield/atom <compute_efield_atom>`
     - PACKAGES/dielectric
     - no
   * - :ref:`DIFFRACTION <PKG-DIFFRACTION>`
     - virtual x-ray and electron diffraction
     - :doc:`compute xrd <compute_xrd>`
     - PACKAGES/diffraction
     - no
   * - :ref:`DIPOLE <PKG-DIPOLE>`
     - point dipole particles
     - :doc:`pair_style lj/.../dipole <pair_dipole>`
     - dipole
     - no
   * - :ref:`DPD-BASIC <PKG-DPD-BASIC>`
     - basic DPD models
     - :doc:`pair_styles dpd <pair_dpd>` :doc:`dpd/ext <pair_dpd_ext>`
     - PACKAGES/dpd-basic
     - no
   * - :ref:`DPD-MESO <PKG-DPD-MESO>`
     - mesoscale DPD models
     - :doc:`pair_style edpd <pair_mesodpd>`
     - PACKAGES/dpd-meso
     - no
   * - :ref:`DPD-REACT <PKG-DPD-REACT>`
     - reactive dissipative particle dynamics
     - src/DPD-REACT/README
     - PACKAGES/dpd-react
     - no
   * - :ref:`DPD-SMOOTH <PKG-DPD-SMOOTH>`
     - smoothed dissipative particle dynamics
     - src/DPD-SMOOTH/README
     - PACKAGES/dpd-smooth
     - no
   * - :ref:`DRUDE <PKG-DRUDE>`
     - Drude oscillators
     - :doc:`Howto drude <Howto_drude>`
     - PACKAGES/drude
     - no
   * - :ref:`EFF <PKG-EFF>`
     - electron force field
     - :doc:`pair_style eff/cut <pair_eff>`
     - PACKAGES/eff
     - no
   * - :ref:`ELECTRODE <PKG-ELECTRODE>`
     - electrode charges to match potential
     - :doc:`fix electrode/conp <fix_electrode>`
     - PACKAGES/electrode
     - no
   * - :ref:`EXTRA-COMPUTE <PKG-EXTRA-COMPUTE>`
     - additional compute styles
     - :doc:`compute <compute>`
     - n/a
     - no
   * - :ref:`EXTRA-DUMP <PKG-EXTRA-DUMP>`
     - additional dump styles
     - :doc:`dump <dump>`
     - n/a
     - no
   * - :ref:`EXTRA-FIX <PKG-EXTRA-FIX>`
     - additional fix styles
     - :doc:`fix <fix>`
     - n/a
     - no
   * - :ref:`EXTRA-MOLECULE <PKG-EXTRA-MOLECULE>`
     - additional molecular styles
     - :doc:`molecular styles <Commands_bond>`
     - n/a
     - no
   * - :ref:`EXTRA-PAIR <PKG-EXTRA-PAIR>`
     - additional pair styles
     - :doc:`pair_style <pair_style>`
     - n/a
     - no
   * - :ref:`FEP <PKG-FEP>`
     - free energy perturbation
     - :doc:`compute fep <compute_fep>`
     - PACKAGES/fep
     - no
   * - :ref:`GPU <PKG-GPU>`
     - GPU-enabled styles
     - :doc:`Section gpu <Speed_gpu>`
     - `Benchmarks <https://www.lammps.org/bench.html>`_
     - int
   * - :ref:`GRANULAR <PKG-GRANULAR>`
     - granular systems
     - :doc:`Howto granular <Howto_granular>`
     - pour
     - no
   * - :ref:`H5MD <PKG-H5MD>`
     - dump output via HDF5
     - :doc:`dump h5md <dump_h5md>`
     - n/a
     - ext
   * - :ref:`INTEL <PKG-INTEL>`
     - optimized Intel CPU and KNL styles
     - :doc:`Speed intel <Speed_intel>`
     - `Benchmarks <https://www.lammps.org/bench.html>`_
     - no
   * - :ref:`INTERLAYER <PKG-INTERLAYER>`
     - Inter-layer pair potentials
     - :doc:`several pair styles <Commands_pair>`
     - PACKAGES/interlayer
     - no
   * - :ref:`KIM <PKG-KIM>`
     - OpenKIM wrapper
     - :doc:`pair_style kim <pair_kim>`
     - kim
     - ext
   * - :ref:`KOKKOS <PKG-KOKKOS>`
     - Kokkos-enabled styles
     - :doc:`Speed kokkos <Speed_kokkos>`
     - `Benchmarks <https://www.lammps.org/bench.html>`_
     - no
   * - :ref:`KSPACE <PKG-KSPACE>`
     - long-range Coulombic solvers
     - :doc:`kspace_style <kspace_style>`
     - peptide
     - no
   * - :ref:`LATBOLTZ <PKG-LATBOLTZ>`
     - Lattice Boltzmann fluid
     - :doc:`fix lb/fluid <fix_lb_fluid>`
     - PACKAGES/latboltz
     - no
   * - :ref:`LEPTON <PKG-LEPTON>`
     - evaluate strings as potential function
     - :doc:`pair_style lepton <pair_lepton>`
     - PACKAGES/lepton
     - int
   * - :ref:`MACHDYN <PKG-MACHDYN>`
     - smoothed Mach dynamics
     - `SMD User Guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_
     - PACKAGES/machdyn
     - ext
   * - :ref:`MANIFOLD <PKG-MANIFOLD>`
     - motion on 2d surfaces
     - :doc:`fix manifoldforce <fix_manifoldforce>`
     - PACKAGES/manifold
     - no
   * - :ref:`MANYBODY <PKG-MANYBODY>`
     - many-body potentials
     - :doc:`pair_style tersoff <pair_tersoff>`
     - shear
     - no
   * - :ref:`MC <PKG-MC>`
     - Monte Carlo options
     - :doc:`fix gcmc <fix_gcmc>`
     - n/a
     - no
   * - :ref:`MDI <PKG-MDI>`
     - client-server code coupling
     - :doc:`MDI Howto <Howto_mdi>`
     - PACKAGES/mdi
     - ext
   * - :ref:`MEAM <PKG-MEAM>`
     - modified EAM potential (C++)
     - :doc:`pair_style meam <pair_meam>`
     - meam
     - no
   * - :ref:`MESONT <PKG-MESONT>`
     - mesoscopic tubular potential model
     - pair styles :doc:`mesocnt <pair_mesocnt>`
     - PACKAGES/mesont
     - no
   * - :ref:`MGPT <PKG-MGPT>`
     - fast MGPT multi-ion potentials
     - :doc:`pair_style mgpt <pair_mgpt>`
     - PACKAGES/mgpt
     - no
   * - :ref:`MISC <PKG-MISC>`
     - miscellaneous single-file commands
     - n/a
     - no
     - no
   * - :ref:`ML-HDNNP <PKG-ML-HDNNP>`
     - High-dimensional neural network potentials
     - :doc:`pair_style hdnnp <pair_hdnnp>`
     - PACKAGES/hdnnp
     - ext
   * - :ref:`ML-IAP <PKG-ML-IAP>`
     - multiple machine learning potentials
     - :doc:`pair_style mliap <pair_mliap>`
     - mliap
     - no
   * - :ref:`ML-PACE <PKG-ML-PACE>`
     - Atomic Cluster Expansion potential
     - :doc:`pair pace <pair_pace>`
     - PACKAGES/pace
     - ext
   * - :ref:`ML-POD <PKG-ML-POD>`
     - Proper orthogonal decomposition potentials
     - :doc:`pair pod <pair_pod>`
     - pod
     - ext
   * - :ref:`ML-QUIP <PKG-ML-QUIP>`
     - QUIP/libatoms interface
     - :doc:`pair_style quip <pair_quip>`
     - PACKAGES/quip
     - ext
   * - :ref:`ML-RANN <PKG-ML-RANN>`
     - Pair style for RANN potentials
     - :doc:`pair rann <pair_rann>`
     - PACKAGES/rann
     - no
   * - :ref:`ML-SNAP <PKG-ML-SNAP>`
     - quantum-fitted potential
     - :doc:`pair_style snap <pair_snap>`
     - snap
     - no
   * - :ref:`MOFFF <PKG-MOFFF>`
     - styles for `MOF-FF <MOFplus_>`_ force field
     - :doc:`pair_style buck6d/coul/gauss <pair_buck6d_coul_gauss>`
     - PACKAGES/mofff
     - no
   * - :ref:`MOLECULE <PKG-MOLECULE>`
     - molecular system force fields
     - :doc:`Howto bioFF <Howto_bioFF>`
     - peptide
     - no
   * - :ref:`MOLFILE <PKG-MOLFILE>`
     - `VMD <VMD_>`_ molfile plug-ins
     - :doc:`dump molfile <dump_molfile>`
     - n/a
     - ext
   * - :ref:`NETCDF <PKG-NETCDF>`
     - dump output via NetCDF
     - :doc:`dump netcdf <dump_netcdf>`
     - n/a
     - ext
   * - :ref:`OPENMP <PKG-OPENMP>`
     - OpenMP-enabled styles
     - :doc:`Speed omp <Speed_omp>`
     - `Benchmarks <https://www.lammps.org/bench.html>`_
     - no
   * - :ref:`OPT <PKG-OPT>`
     - optimized pair styles
     - :doc:`Speed opt <Speed_opt>`
     - `Benchmarks <https://www.lammps.org/bench.html>`_
     - no
   * - :ref:`ORIENT <PKG-ORIENT>`
     - fixes for orientation depended forces
     - :doc:`fix orient/* <fix_orient>`
     - PACKAGES/orient_eco
     - no
   * - :ref:`PERI <PKG-PERI>`
     - Peridynamics models
     - :doc:`pair_style peri <pair_peri>`
     - peri
     - no
   * - :ref:`PHONON <PKG-PHONON>`
     - phonon dynamical matrix
     - :doc:`fix phonon <fix_phonon>`
     - PACKAGES/phonon
     - no
   * - :ref:`PLUGIN <PKG-PLUGIN>`
     - Plugin loader command
     - :doc:`plugin <plugin>`
     - plugins
     - no
   * - :ref:`PLUMED <PKG-PLUMED>`
     - `PLUMED free energy library <https://www.plumed.org>`_
     - :doc:`fix plumed <fix_plumed>`
     - PACKAGES/plumed
     - ext
   * - :ref:`POEMS <PKG-POEMS>`
     - coupled rigid body motion
     - :doc:`fix poems <fix_poems>`
     - rigid
     - int
   * - :ref:`PTM <PKG-PTM>`
     - Polyhedral Template Matching
     - :doc:`compute ptm/atom <compute_ptm_atom>`
     - n/a
     - no
   * - :ref:`PYTHON <PKG-PYTHON>`
     - embed Python code in an input script
     - :doc:`python <python>`
     - python
     - sys
   * - :ref:`QEQ <PKG-QEQ>`
     - QEq charge equilibration
     - :doc:`fix qeq <fix_qeq>`
     - qeq
     - no
   * - :ref:`QMMM <PKG-QMMM>`
     - QM/MM coupling
     - :doc:`fix qmmm <fix_qmmm>`
     - PACKAGES/qmmm
     - ext
   * - :ref:`QTB <PKG-QTB>`
     - quantum nuclear effects
     - :doc:`fix qtb <fix_qtb>` :doc:`fix qbmsst <fix_qbmsst>`
     - qtb
     - no
   * - :ref:`REACTION <PKG-REACTION>`
     - chemical reactions in classical MD
     - :doc:`fix bond/react <fix_bond_react>`
     - PACKAGES/reaction
     - no
   * - :ref:`REAXFF <PKG-REAXFF>`
     - ReaxFF potential (C/C++)
     - :doc:`pair_style reaxff <pair_reaxff>`
     - reax
     - no
   * - :ref:`REPLICA <PKG-REPLICA>`
     - multi-replica methods
     - :doc:`Howto replica <Howto_replica>`
     - tad
     - no
   * - :ref:`RIGID <PKG-RIGID>`
     - rigid bodies and constraints
     - :doc:`fix rigid <fix_rigid>`
     - rigid
     - no
   * - :ref:`SCAFACOS <PKG-SCAFACOS>`
     - wrapper for ScaFaCoS Kspace solver
     - :doc:`kspace_style scafacos <kspace_style>`
     - PACKAGES/scafacos
     - ext
   * - :ref:`SHOCK <PKG-SHOCK>`
     - shock loading methods
     - :doc:`fix msst <fix_msst>`
     - n/a
     - no
   * - :ref:`SMTBQ <PKG-SMTBQ>`
     - second moment tight binding potentials
     - pair styles :doc:`smtbq <pair_smtbq>`, :doc:`smatb <pair_smatb>`
     - PACKAGES/smtbq
     - no
   * - :ref:`SPH <PKG-SPH>`
     - smoothed particle hydrodynamics
     - `SPH User Guide <PDF/SPH_LAMMPS_userguide.pdf>`_
     - PACKAGES/sph
     - no
   * - :ref:`SPIN <PKG-SPIN>`
     - magnetic atomic spin dynamics
     - :doc:`Howto spins <Howto_spins>`
     - SPIN
     - no
   * - :ref:`SRD <PKG-SRD>`
     - stochastic rotation dynamics
     - :doc:`fix srd <fix_srd>`
     - srd
     - no
   * - :ref:`TALLY <PKG-TALLY>`
     - pairwise tally computes
     - :doc:`compute XXX/tally <compute_tally>`
     - PACKAGES/tally
     - no
   * - :ref:`UEF <PKG-UEF>`
     - extensional flow
     - :doc:`fix nvt/uef <fix_nh_uef>`
     - PACKAGES/uef
     - no
   * - :ref:`VORONOI <PKG-VORONOI>`
     - Voronoi tesselation
     - :doc:`compute voronoi/atom <compute_voronoi_atom>`
     - n/a
     - ext
   * - :ref:`VTK <PKG-VTK>`
     - dump output via VTK
     - :doc:`compute vtk <dump_vtk>`
     - n/a
     - ext
   * - :ref:`YAFF <PKG-YAFF>`
     - additional styles implemented in YAFF
     - :doc:`angle_style cross <angle_cross>`
     - PACKAGES/yaff
     - no

.. _MOFplus: https://www.mofplus.org/content/show/MOF-FF
.. _PLUMED: https://www.plumed.org
.. _VMD: https://www.ks.uiuc.edu/Research/vmd/
