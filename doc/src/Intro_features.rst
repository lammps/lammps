LAMMPS features
---------------

LAMMPS is a classical molecular dynamics (MD) code with these general
classes of functionality:

1. :ref:`General features <general>`
2. :ref:`Particle and model types <particle>`
3. :ref:`Interatomic potentials (force fields) <ff>`
4. :ref:`Atom creation <create>`
5. :ref:`Ensembles, constraints, and boundary conditions <ensemble>`
6. :ref:`Integrators <integrate>`
7. :ref:`Diagnostics <diag>`
8. :ref:`Output <output>`
9. :ref:`Multi-replica models <replica1>`
10. :ref:`Pre- and post-processing <prepost>`
11. :ref:`Specialized features (beyond MD itself) <special>`

----------

.. _general:

General features
^^^^^^^^^^^^^^^^

* runs on a single processor or in parallel
* distributed memory message-passing parallelism (MPI)
* shared memory multi-threading parallelism (OpenMP)
* spatial decomposition of simulation domain for MPI parallelism
* particle decomposition inside spatial decomposition for OpenMP and GPU parallelism
* GPLv2 licensed open-source distribution
* highly portable C++-11
* modular code with most functionality in optional packages
* only depends on MPI library for basic parallel functionality, MPI stub for serial compilation
* other libraries are optional and only required for specific packages
* GPU (CUDA, OpenCL, HIP, SYCL), Intel Xeon Phi, and OpenMP support for many code features
* easy to extend with new features and functionality
* runs from an input script
* syntax for defining and using variables and formulas
* syntax for looping over runs and breaking out of loops
* run one or multiple simulations simultaneously (in parallel) from one script
* build as library, invoke LAMMPS through library interface (from C, C++, Fortran) or provided Python wrapper or SWIG based wrappers
* couple with other codes: LAMMPS calls other code, other code calls LAMMPS, umbrella code calls both, MDI coupling interface
* call out to Python for computing forces, time integration, or other tasks
* plugin interface for loading external features at runtime
* large integrated collection of tests

.. _particle:

Particle and model types
^^^^^^^^^^^^^^^^^^^^^^^^

(See :doc:`atom style <atom_style>` command)

* atoms
* coarse-grained particles (e.g. bead-spring polymers)
* united-atom polymers or organic molecules
* all-atom polymers, organic molecules, proteins, DNA
* metals
* metal oxides
* granular materials
* coarse-grained mesoscale models
* finite-size spherical and ellipsoidal particles
* finite-size line segment (2d) and triangle (3d) particles
* finite-size rounded polygons (2d) and polyhedra (3d) particles
* point dipole particles
* particles with magnetic spin
* rigid collections of n particles
* hybrid combinations of these

.. _ff:

Interatomic potentials (force fields)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(See :doc:`pair style <pair_style>`, :doc:`bond style <bond_style>`,
:doc:`angle style <angle_style>`, :doc:`dihedral style <dihedral_style>`,
:doc:`improper style <improper_style>`, :doc:`kspace style <kspace_style>`
commands)

* pairwise potentials: Lennard-Jones, Buckingham, Morse, Born-Mayer-Huggins, Yukawa, soft, Class II (COMPASS), hydrogen bond, harmonic, gaussian, tabulated, scripted
* charged pairwise potentials: Coulombic, point-dipole
* many-body potentials: EAM, Finnis/Sinclair, MEAM, MEAM+SW, EIM, EDIP, ADP, Stillinger-Weber, Tersoff, REBO, AIREBO, ReaxFF, COMB, Streitz-Mintmire, 3-body polymorphic, BOP, Vashishta
* machine learning potentials: ACE, AGNI, GAP, Behler-Parrinello (N2P2), POD, RANN
* interfaces to ML potentials distributed by external groups: ANI, ChIMES, DeepPot, HIPNN, MTP
* long-range interactions for charge, point-dipoles, and LJ dispersion:  Ewald, Wolf, PPPM (similar to particle-mesh Ewald), MSM, ScaFaCoS
* polarization models: :doc:`QEq <fix_qeq>`, :doc:`core/shell model <Howto_coreshell>`, :doc:`Drude dipole model <Howto_drude>`
* charge equilibration (QEq via dynamic, point, shielded, Slater methods)
* coarse-grained potentials: DPD, GayBerne, REsquared, colloidal, DLVO, oxDNA / oxRNA, SPICA
* mesoscopic potentials: granular, Peridynamics, SPH, mesoscopic tubular potential (MESONT)
* semi-empirical potentials: multi-ion generalized pseudopotential theory (MGPT), second moment tight binding + QEq (SMTB-Q), density functional tight-binding (LATTE)
* electron force field (eFF, AWPMD)
* bond potentials: harmonic, FENE, Morse, nonlinear, Class II (COMPASS), quartic (breakable), tabulated, scripted
* angle potentials: harmonic, CHARMM, cosine, cosine/squared, cosine/periodic, Class II (COMPASS), tabulated, scripted
* dihedral potentials: harmonic, CHARMM, multi-harmonic, helix, Class II (COMPASS), OPLS, tabulated, scripted
* improper potentials: harmonic, cvff, umbrella, Class II (COMPASS), tabulated
* polymer potentials: all-atom, united-atom, bead-spring, breakable
* water potentials: TIP3P, TIP4P, SPC, SPC/E and variants
* interlayer potentials for graphene and analogues, hetero-junctions
* metal-organic framework potentials (QuickFF, MO-FF)
* implicit solvent potentials: hydrodynamic lubrication, Debye
* force-field compatibility with CHARMM, AMBER, DREIDING, OPLS, GROMACS, Class II (COMPASS), UFF, ClayFF, DREIDING, AMOEBA, INTERFACE
* access to the `OpenKIM Repository <https://openkim.org>`_ of potentials via the :doc:`kim command <kim_commands>`
* hybrid potentials: multiple pair, bond, angle, dihedral, improper potentials can be used in one simulation
* overlaid potentials: superposition of multiple pair potentials (including many-body) with optional scale factor

.. _create:

Atom creation
^^^^^^^^^^^^^

(See :doc:`read_data <read_data>`, :doc:`lattice <lattice>`,
:doc:`create_atoms <create_atoms>`, :doc:`delete_atoms <delete_atoms>`,
:doc:`displace_atoms <displace_atoms>`, :doc:`replicate <replicate>` commands)

* read in atom coordinates from files
* create atoms on one or more lattices (e.g. grain boundaries)
* delete geometric or logical groups of atoms (e.g. voids)
* replicate existing atoms multiple times
* displace atoms

.. _ensemble:

Ensembles, constraints, and boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(See :doc:`fix <fix>` command)

* 2d or 3d systems
* orthogonal or non-orthogonal (triclinic symmetry) simulation domains
* constant NVE, NVT, NPT, NPH, Parrinello/Rahman integrators
* thermostatting options for groups and geometric regions of atoms
* pressure control via Nose/Hoover or Berendsen barostatting in 1 to 3 dimensions
* simulation box deformation (tensile and shear)
* harmonic (umbrella) constraint forces
* rigid body constraints
* SHAKE / RATTLE bond and angle constraints
* motion constraints to manifold surfaces
* Monte Carlo bond breaking, formation, swapping, template based reaction modeling
* atom/molecule insertion and deletion
* walls of various kinds, static and moving
* non-equilibrium molecular dynamics (NEMD)
* variety of additional boundary conditions and constraints

.. _integrate:

Integrators
^^^^^^^^^^^

(See :doc:`run <run>`, :doc:`run_style <run_style>`, :doc:`minimize <minimize>` commands)

* velocity-Verlet integrator
* Brownian dynamics
* rigid body integration
* energy minimization via conjugate gradient, steepest descent relaxation, or damped dynamics (FIRE, Quickmin)
* rRESPA hierarchical timestepping
* fixed or adaptive time step
* rerun command for post-processing of dump files

.. _diag:

Diagnostics
^^^^^^^^^^^

* see various flavors of the :doc:`fix <fix>` and :doc:`compute <compute>` commands
* introspection command for system, simulation, and compile time settings and configurations

.. _output:

Output
^^^^^^

(:doc:`dump <dump>`, :doc:`restart <restart>` commands)

* log file of thermodynamic info
* text dump files of atom coordinates, velocities, other per-atom quantities
* dump output on fixed and variable intervals, based timestep or simulated time
* binary restart files
* parallel I/O of dump and restart files
* per-atom quantities (energy, stress, centro-symmetry parameter, CNA, etc.)
* user-defined system-wide (log file) or per-atom (dump file) calculations
* custom partitioning (chunks) for binning, and static or dynamic grouping of atoms for analysis
* spatial, time, and per-chunk averaging of per-atom quantities
* time averaging and histogramming of system-wide quantities
* atom snapshots in native, XYZ, XTC, DCD, CFG, NetCDF, HDF5, ADIOS2, YAML formats
* on-the-fly compression of output and decompression of read in files

.. _replica1:

Multi-replica models
^^^^^^^^^^^^^^^^^^^^

* :doc:`nudged elastic band <neb>`
* :doc:`hyperdynamics <hyper>`
* :doc:`parallel replica dynamics <prd>`
* :doc:`temperature accelerated dynamics <tad>`
* :doc:`parallel tempering <temper>`
* path-integral MD: :doc:`first variant <fix_pimd>`, :doc:`second variant <fix_ipi>`
* multi-walker collective variables with :doc:`Colvars <fix_colvars>` and :doc:`Plumed <fix_plumed>`

.. _prepost:

Pre- and post-processing
^^^^^^^^^^^^^^^^^^^^^^^^

* A handful of pre- and post-processing tools are packaged with LAMMPS,
  some of which can convert input and output files to/from formats used
  by other codes; see the :doc:`Tools <Tools>` page.
* Our group has also written and released a separate toolkit called
  `Pizza.py <pizza_>`_ which provides tools for doing setup, analysis,
  plotting, and visualization for LAMMPS simulations.  Pizza.py is
  written in `Python <python_>`_ and is available for download from `the Pizza.py WWW site <pizza_>`_.

.. _pizza: https://lammps.github.io/pizza

.. _python: https://www.python.org

.. _special:

Specialized features
^^^^^^^^^^^^^^^^^^^^

LAMMPS can be built with optional packages which implement a variety of
additional capabilities.  See the :doc:`Optional Packages <Packages>`
page for details.

These are LAMMPS capabilities which you may not think of as typical
classical MD options:

* :doc:`static <balance>` and :doc:`dynamic load-balancing <fix_balance>`, optional with recursive bisectioning decomposition
* :doc:`generalized aspherical particles <Howto_body>`
* :doc:`stochastic rotation dynamics (SRD) <fix_srd>`
* :doc:`real-time visualization and interactive MD <fix_imd>`, :doc:`built-in renderer for images and movies <dump_image>`
* calculate :doc:`virtual diffraction patterns <compute_xrd>`
* calculate :doc:`finite temperature phonon dispersion <fix_phonon>` and the :doc:`dynamical matrix of minimized structures <dynamical_matrix>`
* :doc:`atom-to-continuum coupling <fix_atc>` with finite elements
* coupled rigid body integration via the :doc:`POEMS <fix_poems>` library
* :doc:`QM/MM coupling <fix_qmmm>`
* Monte Carlo via :doc:`GCMC <fix_gcmc>` and :doc:`tfMC <fix_tfmc>` and :doc:`atom swapping <fix_atom_swap>`
* :doc:`path-integral molecular dynamics (PIMD) <fix_ipi>` and :doc:`this as well <fix_pimd>`
* :doc:`Direct Simulation Monte Carlo <pair_dsmc>` for low-density fluids
* :doc:`Peridynamics modeling <pair_peri>`
* :doc:`Lattice Boltzmann fluid <fix_lb_fluid>`
* :doc:`targeted <fix_tmd>` and :doc:`steered <fix_smd>` molecular dynamics
* :doc:`two-temperature electron model <fix_ttm>`
