LAMMPS features
===============

LAMMPS is a classical molecular dynamics (MD) code with these general
classes of functionality:

* :ref:`General features <general>`
* :ref:`Particle and model types <particle>`
* :ref:`Interatomic potentials (force fields) <ff>`
* :ref:`Atom creation <create>`
* :ref:`Ensembles, constraints, and boundary conditions <ensemble>`
* :ref:`Integrators <integrate>`
* :ref:`Diagnostics <diag>`
* :ref:`Output <output>`
* :ref:`Multi-replica models <replica1>`
* :ref:`Pre- and post-processing <prepost>`
* :ref:`Specialized features (beyond MD itself) <special>`


----------


.. _general:

General features
------------------------------

* runs on a single processor or in parallel
* distributed-memory message-passing parallelism (MPI)
* spatial-decomposition of simulation domain for parallelism
* open-source distribution
* highly portable C++
* optional libraries used: MPI and single-processor FFT
* GPU (CUDA and OpenCL), Intel Xeon Phi, and OpenMP support for many code features
* easy to extend with new features and functionality
* runs from an input script
* syntax for defining and using variables and formulas
* syntax for looping over runs and breaking out of loops
* run one or multiple simulations simultaneously (in parallel) from one script
* build as library, invoke LAMMPS through library interface or provided Python wrapper
* couple with other codes: LAMMPS calls other code, other code calls LAMMPS, umbrella code calls both

.. _particle:

Particle and model types
---------------------------------------

(:doc:`atom style <atom_style>` command)

* atoms
* coarse-grained particles (e.g. bead-spring polymers)
* united-atom polymers or organic molecules
* all-atom polymers, organic molecules, proteins, DNA
* metals
* granular materials
* coarse-grained mesoscale models
* finite-size spherical and ellipsoidal particles
* finite-size  line segment (2d) and triangle (3d) particles
* point dipole particles
* rigid collections of particles
* hybrid combinations of these

.. _ff:

Interatomic potentials (force fields)
----------------------------------------------

(:doc:`pair style <pair_style>`, :doc:`bond style <bond_style>`,
:doc:`angle style <angle_style>`, :doc:`dihedral style <dihedral_style>`,
:doc:`improper style <improper_style>`, :doc:`kspace style <kspace_style>`
commands)

* pairwise potentials: Lennard-Jones, Buckingham, Morse, Born-Mayer-Huggins,     Yukawa, soft, class 2 (COMPASS), hydrogen bond, tabulated
* charged pairwise potentials: Coulombic, point-dipole
* many-body potentials: EAM, Finnis/Sinclair EAM, modified EAM (MEAM),     embedded ion method (EIM), EDIP, ADP, Stillinger-Weber, Tersoff,     REBO, AIREBO, ReaxFF, COMB, SNAP, Streitz-Mintmire, 3-body polymorphic
* long-range interactions for charge, point-dipoles, and LJ dispersion:     Ewald, Wolf, PPPM (similar to particle-mesh Ewald)
* polarization models: :doc:`QEq <fix_qeq>`,     :doc:`core/shell model <Howto_coreshell>`,     :doc:`Drude dipole model <Howto_drude>`
* charge equilibration (QEq via dynamic, point, shielded, Slater methods)
* coarse-grained potentials: DPD, GayBerne, REsquared, colloidal, DLVO
* mesoscopic potentials: granular, Peridynamics, SPH
* electron force field (eFF, AWPMD)
* bond potentials: harmonic, FENE, Morse, nonlinear, class 2,     quartic (breakable)
* angle potentials: harmonic, CHARMM, cosine, cosine/squared, cosine/periodic,     class 2 (COMPASS)
* dihedral potentials: harmonic, CHARMM, multi-harmonic, helix,     class 2 (COMPASS), OPLS
* improper potentials: harmonic, cvff, umbrella, class 2 (COMPASS)
* polymer potentials: all-atom, united-atom, bead-spring, breakable
* water potentials: TIP3P, TIP4P, SPC
* implicit solvent potentials: hydrodynamic lubrication, Debye
* force-field compatibility with common CHARMM, AMBER, DREIDING,     OPLS, GROMACS, COMPASS options
* access to the `OpenKIM Repository <http://openkim.org>`_ of potentials via     :doc:`kim\_init, kim\_interactions, and kim\_query <kim_commands>` commands
* hybrid potentials: multiple pair, bond, angle, dihedral, improper     potentials can be used in one simulation
* overlaid potentials: superposition of multiple pair potentials

.. _create:

Atom creation
--------------------------

(:doc:`read\_data <read_data>`, :doc:`lattice <lattice>`,
:doc:`create\_atoms <create_atoms>`, :doc:`delete\_atoms <delete_atoms>`,
:doc:`displace\_atoms <displace_atoms>`, :doc:`replicate <replicate>` commands)

* read in atom coords from files
* create atoms on one or more lattices (e.g. grain boundaries)
* delete geometric or logical groups of atoms (e.g. voids)
* replicate existing atoms multiple times
* displace atoms

.. _ensemble:

Ensembles, constraints, and boundary conditions
--------------------------------------------------------------

(:doc:`fix <fix>` command)

* 2d or 3d systems
* orthogonal or non-orthogonal (triclinic symmetry) simulation domains
* constant NVE, NVT, NPT, NPH, Parrinello/Rahman integrators
* thermostatting options for groups and geometric regions of atoms
* pressure control via Nose/Hoover or Berendsen barostatting in 1 to 3 dimensions
* simulation box deformation (tensile and shear)
* harmonic (umbrella) constraint forces
* rigid body constraints
* SHAKE bond and angle constraints
* Monte Carlo bond breaking, formation, swapping
* atom/molecule insertion and deletion
* walls of various kinds
* non-equilibrium molecular dynamics (NEMD)
* variety of additional boundary conditions and constraints

.. _integrate:

Integrators
---------------------------

(:doc:`run <run>`, :doc:`run\_style <run_style>`, :doc:`minimize <minimize>` commands)

* velocity-Verlet integrator
* Brownian dynamics
* rigid body integration
* energy minimization via conjugate gradient or steepest descent relaxation
* rRESPA hierarchical timestepping
* rerun command for post-processing of dump files

.. _diag:

Diagnostics
----------------------

* see various flavors of the :doc:`fix <fix>` and :doc:`compute <compute>` commands

.. _output:

Output
-------------------

(:doc:`dump <dump>`, :doc:`restart <restart>` commands)

* log file of thermodynamic info
* text dump files of atom coords, velocities, other per-atom quantities
* binary restart files
* parallel I/O of dump and restart files
* per-atom quantities (energy, stress, centro-symmetry parameter, CNA, etc)
* user-defined system-wide (log file) or per-atom (dump file) calculations
* spatial and time averaging of per-atom quantities
* time averaging of system-wide quantities
* atom snapshots in native, XYZ, XTC, DCD, CFG formats

.. _replica1:

Multi-replica models
-----------------------------------

* :doc:`nudged elastic band <neb>`
* :doc:`parallel replica dynamics <prd>`
* :doc:`temperature accelerated dynamics <tad>`
* :doc:`parallel tempering <temper>`

.. _prepost:

Pre- and post-processing
--------------------------------------

* A handful of pre- and post-processing tools are packaged with LAMMPS,
  some of which can convert input and output files to/from formats used
  by other codes; see the :doc:`Toos <Tools>` doc page.
* Our group has also written and released a separate toolkit called
  `Pizza.py <pizza_>`_ which provides tools for doing setup, analysis,
  plotting, and visualization for LAMMPS simulations.  Pizza.py is
  written in `Python <python_>`_ and is available for download from `the Pizza.py WWW site <pizza_>`_.

.. _pizza: http://www.sandia.gov/~sjplimp/pizza.html



.. _python: http://www.python.org



.. _special:

Specialized features
----------------------------------

LAMMPS can be built with optional packages which implement a variety
of additional capabilities.  See the :doc:`Packages <Packages>` doc
page for details.

These are LAMMPS capabilities which you may not think of as typical
classical MD options:

* :doc:`static <balance>` and :doc:`dynamic load-balancing <fix_balance>`
* :doc:`generalized aspherical particles <Howto_body>`
* :doc:`stochastic rotation dynamics (SRD) <fix_srd>`
* :doc:`real-time visualization and interactive MD <fix_imd>`
* calculate :doc:`virtual diffraction patterns <compute_xrd>`
* :doc:`atom-to-continuum coupling <fix_atc>` with finite elements
* coupled rigid body integration via the :doc:`POEMS <fix_poems>` library
* :doc:`QM/MM coupling <fix_qmmm>`
* Monte Carlo via :doc:`GCMC <fix_gcmc>` and :doc:`tfMC <fix_tfmc>` and :doc:`atom swapping <fix_atom_swap>`
* :doc:`path-integral molecular dynamics (PIMD) <fix_ipi>` and :doc:`this as well <fix_pimd>`
* :doc:`Direct Simulation Monte Carlo <pair_dsmc>` for low-density fluids
* :doc:`Peridynamics mesoscale modeling <pair_peri>`
* :doc:`Lattice Boltzmann fluid <fix_lb_fluid>`
* :doc:`targeted <fix_tmd>` and :doc:`steered <fix_smd>` molecular dynamics
* :doc:`two-temperature electron model <fix_ttm>`


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
