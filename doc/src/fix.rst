.. index:: fix

fix command
===========

Syntax
""""""

.. parsed-literal::

   fix ID group-ID style args

* ID = user-assigned name for the fix
* group-ID = ID of the group of atoms to apply the fix to
* style = one of a long list of possible style names (see below)
* args = arguments used by a particular style

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve
   fix 3 all nvt temp 300.0 300.0 0.01
   fix mine top setforce 0.0 NULL 0.0

Description
"""""""""""

Set a fix that will be applied to a group of atoms.  In LAMMPS, a
"fix" is any operation that is applied to the system during
timestepping or minimization.  Examples include updating of atom
positions and velocities due to time integration, controlling
temperature, applying constraint forces to atoms, enforcing boundary
conditions, computing diagnostics, etc.  There are hundreds of fixes
defined in LAMMPS and new ones can be added; see the
:doc:`Modify <Modify>` page for details.

Fixes perform their operations at different stages of the timestep.
If 2 or more fixes operate at the same stage of the timestep, they are
invoked in the order they were specified in the input script.

The ID of a fix can only contain alphanumeric characters and
underscores.

Fixes can be deleted with the :doc:`unfix <unfix>` command.

.. note::

   The :doc:`unfix <unfix>` command is the only way to turn off a
   fix; simply specifying a new fix with a similar style will not turn
   off the first one.  This is especially important to realize for
   integration fixes.  For example, using a :doc:`fix nve <fix_nve>`
   command for a second run after using a :doc:`fix nvt <fix_nh>` command
   for the first run, will not cancel out the NVT time integration
   invoked by the "fix nvt" command.  Thus two time integrators would be
   in place!

If you specify a new fix with the same ID and style as an existing
fix, the old fix is deleted and the new one is created (presumably
with new settings).  This is the same as if an "unfix" command were
first performed on the old fix, except that the new fix is kept in the
same order relative to the existing fixes as the old one originally
was.  Note that this operation also wipes out any additional changes
made to the old fix via the :doc:`fix_modify <fix_modify>` command.

The :doc:`fix modify <fix_modify>` command allows settings for some
fixes to be reset.  See the page for individual fixes for details.

Some fixes store an internal "state" which is written to binary
restart files via the :doc:`restart <restart>` or
:doc:`write_restart <write_restart>` commands.  This allows the fix to
continue on with its calculations in a restarted simulation.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file.  See the doc pages
for individual fixes for info on which ones can be restarted.

----------

Some fixes calculate one of three styles of quantities: global,
per-atom, or local, which can be used by other commands or output as
described below.  A global quantity is one or more system-wide values,
e.g. the energy of a wall interacting with particles.  A per-atom
quantity is one or more values per atom, e.g. the displacement vector
for each atom since time 0.  Per-atom values are set to 0.0 for atoms
not in the specified fix group.  Local quantities are calculated by
each processor based on the atoms it owns, but there may be zero or
more per atoms.

Note that a single fix can produce either global or per-atom or local
quantities (or none at all), but not both global and per-atom.  It can
produce local quantities in tandem with global or per-atom quantities.
The fix page will explain.

Global, per-atom, and local quantities each come in three kinds: a
single scalar value, a vector of values, or a 2d array of values.  The
doc page for each fix describes the style and kind of values it
produces, e.g. a per-atom vector.  Some fixes produce more than one
kind of a single style, e.g. a global scalar and a global vector.

When a fix quantity is accessed, as in many of the output commands
discussed below, it can be referenced via the following bracket
notation, where ID is the ID of the fix:

+-------------+--------------------------------------------+
| f_ID       | entire scalar, vector, or array             |
+-------------+--------------------------------------------+
| f_ID[I]    | one element of vector, one column of array  |
+-------------+--------------------------------------------+
| f_ID[I][J] | one element of array                        |
+-------------+--------------------------------------------+

In other words, using one bracket reduces the dimension of the
quantity once (vector -> scalar, array -> vector).  Using two brackets
reduces the dimension twice (array -> scalar).  Thus a command that
uses scalar fix values as input can also process elements of a vector
or array.

Note that commands and :doc:`variables <variable>` which use fix
quantities typically do not allow for all kinds, e.g. a command may
require a vector of values, not a scalar.  This means there is no
ambiguity about referring to a fix quantity as f_ID even if it
produces, for example, both a scalar and vector.  The doc pages for
various commands explain the details.

----------

In LAMMPS, the values generated by a fix can be used in several ways:

* Global values can be output via the :doc:`thermo_style custom <thermo_style>` or :doc:`fix ave/time <fix_ave_time>` command.
  Or the values can be referenced in a :doc:`variable equal <variable>` or
  :doc:`variable atom <variable>` command.
* Per-atom values can be output via the :doc:`dump custom <dump>` command.
  Or they can be time-averaged via the :doc:`fix ave/atom <fix_ave_atom>`
  command or reduced by the :doc:`compute reduce <compute_reduce>`
  command.  Or the per-atom values can be referenced in an :doc:`atom-style variable <variable>`.
* Local values can be reduced by the :doc:`compute reduce <compute_reduce>` command, or histogrammed by the :doc:`fix ave/histo <fix_ave_histo>` command.

See the :doc:`Howto output <Howto_output>` page for a summary of
various LAMMPS output options, many of which involve fixes.

The results of fixes that calculate global quantities can be either
"intensive" or "extensive" values.  Intensive means the value is
independent of the number of atoms in the simulation,
e.g. temperature.  Extensive means the value scales with the number of
atoms in the simulation, e.g. total rotational kinetic energy.
:doc:`Thermodynamic output <thermo_style>` will normalize extensive
values by the number of atoms in the system, depending on the
"thermo_modify norm" setting.  It will not normalize intensive values.
If a fix value is accessed in another way, e.g. by a
:doc:`variable <variable>`, you may want to know whether it is an
intensive or extensive value.  See the page for individual fixes
for further info.

----------

Each fix style has its own page which describes its arguments and
what it does, as listed below.  Here is an alphabetic list of fix
styles available in LAMMPS.  They are also listed in more compact form
on the :doc:`Commands fix <Commands_fix>` doc page.

There are also additional accelerated fix styles included in the
LAMMPS distribution for faster performance on CPUs, GPUs, and KNLs.
The individual style names on the :doc:`Commands fix <Commands_fix>` doc
page are followed by one or more of (g,i,k,o,t) to indicate which
accelerated styles exist.

* :doc:`accelerate/cos <fix_accelerate_cos>` - apply cosine-shaped acceleration to atoms
* :doc:`adapt <fix_adapt>` - change a simulation parameter over time
* :doc:`adapt/fep <fix_adapt_fep>` - enhanced version of fix adapt
* :doc:`addforce <fix_addforce>` - add a force to each atom
* :doc:`addtorque <fix_addtorque>` - add a torque to a group of atoms
* :doc:`append/atoms <fix_append_atoms>` - append atoms to a running simulation
* :doc:`atc <fix_atc>` - initiates a coupled MD/FE simulation
* :doc:`atom/swap <fix_atom_swap>` - Monte Carlo atom type swapping
* :doc:`ave/atom <fix_ave_atom>` - compute per-atom time-averaged quantities
* :doc:`ave/chunk <fix_ave_chunk>` - compute per-chunk time-averaged quantities
* :doc:`ave/correlate <fix_ave_correlate>` - compute/output time correlations
* :doc:`ave/correlate/long <fix_ave_correlate_long>` -
* :doc:`ave/histo <fix_ave_histo>` - compute/output time-averaged histograms
* :doc:`ave/histo/weight <fix_ave_histo>` - weighted version of fix ave/histo
* :doc:`ave/time <fix_ave_time>` - compute/output global time-averaged quantities
* :doc:`aveforce <fix_aveforce>` - add an averaged force to each atom
* :doc:`balance <fix_balance>` - perform dynamic load-balancing
* :doc:`brownian <fix_brownian>` - overdamped translational brownian motion
* :doc:`brownian/asphere <fix_brownian>` - overdamped translational and rotational brownian motion for ellipsoids
* :doc:`brownian/sphere <fix_brownian>` - overdamped translational and rotational brownian motion for spheres
* :doc:`bocs <fix_bocs>` - NPT style time integration with pressure correction
* :doc:`bond/break <fix_bond_break>` - break bonds on the fly
* :doc:`bond/create <fix_bond_create>` - create bonds on the fly
* :doc:`bond/create/angle <fix_bond_create>` - create bonds on the fly with angle constraints
* :doc:`bond/react <fix_bond_react>` - apply topology changes to model reactions
* :doc:`bond/swap <fix_bond_swap>` - Monte Carlo bond swapping
* :doc:`box/relax <fix_box_relax>` - relax box size during energy minimization
* :doc:`charge/regulation <fix_charge_regulation>` - Monte Carlo sampling of charge regulation
* :doc:`client/md <fix_client_md>` - MD client for client/server simulations
* :doc:`cmap <fix_cmap>` - enables CMAP cross-terms of the CHARMM force field
* :doc:`colvars <fix_colvars>` - interface to the collective variables "Colvars" library
* :doc:`controller <fix_controller>` - apply control loop feedback mechanism
* :doc:`deform <fix_deform>` - change the simulation box size/shape
* :doc:`deposit <fix_deposit>` - add new atoms above a surface
* :doc:`dpd/energy <fix_dpd_energy>` - constant energy dissipative particle dynamics
* :doc:`drag <fix_drag>` - drag atoms towards a defined coordinate
* :doc:`drude <fix_drude>` - part of Drude oscillator polarization model
* :doc:`drude/transform/direct <fix_drude_transform>` -  part of Drude oscillator polarization model
* :doc:`drude/transform/inverse <fix_drude_transform>` -  part of Drude oscillator polarization model
* :doc:`dt/reset <fix_dt_reset>` - reset the timestep based on velocity, forces
* :doc:`edpd/source <fix_dpd_source>` - add heat source to eDPD simulations
* :doc:`efield <fix_efield>` - impose electric field on system
* :doc:`ehex <fix_ehex>` - enhanced heat exchange algorithm
* :doc:`electron/stopping <fix_electron_stopping>` - electronic stopping power as a friction force
* :doc:`electron/stopping/fit <fix_electron_stopping>` - electronic stopping power as a friction force
* :doc:`enforce2d <fix_enforce2d>` - zero out z-dimension velocity and force
* :doc:`eos/cv <fix_eos_cv>` -
* :doc:`eos/table <fix_eos_table>` -
* :doc:`eos/table/rx <fix_eos_table_rx>` -
* :doc:`evaporate <fix_evaporate>` - remove atoms from simulation periodically
* :doc:`external <fix_external>` - callback to an external driver program
* :doc:`ffl <fix_ffl>` - apply a Fast-Forward Langevin equation thermostat
* :doc:`filter/corotate <fix_filter_corotate>` - implement corotation filter to allow larger timesteps with r-RESPA
* :doc:`flow/gauss <fix_flow_gauss>` - Gaussian dynamics for constant mass flux
* :doc:`freeze <fix_freeze>` - freeze atoms in a granular simulation
* :doc:`gcmc <fix_gcmc>` - grand canonical insertions/deletions
* :doc:`gld <fix_gld>` - generalized Langevin dynamics integrator
* :doc:`gle <fix_gle>` - generalized Langevin equation thermostat
* :doc:`gravity <fix_gravity>` - add gravity to atoms in a granular simulation
* :doc:`grem <fix_grem>` - implements the generalized replica exchange method
* :doc:`halt <fix_halt>` - terminate a dynamics run or minimization
* :doc:`heat <fix_heat>` - add/subtract momentum-conserving heat
* :doc:`hyper/global <fix_hyper_global>` - global hyperdynamics
* :doc:`hyper/local <fix_hyper_local>` - local hyperdynamics
* :doc:`imd <fix_imd>` - implements the "Interactive MD" (IMD) protocol
* :doc:`indent <fix_indent>` - impose force due to an indenter
* :doc:`ipi <fix_ipi>` - enable LAMMPS to run as a client for i-PI path-integral simulations
* :doc:`langevin <fix_langevin>` - Langevin temperature control
* :doc:`langevin/drude <fix_langevin_drude>` - Langevin temperature control of Drude oscillators
* :doc:`langevin/eff <fix_langevin_eff>` - Langevin temperature control for the electron force field model
* :doc:`langevin/spin <fix_langevin_spin>` - Langevin temperature control for a spin or spin-lattice system
* :doc:`latte <fix_latte>` - wrapper on LATTE density-functional tight-binding code
* :doc:`lb/fluid <fix_lb_fluid>` -
* :doc:`lb/momentum <fix_lb_momentum>` -
* :doc:`lb/pc <fix_lb_pc>` -
* :doc:`lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>` -
* :doc:`lb/viscous <fix_lb_viscous>` -
* :doc:`lineforce <fix_lineforce>` - constrain atoms to move in a line
* :doc:`manifoldforce <fix_manifoldforce>` - restrain atoms to a manifold during minimization
* :doc:`mdi/engine <fix_mdi_engine>` - connect LAMMPS to external programs via the MolSSI Driver Interface (MDI)
* :doc:`meso/move <fix_meso_move>` - move mesoscopic SPH/SDPD particles in a prescribed fashion
* :doc:`momentum <fix_momentum>` - zero the linear and/or angular momentum of a group of atoms
* :doc:`momentum/chunk <fix_momentum>` - zero the linear and/or angular momentum of a chunk of atoms
* :doc:`move <fix_move>` - move atoms in a prescribed fashion
* :doc:`mscg <fix_mscg>` - apply MSCG method for force-matching to generate coarse grain models
* :doc:`msst <fix_msst>` - multi-scale shock technique (MSST) integration
* :doc:`mvv/dpd <fix_mvv_dpd>` - DPD using the modified velocity-Verlet integration algorithm
* :doc:`mvv/edpd <fix_mvv_dpd>` - constant energy DPD using the modified velocity-Verlet algorithm
* :doc:`mvv/tdpd <fix_mvv_dpd>` - constant temperature DPD using the modified velocity-Verlet algorithm
* :doc:`neb <fix_neb>` - nudged elastic band (NEB) spring forces
* :doc:`neb/spin <fix_neb_spin>` - nudged elastic band (NEB) spring forces for spins
* :doc:`nph <fix_nh>` - constant NPH time integration via Nose/Hoover
* :doc:`nph/asphere <fix_nph_asphere>` - NPH for aspherical particles
* :doc:`nph/body <fix_nph_body>` - NPH for body particles
* :doc:`nph/eff <fix_nh_eff>` - NPH for  nuclei and electrons in the electron force field model
* :doc:`nph/sphere <fix_nph_sphere>` - NPH for spherical particles
* :doc:`nphug <fix_nphug>` - constant-stress Hugoniostat integration
* :doc:`npt <fix_nh>` - constant NPT time integration via Nose/Hoover
* :doc:`npt/asphere <fix_npt_asphere>` - NPT for aspherical particles
* :doc:`npt/body <fix_npt_body>` - NPT for body particles
* :doc:`npt/cauchy <fix_npt_cauchy>` - NPT with Cauchy stress
* :doc:`npt/eff <fix_nh_eff>` - NPT for  nuclei and electrons in the electron force field model
* :doc:`npt/sphere <fix_npt_sphere>` - NPT for spherical particles
* :doc:`npt/uef <fix_nh_uef>` - NPT style time integration with diagonal flow
* :doc:`numdiff <fix_numdiff>` - compute derivatives of per-atom data from finite differences
* :doc:`nve <fix_nve>` - constant NVE time integration
* :doc:`nve/asphere <fix_nve_asphere>` - NVE for aspherical particles
* :doc:`nve/asphere/noforce <fix_nve_asphere_noforce>` - NVE for aspherical particles without forces
* :doc:`nve/awpmd <fix_nve_awpmd>` - NVE for the Antisymmetrized Wave Packet Molecular Dynamics model
* :doc:`nve/body <fix_nve_body>` - NVE for body particles
* :doc:`nve/dot <fix_nve_dot>` - rigid body constant energy time integrator for coarse grain models
* :doc:`nve/dotc/langevin <fix_nve_dotc_langevin>` - Langevin style rigid body time integrator for coarse grain models
* :doc:`nve/eff <fix_nve_eff>` - NVE for  nuclei and electrons in the electron force field model
* :doc:`nve/limit <fix_nve_limit>` - NVE with limited step length
* :doc:`nve/line <fix_nve_line>` - NVE for line segments
* :doc:`nve/manifold/rattle <fix_nve_manifold_rattle>` -
* :doc:`nve/noforce <fix_nve_noforce>` - NVE without forces (v only)
* :doc:`nve/sphere <fix_nve_sphere>` - NVE for spherical particles
* :doc:`nve/spin <fix_nve_spin>` - NVE for a spin or spin-lattice system
* :doc:`nve/tri <fix_nve_tri>` - NVE for triangles
* :doc:`nvk <fix_nvk>` - constant kinetic energy time integration
* :doc:`nvt <fix_nh>` - NVT time integration via Nose/Hoover
* :doc:`nvt/asphere <fix_nvt_asphere>` - NVT for aspherical particles
* :doc:`nvt/body <fix_nvt_body>` - NVT for body particles
* :doc:`nvt/eff <fix_nh_eff>` - NVE for  nuclei and electrons in the electron force field model
* :doc:`nvt/manifold/rattle <fix_nvt_manifold_rattle>` -
* :doc:`nvt/sllod <fix_nvt_sllod>` - NVT for NEMD with SLLOD equations
* :doc:`nvt/sllod/eff <fix_nvt_sllod_eff>` - NVT for NEMD with SLLOD equations for the electron force field model
* :doc:`nvt/sphere <fix_nvt_sphere>` - NVT for spherical particles
* :doc:`nvt/uef <fix_nh_uef>` - NVT style time integration with diagonal flow
* :doc:`oneway <fix_oneway>` - constrain particles on move in one direction
* :doc:`orient/bcc <fix_orient>` - add grain boundary migration force for BCC
* :doc:`orient/fcc <fix_orient>` - add grain boundary migration force for FCC
* :doc:`orient/eco <fix_orient_eco>` - add generalized grain boundary migration force
* :doc:`pafi <fix_pafi>` - constrained force averages on hyper-planes to compute free energies (PAFI)
* :doc:`pair/tracker <fix_pair_tracker>` - track properties of pairwise interactions
* :doc:`phonon <fix_phonon>` - calculate dynamical matrix from MD simulations
* :doc:`pimd <fix_pimd>` - Feynman path integral molecular dynamics
* :doc:`planeforce <fix_planeforce>` - constrain atoms to move in a plane
* :doc:`plumed <fix_plumed>` - wrapper on PLUMED free energy library
* :doc:`poems <fix_poems>` - constrain clusters of atoms to move as coupled rigid bodies
* :doc:`polarize/bem/gmres <fix_polarize>` -
* :doc:`polarize/bem/icc <fix_polarize>` -
* :doc:`polarize/functional <fix_polarize>` -
* :doc:`pour <fix_pour>` - pour new atoms/molecules into a granular simulation domain
* :doc:`precession/spin <fix_precession_spin>` -
* :doc:`press/berendsen <fix_press_berendsen>` - pressure control by Berendsen barostat
* :doc:`print <fix_print>` - print text and variables during a simulation
* :doc:`propel/self <fix_propel_self>` - model self-propelled particles
* :doc:`property/atom <fix_property_atom>` - add customized per-atom values
* :doc:`python/invoke <fix_python_invoke>` - call a Python function during a simulation
* :doc:`python/move <fix_python_move>` -  call a Python function during a simulation run
* :doc:`qbmsst <fix_qbmsst>` - quantum bath multi-scale shock technique time integrator
* :doc:`qeq/comb <fix_qeq_comb>` - charge equilibration for COMB potential
* :doc:`qeq/dynamic <fix_qeq>` - charge equilibration via dynamic method
* :doc:`qeq/fire <fix_qeq>` - charge equilibration via FIRE minimizer
* :doc:`qeq/point <fix_qeq>` - charge equilibration via point method
* :doc:`qeq/reaxff <fix_qeq_reaxff>` - charge equilibration for ReaxFF potential
* :doc:`qeq/shielded <fix_qeq>` - charge equilibration via shielded method
* :doc:`qeq/slater <fix_qeq>` - charge equilibration via Slater method
* :doc:`qmmm <fix_qmmm>` - functionality to enable a quantum mechanics/molecular mechanics coupling
* :doc:`qtb <fix_qtb>` - implement quantum thermal bath scheme
* :doc:`rattle <fix_shake>` - RATTLE constraints on bonds and/or angles
* :doc:`reaxff/bonds <fix_reaxff_bonds>` - write out ReaxFF bond information
* :doc:`reaxff/species <fix_reaxff_species>` - write out ReaxFF molecule information
* :doc:`recenter <fix_recenter>` - constrain the center-of-mass position of a group of atoms
* :doc:`restrain <fix_restrain>` - constrain a bond, angle, dihedral
* :doc:`rhok <fix_rhok>` - add bias potential for long-range ordered systems
* :doc:`rigid <fix_rigid>` - constrain one or more clusters of atoms to move as a rigid body with NVE integration
* :doc:`rigid/meso <fix_rigid_meso>` - constrain clusters of mesoscopic SPH/SDPD particles to move as a rigid body
* :doc:`rigid/nph <fix_rigid>` - constrain one or more clusters of atoms to move as a rigid body with NPH integration
* :doc:`rigid/nph/small <fix_rigid>` - constrain many small clusters of atoms to move as a rigid body with NPH integration
* :doc:`rigid/npt <fix_rigid>` - constrain one or more clusters of atoms to move as a rigid body with NPT integration
* :doc:`rigid/npt/small <fix_rigid>` - constrain many small clusters of atoms to move as a rigid body with NPT integration
* :doc:`rigid/nve <fix_rigid>` - constrain one or more clusters of atoms to move as a rigid body with alternate NVE integration
* :doc:`rigid/nve/small <fix_rigid>` - constrain many small clusters of atoms to move as a rigid body with alternate NVE integration
* :doc:`rigid/nvt <fix_rigid>` - constrain one or more clusters of atoms to move as a rigid body with NVT integration
* :doc:`rigid/nvt/small <fix_rigid>` - constrain many small clusters of atoms to move as a rigid body with NVT integration
* :doc:`rigid/small <fix_rigid>` - constrain many small clusters of atoms to move as a rigid body with NVE integration
* :doc:`rx <fix_rx>` -
* :doc:`saed/vtk <fix_saed_vtk>` -
* :doc:`setforce <fix_setforce>` - set the force on each atom
* :doc:`setforce/spin <fix_setforce>` - set magnetic precession vectors on each atom
* :doc:`shake <fix_shake>` - SHAKE constraints on bonds and/or angles
* :doc:`shardlow <fix_shardlow>` - integration of DPD equations of motion using the Shardlow splitting
* :doc:`smd <fix_smd>` - applied a steered MD force to a group
* :doc:`smd/adjust_dt <fix_smd_adjust_dt>` -
* :doc:`smd/integrate_tlsph <fix_smd_integrate_tlsph>` -
* :doc:`smd/integrate_ulsph <fix_smd_integrate_ulsph>` -
* :doc:`smd/move_tri_surf <fix_smd_move_triangulated_surface>` -
* :doc:`smd/setvel <fix_smd_setvel>` -
* :doc:`smd/wall_surface <fix_smd_wall_surface>` -
* :doc:`sph <fix_sph>` - time integration for SPH/DPDE particles
* :doc:`sph/stationary <fix_sph_stationary>` -
* :doc:`spring <fix_spring>` - apply harmonic spring force to group of atoms
* :doc:`spring/chunk <fix_spring_chunk>` - apply harmonic spring force to each chunk of atoms
* :doc:`spring/rg <fix_spring_rg>` - spring on radius of gyration of group of atoms
* :doc:`spring/self <fix_spring_self>` - spring from each atom to its origin
* :doc:`srd <fix_srd>` - stochastic rotation dynamics (SRD)
* :doc:`store/force <fix_store_force>` - store force on each atom
* :doc:`store/state <fix_store_state>` - store attributes for each atom
* :doc:`tdpd/source <fix_dpd_source>` -
* :doc:`temp/berendsen <fix_temp_berendsen>` - temperature control by Berendsen thermostat
* :doc:`temp/csld <fix_temp_csvr>` - canonical sampling thermostat with Langevin dynamics
* :doc:`temp/csvr <fix_temp_csvr>` - canonical sampling thermostat with Hamiltonian dynamics
* :doc:`temp/rescale <fix_temp_rescale>` - temperature control by velocity rescaling
* :doc:`temp/rescale/eff <fix_temp_rescale_eff>` - temperature control by velocity rescaling in the electron force field model
* :doc:`tfmc <fix_tfmc>` - perform force-bias Monte Carlo with time-stamped method
* :doc:`tgnvt/drude <fix_tgnh_drude>` - NVT time integration for Drude polarizable model via temperature-grouped Nose-Hoover
* :doc:`tgnpt/drude <fix_tgnh_drude>` - NPT time integration for Drude polarizable model via temperature-grouped Nose-Hoover
* :doc:`thermal/conductivity <fix_thermal_conductivity>` - Muller-Plathe kinetic energy exchange for thermal conductivity calculation
* :doc:`ti/spring <fix_ti_spring>` -
* :doc:`tmd <fix_tmd>` - guide a group of atoms to a new configuration
* :doc:`ttm <fix_ttm>` - two-temperature model for electronic/atomic coupling
* :doc:`ttm/mod <fix_ttm>` - enhanced two-temperature model with additional options
* :doc:`tune/kspace <fix_tune_kspace>` - auto-tune KSpace parameters
* :doc:`vector <fix_vector>` - accumulate a global vector every N timesteps
* :doc:`viscosity <fix_viscosity>` - Muller-Plathe momentum exchange for viscosity calculation
* :doc:`viscous <fix_viscous>` - viscous damping for granular simulations
* :doc:`wall/body/polygon <fix_wall_body_polygon>` -
* :doc:`wall/body/polyhedron <fix_wall_body_polyhedron>` -
* :doc:`wall/colloid <fix_wall>` - Lennard-Jones wall interacting with finite-size particles
* :doc:`wall/ees <fix_wall_ees>` - wall for ellipsoidal particles
* :doc:`wall/gran <fix_wall_gran>` - frictional wall(s) for granular simulations
* :doc:`wall/gran/region <fix_wall_gran_region>` -
* :doc:`wall/harmonic <fix_wall>` - harmonic spring wall
* :doc:`wall/lj1043 <fix_wall>` - Lennard-Jones 10-4-3 wall
* :doc:`wall/lj126 <fix_wall>` - Lennard-Jones 12-6 wall
* :doc:`wall/lj93 <fix_wall>` - Lennard-Jones 9-3 wall
* :doc:`wall/morse <fix_wall>` - Morse potential wall
* :doc:`wall/piston <fix_wall_piston>` - moving reflective piston wall
* :doc:`wall/reflect <fix_wall_reflect>` - reflecting wall(s)
* :doc:`wall/reflect/stochastic <fix_wall_reflect_stochastic>` - reflecting wall(s) with finite temperature
* :doc:`wall/region <fix_wall_region>` - use region surface as wall
* :doc:`wall/region/ees <fix_wall_ees>` - use region surface as wall for ellipsoidal particles
* :doc:`wall/srd <fix_wall_srd>` - slip/no-slip wall for SRD particles
* :doc:`widom <fix_widom>` - Widom insertions of atoms or molecules

Restrictions
""""""""""""

Some fix styles are part of specific packages.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.  The doc pages for
individual fixes tell if it is part of a package.

Related commands
""""""""""""""""

:doc:`unfix <unfix>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

none
