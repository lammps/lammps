Fix styles
==========

In LAMMPS, a "fix" is any operation that is computed during timestepping
that alters some property of the system.  Essentially everything that
happens during a simulation besides force computation, neighbor list
construction, and output, is a "fix".  This includes time integration
(update of coordinates and velocities), force constraints or boundary
conditions (SHAKE or walls), and diagnostics (compute a diffusion
coefficient).  New styles can be created to add new options to LAMMPS.

The file ``src/fix_setforce.cpp`` is a simple example of setting forces
on atoms to prescribed values.  There are dozens of fix options already
in LAMMPS; choose one as a template that is similar to what you want to
implement.  There also is a detailed discussion of :doc:`how to write
new fix styles <Developer_write_fix>` in LAMMPS.

Here is a brief description of methods you can define in your new
derived class.  See ``src/fix.h`` for additional details.

+-------------------------------+-------------------------------------------------------------+
| Required                      | "pure" method that *must* be overridden in a derived class  |
+===============================+=============================================================+
| setmask                       | determines when the fix is called during the timestep       |
+-------------------------------+-------------------------------------------------------------+

+-------------------------------+----------------------------------------------------------------------------------+
| Optional                      | methods that have a default or empty implementation                              |
+===============================+==================================================================================+
| post_constructor              | perform tasks that cannot be run in the constructor                              |
+-------------------------------+----------------------------------------------------------------------------------+
| init                          | initialization before a run                                                      |
+-------------------------------+----------------------------------------------------------------------------------+
| init_list                     | store pointer to neighbor list; called by neighbor list code                     |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_pre_exchange            | called before atom exchange in setup                                             |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_pre_neighbor            | called before neighbor list build in setup                                       |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_post_neighbor           | called after neighbor list build in setup                                        |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_pre_force               | called before force computation in setup                                         |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_pre_reverse             | called before reverse communication in setup                                     |
+-------------------------------+----------------------------------------------------------------------------------+
| setup                         | called immediately before the first timestep and after forces are computed       |
+-------------------------------+----------------------------------------------------------------------------------+
| min_setup                     | like setup, but for minimizations instead of MD runs                             |
+-------------------------------+----------------------------------------------------------------------------------+
| initial_integrate             | called at very beginning of each timestep                                        |
+-------------------------------+----------------------------------------------------------------------------------+
| post_integrate                | called after the first half-step position/velocity update                        |
+-------------------------------+----------------------------------------------------------------------------------+
| pre_exchange                  | called before atom exchange on re-neighboring steps                              |
+-------------------------------+----------------------------------------------------------------------------------+
| pre_neighbor                  | called before neighbor list build                                                |
+-------------------------------+----------------------------------------------------------------------------------+
| post_neighbor                 | called after neighbor list build                                                 |
+-------------------------------+----------------------------------------------------------------------------------+
| pre_force                     | called before pair & molecular forces are computed                               |
+-------------------------------+----------------------------------------------------------------------------------+
| pre_reverse                   | called before reverse communication of forces                                    |
+-------------------------------+----------------------------------------------------------------------------------+
| post_force                    | called after pair & molecular forces are computed and communicated               |
+-------------------------------+----------------------------------------------------------------------------------+
| final_integrate               | called at end of each timestep                                                   |
+-------------------------------+----------------------------------------------------------------------------------+
| end_of_step                   | called at very end of timestep                                                   |
+-------------------------------+----------------------------------------------------------------------------------+
| post_run                      | called at the end of each run command                                            |
+-------------------------------+----------------------------------------------------------------------------------+
| write_restart                 | dumps fix info to restart file                                                   |
+-------------------------------+----------------------------------------------------------------------------------+
| restart                       | uses info from restart file to re-initialize the fix                             |
+-------------------------------+----------------------------------------------------------------------------------+
| grow_arrays                   | allocate memory for atom-based arrays used by fix                                |
+-------------------------------+----------------------------------------------------------------------------------+
| copy_arrays                   | copy atom info when an atom migrates to a new processor                          |
+-------------------------------+----------------------------------------------------------------------------------+
| set_arrays                    | set per-atom values when a new atom is created                                   |
+-------------------------------+----------------------------------------------------------------------------------+
| pack_exchange                 | store atom's data in a buffer                                                    |
+-------------------------------+----------------------------------------------------------------------------------+
| unpack_exchange               | retrieve atom's data from a buffer                                               |
+-------------------------------+----------------------------------------------------------------------------------+
| pack_restart                  | store atom's data for writing to restart file                                    |
+-------------------------------+----------------------------------------------------------------------------------+
| unpack_restart                | retrieve atom's data from a restart file buffer                                  |
+-------------------------------+----------------------------------------------------------------------------------+
| size_restart                  | size of atom's data                                                              |
+-------------------------------+----------------------------------------------------------------------------------+
| maxsize_restart               | max size of atom's data                                                          |
+-------------------------------+----------------------------------------------------------------------------------+
| setup_pre_force_respa         | same as setup_pre_force, but for rRESPA                                          |
+-------------------------------+----------------------------------------------------------------------------------+
| initial_integrate_respa       | same as initial_integrate, but for rRESPA                                        |
+-------------------------------+----------------------------------------------------------------------------------+
| post_integrate_respa          | called after the first half integration step is done in rRESPA                   |
+-------------------------------+----------------------------------------------------------------------------------+
| pre_force_respa               | same as pre_force, but for rRESPA                                                |
+-------------------------------+----------------------------------------------------------------------------------+
| post_force_respa              | same as post_force, but for rRESPA                                               |
+-------------------------------+----------------------------------------------------------------------------------+
| final_integrate_respa         | same as final_integrate, but for rRESPA                                          |
+-------------------------------+----------------------------------------------------------------------------------+
| min_pre_exchange              | called before atom exchange during minimization                                  |
+-------------------------------+----------------------------------------------------------------------------------+
| min_pre_neighbor              | called before neighbor list build during minimization                            |
+-------------------------------+----------------------------------------------------------------------------------+
| min_post_neighbor             | called after neighbor list build during minimization                             |
+-------------------------------+----------------------------------------------------------------------------------+
| min_pre_force                 | called before pair & molecular forces are computed in minimizer                  |
+-------------------------------+----------------------------------------------------------------------------------+
| min_pre_reverse               | called before reverse communication during minimization                          |
+-------------------------------+----------------------------------------------------------------------------------+
| min_post_force                | called after pair & molecular forces are computed and communicated in minimizer  |
+-------------------------------+----------------------------------------------------------------------------------+
| min_energy                    | return an extra energy contribution during minimization                          |
+-------------------------------+----------------------------------------------------------------------------------+
| min_store                     | store extra data for linesearch based minimization on a LIFO stack               |
+-------------------------------+----------------------------------------------------------------------------------+
| min_pushstore                 | push the minimization LIFO stack one element down                                |
+-------------------------------+----------------------------------------------------------------------------------+
| min_popstore                  | pop the minimization LIFO stack one element up                                   |
+-------------------------------+----------------------------------------------------------------------------------+
| min_clearstore                | clear minimization LIFO stack                                                    |
+-------------------------------+----------------------------------------------------------------------------------+
| min_reset_ref                 | reset the reference state for the line search                                    |
+-------------------------------+----------------------------------------------------------------------------------+
| min_step                      | reset or move forward on line search minimization                                |
+-------------------------------+----------------------------------------------------------------------------------+
| min_dof                       | report number of degrees of freedom *added* by this fix in minimization          |
+-------------------------------+----------------------------------------------------------------------------------+
| max_alpha                     | report maximum allowed step size during linesearch minimization                  |
+-------------------------------+----------------------------------------------------------------------------------+
| pack_forward_comm             | pack a buffer to communicate a per-atom quantity                                 |
+-------------------------------+----------------------------------------------------------------------------------+
| unpack_forward_comm           | unpack a buffer to communicate a per-atom quantity                               |
+-------------------------------+----------------------------------------------------------------------------------+
| pack_reverse_comm_size        | return the dynamic size of the reverse communication buffer                      |
+-------------------------------+----------------------------------------------------------------------------------+
| pack_reverse_comm             | pack a buffer to reverse communicate a per-atom quantity                         |
+-------------------------------+----------------------------------------------------------------------------------+
| unpack_reverse_comm           | unpack a buffer to reverse communicate a per-atom quantity                       |
+-------------------------------+----------------------------------------------------------------------------------+
| dof                           | report number of degrees of freedom *removed* by this fix during MD              |
+-------------------------------+----------------------------------------------------------------------------------+
| compute_scalar                | return a global scalar property that the fix computes                            |
+-------------------------------+----------------------------------------------------------------------------------+
| compute_vector                | return a component of a vector property that the fix computes                    |
+-------------------------------+----------------------------------------------------------------------------------+
| compute_array                 | return a component of an array property that the fix computes                    |
+-------------------------------+----------------------------------------------------------------------------------+
| get_thermo_colname            | return custom thermo column label for a vector or array component                |
+-------------------------------+----------------------------------------------------------------------------------+
| deform                        | called when the box size is changed                                              |
+-------------------------------+----------------------------------------------------------------------------------+
| reset_target                  | called when a change of the target temperature is requested during a run         |
+-------------------------------+----------------------------------------------------------------------------------+
| reset_dt                      | called when a change of the time step is requested during a run                  |
+-------------------------------+----------------------------------------------------------------------------------+
| modify_param                  | called when a fix_modify request is executed                                     |
+-------------------------------+----------------------------------------------------------------------------------+
| memory_usage                  | report memory used by fix                                                        |
+-------------------------------+----------------------------------------------------------------------------------+
| image                         | pass lists of graphics objects to :doc:`dump image fix <dump_image>`             |
+-------------------------------+----------------------------------------------------------------------------------+

Typically, only a small fraction of these methods are defined for a
particular fix.  The ``setmask()`` method is mandatory, as it determines
*when* the fix will be invoked during :doc:`the evolution of a timestep
<Developer_flow>`.  Fixes that perform time integration (\ *nve*, *nvt*,
*npt*\ ) implement ``initial_integrate()`` and ``final_integrate()`` to
perform velocity Verlet time stepping updates.  Fixes that apply forces
implement ``post_force()``, i.e. after the forces on atoms have been
computed and collected to the local atoms.

Fixes that perform diagnostics typically implement ``end_of_step()``.
For such a fix, one of your fix arguments must set the variable "nevery"
which is used to determine when to call the fix and you **must** set
this variable in the constructor of your fix.  By convention, this is
the first argument the fix defines (after the fix-ID, group-ID, and fix
style).

If the fix needs to store information for each atom that persists from
timestep to timestep, it can manage that memory and migrate the info
with the atoms as they move from processors to processor by implementing
the ``grow_arrays()``, ``copy_arrays()``, ``pack_exchange()``, and
``unpack_exchange()`` methods.  Similarly, the ``pack_restart()`` and
``unpack_restart()`` methods can be implemented to store information
about the fix in binary restart files.  If you wish an integrator or
force constraint fix to work with rRESPA (see the :doc:`run_style
<run_style>` command), the ``initial_integrate_respa()``,
``post_force_respa()``, and ``final_integrate_respa()`` methods can be
implemented.  The ``compute_scalar()`` and ``compute_vector()`` methods
enable a fix to contribute values to thermodynamic output, as printed
quantities and/or to be summed to the potential energy of the system.
The corresponding flags (``scalar_flag``, ``vector_flag``, etc.) must be
set in the constructor to tell LAMMPS which of these methods are
implemented.
