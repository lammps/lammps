Fix styles
==========

In LAMMPS, a "fix" is any operation that is computed during
timestepping that alters some property of the system.  Essentially
everything that happens during a simulation besides force computation,
neighbor list construction, and output, is a "fix".  This includes
time integration (update of coordinates and velocities), force
constraints or boundary conditions (SHAKE or walls), and diagnostics
(compute a diffusion coefficient).  New styles can be created to add
new options to LAMMPS.

Fix_setforce.cpp is a simple example of setting forces on atoms to
prescribed values.  There are dozens of fix options already in LAMMPS;
choose one as a template that is similar to what you want to
implement.

Here is a brief description of methods you can define in your new
derived class.  See fix.h for details.

+---------------------------+--------------------------------------------------------------------------------------------+
| setmask                   | determines when the fix is called during the timestep (required)                           |
+---------------------------+--------------------------------------------------------------------------------------------+
| init                      | initialization before a run (optional)                                                     |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup_pre_exchange        | called before atom exchange in setup (optional)                                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup_pre_force           | called before force computation in setup (optional)                                        |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup                     | called immediately before the first timestep and after forces are computed (optional)      |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_setup_pre_force       | like setup_pre_force, but for minimizations instead of MD runs (optional)                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_setup                 | like setup, but for minimizations instead of MD runs (optional)                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| initial_integrate         | called at very beginning of each timestep (optional)                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre_exchange              | called before atom exchange on re-neighboring steps (optional)                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre_neighbor              | called before neighbor list build (optional)                                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre_force                 | called before pair & molecular forces are computed (optional)                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| post_force                | called after pair & molecular forces are computed and communicated (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| final_integrate           | called at end of each timestep (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| end_of_step               | called at very end of timestep (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| write_restart             | dumps fix info to restart file (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| restart                   | uses info from restart file to re-initialize the fix (optional)                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| grow_arrays               | allocate memory for atom-based arrays used by fix (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| copy_arrays               | copy atom info when an atom migrates to a new processor (optional)                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack_exchange             | store atom's data in a buffer (optional)                                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack_exchange           | retrieve atom's data from a buffer (optional)                                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack_restart              | store atom's data for writing to restart file (optional)                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack_restart            | retrieve atom's data from a restart file buffer (optional)                                 |
+---------------------------+--------------------------------------------------------------------------------------------+
| size_restart              | size of atom's data (optional)                                                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| maxsize_restart           | max size of atom's data (optional)                                                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup_pre_force_respa     | same as setup_pre_force, but for rRESPA (optional)                                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| initial_integrate_respa   | same as initial_integrate, but for rRESPA (optional)                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| post_integrate_respa      | called after the first half integration step is done in rRESPA (optional)                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre_force_respa           | same as pre_force, but for rRESPA (optional)                                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| post_force_respa          | same as post_force, but for rRESPA (optional)                                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| final_integrate_respa     | same as final_integrate, but for rRESPA (optional)                                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_pre_force             | called after pair & molecular forces are computed in minimizer (optional)                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_post_force            | called after pair & molecular forces are computed and communicated in minimizer (optional) |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_store                 | store extra data for linesearch based minimization on a LIFO stack (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_pushstore             | push the minimization LIFO stack one element down (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_popstore              | pop the minimization LIFO stack one element up (optional)                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_clearstore            | clear minimization LIFO stack (optional)                                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_step                  | reset or move forward on line search minimization (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| min_dof                   | report number of degrees of freedom *added* by this fix in minimization (optional)         |
+---------------------------+--------------------------------------------------------------------------------------------+
| max_alpha                 | report maximum allowed step size during linesearch minimization (optional)                 |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack_comm                 | pack a buffer to communicate a per-atom quantity (optional)                                |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack_comm               | unpack a buffer to communicate a per-atom quantity (optional)                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack_reverse_comm         | pack a buffer to reverse communicate a per-atom quantity (optional)                        |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack_reverse_comm       | unpack a buffer to reverse communicate a per-atom quantity (optional)                      |
+---------------------------+--------------------------------------------------------------------------------------------+
| dof                       | report number of degrees of freedom *removed* by this fix during MD (optional)             |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute_scalar            | return a global scalar property that the fix computes (optional)                           |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute_vector            | return a component of a vector property that the fix computes (optional)                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute_array             | return a component of an array property that the fix computes (optional)                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| deform                    | called when the box size is changed (optional)                                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| reset_target              | called when a change of the target temperature is requested during a run (optional)        |
+---------------------------+--------------------------------------------------------------------------------------------+
| reset_dt                  | is called when a change of the time step is requested during a run (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| modify_param              | called when a fix_modify request is executed (optional)                                    |
+---------------------------+--------------------------------------------------------------------------------------------+
| memory_usage              | report memory used by fix (optional)                                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| thermo                    | compute quantities for thermodynamic output (optional)                                     |
+---------------------------+--------------------------------------------------------------------------------------------+

Typically, only a small fraction of these methods are defined for a
particular fix.  Setmask is mandatory, as it determines when the fix
will be invoked during the timestep.  Fixes that perform time
integration (\ *nve*, *nvt*, *npt*\ ) implement initial_integrate() and
final_integrate() to perform velocity Verlet updates.  Fixes that
constrain forces implement post_force().

Fixes that perform diagnostics typically implement end_of_step().  For
an end_of_step fix, one of your fix arguments must be the variable
"nevery" which is used to determine when to call the fix and you must
set this variable in the constructor of your fix.  By convention, this
is the first argument the fix defines (after the ID, group-ID, style).

If the fix needs to store information for each atom that persists from
timestep to timestep, it can manage that memory and migrate the info
with the atoms as they move from processors to processor by
implementing the grow_arrays, copy_arrays, pack_exchange, and
unpack_exchange methods.  Similarly, the pack_restart and
unpack_restart methods can be implemented to store information about
the fix in restart files.  If you wish an integrator or force
constraint fix to work with rRESPA (see the :doc:`run_style <run_style>`
command), the initial_integrate, post_force_integrate, and
final_integrate_respa methods can be implemented.  The thermo method
enables a fix to contribute values to thermodynamic output, as printed
quantities and/or to be summed to the potential energy of the system.
