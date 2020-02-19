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

Fix\_setforce.cpp is a simple example of setting forces on atoms to
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
| setup\_pre\_exchange      | called before atom exchange in setup (optional)                                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup\_pre\_force         | called before force computation in setup (optional)                                        |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup                     | called immediately before the 1st timestep and after forces are computed (optional)        |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_setup\_pre\_force    | like setup\_pre\_force, but for minimizations instead of MD runs (optional)                |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_setup                | like setup, but for minimizations instead of MD runs (optional)                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| initial\_integrate        | called at very beginning of each timestep (optional)                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre\_exchange             | called before atom exchange on re-neighboring steps (optional)                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre\_neighbor             | called before neighbor list build (optional)                                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre\_force                | called before pair & molecular forces are computed (optional)                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| post\_force               | called after pair & molecular forces are computed and communicated (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| final\_integrate          | called at end of each timestep (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| end\_of\_step             | called at very end of timestep (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| write\_restart            | dumps fix info to restart file (optional)                                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| restart                   | uses info from restart file to re-initialize the fix (optional)                            |
+---------------------------+--------------------------------------------------------------------------------------------+
| grow\_arrays              | allocate memory for atom-based arrays used by fix (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| copy\_arrays              | copy atom info when an atom migrates to a new processor (optional)                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack\_exchange            | store atom's data in a buffer (optional)                                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack\_exchange          | retrieve atom's data from a buffer (optional)                                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack\_restart             | store atom's data for writing to restart file (optional)                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack\_restart           | retrieve atom's data from a restart file buffer (optional)                                 |
+---------------------------+--------------------------------------------------------------------------------------------+
| size\_restart             | size of atom's data (optional)                                                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| maxsize\_restart          | max size of atom's data (optional)                                                         |
+---------------------------+--------------------------------------------------------------------------------------------+
| setup\_pre\_force\_respa  | same as setup\_pre\_force, but for rRESPA (optional)                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| initial\_integrate\_respa | same as initial\_integrate, but for rRESPA (optional)                                      |
+---------------------------+--------------------------------------------------------------------------------------------+
| post\_integrate\_respa    | called after the first half integration step is done in rRESPA (optional)                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| pre\_force\_respa         | same as pre\_force, but for rRESPA (optional)                                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| post\_force\_respa        | same as post\_force, but for rRESPA (optional)                                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| final\_integrate\_respa   | same as final\_integrate, but for rRESPA (optional)                                        |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_pre\_force           | called after pair & molecular forces are computed in minimizer (optional)                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_post\_force          | called after pair & molecular forces are computed and communicated in minimizer (optional) |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_store                | store extra data for linesearch based minimization on a LIFO stack (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_pushstore            | push the minimization LIFO stack one element down (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_popstore             | pop the minimization LIFO stack one element up (optional)                                  |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_clearstore           | clear minimization LIFO stack (optional)                                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_step                 | reset or move forward on line search minimization (optional)                               |
+---------------------------+--------------------------------------------------------------------------------------------+
| min\_dof                  | report number of degrees of freedom *added* by this fix in minimization (optional)         |
+---------------------------+--------------------------------------------------------------------------------------------+
| max\_alpha                | report maximum allowed step size during linesearch minimization (optional)                 |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack\_comm                | pack a buffer to communicate a per-atom quantity (optional)                                |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack\_comm              | unpack a buffer to communicate a per-atom quantity (optional)                              |
+---------------------------+--------------------------------------------------------------------------------------------+
| pack\_reverse\_comm       | pack a buffer to reverse communicate a per-atom quantity (optional)                        |
+---------------------------+--------------------------------------------------------------------------------------------+
| unpack\_reverse\_comm     | unpack a buffer to reverse communicate a per-atom quantity (optional)                      |
+---------------------------+--------------------------------------------------------------------------------------------+
| dof                       | report number of degrees of freedom *removed* by this fix during MD (optional)             |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute\_scalar           | return a global scalar property that the fix computes (optional)                           |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute\_vector           | return a component of a vector property that the fix computes (optional)                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| compute\_array            | return a component of an array property that the fix computes (optional)                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| deform                    | called when the box size is changed (optional)                                             |
+---------------------------+--------------------------------------------------------------------------------------------+
| reset\_target             | called when a change of the target temperature is requested during a run (optional)        |
+---------------------------+--------------------------------------------------------------------------------------------+
| reset\_dt                 | is called when a change of the time step is requested during a run (optional)              |
+---------------------------+--------------------------------------------------------------------------------------------+
| modify\_param             | called when a fix\_modify request is executed (optional)                                   |
+---------------------------+--------------------------------------------------------------------------------------------+
| memory\_usage             | report memory used by fix (optional)                                                       |
+---------------------------+--------------------------------------------------------------------------------------------+
| thermo                    | compute quantities for thermodynamic output (optional)                                     |
+---------------------------+--------------------------------------------------------------------------------------------+

Typically, only a small fraction of these methods are defined for a
particular fix.  Setmask is mandatory, as it determines when the fix
will be invoked during the timestep.  Fixes that perform time
integration (\ *nve*\ , *nvt*\ , *npt*\ ) implement initial\_integrate() and
final\_integrate() to perform velocity Verlet updates.  Fixes that
constrain forces implement post\_force().

Fixes that perform diagnostics typically implement end\_of\_step().  For
an end\_of\_step fix, one of your fix arguments must be the variable
"nevery" which is used to determine when to call the fix and you must
set this variable in the constructor of your fix.  By convention, this
is the first argument the fix defines (after the ID, group-ID, style).

If the fix needs to store information for each atom that persists from
timestep to timestep, it can manage that memory and migrate the info
with the atoms as they move from processors to processor by
implementing the grow\_arrays, copy\_arrays, pack\_exchange, and
unpack\_exchange methods.  Similarly, the pack\_restart and
unpack\_restart methods can be implemented to store information about
the fix in restart files.  If you wish an integrator or force
constraint fix to work with rRESPA (see the :doc:`run_style <run_style>`
command), the initial\_integrate, post\_force\_integrate, and
final\_integrate\_respa methods can be implemented.  The thermo method
enables a fix to contribute values to thermodynamic output, as printed
quantities and/or to be summed to the potential energy of the system.
