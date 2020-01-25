Warning messages
================

This is an alphabetic list of the WARNING messages LAMMPS prints out
and the reason why.  If the explanation here is not sufficient, the
documentation for the offending command may help.  Warning messages
also list the source file and line number where the warning was
generated.  For example, a message like this:


.. parsed-literal::

   WARNING: Bond atom missing in box size check (domain.cpp:187)

means that line #187 in the file src/domain.cpp generated the error.
Looking in the source code may help you figure out what went wrong.

Note that warning messages from :doc:`user-contributed packages <Packages_user>` are not listed here.  If such a warning
occurs and is not self-explanatory, you'll need to look in the source
code or contact the author of the package.

Doc page with :doc:`ERROR messages <Errors_messages>`


----------




*Adjusting Coulombic cutoff for MSM, new cutoff = %g*
   The adjust/cutoff command is turned on and the Coulombic cutoff has been
   adjusted to match the user-specified accuracy.

*Angle atoms missing at step %ld*
   One or more of 3 atoms needed to compute a particular angle are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the angle has blown apart and an atom is
   too far away.

*Angle style in data file differs from currently defined angle style*
   Self-explanatory.

*Angles are defined but no angle style is set*
   The topology contains angles, but there are no angle forces computed
   since there was no angle\_style command.

*Atom style in data file differs from currently defined atom style*
   Self-explanatory.

*Bond atom missing in box size check*
   The 2nd atoms needed to compute a particular bond is missing on this
   processor.  Typically this is because the pairwise cutoff is set too
   short or the bond has blown apart and an atom is too far away.

*Bond atom missing in image check*
   The 2nd atom in a particular bond is missing on this processor.
   Typically this is because the pairwise cutoff is set too short or the
   bond has blown apart and an atom is too far away.

*Bond atoms missing at step %ld*
   The 2nd atom needed to compute a particular bond is missing on this
   processor.  Typically this is because the pairwise cutoff is set too
   short or the bond has blown apart and an atom is too far away.

*Bond style in data file differs from currently defined bond style*
   Self-explanatory.

*Bonds are defined but no bond style is set*
   The topology contains bonds, but there are no bond forces computed
   since there was no bond\_style command.

*Bond/angle/dihedral extent > half of periodic box length*
   This is a restriction because LAMMPS can be confused about which image
   of an atom in the bonded interaction is the correct one to use.
   "Extent" in this context means the maximum end-to-end length of the
   bond/angle/dihedral.  LAMMPS computes this by taking the maximum bond
   length, multiplying by the number of bonds in the interaction (e.g. 3
   for a dihedral) and adding a small amount of stretch.

*Bond/react: Atom affected by reaction too close to template edge*
   This means an atom which changes type or connectivity during the
   reaction is too close to an 'edge' atom defined in the superimpose
   file. This could cause incorrect assignment of bonds, angle, etc.
   Generally, this means you must include more atoms in your templates,
   such that there are at least two atoms between each atom involved in
   the reaction and an edge atom.

*Both groups in compute group/group have a net charge; the Kspace boundary correction to energy will be non-zero*
   Self-explanatory.

*Calling write\_dump before a full system init.*
   The write\_dump command is used before the system has been fully
   initialized as part of a 'run' or 'minimize' command. Not all dump
   styles and features are fully supported at this point and thus the
   command may fail or produce incomplete or incorrect output. Insert
   a "run 0" command, if a full system init is required.

*Cannot count rigid body degrees-of-freedom before bodies are fully initialized*
   This means the temperature associated with the rigid bodies may be
   incorrect on this timestep.

*Cannot count rigid body degrees-of-freedom before bodies are initialized*
   This means the temperature associated with the rigid bodies may be
   incorrect on this timestep.

*Cannot include log terms without 1/r terms; setting flagHI to 1*
   Self-explanatory.

*Cannot include log terms without 1/r terms; setting flagHI to 1.*
   Self-explanatory.

*Charges are set, but coulombic solver is not used*
   Self-explanatory.

*Charges did not converge at step %ld: %lg*
   Self-explanatory.

*Communication cutoff is 0.0. No ghost atoms will be generated. Atoms may get lost*
   The communication cutoff defaults to the maximum of what is inferred from
   pair and bond styles (will be zero, if none are defined) and what is specified
   via :doc:`comm\_modify cutoff <comm_modify>` (defaults to 0.0).  If this results
   to 0.0, no ghost atoms will be generated and LAMMPS may lose atoms or use
   incorrect periodic images of atoms in interaction lists.  To avoid, either use
   :doc:`pair style zero <pair_zero>` with a suitable cutoff or use :doc:`comm\_modify cutoff <comm_modify>`.

*Communication cutoff is too small for SNAP micro load balancing, increased to %lf*
   Self-explanatory.

*Compute cna/atom cutoff may be too large to find ghost atom neighbors*
   The neighbor cutoff used may not encompass enough ghost atoms
   to perform this operation correctly.

*Computing temperature of portions of rigid bodies*
   The group defined by the temperature compute does not encompass all
   the atoms in one or more rigid bodies, so the change in
   degrees-of-freedom for the atoms in those partial rigid bodies will
   not be accounted for.

*Create\_bonds max distance > minimum neighbor cutoff*
   This means atom pairs for some atom types may not be in the neighbor
   list and thus no bond can be created between them.

*Delete\_atoms cutoff > minimum neighbor cutoff*
   This means atom pairs for some atom types may not be in the neighbor
   list and thus an atom in that pair cannot be deleted.

*Dihedral atoms missing at step %ld*
   One or more of 4 atoms needed to compute a particular dihedral are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the dihedral has blown apart and an atom is
   too far away.

*Dihedral problem*
   Conformation of the 4 listed dihedral atoms is extreme; you may want
   to check your simulation geometry.

*Dihedral problem: %d %ld %d %d %d %d*
   Conformation of the 4 listed dihedral atoms is extreme; you may want
   to check your simulation geometry.

*Dihedral style in data file differs from currently defined dihedral style*
   Self-explanatory.

*Dihedrals are defined but no dihedral style is set*
   The topology contains dihedrals, but there are no dihedral forces computed
   since there was no dihedral\_style command.

*Dump dcd/xtc timestamp may be wrong with fix dt/reset*
   If the fix changes the timestep, the dump dcd file will not
   reflect the change.

*Energy due to X extra global DOFs will be included in minimizer energies*
   When using fixes like box/relax, the potential energy used by the minimizer
   is augmented by an additional energy provided by the fix. Thus the printed
   converged energy may be different from the total potential energy.

*Estimated error in splitting of dispersion coeffs is %g*
   Error is greater than 0.0001 percent.

*Ewald/disp Newton solver failed, using old method to estimate g\_ewald*
   Self-explanatory. Choosing a different cutoff value may help.

*FENE bond too long*
   A FENE bond has stretched dangerously far.  It's interaction strength
   will be truncated to attempt to prevent the bond from blowing up.

*FENE bond too long: %ld %d %d %g*
   A FENE bond has stretched dangerously far.  It's interaction strength
   will be truncated to attempt to prevent the bond from blowing up.

*FENE bond too long: %ld %g*
   A FENE bond has stretched dangerously far.  It's interaction strength
   will be truncated to attempt to prevent the bond from blowing up.

*Fix SRD walls overlap but fix srd overlap not set*
   You likely want to set this in your input script.

* Fix bond/create is used multiple times or with fix bond/break - may not work as expected*
   When using fix bond/create multiple times or in combination with
   fix bond/break, the individual fix instances do not share information
   about changes they made at the same time step and thus it may result
   in unexpected behavior.

*Fix bond/swap will ignore defined angles*
   See the doc page for fix bond/swap for more info on this
   restriction.

*Fix deposit near setting < possible overlap separation %g*
   This test is performed for finite size particles with a diameter, not
   for point particles.  The near setting is smaller than the particle
   diameter which can lead to overlaps.

*Fix evaporate may delete atom with non-zero molecule ID*
   This is probably an error, since you should not delete only one atom
   of a molecule.

*Fix gcmc using full\_energy option*
   Fix gcmc has automatically turned on the full\_energy option since it
   is required for systems like the one specified by the user. User input
   included one or more of the following: kspace, triclinic, a hybrid
   pair style, an eam pair style, or no "single" function for the pair
   style.

*Fix langevin gjf using random gaussians is not implemented with kokkos*
This will most likely cause errors in kinetic fluctuations.

*Fix property/atom mol or charge w/out ghost communication*
   A model typically needs these properties defined for ghost atoms.

*Fix qeq CG convergence failed (%g) after %d iterations at %ld step*
   Self-explanatory.

*Fix qeq has non-zero lower Taper radius cutoff*
   Absolute value must be <= 0.01.

*Fix qeq has very low Taper radius cutoff*
   Value should typically be >= 5.0.

*Fix qeq/dynamic tolerance may be too small for damped dynamics*
   Self-explanatory.

*Fix qeq/fire tolerance may be too small for damped fires*
   Self-explanatory.

*Fix rattle should come after all other integration fixes*
   This fix is designed to work after all other integration fixes change
   atom positions.  Thus it should be the last integration fix specified.
   If not, it will not satisfy the desired constraints as well as it
   otherwise would.

*Fix recenter should come after all other integration fixes*
   Other fixes may change the position of the center-of-mass, so
   fix recenter should come last.

*Fix srd SRD moves may trigger frequent reneighboring*
   This is because the SRD particles may move long distances.

*Fix srd grid size > 1/4 of big particle diameter*
   This may cause accuracy problems.

*Fix srd particle moved outside valid domain*
   This may indicate a problem with your simulation parameters.

*Fix srd particles may move > big particle diameter*
   This may cause accuracy problems.

*Fix srd viscosity < 0.0 due to low SRD density*
   This may cause accuracy problems.

*Fixes cannot send data in Kokkos communication, switching to classic communication*
   This is current restriction with Kokkos.

*For better accuracy use 'pair\_modify table 0'*
   The user-specified force accuracy cannot be achieved unless the table
   feature is disabled by using 'pair\_modify table 0'.

*Geometric mixing assumed for 1/r\^6 coefficients*
   Self-explanatory.

*Group for fix\_modify temp != fix group*
   The fix\_modify command is specifying a temperature computation that
   computes a temperature on a different group of atoms than the fix
   itself operates on.  This is probably not what you want to do.

*H matrix size has been exceeded: m\_fill=%d H.m=%d\n*
   This is the size of the matrix.

*Ignoring unknown or incorrect info command flag*
   Self-explanatory.  An unknown argument was given to the info command.
   Compare your input with the documentation.

*Improper atoms missing at step %ld*
   One or more of 4 atoms needed to compute a particular improper are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the improper has blown apart and an atom is
   too far away.

*Improper problem: %d %ld %d %d %d %d*
   Conformation of the 4 listed improper atoms is extreme; you may want
   to check your simulation geometry.

*Improper style in data file differs from currently defined improper style*
   Self-explanatory.

*Impropers are defined but no improper style is set*
   The topology contains impropers, but there are no improper forces computed
   since there was no improper\_style command.

*Inconsistent image flags*
   The image flags for a pair on bonded atoms appear to be inconsistent.
   Inconsistent means that when the coordinates of the two atoms are
   unwrapped using the image flags, the two atoms are far apart.
   Specifically they are further apart than half a periodic box length.
   Or they are more than a box length apart in a non-periodic dimension.
   This is usually due to the initial data file not having correct image
   flags for the 2 atoms in a bond that straddles a periodic boundary.
   They should be different by 1 in that case.  This is a warning because
   inconsistent image flags will not cause problems for dynamics or most
   LAMMPS simulations.  However they can cause problems when such atoms
   are used with the fix rigid or replicate commands.  Note that if you
   have an infinite periodic crystal with bonds then it is impossible to
   have fully consistent image flags, since some bonds will cross
   periodic boundaries and connect two atoms with the same image
   flag.

*Increasing communication cutoff for GPU style*
   The pair style has increased the communication cutoff to be consistent with
   the communication cutoff requirements for this pair style when run on the GPU.

*KIM Model does not provide 'energy'; Potential energy will be zero*
   Self-explanatory.

*KIM Model does not provide 'forces'; Forces will be zero*
   Self-explanatory.

*KIM Model does not provide 'particleEnergy'; energy per atom will be zero*
   Self-explanatory.

*KIM Model does not provide 'particleVirial'; virial per atom will be zero*
   Self-explanatory.

*Kspace\_modify slab param < 2.0 may cause unphysical behavior*
   The kspace\_modify slab parameter should be larger to insure periodic
   grids padded with empty space do not overlap.

*Less insertions than requested*
   The fix pour command was unsuccessful at finding open space
   for as many particles as it tried to insert.

*Library error in lammps\_gather\_atoms*
   This library function cannot be used if atom IDs are not defined
   or are not consecutively numbered.

*Library error in lammps\_scatter\_atoms*
   This library function cannot be used if atom IDs are not defined or
   are not consecutively numbered, or if no atom map is defined.  See the
   atom\_modify command for details about atom maps.

*Likewise 1-2 special neighbor interactions != 1.0*
   The topology contains bonds, but there is no bond style defined
   and a 1-2 special neighbor scaling factor was not 1.0. This
   means that pair style interactions may have scaled or missing
   pairs in the neighbor list in expectation of interactions for
   those pairs being computed from the bond style.

*Likewise 1-3 special neighbor interactions != 1.0*
   The topology contains angles, but there is no angle style defined
   and a 1-3 special neighbor scaling factor was not 1.0. This
   means that pair style interactions may have scaled or missing
   pairs in the neighbor list in expectation of interactions for
   those pairs being computed from the angle style.

*Likewise 1-4 special neighbor interactions != 1.0*
   The topology contains dihedrals, but there is no dihedral style defined
   and a 1-4 special neighbor scaling factor was not 1.0. This
   means that pair style interactions may have scaled or missing
   pairs in the neighbor list in expectation of interactions for
   those pairs being computed from the dihedral style.

*Lost atoms via change\_box: original %ld current %ld*
   The command options you have used caused atoms to be lost.

*Lost atoms via displace\_atoms: original %ld current %ld*
   The command options you have used caused atoms to be lost.

*Lost atoms: original %ld current %ld*
   Lost atoms are checked for each time thermo output is done.  See the
   thermo\_modify lost command for options.  Lost atoms usually indicate
   bad dynamics, e.g. atoms have been blown far out of the simulation
   box, or moved further than one processor's sub-domain away before
   reneighboring.

*MSM mesh too small, increasing to 2 points in each direction*
   Self-explanatory.

*Mismatch between velocity and compute groups*
   The temperature computation used by the velocity command will not be
   on the same group of atoms that velocities are being set for.

*Mixing forced for lj coefficients*
   Self-explanatory.

*Molecule attributes do not match system attributes*
   An attribute is specified (e.g. diameter, charge) that is
   not defined for the specified atom style.

*Molecule has bond topology but no special bond settings*
   This means the bonded atoms will not be excluded in pair-wise
   interactions.

*Molecule template for create\_atoms has multiple molecules*
   The create\_atoms command will only create molecules of a single type,
   i.e. the first molecule in the template.

*Molecule template for fix gcmc has multiple molecules*
   The fix gcmc command will only create molecules of a single type,
   i.e. the first molecule in the template.

*Molecule template for fix shake has multiple molecules*
   The fix shake command will only recognize molecules of a single
   type, i.e. the first molecule in the template.

*More than one compute centro/atom*
   It is not efficient to use compute centro/atom more than once.

*More than one compute cluster/atom*
   It is not efficient to use compute cluster/atom  more than once.

*More than one compute cna/atom defined*
   It is not efficient to use compute cna/atom  more than once.

*More than one compute contact/atom*
   It is not efficient to use compute contact/atom more than once.

*More than one compute coord/atom*
   It is not efficient to use compute coord/atom more than once.

*More than one compute damage/atom*
   It is not efficient to use compute ke/atom more than once.

*More than one compute dilatation/atom*
   Self-explanatory.

*More than one compute erotate/sphere/atom*
   It is not efficient to use compute erorate/sphere/atom more than once.

*More than one compute hexorder/atom*
   It is not efficient to use compute hexorder/atom more than once.

*More than one compute ke/atom*
   It is not efficient to use compute ke/atom more than once.

*More than one compute orientorder/atom*
   It is not efficient to use compute orientorder/atom more than once.

*More than one compute plasticity/atom*
   Self-explanatory.

*More than one compute sna/atom*
   Self-explanatory.

*More than one compute snad/atom*
   Self-explanatory.

*More than one compute snav/atom*
   Self-explanatory.

*More than one fix poems*
   It is not efficient to use fix poems more than once.

*More than one fix rigid*
   It is not efficient to use fix rigid more than once.

*Neighbor exclusions used with KSpace solver may give inconsistent Coulombic energies*
   This is because excluding specific pair interactions also excludes
   them from long-range interactions which may not be the desired effect.
   The special\_bonds command handles this consistently by insuring
   excluded (or weighted) 1-2, 1-3, 1-4 interactions are treated
   consistently by both the short-range pair style and the long-range
   solver.  This is not done for exclusions of charged atom pairs via the
   neigh\_modify exclude command.

*New thermo\_style command, previous thermo\_modify settings will be lost*
   If a thermo\_style command is used after a thermo\_modify command, the
   settings changed by the thermo\_modify command will be reset to their
   default values.  This is because the thermo\_modify command acts on
   the currently defined thermo style, and a thermo\_style command creates
   a new style.

*No Kspace calculation with verlet/split*
   The 2nd partition performs a kspace calculation so the kspace\_style
   command must be used.

*No automatic unit conversion to XTC file format conventions possible for units lj*
   This means no scaling will be performed.

*No fixes defined, atoms won't move*
   If you are not using a fix like nve, nvt, npt then atom velocities and
   coordinates will not be updated during timestepping.

*No joints between rigid bodies, use fix rigid instead*
   The bodies defined by fix poems are not connected by joints.  POEMS
   will integrate the body motion, but it would be more efficient to use
   fix rigid.

*Not using real units with pair reax*
   This is most likely an error, unless you have created your own ReaxFF
   parameter file in a different set of units.

*Number of MSM mesh points changed to be a multiple of 2*
   MSM requires that the number of grid points in each direction be a multiple
   of two and the number of grid points in one or more directions have been
   adjusted to meet this requirement.

*OMP\_NUM\_THREADS environment is not set.*
   This environment variable must be set appropriately to use the
   USER-OMP package.

*One or more atoms are time integrated more than once*
   This is probably an error since you typically do not want to
   advance the positions or velocities of an atom more than once
   per timestep.

*One or more chunks do not contain all atoms in molecule*
   This may not be what you intended.

*One or more dynamic groups may not be updated at correct point in timestep*
   If there are other fixes that act immediately after the initial stage
   of time integration within a timestep (i.e. after atoms move), then
   the command that sets up the dynamic group should appear after those
   fixes.  This will insure that dynamic group assignments are made
   after all atoms have moved.

*One or more respa levels compute no forces*
   This is computationally inefficient.

*Pair COMB charge %.10f with force %.10f hit max barrier*
   Something is possibly wrong with your model.

*Pair COMB charge %.10f with force %.10f hit min barrier*
   Something is possibly wrong with your model.

*Pair brownian needs newton pair on for momentum conservation*
   Self-explanatory.

*Pair dpd needs newton pair on for momentum conservation*
   Self-explanatory.

*Pair dsmc: num\_of\_collisions > number\_of\_A*
   Collision model in DSMC is breaking down.

*Pair dsmc: num\_of\_collisions > number\_of\_B*
   Collision model in DSMC is breaking down.

*Pair style in data file differs from currently defined pair style*
   Self-explanatory.

*Pair style restartinfo set but has no restart support*
   This pair style has a bug, where it does not support reading and
   writing information to a restart file, but does not set the member
   variable "restartinfo" to 0 as required in that case.

*Particle deposition was unsuccessful*
   The fix deposit command was not able to insert as many atoms as
   needed.  The requested volume fraction may be too high, or other atoms
   may be in the insertion region.

*Proc sub-domain size < neighbor skin, could lead to lost atoms*
   The decomposition of the physical domain (likely due to load
   balancing) has led to a processor's sub-domain being smaller than the
   neighbor skin in one or more dimensions.  Since reneighboring is
   triggered by atoms moving the skin distance, this may lead to lost
   atoms, if an atom moves all the way across a neighboring processor's
   sub-domain before reneighboring is triggered.

*Reducing PPPM order b/c stencil extends beyond nearest neighbor processor*
   This may lead to a larger grid than desired.  See the kspace\_modify overlap
   command to prevent changing of the PPPM order.

*Reducing PPPMDisp Coulomb order b/c stencil extends beyond neighbor processor*
   This may lead to a larger grid than desired.  See the kspace\_modify overlap
   command to prevent changing of the PPPM order.

*Reducing PPPMDisp dispersion order b/c stencil extends beyond neighbor processor*
   This may lead to a larger grid than desired.  See the kspace\_modify overlap
   command to prevent changing of the PPPM order.

*Replacing a fix, but new group != old group*
   The ID and style of a fix match for a fix you are changing with a fix
   command, but the new group you are specifying does not match the old
   group.

*Replicating in a non-periodic dimension*
   The parameters for a replicate command will cause a non-periodic
   dimension to be replicated; this may cause unwanted behavior.

*Resetting reneighboring criteria during PRD*
   A PRD simulation requires that neigh\_modify settings be delay = 0,
   every = 1, check = yes.  Since these settings were not in place,
   LAMMPS changed them and will restore them to their original values
   after the PRD simulation.

*Resetting reneighboring criteria during TAD*
   A TAD simulation requires that neigh\_modify settings be delay = 0,
   every = 1, check = yes.  Since these settings were not in place,
   LAMMPS changed them and will restore them to their original values
   after the PRD simulation.

*Resetting reneighboring criteria during minimization*
   Minimization requires that neigh\_modify settings be delay = 0, every =
   1, check = yes.  Since these settings were not in place, LAMMPS
   changed them and will restore them to their original values after the
   minimization.

*Restart file used different # of processors*
   The restart file was written out by a LAMMPS simulation running on a
   different number of processors.  Due to round-off, the trajectories of
   your restarted simulation may diverge a little more quickly than if
   you ran on the same # of processors.

*Restart file used different 3d processor grid*
   The restart file was written out by a LAMMPS simulation running on a
   different 3d grid of processors.  Due to round-off, the trajectories
   of your restarted simulation may diverge a little more quickly than if
   you ran on the same # of processors.

*Restart file used different boundary settings, using restart file values*
   Your input script cannot change these restart file settings.

*Restart file used different newton bond setting, using restart file value*
   The restart file value will override the setting in the input script.

*Restart file used different newton pair setting, using input script value*
   The input script value will override the setting in the restart file.

*Restrain problem: %d %ld %d %d %d %d*
   Conformation of the 4 listed dihedral atoms is extreme; you may want
   to check your simulation geometry.

*Running PRD with only one replica*
   This is allowed, but you will get no parallel speed-up.

*SRD bin shifting turned on due to small lamda*
   This is done to try to preserve accuracy.

*SRD bin size for fix srd differs from user request*
   Fix SRD had to adjust the bin size to fit the simulation box.  See the
   cubic keyword if you want this message to be an error vs warning.

*SRD bins for fix srd are not cubic enough*
   The bin shape is not within tolerance of cubic.  See the cubic
   keyword if you want this message to be an error vs warning.

*SRD particle %d started inside big particle %d on step %ld bounce %d*
   See the inside keyword if you want this message to be an error vs
   warning.

*SRD particle %d started inside wall %d on step %ld bounce %d*
   See the inside keyword if you want this message to be an error vs
   warning.

*Shake determinant < 0.0*
   The determinant of the quadratic equation being solved for a single
   cluster specified by the fix shake command is numerically suspect.  LAMMPS
   will set it to 0.0 and continue.

*Shell command '%s' failed with error '%s'*
   Self-explanatory.

*Shell command returned with non-zero status*
   This may indicate the shell command did not operate as expected.

*Should not allow rigid bodies to bounce off reflecting walls*
   LAMMPS allows this, but their dynamics are not computed correctly.

*Should not use fix nve/limit with fix shake or fix rattle*
   This will lead to invalid constraint forces in the SHAKE/RATTLE
   computation.

*Simulations might be very slow because of large number of structure factors*
   Self-explanatory.

*Slab correction not needed for MSM*
   Slab correction is intended to be used with Ewald or PPPM and is not needed by MSM.

*Specifying an 'subset' value of '0' is equivalent to no 'subset' keyword*
   Self-explanatory.

*System is not charge neutral, net charge = %g*
   The total charge on all atoms on the system is not 0.0.
   For some KSpace solvers this is only a warning.

*Table inner cutoff >= outer cutoff*
   You specified an inner cutoff for a Coulombic table that is longer
   than the global cutoff.  Probably not what you wanted.

*Temperature for MSST is not for group all*
   User-assigned temperature to MSST fix does not compute temperature for
   all atoms.  Since MSST computes a global pressure, the kinetic energy
   contribution from the temperature is assumed to also be for all atoms.
   Thus the pressure used by MSST could be inaccurate.

*Temperature for NPT is not for group all*
   User-assigned temperature to NPT fix does not compute temperature for
   all atoms.  Since NPT computes a global pressure, the kinetic energy
   contribution from the temperature is assumed to also be for all atoms.
   Thus the pressure used by NPT could be inaccurate.

*Temperature for fix modify is not for group all*
   The temperature compute is being used with a pressure calculation
   which does operate on group all, so this may be inconsistent.

*Temperature for thermo pressure is not for group all*
   User-assigned temperature to thermo via the thermo\_modify command does
   not compute temperature for all atoms.  Since thermo computes a global
   pressure, the kinetic energy contribution from the temperature is
   assumed to also be for all atoms.  Thus the pressure printed by thermo
   could be inaccurate.

*The fix ave/spatial command has been replaced by the more flexible fix ave/chunk and compute chunk/atom commands -- fix ave/spatial will be removed in the summer of 2015*
   Self-explanatory.

*The minimizer does not re-orient dipoles when using fix efield*
   This means that only the atom coordinates will be minimized,
   not the orientation of the dipoles.

*Too many common neighbors in CNA %d times*
   More than the maximum # of neighbors was found multiple times.  This
   was unexpected.

*Too many inner timesteps in fix ttm*
   Self-explanatory.

*Too many neighbors in CNA for %d atoms*
   More than the maximum # of neighbors was found multiple times.  This
   was unexpected.

*Triclinic box skew is large*
   The displacement in a skewed direction is normally required to be less
   than half the box length in that dimension.  E.g. the xy tilt must be
   between -half and +half of the x box length.  You have relaxed the
   constraint using the box tilt command, but the warning means that a
   LAMMPS simulation may be inefficient as a result.

*Use special bonds = 0,1,1 with bond style fene*
   Most FENE models need this setting for the special\_bonds command.

*Use special bonds = 0,1,1 with bond style fene/expand*
   Most FENE models need this setting for the special\_bonds command.

*Using a many-body potential with bonds/angles/dihedrals and special\_bond exclusions*
   This is likely not what you want to do.  The exclusion settings will
   eliminate neighbors in the neighbor list, which the many-body potential
   needs to calculated its terms correctly.

*Using compute temp/deform with inconsistent fix deform remap option*
   Fix nvt/sllod assumes deforming atoms have a velocity profile provided
   by "remap v" or "remap none" as a fix deform option.

*Using compute temp/deform with no fix deform defined*
   This is probably an error, since it makes little sense to use
   compute temp/deform in this case.

*Using fix srd with box deformation but no SRD thermostat*
   The deformation will heat the SRD particles so this can
   be dangerous.

*Using kspace solver on system with no charge*
   Self-explanatory.

*Using largest cut-off for lj/long/dipole/long long long*
   Self-explanatory.

*Using largest cutoff for buck/long/coul/long*
   Self-explanatory.

*Using largest cutoff for lj/long/coul/long*
   Self-explanatory.

*Using largest cutoff for pair\_style lj/long/tip4p/long*
   Self-explanatory.

*Using package gpu without any pair style defined*
   Self-explanatory.

*Using pair potential shift with pair\_modify compute no*
   The shift effects will thus not be computed.

*Using pair tail corrections with nonperiodic system*
   This is probably a bogus thing to do, since tail corrections are
   computed by integrating the density of a periodic system out to
   infinity.

*Using pair tail corrections with pair\_modify compute no*
   The tail corrections will thus not be computed.

*pair style reax is now deprecated and will soon be retired. Users should switch to pair\_style reax/c*
   Self-explanatory.




.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
