# Run all the examples

# -------------------------------------------------
# -------------------------------------------------

# Example 1 = run ab initio MD (AIMD)

# ---

# Run without MDI

lmp_mpi -log log.aimd.alone.1 < in.aimd.alone

# ---

# Run with TCP: 1 proc each

lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.aimd.driver.tcp.1 -in in.aimd.driver &

lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.aimd.engine.tcp.1 -in in.aimd.engine

# ---

# Run with TCP: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.aimd.driver.tcp.3 -in in.aimd.driver &

mpirun -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.aimd.engine.tcp.4 -in in.aimd.engine

# ---

# Run with MPI: 1 proc each

mpirun -np 1 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.aimd.driver.mpi.1 -in in.aimd.driver : -np 1 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.aimd.engine.mpi.1 -in in.aimd.engine

# ---

# Run with MPI: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.aimd.driver.mpi.3 -in in.aimd.driver : -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.aimd.engine.mpi.4 -in in.aimd.engine

# ---

# Run in plugin mode: 1 proc

lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.aimd.driver.plugin.1 -in in.aimd.driver.plugin
mv log.aimd.engine.plugin log.aimd.engine.plugin.1

# ---

# Run in plugin mode: 3 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.aimd.driver.plugin.3 -in in.aimd.driver.plugin
mv log.aimd.engine.plugin log.aimd.engine.plugin.3

# -------------------------------------------------
# -------------------------------------------------

# Example 2 = run LAMMPS, compute QM forces on snapshots from a long run

# ---

# Run without MDI

lmp_mpi -log log.snapshot.alone.1 < in.snapshot.alone
mv dump.snapshot.alone dump.snapshot.alone.1

# ---

# Run with TCP: 1 proc each

lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.snapshot.driver.tcp.1 -in in.snapshot.driver &

lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.snapshot.engine.tcp.1 -in in.snapshot.engine
mv dump.snapshot.driver dump.snapshot.driver.tcp.1

# ---

# Run with TCP: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.snapshot.driver.tcp.3 -in in.snapshot.driver &

mpirun -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.snapshot.engine.tcp.4 -in in.snapshot.engine
mv dump.snapshot.driver dump.snapshot.driver.tcp.4

# ---

# Run with MPI: 1 proc each

mpirun -np 1 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.snapshot.driver.mpi.1 -in in.snapshot.driver : -np 1 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.snapshot.engine.mpi.1 -in in.snapshot.engine
mv dump.snapshot.driver dump.snapshot.driver.mpi.1

# ---

# Run with MPI: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.snapshot.driver.mpi.3 -in in.snapshot.driver : -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.snapshot.engine.mpi.3 -in in.snapshot.engine
mv dump.snapshot.driver dump.snapshot.driver.mpi.4

# ---

# Run in plugin mode: 1 proc

lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.snapshot.driver.plugin.1 -in in.snapshot.driver.plugin
mv log.snapshot.engine.plugin log.snapshot.engine.plugin.1
mv dump.snapshot.driver.plugin dump.snapshot.driver.plugin.1

# ---

# Run in plugin mode: 3 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.snapshot.driver.plugin.3 -in in.snapshot.driver.plugin
mv log.snapshot.engine.plugin log.snapshot.engine.plugin.3
mv dump.snapshot.driver.plugin dump.snapshot.driver.plugin.3

# -------------------------------------------------
# -------------------------------------------------

# Example 3 = run LAMMPS, compute QM forces on series of independent systems

# ---

# Run without MDI

lmp_mpi -log log.series.alone.1 < in.series.alone
mv dump.series.alone dump.series.alone.1

# ---

# Run with TCP: 1 proc each

lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.series.driver.tcp.1 -in in.series.driver &

lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.series.engine.tcp.1 -in in.series.engine
mv dump.series.driver dump.series.driver.tcp.1

# ---

# Run with TCP: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method TCP -port 8021" -log log.series.driver.tcp.3 -in in.series.driver &

mpirun -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method TCP -port 8021 -hostname localhost" -log log.series.engine.tcp.4 -in in.series.engine
mv dump.series.driver dump.series.driver.tcp.4

# ---

# Run with MPI: 1 proc each

mpirun -np 1 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.series.driver.mpi.1 -in in.series.driver : -np 1 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.series.engine.mpi.1 -in in.series.engine
mv dump.series.driver dump.series.driver.mpi.1

# ---

# Run with MPI: 3 procs + 4 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method MPI" -log log.series.driver.mpi.3 -in in.series.driver : -np 4 lmp_mpi -mdi "-name LMP2 -role ENGINE -method MPI" -log log.series.engine.mpi.4 -in in.series.engine
mv dump.series.driver dump.series.driver.mpi.4

# ---

# Run in plugin mode: 1 proc

lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.series.driver.plugin.1 -in in.series.driver.plugin
mv log.series.engine.plugin log.series.engine.plugin.1
mv dump.series.driver.plugin dump.series.driver.plugin.1

# ---

# Run in plugin mode: 3 procs

mpirun -np 3 lmp_mpi -mdi "-name LMP1 -role DRIVER -method LINK -plugin_path /home/sjplimp/lammps/git/src" -log log.series.driver.plugin.3 -in in.series.driver.plugin
mv log.series.engine.plugin log.series.engine.plugin.3
mv dump.series.driver.plugin dump.series.driver.plugin.3

# -------------------------------------------------
# -------------------------------------------------

# Example 4 = Python driver runs a sequence of unrelated LAMMPS calculations

# ---

# Run with TCP: 1 proc each

python3 sequence_driver.py -mdi "-role DRIVER -name sequence -method TCP -port 8021" &

lmp_mpi -mdi "-role ENGINE -name LMP -method TCP -port 8021 -hostname localhost" -log log.sequence.engine.tcp.1 -in in.sequence.python

# ---

# Run with TCP: 2 proc + 4 procs

mpirun -np 2 python3 sequence_driver.py -mdi "-role DRIVER -name sequence -method TCP -port 8021" &

mpirun -np 4 lmp_mpi -mdi "-role ENGINE -name LMP -method TCP -port 8021 -hostname localhost" -log log.sequence.engine.tcp.4 -in in.sequence.python

# ---

# Run with MPI: 1 proc each

mpirun -np 1 python3 sequence_driver.py -mdi "-role DRIVER -name sequence -method MPI" : -np 1 lmp_mpi -mdi "-role ENGINE -name LAMMPS -method MPI" -log log.sequence.engine.mpi.1 -in in.sequence.python

# ---

# Run with MPI: 2 procs + 4 procs

mpirun -np 2 python3 sequence_driver.py -mdi "-role DRIVER -name sequence -method MPI" : -np 4 lmp_mpi -mdi "-role ENGINE -name LMP -method MPI" -log log.sequence.engine.mpi.4 -in in.sequence.python

# ---

# Run in plugin mode: 1 proc

python3 sequence_driver.py -plugin lammps -mdi "-role DRIVER -name sequence -method LINK -plugin_path /home/sjplimp/lammps/git/src" -plugin_args "-log log.sequence.engine.plugin.1 -in in.sequence".python

# ---

# Run in plugin mode: 3 procs

mpirun -np 3 python3 sequence_driver.py -plugin lammps -mdi "-role DRIVER -name sequence -method LINK -plugin_path /home/sjplimp/lammps/git/src" -plugin_args "-log log.sequence.engine.plugin.3 -in in.sequence".python

# -------------------------------------------------
# -------------------------------------------------

# Example 5 = run AIMD with Python driver code and 2 LAMMPS instances as engines

# ---

# Run with TCP: 1 proc each

python3 aimd_driver.py -mdi "-role DRIVER -name aimdpy -method TCP -port 8021" > log.aimdpy.driver.tcp.1 &

lmp_mpi -mdi "-role ENGINE -name MM -method TCP -port 8021 -hostname localhost" -log log.aimdpy.mm.tcp.1 -in in.aimdpy.mm &

lmp_mpi -mdi "-role ENGINE -name QM -method TCP -port 8021 -hostname localhost" -log log.aimdpy.qm.tcp.1 -in in.aimdpy.qm

# ---

# Run with TCP: 2 procs + 2 procs + 3 procs

mpirun -np 2 python3 aimd_driver.py -mdi "-role DRIVER -name aimdpy -method TCP -port 8021" > log.aimdpy.driver.tcp.2 &

mpirun -np 2 lmp_mpi -mdi "-role ENGINE -name MM -method TCP -port 8021 -hostname localhost" -log log.aimdpy.mm.tcp.2 -in in.aimdpy.mm &

mpirun -np 3 lmp_mpi -mdi "-role ENGINE -name QM -method TCP -port 8021 -hostname localhost" -log log.aimdpy.qm.tcp.3 -in in.aimdpy.qm

# ---

# Run with MPI: 1 proc each

mpirun -np 1 python3 aimd_driver.py -mdi "-role DRIVER -name aimdpy -method MPI" : -np 1 lmp_mpi -mdi "-role ENGINE -name MM -method MPI" -log log.aimdpy.mm.mpi.1 -in in.aimdpy.mm : -np 1 lmp_mpi -mdi "-role ENGINE -name QM -method MPI" -log log.aimdpy.qm.mpi.1 -in in.aimdpy.qm > log.aimdpy.driver.mpi.1

# ---

# Run with MPI: 2 procs + 2 procs + 3 procs

mpirun -np 2 python3 aimd_driver.py -mdi "-role DRIVER -name aimdpy -method MPI" : -np 2 lmp_mpi -mdi "-role ENGINE -name MM -method MPI" -log log.aimdpy.mm.mpi.2 -in in.aimdpy.mm : -np 3 lmp_mpi -mdi "-role ENGINE -name QM -method MPI" -log log.aimdpy.qm.mpi.3 -in in.aimdpy.qm > log.aimdpy.driver.mpi.2
