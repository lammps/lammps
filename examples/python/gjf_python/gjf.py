"""Made by Charlie Sievers Ph.D. Candidate, UC Davis, Donadio Lab 2019"""

from mpi4py import MPI
from lammps import lammps
import lammps_tools as lt
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" LAMMPS  VARIABLES """

# new file or restart
run_no = 0

# data files
infile = "argon.lmp"
restart_file = "final_restart.{}".format(run_no)
ff_file = "ff-argon.lmp"
outfile = "output.dat"

# write final_restart
write_final_restart = False

# random numbers
seed0 = 2357
seed1 = 26588
seed2 = 10669

# MD Parameters
# number of steps
nsteps = 50000
# timestep
# dt = 0.001
# starting simulation temp
temp_start = 10
# final simulation temp
temp_final = 10
# relaxation time
trel = 1
# trajectory frequency
ntraj = 0

# Ensemble 0 = GJF u, 1 = GJF v, 2 = Nose-Hoover, 3 = Langevin, 4 = BDP (Currently all NVT)
ensemble = 0

# Output Parameters
nthermo = 200
nout = int(nsteps / nthermo)  # Important

# output to screen and log file?
lammps_output = False
# Lammps Thermo
thermo = False

python_output = True

# Write output to file?
write_output = False

if write_output is True:
    data = open("{}".format(outfile), "w")

if python_output is True:
    if rank == 0:
        print("dt, temp, ke, fke, pe, fpe")

for j in range(20):

    # timestep
    dt = 0.005*(j+1)

    if lammps_output is True:
        lmp = lammps()
    else:
        lmp = lammps(cmdargs=["-screen", "none", "-log", "none"])

    lmp.command("atom_style full")
    lmp.command("units metal")
    lmp.command("processors * * *")
    lmp.command("neighbor 1 bin")
    lmp.command("boundary p p p")

    if run_no is 0:
        lmp.command("read_data {}".format(infile))
    else:
        lmp.command("read_restart final_restart".format(run_no-1))

    if thermo is True:
        lmp.command("thermo_style custom time temp  pe ke press vol  cpu")
        lmp.command("thermo {}".format(nthermo))
        lmp.command("thermo_modify flush yes")

    lmp.file("{}".format(ff_file))
    lmp.command("timestep {}".format(dt))

    # get_per_atom_compute example with dim of two and within a group
    # lmp.command("region rand block 5 20 5 20 5 20")
    # lmp.command("group rand region rand")
    # lmp.command("compute x rand property/atom x y")
    # test = get_per_atom_compute(comm, lmp, "x", 2, group="rand")

    lmp.command("compute ke all ke/atom")

    lmp.command("compute pe all pe")

    if ntraj != 0:
        lmp.command("dump 1 all dcd {} trajectory.dcd".format(ntraj))
        lmp.command("dump_modify 1 unwrap yes")

    if run_no == 0:
        lmp.command("velocity all create {} {} mom yes dist gaussian".format(temp_start, seed0))
    lmp.command("fix nve all nve")

    if ensemble == 0:
        # gjf u
        lmp.command("fix lang all langevin {} {} {} {} gjf yes halfstep yes".format(
            temp_start, temp_final, trel, seed1))
    elif ensemble == 1:
        # gjf v
        lmp.command("fix lang all langevin {} {} {} {} gjf yes".format(
            temp_start, temp_final, trel, seed1))
    elif ensemble == 2:
        # NH
        lmp.command("fix nvt all nvt temp {} {} {}".format(
            temp_start, temp_final, trel))
    elif ensemble == 3:
        # lang
        lmp.command("fix lang all langevin {} {} {} {} tally yes zero yes".format(
            temp_start, temp_final, trel, seed1))
    elif ensemble == 4:
        # BDP
        lmp.command("fix stoch all temp/csvr {} {} {} {}".format(
            temp_start, temp_final, trel, seed1))

    natoms = lmp.extract_global("natoms", 0)
    nlocal = lmp.extract_global("nlocal", 0)
    ke_sum = lt.get_per_atom_compute(comm, lmp, "ke")
    ke_2 = ke_sum**2
    pe_sum = 0
    pe_2 = 0
    temp_sum = 0

    for i in range(nout):
        nlocal = lmp.extract_global("nlocal", 0)
        lmp.command("run {} pre no post no".format(nthermo))
        temp = lmp.extract_compute("thermo_temp", 0, 0)
        ke = lt.get_per_atom_compute(comm, lmp, "ke")
        pe = lmp.extract_compute("pe", 0, 0)
        ke_sum += ke
        ke_2 += ke**2
        pe_sum += pe
        pe_2 += pe**2
        temp_sum += temp

        if python_output is True:
            if rank == 0:
                print("Time: {:.6f}, Temp: {:.6f}, KE: {:.6f}, PE: {:.6f}".format(
                    i*nthermo*dt, temp, ke.sum(), pe))

    if write_final_restart is True:
        lmp.command("write_restart {}".format(restart_file))

    if rank == 0:
        ke = ke_sum.sum() / (nout + 1)
        fke = (np.sqrt((ke_2 - ke_sum ** 2 / (nout + 1)) / (nout + 1))).sum()
        pe = pe_sum / nout
        fpe = np.sqrt((pe_2 - pe_sum ** 2 / nout) / nout)
        temp = temp_sum / nout

        if python_output is True:
            print(dt, temp, ke, fke, pe, fpe)

        if write_output is True:
            data.write("{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}\n".format(
                dt, temp, ke, fke, pe, fpe))
            data.flush()

if write_output is True:
    data.close()
