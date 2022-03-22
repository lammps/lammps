# MDI driver to perform a series of independent calculations
# using LAMMPS as an engine

# Syntax: python3 series_driver.py switch arg switch arg ...
#   possible switches:
#   -mdi "-role DRIVER ..."
#     required switch
#   -n 10
#     number of calculations to perform, default = 1
#   -mode eval/run/min
#     style of calculations: single snapshot evals, dynamics, minimization
#     default = eval
#   -size Nx Ny Nz
#     cubic lattice, default = 2 2 2
#   -rho 0.75 0.1
#     reduced density and random variation thereof, default = 0.75 0.1
#   -delta 0.1
#     randomly perturb atoms initially by this distance, default 0.0
#   -nsteps 100
#     number of timesteps in dynamics runs, default = 100
#   -temp 1.0
#     initial temperature in dynamics runs, default = 1.0
#   -tol 0.001
#     tolerance for minimizations, default = 0.001
#   -seed 12345
#     random number seed > 0, default = 12345

import sys,math,random
import mdi
import numpy as np
from mpi4py import MPI

# error message

def error(txt=None):
  if txt: raise Exception(txt)
  raise Exception("Syntax: python3 series_driver.py switch arg switch arg ...")

# send a LAMMPS input script command to MDI engine

def send_command(cmd):
  mdi.MDI_Send_Command("NBYTES",mdicomm)
  mdi.MDI_Send(len(cmd),1,mdi.MDI_INT,mdicomm)
  mdi.MDI_Send_Command("COMMAND",mdicomm)
  mdi.MDI_Send(cmd,len(cmd)+1,mdi.MDI_CHAR,mdicomm)

# parse command-line args

args = sys.argv[1:]
narg = len(args)

mdiarg = 0
ncalc = 1
mode = "eval"
nx = ny = nz = 2
rho = 0.75
rhodelta = 0.1
delta = 0.0
nsteps = 100
tinitial = 1.0
tol = 0.001
seed = 12345

iarg = 0
while iarg < narg:
  if args[iarg] == "-mdi": 
    if iarg+2 > narg: error()
    mdiarg = iarg + 1
    iarg += 2
  elif args[iarg] == "-n":
    if iarg+2 > narg: error()
    ncalc = int(args[iarg+1])
    iarg += 2
  elif args[iarg] == "-mode":
    if iarg+2 > narg: error()
    mode = args[iarg+1]
    if mode != "eval" and mode != "run" and mode != "min": error()
    iarg += 2
  elif args[iarg] == "-size":
    if iarg+4 > narg: error()
    nx = int(args[iarg+1])
    ny = int(args[iarg+2])
    nz = int(args[iarg+3])
    if nx <= 0 or ny <= 0 or nz <= 0: error()
    iarg += 4
  elif args[iarg] == "-rho":
    if iarg+3 > narg: error()
    rho = float(args[iarg+1])
    rhodelta = float(args[iarg+2])
    if rho-rhodelta <= 0.0: error()
    iarg += 4
  elif args[iarg] == "-delta":
    if iarg+2 > narg: error()
    delta = float(args[iarg+1])
    if delta < 0.0: error()
    iarg += 2
  elif args[iarg] == "-nsteps":
    if iarg+2 > narg: error()
    nsteps = int(args[iarg+1])
    if nsteps < 0: error()
    iarg += 2 
  elif args[iarg] == "-temp":
    if iarg+2 > narg: error()
    tinitial = float(args[iarg+1])
    if tinitial < 0.0: error()
    iarg += 2
  elif args[iarg] == "-tol":
    if iarg+2 > narg: error()
    tol = float(args[iarg+1]) 
    if tol < 0.0: error()
    iarg += 2
  elif args[iarg] == "-seed":
    if iarg+2 > narg: error()
    seed = int(args[iarg+1]) 
    if seed <= 0: error()
    iarg += 2
  else: error()

if not mdiarg: error()

# initialize MDI Library

mdi.MDI_Init(args[mdiarg])

# MPI communicator for just the driver

world = mdi.MDI_MPI_get_world_comm()
me = world.Get_rank()
nprocs = world.Get_size()

# connect to engine

mdicomm = mdi.MDI_Accept_Communicator()

# allocate vectors for per-atom types, coords, vels, forces

natoms = nx * ny * nz
atypes = np.zeros(natoms,dtype=np.int)
coords = np.zeros(3*natoms,dtype=np.float64)
vels = np.zeros(3*natoms,dtype=np.float64)
forces = np.zeros(3*natoms,dtype=np.float64)

atypes[:] = 1

# initialize RN generator

random.seed(seed)

# loop over sequence of calculations

for icalc in range(ncalc):

  # define simulation box

  onerho = rho + (random.random()-0.5)*rhodelta;
  sigma = pow(1.0/onerho,1.0/3.0)
  
  xlo = ylo = zlo = 0.0
  xhi = nx * sigma
  yhi = ny * sigma
  zhi = nz * sigma
  
  # send simulation box to engine

  vec = [xhi-xlo,0.0,0.0] + [0.0,yhi-ylo,0.0] + [0.0,0.0,zhi-zlo]
  mdi.MDI_Send_command(">CELL",mdicomm)
  mdi.MDI_Send(vec,9,mdi.MDI_DOUBLE,mdicomm)
  
  # create atoms on perfect lattice

  m = 0
  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        coords[m] = i * sigma
        coords[m+1] = j * sigma
        coords[m+2] = k * sigma
        m += 3
        
  # perturb lattice

  for m in range(3*natoms):
    coords[m] += 2.0*random.random()*delta - delta

  # define initial velocities

  for m in range(3*natoms):
    vels[m] = random.random() - 0.5

  tcurrent = 0.0
  for m in range(3*natoms):
    tcurrent += vels[m]*vels[m]
  tcurrent /= 3*(natoms-1)
    
  factor = math.sqrt(tinitial/tcurrent)

  for m in range(3*natoms):
    vels[m] *= factor

  # send atoms and their properties to engine
  
  mdi.MDI_Send_command(">NATOMS",mdicomm)
  mdi.MDI_Send(natoms,1,mdi.MDI_INT,mdicomm)
  mdi.MDI_Send_command(">TYPES",mdicomm)
  mdi.MDI_Send(atypes,natoms,mdi.MDI_INT,mdicomm)
  mdi.MDI_Send_command(">COORDS",mdicomm)
  mdi.MDI_Send(coords,3*natoms,mdi.MDI_DOUBLE,mdicomm)
  mdi.MDI_Send_command(">VELOCITIES",mdicomm)
  mdi.MDI_Send(vels,3*natoms,mdi.MDI_DOUBLE,mdicomm)

  # eval or run or minimize

  if mode == "eval":
    pass
  elif mode == "run":
    print("PRE MD")
    mdi.MDI_Send_command("@INIT_MD",mdicomm)
    mdi.MDI_Send_command(">NITERATE",mdicomm)
    mdi.MDI_Send(nsteps,1,mdi.MDI_INT,mdicomm)
    mdi.MDI_Send_command("@DEFAULT",mdicomm)
    print("POST MD")
  elif mode == "min":
    mdi.MDI_Send_command(">TOLERANCE",mdicomm)
    params = [1.0e-4,1.0e-4,1000.0,1000.0]
    mdi.MDI_Send(params,4,mdi.MDI_DOUBLE,mdicomm)

  # request potential energy

  print("PRE PE")

  mdi.MDI_Send_command("<PE",mdicomm)
  pe = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,mdicomm)
  pe = world.bcast(pe,root=0)

  # request virial tensor

  mdi.MDI_Send_command("<STRESS",mdicomm)
  virial = mdi.MDI_Recv(6,mdi.MDI_DOUBLE,mdicomm)
  virial = world.bcast(virial,root=0)

  # request forces
  
  mdi.MDI_Send_command("<FORCES",mdicomm)
  mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,mdicomm,buf=forces)
  world.Bcast(forces,root=0)

  # final output from each calculation
  # pressure = just virial component, no kinetic component
  
  aveeng = pe/natoms
  pressure = (virial[0] + virial[1] + virial[2]) / 3.0

  m = 0
  fx = fy = fz = 0.0
  for i in range(natoms):
    fx += forces[m]
    fy += forces[m+1]
    fz += forces[m+2]
    m += 3

  fx /= natoms
  fy /= natoms
  fz /= natoms

  line = "Calc %d: eng %7.5g pressure %7.5g aveForce %7.5g %7.5g %7.5g" % \
         (icalc+1,aveeng,pressure,fx,fy,fz)
  if me == 0: print(line)
  
# send EXIT command to engine

mdi.MDI_Send_command("EXIT",mdicomm)
