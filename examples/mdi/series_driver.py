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

  # delete all atoms so can run a new calculation

  send_command("delete_atoms group all")

  # define simulation box

  onerho = rho + (random.random()-0.5)*rhodelta;
  sigma = pow(1.0/onerho,1.0/3.0)
  
  xlo = ylo = zlo = 0.0
  xhi = nx * sigma
  yhi = ny * sigma
  zhi = nz * sigma
  
  # send simulation box to engine

  vec = [xlo,ylo,zlo,xhi,yhi,zhi,0.0,0.0,0.0]
  mdi.MDI_Send_command("RESET_BOX",mdicomm)
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
  
  mdi.MDI_Send_command("CREATE_ATOM",mdicomm)
  mdi.MDI_Send(natoms,1,mdi.MDI_INT,mdicomm)
  mdi.MDI_Send_command("CREATE_TYPE",mdicomm)
  mdi.MDI_Send(atypes,natoms,mdi.MDI_INT,mdicomm)
  mdi.MDI_Send_command("CREATE_X",mdicomm)
  mdi.MDI_Send(coords,3*natoms,mdi.MDI_DOUBLE,mdicomm)
  mdi.MDI_Send_command("CREATE_V",mdicomm)
  mdi.MDI_Send(vels,3*natoms,mdi.MDI_DOUBLE,mdicomm)
  mdi.MDI_Send_command("CREATE_GO",mdicomm)
  
  # eval or run or minimize

  if mode == "eval":
    mdi.MDI_Send_command("EVAL",mdicomm)
  elif mode == "run":
    send_command("run %d" % nsteps)
  elif mode == "min":
    send_command("minimize %g %g 1000 1000" % (tol,tol))

  # request energy

  mdi.MDI_Send_command("<ENERGY",mdicomm)
  energy = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,mdicomm)
  energy = world.bcast(energy,root=0)

  # request pressure tensor

  mdi.MDI_Send_command("<PTENSOR",mdicomm)
  ptensor = mdi.MDI_Recv(6,mdi.MDI_DOUBLE,mdicomm)
  ptensor = world.bcast(ptensor,root=0)

  # request forces
  
  mdi.MDI_Send_command("<FORCES",mdicomm)
  mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,mdicomm,buf=forces)
  world.Bcast(forces,root=0)

  # final output from each calculation

  aveeng = energy/natoms
  pressure = (ptensor[0] + ptensor[1] + ptensor[2]) / 3.0

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

  line = "Calc %d: eng %7.5g press %7.5g aveForce %7.5g %7.5g %7.5g" % \
         (icalc+1,aveeng,pressure,fx,fy,fz)
  if me == 0: print(line)
  
# send EXIT command to engine

mdi.MDI_Send_command("EXIT",mdicomm)
