# MDI driver to perform an AIMD simulation
# using one instance of LAMMPS as the MD timestepper
# using second instance of LAMMPS as a QM surrogate to compute forces

# NOTE: this script is derived from the MDI_AIMD_Driver.cpp code
#       included in the MDI distribution
#  it alters the timestepping to match a velocity Verlet algorithm
#    forces are computed once before timestepping beings
#  both the @COORDS and @FORCES nodes are triggered in the MM code
#    as the appropriate places to extract COORDS and provide FORCES

# Syntax: python3 aimd_driver.py switch arg switch arg ...
#   possible switches:
#   -mdi "-role DRIVER ..."
#     required switch
#     example for stand-alone mode:
#       -mdi "-role DRIVER -name sequence -method TCP -port 8021"
#     example for plugin mode:
#       -mdi "-role DRIVER -name sequemce -method LINK
#             -plugin_path /home/sjplimp/lammps/src/"
#   -plugin name
#     name of plugin library, only when using plugin mode
#   -plugin_args arglist
#     args to add when launching plugin library, only when using plugin mode
#     enclose arglist in quotes if multiple words
#   -nsteps 10
#     number of timesteps, default = 10

import sys
import mdi
import numpy as np
from mpi4py import MPI

# error message

def error(txt=None):
  if txt: raise Exception(txt)
  raise Exception("Syntax: python3 aimd_driver.py switch arg switch arg ...")

# run an AIMD simulation

def perform_aimd(world,mm_comm,qm_comm):

  me = world.Get_rank()

  # receive number of atoms from the MM engine

  mdi.MDI_Send_command("<NATOMS",mm_comm)
  natoms = mdi.MDI_Recv(1,mdi.MDI_INT,mm_comm)
  natoms = world.bcast(natoms,root=0)

  # allocate arrays for coordinates and forces

  coords = np.zeros(3*natoms,dtype=np.float64)
  forces = np.zeros(3*natoms,dtype=np.float64)

  # MM engine initializes a new MD simulation
  
  mdi.MDI_Send_command("@INIT_MD",mm_comm)
  
  # -----------------
  # compute initial forces for Verlet timestepping
  #   and initial energy for output on step 0
  # -----------------
  
  # MM engine proceeds to @FORCES node in setup()
  
  mdi.MDI_Send_command("@FORCES",mm_comm)

  # get coords from MM engine
    
  mdi.MDI_Send_command("<COORDS",mm_comm)
  mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,mm_comm,buf=coords)
  world.Bcast(coords,root=0)

  # send coords to QM engine
  
  mdi.MDI_Send_command(">COORDS",qm_comm)
  mdi.MDI_Send(coords,3*natoms,mdi.MDI_DOUBLE,qm_comm)

  # get QM potential energy
    
  mdi.MDI_Send_command("<PE",qm_comm)
  qm_pe = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,qm_comm)
  qm_pe = world.bcast(qm_pe,root=0)
  
  # get forces from QM engine
    
  mdi.MDI_Send_command("<FORCES",qm_comm)
  mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,qm_comm,buf=forces)
  world.Bcast(forces,root=0)

  # send forces to MM engine
  
  mdi.MDI_Send_command(">FORCES",mm_comm)
  mdi.MDI_Send(forces,3*natoms,mdi.MDI_DOUBLE,mm_comm)

  # get MM kinetic energy
  
  mdi.MDI_Send_command("<KE",mm_comm)
  mm_ke = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,mm_comm)
  mm_ke = world.bcast(mm_ke,root=0)

  # output by driver
  # normalize energies by atom count
    
  if me == 0:
    print("Step %d: MM energy %g, QM energy %g, Total energy %g" % \
          (0,mm_ke/natoms,qm_pe/natoms,(mm_ke+qm_pe)/natoms))

  # -----------------
  # timestepping loop
  # -----------------

  for istep in range(nsteps):

    # MM engine proceeds to @FORCES node
    
    mdi.MDI_Send_command("@FORCES",mm_comm)

    # get coords from MM engine
    
    mdi.MDI_Send_command("<COORDS",mm_comm)
    mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,mm_comm,buf=coords)
    world.Bcast(coords,root=0)

    # send coords to QM engine
    
    mdi.MDI_Send_command(">COORDS",qm_comm)
    mdi.MDI_Send(coords,3*natoms,mdi.MDI_DOUBLE,qm_comm)

    # get QM potential energy
    
    mdi.MDI_Send_command("<PE",qm_comm)
    qm_pe = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,qm_comm)
    qm_pe = world.bcast(qm_pe,root=0)

    # get forces from QM engine
    
    mdi.MDI_Send_command("<FORCES",qm_comm)
    mdi.MDI_Recv(3*natoms,mdi.MDI_DOUBLE,qm_comm,buf=forces)
    world.Bcast(forces,root=0)

    # send forces to MM engine
    
    mdi.MDI_Send_command(">FORCES",mm_comm)
    mdi.MDI_Send(forces,3*natoms,mdi.MDI_DOUBLE,mm_comm)

    # MM engine proceeds to @ENDSTEP node
    # so that KE will be for fully updated velocity
    
    mdi.MDI_Send_command("@ENDSTEP",mm_comm)

    # get MM kinetic energy
    
    mdi.MDI_Send_command("<KE",mm_comm)
    mm_ke = mdi.MDI_Recv(1,mdi.MDI_DOUBLE,mm_comm)
    mm_ke = world.bcast(mm_ke,root=0)

    # output by driver
    # normalize energies by atom count
    
    if me == 0:
      print("Step %d: MM energy %g, QM energy %g, Total energy %g" % \
            (istep+1,mm_ke/natoms,qm_pe/natoms,(mm_ke+qm_pe)/natoms))

  #  send EXIT to each engine

  mdi.MDI_Send_command("EXIT",mm_comm)
  mdi.MDI_Send_command("EXIT",qm_comm)

# ------------------------
# main program
# ------------------------

args = sys.argv[1:]
narg = len(args)

# defaults for command-line args

mdiarg = ""
plugin = ""
plugin_args = ""

nsteps = 10

# parse command-line args

iarg = 0
while iarg < narg:
  if args[iarg] == "-mdi": 
    if iarg+2 > narg: error()
    mdiarg = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-plugin": 
    if iarg+2 > narg: error()
    plugin = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-plugin_args": 
    if iarg+2 > narg: error()
    plugin_args = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-nsteps":
    if iarg+2 > narg: error()
    nsteps = int(args[iarg+1])
    if nsteps < 0: error()
    iarg += 2 
  else: error()

if not mdiarg: error()

# LAMMPS engines are stand-alone codes
# world = MPI communicator for just this driver
# invoke perform_tasks() directly

if not plugin:
  mdi.MDI_Init(mdiarg)
  world = mdi.MDI_MPI_get_world_comm()
  
  # connect to 2 engines, determine which is MM vs QM

  mdicomm1 = mdi.MDI_Accept_Communicator()
  mdicomm2 = mdi.MDI_Accept_Communicator()

  mdi.MDI_Send_command("<NAME",mdicomm1)
  name1 = mdi.MDI_Recv(mdi.MDI_NAME_LENGTH,mdi.MDI_CHAR,mdicomm1)
  name1 = world.bcast(name1,root=0)
  mdi.MDI_Send_command("<NAME",mdicomm2)
  name2 = mdi.MDI_Recv(mdi.MDI_NAME_LENGTH,mdi.MDI_CHAR,mdicomm2)
  name2 = world.bcast(name2,root=0)

  if name1 == "MM" and name2 == "QM":
    mm_comm = mdicomm1
    qm_comm = mdicomm2
  elif name1 == "QM" and name2 == "MM":
    mm_comm = mdicomm2
    qm_comm = mdicomm1
  else: error("Two engines have invalid names")
    
  perform_aimd(world,mm_comm,qm_comm)

# LAMMPS engines are plugin libraries
# launch plugins
# NOTE: need to run driver on 2 or more procs
#       partition driver into 2 MPI comms
#       launch one plugin on each partition
#       each with their own callback function

if plugin:
  error("Cannot yet run in plugin mode")
  mdi.MDI_Init(mdiarg)
  world = MPI.COMM_WORLD
  plugin_args += " -mdi \"-role ENGINE -name lammps -method LINK\""
  mdi.MDI_Launch_plugin(plugin,plugin_args,world,perform_tasks,None)
