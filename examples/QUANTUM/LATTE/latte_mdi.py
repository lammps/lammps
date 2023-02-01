# MDI wrapper on LATTE code

import sys,time
from ctypes import *

import numpy as np
from mpi4py import MPI
import MDI_Library as mdi

# conversions of atomic number to element symbol

atomic_number_to_symbol = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 17: 'Cl'}

# --------------------------------------------
# global data
# --------------------------------------------

world = 0
me = nprocs = 0

exitflag = False

AIMD = 0
QMMM = 1
mode = AIMD

# LATTE library

libname = "liblatte.so"
liblatte = None

# QM inputs

flag_qm_natoms = flag_mm_natoms = 0
flag_box = flag_box_displ = 0
flag_qm_elements = flag_qm_types = 0
flag_qm_coords = flag_qm_potential = 0
flag_mm_elements = 0
flag_mm_coords = flag_mm_charges = 0

box = np.empty(9)
box_displ = np.empty(3)

qm_natoms = 0
qm_elements = None
qm_coords = None
qm_potential = None

mm_natoms = 0
mm_coords = None
mm_charges = None
mm_elements = None

# QM outputs

qm_pe = 0.0
qm_stress = np.empty(9)
qm_forces = None
qm_charges = None

mm_forces = None

# --------------------------------------------
# print error message and halt
# --------------------------------------------

def error(txt):
  if me == 0: print("ERROR:",txt)
  world.Abort()
  
# --------------------------------------------
# process non-MDI command-line options to LATTE
# --------------------------------------------

def options(other_options):
  if len(other_options) != 0:
    error("No args currently used by LATTE wrapper")

# --------------------------------------------
# operate as an engine
# --------------------------------------------

def mdi_engine(other_options):
  global world

  # get the MPI intra-communicator for this engine

  world = mdi.MDI_MPI_get_world_comm()
  me = world.Get_rank()
  nprocs = world.Get_size()

  # process non-MDI command line args

  options(other_options)

  # confirm LATTE is being run as an engine
    
  role = mdi.MDI_Get_Role()
  if not role == mdi.MDI_ENGINE:
    error("Must run LATTE as an MDI engine")

  # supported MDI commands

  mdi.MDI_Register_Node("@DEFAULT")
  mdi.MDI_Register_Command("@DEFAULT","EXIT")

  # driver --> engine
  
  mdi.MDI_Register_Command("@DEFAULT",">NATOMS")
  mdi.MDI_Register_Command("@DEFAULT",">CELL")
  mdi.MDI_Register_Command("@DEFAULT",">CELL_DISPL")
  mdi.MDI_Register_Command("@DEFAULT",">ELEMENTS")
  mdi.MDI_Register_Command("@DEFAULT",">TYPES")
  mdi.MDI_Register_Command("@DEFAULT",">COORDS")
  mdi.MDI_Register_Command("@DEFAULT",">POTENTIAL_AT_NUCLEI")
  mdi.MDI_Register_Command("@DEFAULT",">NLATTICE")
  mdi.MDI_Register_Command("@DEFAULT",">LATTICE_ELEMENTS")
  mdi.MDI_Register_Command("@DEFAULT",">CLATTICE")
  mdi.MDI_Register_Command("@DEFAULT",">LATTICE")

  # engine --> driver

  mdi.MDI_Register_Command("@DEFAULT","<PE")
  mdi.MDI_Register_Command("@DEFAULT","<FORCES")
  mdi.MDI_Register_Command("@DEFAULT",">LATTICE_FORCES")
  mdi.MDI_Register_Command("@DEFAULT","<STRESS")
  mdi.MDI_Register_Command("@DEFAULT","<CHARGES")

  # load LATTE lib and set ctypes signatures for function calls

  latte_load()

  # one-time operation to establish a connection with the driver

  mdicomm = mdi.MDI_Accept_Communicator()

  # set callback to execute_command
  
  mdi.MDI_Set_Execute_Command_Func(execute_command,None)

  # infinite loop to exchange messages with driver
  # until EXIT command received

  while not exitflag:
    command = mdi.MDI_Recv_Command(mdicomm)
    command = world.bcast(command,root=0)
    execute_command(command,mdicomm,None)
    
# --------------------------------------------
# called in loop by mdi_engine()
# called internally from MDI_Recv_command() until an EXIT command is received
# command = command name
# mdicomm = MDI communicator for all MDI commands
#           can also use world for MPI commands
# object_ptr = ptr to data if necessary
# --------------------------------------------

def execute_command(command,mdicomm,object_ptr):
  global exitflag

  # driver says done
  
  if command == "EXIT":
    exitflag = True

  # MDI commands which setup quantum calculation
  # NOTE: not sure if it makes sense for LATTE wrapper to support >TYPES
  
  elif command == ">NATOMS":
    receive_qm_natoms(mdicomm)

  elif command == ">CELL":
    receive_box(mdicomm)
    
  elif command == ">CELL_DISPL":
    receive_box_displ(mdicomm)
  
  elif command == ">ELEMENTS":
    receive_qm_elements(mdicomm)
    
  elif command == ">COORDS":
    receive_qm_coords(mdicomm)

  elif command == ">POTENTIAL_AT_NUCLEI":
    receive_qm_potential(mdicomm)

  elif command == ">NLATTICE":
    receive_mm_natoms(mdicomm)
  
  elif command == ">LATTICE_ELEMENTS":
    receive_mm_elements(mdicomm)

  elif command == ">CLATTICE":
    receive_mm_coords(mdicomm)
  
  elif command == ">LATTICE":
    receive_mm_charges(mdicomm)

  # MDI commands which retreive quantum results
  # each may also trigger the quantum calculation
  
  elif command == "<PE":
    evaluate()
    ierr = mdi.MDI_Send(qm_pe,1,mdi.MDI_DOUBLE,mdicomm)
    if ierr: error("MDI: <PE data")
    
  elif command == "<FORCES":
    evaluate()
    ierr = mdi.MDI_Send(qm_forces,3*qm_natoms,mdi.MDI_DOUBLE,mdicomm)
    if ierr: error("MDI: <FORCES data")
    
  elif command == "<LATTICE_FORCES":
    evaluate()
    ierr = mdi.MDI_Send(mm_forces,3*mm_natoms,mdi.MDI_DOUBLE,mdicomm)
    if ierr: error("MDI: <LATTICE_FORCES data")
    
  elif command == "<STRESS":
    evaluate()
    ierr = mdi.MDI_Send(qm_stress,1,mdi.MDI_DOUBLE,mdicomm)
    if ierr: error("MDI: <STRESS data")

  elif command == "<CHARGES":
    evaluate()
    ierr = mdi.MDI_Send(qm_charges,qm_natoms,mdi.MDI_DOUBLE,mdicomm)
    if ierr: error("MDI: <CHARGES data")

  # unrecognized command
  
  else: error("Unrecognized MDI command")

  return 0

# --------------------------------------------
# receive count of QM atoms from driver
# --------------------------------------------

def receive_qm_natoms(mdicomm):
  global flag_qm_natoms,qm_natoms
  flag_qm_natoms = 1
  
  qm_natoms = mdi.MDI_Recv(1,mdi.MDI_INT,mdicomm)
  qm_natoms = world.bcast(qm_natoms,root=0)
  allocate("qm")
    
# --------------------------------------------
# receive 3 simulation box edge vectors from driver
# --------------------------------------------

def receive_box(mdicomm):
  global flag_box
  flag_box = 1
  
  ierr = mdi.MDI_Recv(9,mdi.MDI_DOUBLE,mdicomm,buf=box)
  if ierr: error("MDI: >CELL data")
  world.Bcast(box,root=0)

# --------------------------------------------
# receive simulation box displacement vector from driver
# --------------------------------------------

def receive_box_displ(mdicomm):
  global flag_box_displ
  flag_box_displ = 1
  
  ierr = mdi.MDI_Recv(3,mdi.MDI_DOUBLE,mdicomm,buf=box_displ)
  if ierr: error("MDI: >CELL_DISPL data")
  world.Bcast(box_displ,root=0)

# --------------------------------------------
# receive QM atom coords from driver
# --------------------------------------------

def receive_qm_coords(mdicomm):
  global flag_qm_coords
  flag_qm_coords = 1

  if not qm_natoms: error("Cannot MDI >COORDS if # of QM atoms = 0")
  
  ierr = mdi.MDI_Recv(3*qm_natoms,mdi.MDI_DOUBLE,mdicomm,buf=qm_coords)
  if ierr: error("MDI: >COORDS data")
  world.Bcast(qm_coords,root=0)

# --------------------------------------------
# receive Coulomb potential at QM nuclei from driver
# --------------------------------------------

def receive_qm_potential(mdicomm):
  global flag_qm_potential
  flag_qm_potential = 1

  if not qm_natoms: error("Cannot MDI >POTENTIAL_AT_NUCLEI if # of QM atoms = 0")

  ierr = mdi.MDI_Recv(qm_natoms,mdi.MDI_DOUBLE,mdicomm,buf=qm_potential)
  if ierr: error("MDI: >POTENTIAL_AT_NUCLEI data")
  world.Bcast(qm_potential,root=0)

# --------------------------------------------
# receive QM atomic numbers from driver
# --------------------------------------------

def receive_qm_elements(mdicomm):
  global flag_qm_elements
  flag_qm_elements = 1
   
  if not qm_natoms: error("Cannot MDI >ELEMENTS if # of QM atoms = 0")

  ierr = mdi.MDI_Recv(qm_natoms,mdi.MDI_INT,mdicomm,buf=qm_elements)
  if ierr: error("MDI: >ELEMENTS data")
  world.Bcast(qm_elements,root=0)

# --------------------------------------------
# receive QM atom types from driver
# --------------------------------------------

def receive_types(mdicomm):
  global flag_qm_types
  flag_qm_types = 1
  
  if not qm_natoms: error("Cannot MDI >TYPES if # of QM atoms = 0")

  ierr = mdi.MDI_Recv(qm_natoms,mdi.MDI_INT,mdicomm,buf=qm_types)
  if ierr: error("MDI: >TYPES data")
  world.Bcast(qm_types,root=0)

# --------------------------------------------
# receive count of MM atoms from driver
# --------------------------------------------

def receive_mm_natoms(mdicomm):
  global flag_mm_natoms,mm_natoms
  flag_mm_natoms = 1
  
  mm_natoms = mdi.MDI_Recv(1,mdi.MDI_INT,mdicomm)
  mm_natoms = world.bcast(mm_natoms,root=0)
  allocate("mm")

# --------------------------------------------
# receive MM atomic numbers from driver
# --------------------------------------------

def receive_mm_elements(mdicomm):
  global flag_mm_elements
  flag_mm_elements = 1
   
  if not mm_natoms: error("Cannot MDI >LATTICE_ELEMENTS if # of MM atoms = 0")

  ierr = mdi.MDI_Recv(mm_natoms,mdi.MDI_INT,mdicomm,buf=mm_elements)
  if ierr: error("MDI: >LATTICE_ELEMENTS data")
  world.Bcast(mm_elements,root=0)

# --------------------------------------------
# receive MM atom coords from driver
# --------------------------------------------

def receive_mm_coords(mdicomm):
  global flag_mm_coords
  flag_mm_coords = 1

  if not mm_natoms: error("Cannot MDI >CLATTICE if # of MM atoms = 0")
  
  ierr = mdi.MDI_Recv(3*mm_natoms,mdi.MDI_DOUBLE,mdicomm,buf=mm_coords)
  if ierr: error("MDI: >CLATTICE data")
  world.Bcast(mm_coords,root=0)

# --------------------------------------------
# receive charge on MM atoms from driver
# --------------------------------------------

def receive_mm_charges(mdicomm):
  global flag_mm_charges
  flag_mm_charges = 1

  if not mm_natoms: error("Cannot MDI >LATTICE if # of MM atoms = 0")

  ierr = mdi.MDI_Recv(mm_natoms,mdi.MDI_DOUBLE,mdicomm,buf=mm_charges)
  if ierr: error("MDI: >LATTICE data")
  world.Bcast(mm_charges,root=0)

# --------------------------------------------
# allocate persistent data for QM or MM atoms
# called when qm_natoms or mm_natoms is reset by MDI driver
# --------------------------------------------

def allocate(which):
  global qm_types,qm_elements,qm_coords,qm_potential,qm_forces,qm_charges
  global mm_coords,mm_charges,mm_forces
  
  if which == "qm":
    n = qm_natoms
    qm_types = np.empty(n,dtype=np.int32)
    qm_elements = np.empty(n,dtype=np.int32)
    qm_coords = np.empty((n,3))
    qm_potential = np.empty(n)
    qm_forces = np.empty((n,3))
    qm_charges = np.empty(n)

  if which == "mm":
    n = mm_natoms
    mm_elements = np.empty(n,dtype=np.int32)
    mm_coords = np.empty((n,3))
    mm_charges = np.empty(n)
    mm_forces = np.empty((n,3))

# --------------------------------------------
# perform a quantum calculation via LATTE
# NOTE: ignore change of box size each step for now, e.g. MD NPT dynamics
# NOTE: assume periodic for now, worry about non-periodic later
# --------------------------------------------
  
def evaluate():
  global mode
  global flag_qm_natoms,flag_mm_natoms
  global flag_box,flag_box_displ
  global flag_qm_elements,flag_qm_types,flag_qm_coords,flag_qm_potential
  global flag_mm_elements,flag_mm_coords,flag_mm_charges
  global qm_pe,qm_stress,qm_forces,qm_charges
  global mm_forces
  global dm_previous

  # just return if the QM system was already evaluated
  # happens when multiple results are requested by driver
  
  any_flag = flag_qm_natoms + flag_mm_natoms + flag_box + flag_box_displ + \
    flag_qm_elements + flag_qm_coords + flag_qm_potential  + \
    flag_mm_elements + flag_mm_coords + flag_mm_charges
  if not any_flag: return

  # if any of these MDI commands received from LAMMPS,
  #   treat it as a brand new system

  new_system = 0
  if flag_qm_natoms or flag_mm_natoms: new_system = 1
  if flag_qm_elements or flag_mm_elements: new_system = 1
  if new_system:
    if flag_mm_natoms or flag_qm_potential: mode = QMMM
    else: mode = AIMD

  # if new system, error check that all needed MDI calls have been made

  if new_system:
    if not flag_qm_natoms: error("QM atom count not specified")
    if not flag_qm_elements or not flag_qm_coords:
      error("QM atom properties not fully specified")
    if mode == QMMM and not flag_qm_potential:
      error("QM atom properties not fully specified")

  # hardwire these unsupported flags for now

  coulombflag = 0
  neighflag = 0
  pbcflag = 1      # NOTE: pass this in as latte_mdi.py command-line arg
  thermo_virial = 1
  eflag_atom = 1
  vflag_global = 1
  vflag_atom = 0
  
  flags_latte = 6*[0]
  flags_latte[0] = pbcflag;   # 1 for fully periodic, 0 for fully non-periodic
  flags_latte[1] = coulombflag;  # 1 for LAMMPS computes Coulombics, 0 for LATTE
  flags_latte[2] = eflag_atom;   # 1 to return per-atom energies, 0 for no
  flags_latte[3] = vflag_global and thermo_virial; # 1 to return global/per-atom
  flags_latte[4] = vflag_atom and thermo_virial;   #   virial, 0 for no
  flags_latte[5] = neighflag;    # 1 to pass neighbor list to LATTE, 0 for no

  boxlo = [0.0,0.0,0.0]         # NOTE: does this matter for LATTE ?
  boxhi = [box[0],box[4],box[8]]
  xy = box[3]
  xz = box[6]
  yz = box[7]
  maxiter = -1

  # QMMM with QM and MM atoms
  # NOTE: need qm_velocity and timestep and mass and types ?
  # all of these are addresses of scalars for Fortran ?

  if mode == QMMM:
    error("QMMM not yet supported with LATTE")
    
  # AIMD with only QM atoms
    
  elif mode == AIMD:

    c_flags_latte = (c_int*6)(*flags_latte)
    c_qm_natoms = c_int(qm_natoms)
    qm_ntypes = 2
    c_qm_ntypes = c_int(qm_ntypes)
    c_xy = c_double(xy)
    c_xz = c_double(xz)
    c_yz = c_double(yz)
    c_maxiter = c_int(maxiter)
    c_qm_pe = c_double(qm_pe)
    c_new_system = c_int(new_system)

    qm_types = np.empty(qm_natoms,dtype=np.int32)
    for i in range(0,qm_natoms,3): qm_types[i] = 1
    for i in range(1,qm_natoms,3): qm_types[i] = 2
    for i in range(2,qm_natoms,3): qm_types[i] = 2
    qm_mass = [15.995, 1.008]
    c_qm_mass = (c_double*qm_ntypes)(*qm_mass)

    c_boxlo = (c_double*3)(*boxlo)
    c_boxhi = (c_double*3)(*boxhi)
    
    qm_velocity = np.empty((qm_natoms,3))
    qm_velocity.fill(0.0)

    timestep = 0.00025
    c_timestep = c_double(timestep)
    
    latte_error = 0
    c_latte_error = c_bool(latte_error)

    print("Calling LATTE ...")
    time1 = time.time()

    print("flags_latte",c_flags_latte[0:6])
    print("qm_natoms",c_qm_natoms.value)
    print("qm_coords",qm_coords)  
    print("qm_types",qm_types)
    print("qm_ntypes",c_qm_ntypes.value)
    print("qm_mass",c_qm_mass[0:2])
    print("boxlo",c_boxlo[0:3])
    print("boxhi",c_boxhi[0:3]) 
    print("xy",c_xy.value)
    print("xz",c_xz.value)
    print("yz",c_yz.value)
    print("maxiter",c_maxiter.value)
    print("timestep",c_timestep.value)
    print("new_system",c_new_system.value)
 
    liblatte.\
      latte(c_flags_latte,byref(c_qm_natoms),qm_coords,
            qm_types,byref(c_qm_ntypes),c_qm_mass,
            c_boxlo,c_boxhi,byref(c_xy),byref(c_xz),byref(c_yz),qm_forces,
            byref(c_maxiter),byref(c_qm_pe),
            qm_velocity,byref(c_timestep),qm_stress,
            byref(c_new_system),byref(c_latte_error))
    # NOTE: check latte_error return?
    latte_error = c_latte_error.value
    qm_pe = c_qm_pe.value

    time2 = time.time()
    print("DONE LATTE",latte_error,time2-time1)
    print("PE",qm_pe)
    print("FORCE",qm_forces)

  # clear flags for all MDI commands for next QM evaluation

  flag_qm_natoms = flag_mm_natoms = 0
  flag_box = flag_box_displ = 0
  flag_qm_elements = 0
  flag_qm_coords = flag_qm_potential = 0
  flag_mm_elements = 0
  flag_mm_coords = flag_mm_charges = 0

# --------------------------------------------
# load LATTE library
# set ctypes signatures for single function calls to LATTE lib
# --------------------------------------------

def latte_load():
  global liblatte

  liblatte = CDLL(libname,RTLD_GLOBAL)

  nparray = np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags="C_CONTIGUOUS")
  npvector_double = np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags="C_CONTIGUOUS")
  npvector_int = np.ctypeslib.ndpointer(dtype=np.int32,ndim=1,flags="C_CONTIGUOUS")

  liblatte.latte_abiversion.restype = None
  liblatte.latte_abiversion.argtypes = None

  liblatte.latte.restype = None
  liblatte.latte.argtypes = \
    [POINTER(c_int), POINTER(c_int), nparray, npvector_int, POINTER(c_int), POINTER(c_double),
     POINTER(c_double), POINTER(c_double), POINTER(c_double),
     POINTER(c_double), POINTER(c_double), 
     nparray, POINTER(c_int), POINTER(c_double), nparray,
     POINTER(c_double), npvector_double, POINTER(c_int), POINTER(c_bool)]

# --------------------------------------------
# function called by MDI driver
# only when it invokes pyscf_mdi.py as a plugin
# --------------------------------------------

def MDI_Plugin_init_latte_mdi(plugin_state):
  
  # other_options = all non-MDI args
  # -mdi arg is processed and stripped internally by MDI
  
  other_options = []

  mdi.MDI_Set_plugin_state(plugin_state)
  narg = mdi.MDI_Plugin_get_argc()

  for iarg in range(narg):
    arg = mdi.MDI_Plugin_get_arg(iarg)
    other_options.append(arg)

  # start running as an MDI engine
  
  mdi_engine(other_options)

# --------------------------------------------
# main program
# invoked when MDI driver and pyscf_mdi.py
#   are run as independent executables
# --------------------------------------------

if __name__== "__main__":
  
  # mdi_option = single arg in quotes that follows -mdi
  # other_options = all non-MDI args

  mdi_option = ""
  other_options = []

  narg = len(sys.argv)
  args = sys.argv
  
  iarg = 1
  while iarg < narg:
    arg = args[iarg]
    if arg == "-mdi" or arg == "--mdi":
      if narg > iarg+1: mdi_option = sys.argv[iarg+1]
      else: error("LATTE -mdi argument not provided")
      iarg += 1
    else: other_options.append(arg)
    iarg += 1

  if not mdi_option: error("LATTE -mdi option not provided")

  # call MDI_Init with just -mdi option
  
  mdi.MDI_Init(mdi_option)

  # start running as an MDI engine
  
  mdi_engine(other_options)
