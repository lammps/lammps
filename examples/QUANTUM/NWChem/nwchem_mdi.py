# MDI wrapper on NWChem PWDFT code

# native PWDFT units are Bohr and Hartree
# but box and atom coord inputs in *.nw file are in Angstroms

import sys,os,time
from ctypes import *

import numpy as np
from mpi4py import MPI
import mdi as mdi

# --------------------------------------------

ELEMENTS = [
  'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
  'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
  'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
  'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
  'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
  'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
  'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
  'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
]

# atomic_number_to_symbol converts atomic number to element symbol

atomic_number_to_symbol = {}
for i,symbol in enumerate(ELEMENTS):
  atomic_number_to_symbol[i+1] = symbol

# --------------------------------------------
# global data
# --------------------------------------------

world = 0
me = nprocs = 0

if MPI._sizeof(MPI.Comm) == sizeof(c_int):
  MPI_Comm = c_int
else:
  MPI_Comm = c_void_p

exitflag = False

AIMD = 0
QMMM = 1
mode = AIMD

# NWChem PWDFT library

libname = "libpwdft.so"
libpwdft = None

# QM inputs

nw_template = ""
nw_infile = ""
nw_outfile = ""

flag_qm_natoms = flag_mm_natoms = 0
flag_box = flag_box_displ = 0
flag_qm_elements = 0
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

def error(txt,mpiexists=1):
  if me == 0: print("ERROR:",txt)
  if mpiexists: world.Abort()
  sys.exit()
  
# --------------------------------------------
# process non-MDI options to PWDFT
# if this script is executed independently:
#   args = command-line args
# if this script is invoked as a plugin library:
#   args = passed via MDI
# --------------------------------------------

def options(other_options):
  global nw_template,nw_infile,nw_outfile
  if len(other_options) != 3:
    error("Invalid args to NWChem wrapper: template_file infile outfile")
  nw_template = other_options[0]
  nw_infile = other_options[1]
  nw_outfile = other_options[2]

# --------------------------------------------
# operate as an engine
# --------------------------------------------

def mdi_engine(other_options):
  global world,me,nprocs,MPI_Comm,libpwdft

  # get the MPI intra-communicator for this engine

  world = mdi.MDI_MPI_get_world_comm()
  me = world.Get_rank()
  nprocs = world.Get_size()

  # process non-MDI command line args

  options(other_options)

  # confirm PWDFT is being run as an engine
    
  role = mdi.MDI_Get_Role()
  if not role == mdi.MDI_ENGINE:
    error("Must run NWChem as an MDI engine")

  # supported MDI commands

  mdi.MDI_Register_Node("@DEFAULT")
  mdi.MDI_Register_Command("@DEFAULT","EXIT")

  # driver --> engine
  
  mdi.MDI_Register_Command("@DEFAULT",">NATOMS")
  mdi.MDI_Register_Command("@DEFAULT",">CELL")
  mdi.MDI_Register_Command("@DEFAULT",">CELL_DISPL")
  mdi.MDI_Register_Command("@DEFAULT",">ELEMENTS")
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
  #mdi.MDI_Register_Command("@DEFAULT","<STRESS")
  mdi.MDI_Register_Command("@DEFAULT","<CHARGES")

  # load PWDFT lib and set ctypes signatures for function calls

  pwdft_load()
  
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
# called internally from MDI_Recv_command() until EXIT received
# command = name of MDI command
# mdicomm = MDI communicator for all MDI commands
# object_ptr = ptr to data if necessary
# --------------------------------------------

def execute_command(command,mdicomm,object_ptr):
  global exitflag

  # driver says done
  
  if command == "EXIT":
    exitflag = True

  # MDI commands which setup quantum calculation
  
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
  global qm_elements,qm_coords,qm_potential,qm_forces,qm_charges
  global mm_elements,mm_coords,mm_charges,mm_forces
  
  if which == "qm":
    n = qm_natoms
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
# perform a quantum calculation via NWChem PWDFT
# --------------------------------------------
  
def evaluate():
  global mode
  global flag_qm_natoms,flag_mm_natoms
  global flag_box,flag_box_displ
  global flag_qm_elements,flag_qm_coords,flag_qm_potential
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

  # setup new system within PWDFT
  # remove nwchemex.movecs file if it exists
  
  if new_system:
    if me == 0:
      if os.path.exists("nwchemex.movecs"): os.remove("nwchemex.movecs")
    pwdft_initialize()

  # QMMM with QM and MM atoms

  world_ptr = MPI._addressof(world)
  c_world = MPI_Comm.from_address(world_ptr)
  c_qm_pe = c_double(qm_pe)
  c_nw_outfile = nw_outfile.encode()
  
  if mode == QMMM:
    #print("QMMM minimizer")
    #print("COORDS",qm_coords)
    #print("POTENTIAL",qm_potential)
    time1 = time.time()
    nwerr = libpwdft.\
      c_lammps_pspw_qmmm_minimizer_filename(c_world,qm_coords,qm_potential,
                                            qm_forces,qm_charges,byref(c_qm_pe),
                                            False,True,c_nw_outfile)
    # NOTE: check nwerr return?
    qm_pe = c_qm_pe.value
    time2 = time.time()
    if me == 0: print("Time for PWDFT qmmm:",time2-time1)
    #print("PE",qm_pe)
    #print("FORCE",qm_forces)
    #print("CHARGES",qm_charges)
    
  # AIMD with only QM atoms
    
  elif mode == AIMD:
    #print("AIMD minimizer")
    #print("COORDS",qm_coords)
    time1 = time.time()
    nwerr = libpwdft.\
      c_lammps_pspw_aimd_minimizer_filename(c_world,qm_coords,qm_forces,
                                            byref(c_qm_pe),c_nw_outfile) 
    # NOTE: check nwerr return?
    qm_pe = c_qm_pe.value
    time2 = time.time()
    if me == 0: print("Time for PWDFT aimd:",time2-time1)
    #print("PE",qm_pe)
    #print("FORCE",qm_forces)

  # clear flags for all MDI commands for next QM evaluation

  flag_qm_natoms = flag_mm_natoms = 0
  flag_box = flag_box_displ = 0
  flag_qm_elements = 0
  flag_qm_coords = flag_qm_potential = 0
  flag_mm_elements = 0
  flag_mm_coords = flag_mm_charges = 0
  
# --------------------------------------------
# load PWDFT library
# set ctypes signatures for 3 function calls to PWDFT lib
# --------------------------------------------

def pwdft_load():
  global libpwdft
  
  libpwdft = CDLL(libname,RTLD_GLOBAL)

  libpwdft.c_lammps_pspw_input_filename.restype = None
  libpwdft.c_lammps_pspw_input_filename.argtypes = \
    [MPI_Comm, c_char_p, c_char_p]

  nparray = np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags="C_CONTIGUOUS")
  npvector = np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags="C_CONTIGUOUS")
  
  libpwdft.c_lammps_pspw_qmmm_minimizer_filename.restype = c_int
  libpwdft.c_lammps_pspw_qmmm_minimizer_filename.argtypes = \
    [MPI_Comm, nparray, npvector, nparray, npvector, POINTER(c_double),
     c_bool, c_bool, c_char_p]

  libpwdft.c_lammps_pspw_aimd_minimizer_filename.restype = c_int
  libpwdft.c_lammps_pspw_aimd_minimizer_filename.argtypes = \
    [MPI_Comm, nparray, nparray, POINTER(c_double), c_char_p]

# --------------------------------------------
# create PWDFT input file with box and list of atoms
# invoke PWDFT function to read it
# --------------------------------------------

def pwdft_initialize():

  # PWDFT input file:
  # box, qm_coords must be converted to Angstroms

  bohr_to_angstrom = mdi.MDI_Conversion_factor("bohr","angstrom")

  box_A = box * bohr_to_angstrom
  qm_coords_A = qm_coords * bohr_to_angstrom

  # proc 0 reads template file, writes PWDFT input file

  if me == 0:
    lines = open(nw_template,'r').readlines()

    fp = open(nw_infile,'w')

    for line in lines:
      word = line.strip()
      if word == "GEOMINSERT":
        print("geometry noautosym noautoz nocenter",file=fp);
        print("system crystal cartesian",file=fp)
        print("lattice_vectors",file=fp)
        print("%20.16g %20.16g %20.16g" % (box_A[0],box_A[1],box_A[2]),file=fp)
        print("%20.16g %20.16g %20.16g" % (box_A[3],box_A[4],box_A[5]),file=fp)
        print("%20.16g %20.16g %20.16g" % (box_A[6],box_A[7],box_A[8]),file=fp)
        print("end\n",file=fp)

        for i in range(qm_natoms):
          symbol = atomic_number_to_symbol[qm_elements[i]]
          print("%s %20.16g %20.16g %20.16g" %
                (symbol,qm_coords_A[i][0],qm_coords_A[i][1],qm_coords_A[i][2]),
                file=fp)
        print("end\n",file=fp)

      else: print(line,file=fp,end="")

    fp.close()

  # all procs call pspw_input_filename() which processes input file
  # performs initial QM calculation within PWDFT

  world_ptr = MPI._addressof(world)
  c_world = MPI_Comm.from_address(world_ptr)
  infile = nw_infile.encode()
  outfile = nw_outfile.encode()

  #print("INPUT filename")
  time1 = time.time()
  nwerr = libpwdft.c_lammps_pspw_input_filename(c_world,infile,outfile)
  time2 = time.time()
  if me == 0: print("Time for PWDFT setup:",time2-time1)

# --------------------------------------------
# function called by MDI driver
# only when it invokes pyscf_mdi.py as a plugin
# --------------------------------------------

def MDI_Plugin_init_nwchem_mdi(plugin_state):
  
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
      else: error("NWChem -mdi argument not provided",0)
      iarg += 1
    else: other_options.append(arg)
    iarg += 1

  if not mdi_option: error("NWChem -mdi option not provided",0)

  # call MDI_Init with just -mdi option
  
  mdi.MDI_Init(mdi_option)

  # start running as an MDI engine
  
  mdi_engine(other_options)
