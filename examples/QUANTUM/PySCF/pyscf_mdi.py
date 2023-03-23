# MDI wrapper on PySCF quantum code

# native PySCF units are Bohr and Hartree
# but box and atom coord inputs are passed in Angstroms

import sys,time

import numpy as np
from mpi4py import MPI
import mdi

from pyscf.gto import Mole
from pyscf.pbc.gto import Cell
from pyscf import qmmm
from pyscf.dft import RKS as RKS_nonpbc
from pyscf.pbc.dft import RKS as RKS_pbc

# --------------------------------------------

# atomic_number_to_radius converts atomic number to radius (Angstroms)
# Chemistry - A European Journal, (2009), 186-197, 15(1)

atomic_number_to_radius = {1: 0.32, 6: 0.75, 7: 0.71, 8: 0.63, 17: 0.99}

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
# PySCF settings
# these are default values
# options() may override them
# --------------------------------------------

periodic = 1
xcstr = "wb97x"
basis = "6-31+G**"

# --------------------------------------------
# global data
# --------------------------------------------

world = 0
me = nprocs = 0
exitflag = False

AIMD = 0
QMMM = 1
mode = AIMD

# QM inputs

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
mm_radii = None

# QM outputs

qm_pe = 0.0
qm_stress = np.empty(9)
qm_forces = None
qm_charges = None

mm_forces = None

# PySCF internal data to persist state from one step to next

dm_previous_exists = 0
dm_previous = None

# --------------------------------------------
# print error message and halt
# --------------------------------------------

def error(txt,mpiexists=1):
  if me == 0: print("ERROR:",txt)
  if mpiexists: world.Abort()
  sys.exit()

# --------------------------------------------
# process non-MDI options to PySCF
# if this script is executed independently:
#   args = command-line args
# if this script is invoked as a plugin library:
#   args = passed via MDI
# --------------------------------------------

def options(args):
  global periodic,xcstr,basis

  narg = len(args)
  iarg = 0
  while iarg < narg:
    if args[iarg] == "-pbc":
      if iarg+1 > narg: error("Invalid PySCF command line args")
      if args[iarg+1] == "yes": periodic = 1
      elif args[iarg+1] == "no": periodic = 0
      else: error("Invalid PySCF command line args")
      iarg += 2
    elif args[iarg] == "-xcstr":
      if iarg+1 > narg: error("Invalid PySCF command line args")
      xcstr = args[iarg+1]
      iarg += 2
    elif args[iarg] == "-basis":
      if iarg+1 > narg: error("Invalid PySCF command line args")
      basis = args[iarg+1]
      iarg += 2
    else: error("Invalid PySCF command line args")

# --------------------------------------------
# operate as an engine
# --------------------------------------------

def mdi_engine(other_options):
  global world

  # get the MPI intra-communicator for this engine

  world = mdi.MDI_MPI_get_world_comm()
  me = world.Get_rank()
  nprocs = world.Get_size()

  # PySCF can only be invoked on a single MPI task

  if nprocs > 1:
    error("PySCF can only run on a single MPI task")

  # process non-MDI command line args

  options(other_options)

  # confirm PySCF is being run as an engine

  role = mdi.MDI_Get_Role()
  if not role == mdi.MDI_ENGINE:
    error("Must run PySCF as an MDI engine")

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
  mdi.MDI_Register_Command("@DEFAULT","<LATTICE_FORCES")
  #mdi.MDI_Register_Command("@DEFAULT","<STRESS")
  mdi.MDI_Register_Command("@DEFAULT","<CHARGES")

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

  else: error("Unrecognized MDI command: %s" % command)

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
# perform a quantum calculation via PySCF
# --------------------------------------------

def evaluate():
  global mode;
  global flag_qm_natoms,flag_mm_natoms
  global flag_box,flag_box_displ
  global flag_qm_elements,flag_qm_coords,flag_qm_potential
  global flag_mm_elements,flag_mm_coords,flag_mm_charges
  global qm_pe,qm_stress,qm_forces,qm_charges
  global mm_forces
  global dm_previous_exists,dm_previous

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
    dm_previous_exists = 0

  # if new system, error check that all needed MDI calls have been made

  if new_system:
    if periodic and not flag_box:
      error("Simulation box not specified for periodic system")
    if not flag_qm_natoms: error("QM atom count not specified")
    if not flag_qm_elements or not flag_qm_coords:
      error("QM atom properties not fully specified")
    if flag_mm_natoms:
      if not flag_mm_elements or not flag_mm_coords or not flag_mm_charges:
        error("MM atom properties not fully specified")

  # PySCF inputs for QM and MM atoms
  # box, qm_coords, mm_coords must be converted to Angstroms

  bohr_to_angstrom = mdi.MDI_Conversion_factor("bohr","angstrom")

  box_A = box * bohr_to_angstrom
  qm_coords_A = qm_coords * bohr_to_angstrom
  mm_coords_A = mm_coords * bohr_to_angstrom

  qm_symbols = [atomic_number_to_symbol[anum] for anum in qm_elements]
  mm_radii = [atomic_number_to_radius[anum] for anum in mm_elements]

  #pos2str = lambda pos: " ".join([str(x) for x in qm_coords_A])
  #atom_str = [f"{a} {pos2str(pos)}\n" for a,pos in zip(qm_symbols,qm_coords_A)]

  clines = ["%s %20.16g %20.16g %20.16g" % (symbol,xyz[0],xyz[1],xyz[2])
            for symbol,xyz in zip(qm_symbols,qm_coords_A)]
  atom_str = "\n".join(clines)

  if periodic:
    edge_vec = "%20.16g %20.16g %20.16g"
    box_str = "%s\n%s\n%s" % (edge_vec,edge_vec,edge_vec)
    box_str = box_str % \
      (box_A[0],box_A[1],box_A[2],box_A[3],box_A[4],box_A[5],box_A[6],box_A[7],box_A[8])
    #print("BOX STR:",box_str)

  #print("ATOM STR:",atom_str)
  #print("QM SYMB:",qm_symbols)
  #print("QM COORDS:",qm_coords)
  #print("MM COORDS:",mm_coords)
  #print("MM CHARGES:",mm_charges)
  #print("MM RADII:",mm_radii)

  # build PySCF system
  # use Cell for periodic, Mole for non-periodic

  if periodic:
    cell = Cell()
    #cell.verbose = 4      # uncomment to see more PySCF output
    cell.atom = atom_str
    cell.a = box_str
    cell.basis = basis
    cell.build()
  else:
    mol = Mole()
    #mol.verbose = 4       # uncomment to see more PySCF output
    mol.atom = atom_str
    mol.basis = basis
    #mol.max_memory = 10000
    mol.build()

  # QMMM with QM and MM atoms
  # mf = mean-field object
  # qm_pe = QM energy + QM/MM energy
  #   QM energy = QM_nuclear/QM_nuclear + electron/QM_nuclear + electron/electron
  #   QM/MM energy = QM_nuclear/MM_charges + electron/MM_charges
  # qm_forces = QM forces = same 3 terms
  # mm_forces = QM/MM forces = same 2 terms
  # dm = molecular orbitals (wave functions) for system

  if mode == QMMM:
    if periodic: mf = RKS_pbc(cell,xc=xcstr)
    else: mf = RKS_nonpbc(mol,xc=xcstr)
    mf = qmmm.mm_charge(mf,mm_coords,mm_charges,mm_radii)

    if dm_previous_exists:
      qm_pe = mf.kernel(dm0=dm_previous)
    else:
      qm_pe = mf.kernel()

    mf_grad = mf.nuc_grad_method()
    qm_forces = -mf_grad.kernel()
    dm = mf.make_rdm1()
    mm_forces = -(mf_grad.grad_nuc_mm() + mf_grad.contract_hcore_mm(dm))

    dm_previous_exists = 1
    dm_previous = dm

    #print("QM FORCES:",qm_forces)

  # AIMD with only QM atoms

  elif mode == AIMD:
    pass

  # clear flags for all MDI commands for next QM evaluation

  flag_qm_natoms = flag_mm_natoms = 0
  flag_box = flag_box_displ = 0
  flag_qm_elements = 0
  flag_qm_coords = flag_qm_potential = 0
  flag_mm_elements = 0
  flag_mm_coords = flag_mm_charges = 0

# --------------------------------------------
# function called by MDI driver
# only when it invokes pyscf_mdi.py as a plugin
# --------------------------------------------

def MDI_Plugin_init_pyscf_mdi(plugin_state):

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

  # mdi_index = index in sys.argv of -mdi
  # mdi_option = single arg in quotes that follows -mdi
  # other_options = all non-MDI args

  mdi_index = -1
  mdi_option = ""
  other_options = []

  narg = len(sys.argv)
  args = sys.argv

  iarg = 1
  while iarg < narg:
    arg = args[iarg]
    if arg == "-mdi" or arg == "--mdi":
      mdi_index = iarg
      if narg > iarg+1: mdi_option = sys.argv[iarg+1]
      else: error("PySCF -mdi argument not provided",0)
      iarg += 1
    else: other_options.append(arg)
    iarg += 1

  if not mdi_option: error("PySCF -mdi option not provided",0)

  # remove -mdi and its string from sys.argv
  # so that PySCF does not try to process it

  sys.argv.pop(mdi_index)
  sys.argv.pop(mdi_index)

  # disable this mode of MDI coupling for now
  # until issue on PySCF side is fixed

  # call MDI_Init with just -mdi option

  mdi.MDI_Init(mdi_option)

  # start running as an MDI engine

  mdi_engine(other_options)
