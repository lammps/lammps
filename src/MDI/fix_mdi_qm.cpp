/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_mdi_qm.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NATIVE, REAL, METAL };    // LAMMPS units which MDI supports

/* ---------------------------------------------------------------------- */

FixMDIQM::FixMDIQM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // check requirements for LAMMPS to work with MDI as an engine

  if (atom->tag_enable == 0) 
    error->all(FLERR, "Cannot use MDI engine without atom IDs");
  if (atom->natoms && atom->tag_consecutive() == 0)
    error->all(FLERR, "MDI engine requires consecutive atom IDs");

  // confirm LAMMPS is being run as a driver

  int role;
  MDI_Get_role(&role);
  if (role != MDI_DRIVER)
    error->all(FLERR, "Must invoke LAMMPS as an MDI driver to use fix mdi/qm");

  // optional args

  virialflag = 0;
  addflag = 1;
  every = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"virial") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix mdi/qm command");
      if (strcmp(arg[iarg],"yes") == 0) virialflag = 1;
      else if (strcmp(arg[iarg],"no") == 0) virialflag = 0;
      else error->all(FLERR,"Illegal fix mdi/qm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"add") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix mdi/qm command");
      if (strcmp(arg[iarg],"yes") == 0) addflag = 1;
      else if (strcmp(arg[iarg],"no") == 0) addflag = 0;
      else error->all(FLERR,"Illegal fix mdi/qm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix mdi/qm command");
      every = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (every <= 0) error->all(FLERR,"Illegal fix mdi/qm command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix mdi/qm command");
  }

  // fix output settings are based on optional keywords

  scalar_flag = 1;
  global_freq = every;
  extscalar = 1;

  if (virialflag) {
    vector_flag = 1;
    size_vector = 6;
    extvector = 1;
  }

  if (addflag) {
    energy_global_flag = 1;
    virial_global_flag = 1;
    thermo_energy = thermo_virial = 1;
  }

  // mdicomm will be one-time initialized in init()
  // cannot be done here for a plugin library, b/c mdi plugin command is later

  mdicomm = MDI_COMM_NULL;

  // peratom storage, both for nlocal and global natoms

  fqm = nullptr;
  maxlocal = 0;

  ibuf1 = ibuf1all = nullptr;
  buf3 = buf3all = nullptr;
  maxbuf = 0;

  // set unit conversion factors

  if (strcmp(update->unit_style, "real") == 0)
    lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0)
    lmpunits = METAL;
  else
    lmpunits = NATIVE;

  unit_conversions();

  nprocs = comm->nprocs;

  // initialize outputs

  qm_energy = 0.0;
  if (virialflag) {
    for (int i = 0; i < 6; i++) {
      qm_virial[i] = 0.0;
      virial[i] = 0.0;
    }
    sumflag = 0;
  }
}

/* ---------------------------------------------------------------------- */

FixMDIQM::~FixMDIQM()
{
  // send exit command to engine if it is a stand-alone code
  // for plugin, this happens in MDIPlugin::plugin_wrapper()

  if (!plugin) {
    int ierr = MDI_Send_command("EXIT", mdicomm);
    if (ierr) error->all(FLERR, "MDI: EXIT command");
  }

  // clean up

  memory->destroy(fqm);

  memory->destroy(ibuf1);
  memory->destroy(ibuf1all);
  memory->destroy(buf3);
  memory->destroy(buf3all);
}

/* ---------------------------------------------------------------------- */

int FixMDIQM::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::init()
{
  // one-time auto-detect whether engine is stand-alone code or plugin library
  // also initializes mdicomm
  // set plugin = 0/1 for engine = stand-alone code vs plugin library

  if (mdicomm == MDI_COMM_NULL) {
    MDI_Get_communicator(&mdicomm, 0);
    if (mdicomm == MDI_COMM_NULL) {
      plugin = 0;
      MDI_Accept_communicator(&mdicomm);
      if (mdicomm == MDI_COMM_NULL) 
        error->all(FLERR, "MDI unable to connect to stand-alone engine");
    } else {
      plugin = 1;
      int method;
      MDI_Get_method(&method, mdicomm);
      if (method != MDI_PLUGIN) 
        error->all(FLERR, "MDI internal error for plugin engine");
    }
  }

  // send natoms, atom types, and simulation box to engine
  // this will trigger setup of a new system
  // subsequent calls in post_force() will be for same system until new init()

  int ierr = MDI_Send_command(">NATOMS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NATOMS command");
  int n = static_cast<int> (atom->natoms);
  ierr = MDI_Send(&n, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NATOMS data");
  
  send_types();
  send_box();
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::post_force(int vflag)
{
  int index, ierr;

  // skip if timestep is not a multiple of every

  if (update->ntimestep % every) return;

  // reallocate peratom storage if necessary, both natoms and nlocal

  reallocate();

  // if simulation box dynamically changes, send current box to MDI engine

  if (domain->box_change_size || domain->box_change_shape)
    send_box();

  // gather all coords, ordered by atomID

  memset(buf3, 0, 3 * atom->natoms * sizeof(double));

  double **x = atom->x;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    index = static_cast<int>(tag[i]) - 1;
    buf3[3 * index + 0] = x[i][0] * lmp2mdi_length;
    buf3[3 * index + 1] = x[i][1] * lmp2mdi_length;
    buf3[3 * index + 2] = x[i][2] * lmp2mdi_length;
  }

  int n = static_cast<int> (atom->natoms);
  MPI_Reduce(buf3, buf3all, 3 * n, MPI_DOUBLE, MPI_SUM, 0, world);

  // send current coords to MDI engine

  ierr = MDI_Send_command(">COORDS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS command");
  ierr = MDI_Send(buf3all, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS data");

  // request potential energy from MDI engine
  // this triggers engine to perform QM calculation
  // qm_energy = fix output for global QM energy

  ierr = MDI_Send_command("<PE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE command");
  ierr = MDI_Recv(&qm_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE data");
  MPI_Bcast(&qm_energy, 1, MPI_DOUBLE, 0, world);
  qm_energy *= mdi2lmp_energy;

  // request forces from MDI engine

  ierr = MDI_Send_command("<FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES command");
  ierr = MDI_Recv(buf3, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES data");
  MPI_Bcast(buf3, 3 * n, MPI_DOUBLE, 0, world);

  // fqm = fix output for peratom QM forces
  // use atomID of local atoms to index into ordered buf3

  for (int i = 0; i < nlocal; i++) {
    index = static_cast<int>(tag[i]) - 1;
    fqm[i][0] = buf3[3 * index + 0] * mdi2lmp_force;
    fqm[i][1] = buf3[3 * index + 1] * mdi2lmp_force;
    fqm[i][2] = buf3[3 * index + 2] * mdi2lmp_force;
  }

  // optionally add forces to owned atoms
  // use atomID of local atoms to index into ordered buf3

  if (addflag) {
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      index = static_cast<int>(tag[i]) - 1;
      f[i][0] += buf3[3 * index + 0] * mdi2lmp_force;
      f[i][1] += buf3[3 * index + 1] * mdi2lmp_force;
      f[i][2] += buf3[3 * index + 2] * mdi2lmp_force;
    }
  }

  // optionally request pressure tensor from MDI engine, convert to virial
  // divide by nprocs so each proc stores a portion
  //   MPI_Allreduce is performed in compute_vector()
  // qm_virial = fix output for global QM virial

  if (virialflag) {
    double ptensor[6];
    ierr = MDI_Send_command("<STRESS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS command");
    ierr = MDI_Recv(ptensor, 6, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS data");
    MPI_Bcast(ptensor, 6, MPI_DOUBLE, 0, world);

    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (int i = 0; i < 6; i++) {
      ptensor[i] *= mdi2lmp_pressure;
      qm_virial[i] = ptensor[i] * volume / force->nktv2p / nprocs;
    }
    sumflag = 0;
  }

  // optionally set fix->virial

  if (virialflag && addflag) {
    for (int i = 0; i < 6; i++)
      virial[i] = qm_virial[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy from MDI engine
------------------------------------------------------------------------- */

double FixMDIQM::compute_scalar()
{
  return qm_energy;
}

/* ----------------------------------------------------------------------
   virial from MDI engine
------------------------------------------------------------------------- */

double FixMDIQM::compute_vector(int n)
{
  // only sum across procs one time

  if (sumflag == 0) {
    MPI_Allreduce(qm_virial, qm_virial_all, 6, MPI_DOUBLE, MPI_SUM, world);
    sumflag = 1;
  }

  return qm_virial_all[n];
}

/* ----------------------------------------------------------------------
   reallocate storage for local and global and atoms if needed
------------------------------------------------------------------------- */

void FixMDIQM::reallocate()
{
  if (atom->nlocal > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(fqm);
    memory->create(fqm, maxlocal, 3, "mdi:fqm");
    array_atom = fqm;
  }

  if (atom->natoms > maxbuf) {
    bigint nsize = atom->natoms * 3;
    if (nsize > MAXSMALLINT)
      error->all(FLERR, "Natoms too large to use with fix mdi/qm");
    
    maxbuf = static_cast<int> (atom->natoms);
    memory->destroy(ibuf1);
    memory->destroy(buf3);
    memory->destroy(buf3all);
    memory->create(ibuf1, maxbuf, "mdi:ibuf1");
    memory->create(ibuf1all, maxbuf, "mdi:ibuf1all");
    memory->create(buf3, 3 * maxbuf, "mdi:buf3");
    memory->create(buf3all, 3 * maxbuf, "mdi:buf3all");
  }
}

/* ----------------------------------------------------------------------
   send numeric atom types to MDI engine
------------------------------------------------------------------------- */

void FixMDIQM::send_types()
{
  memset(ibuf1, 0, atom->natoms * sizeof(int));

  // use local atomID to index into ordered ibuf1

  tagint *tag = atom->tag; 
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  int index;
  for (int i = 0; i < nlocal; i++) {
    index = static_cast<int>(tag[i]) - 1;
    ibuf1[index] = type[i];
  }

  int n = static_cast<int> (atom->natoms);
  MPI_Reduce(ibuf1, ibuf1all, n, MPI_INT, MPI_SUM, 0, world);

  int ierr = MDI_Send(ibuf1all, n, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES data");
}


/* ----------------------------------------------------------------------
   send simulation box size and shape to MDI engine
------------------------------------------------------------------------- */

void FixMDIQM::send_box()
{
  double cell[9];

  int ierr = MDI_Send_command(">CELL_DISPL", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL_DISPL command");
  cell[0] = domain->boxlo[0] * lmp2mdi_length;
  cell[1] = domain->boxlo[1] * lmp2mdi_length;
  cell[2] = domain->boxlo[2] * lmp2mdi_length;
  ierr = MDI_Send(cell, 3, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL_DISPL data");

  ierr = MDI_Send_command(">CELL", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL command");
  cell[0] = domain->boxhi[0] - domain->boxlo[0];
  cell[1] = 0.0;
  cell[2] = 0.0;
  cell[3] = domain->xy;
  cell[4] = domain->boxhi[1] - domain->boxlo[1];
  cell[5] = 0.0;
  cell[6] = domain->xz;
  cell[7] = domain->yz;
  cell[8] = domain->boxhi[2] - domain->boxlo[2];
  ierr = MDI_Send(cell, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL data");
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void FixMDIQM::unit_conversions()
{
  double angstrom_to_bohr, kelvin_to_hartree, ev_to_hartree, second_to_aut;

  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
  MDI_Conversion_Factor("second", "atomic_unit_of_time", &second_to_aut);

  // length units

  mdi2lmp_length = 1.0;
  lmp2mdi_length = 1.0;

  if (lmpunits == REAL || lmpunits == METAL) {
    lmp2mdi_length = angstrom_to_bohr;
    mdi2lmp_length = 1.0 / angstrom_to_bohr;
  }

  // energy units

  mdi2lmp_energy = 1.0;
  lmp2mdi_energy = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_energy = kelvin_to_hartree / force->boltz;
    mdi2lmp_energy = force->boltz / kelvin_to_hartree;
  } else if (lmpunits == METAL) {
    lmp2mdi_energy = ev_to_hartree;
    mdi2lmp_energy = 1.0 / ev_to_hartree;
  }

  // force units = energy/length

  mdi2lmp_force = 1.0;
  lmp2mdi_force = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_force = (kelvin_to_hartree / force->boltz) / angstrom_to_bohr;
    mdi2lmp_force = 1.0 / lmp2mdi_force;
  } else if (lmpunits == METAL) {
    lmp2mdi_force = ev_to_hartree / angstrom_to_bohr;
    mdi2lmp_force = angstrom_to_bohr / ev_to_hartree;
  }

  // pressure or stress units = force/area = energy/volume

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_pressure = (kelvin_to_hartree / force->boltz) /
        (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  } else if (lmpunits == METAL) {
    lmp2mdi_pressure =
        ev_to_hartree / (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  }
}
