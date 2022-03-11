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

#include "error.h"
#include "fix_mdi_aimd.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NATIVE,REAL,METAL};   // LAMMPS units which MDI supports

/* ---------------------------------------------------------------------- */

FixMDIAimd::FixMDIAimd(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix mdi/aimd command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // check requirements for LAMMPS to work with MDI as an engine

  if (atom->tag_enable == 0) 
    error->all(FLERR, "Cannot use MDI engine without atom IDs");

  if (atom->natoms && atom->tag_consecutive() == 0) 
    error->all(FLERR, "MDI engine requires consecutive atom IDs");

  // confirm LAMMPS is being run as a driver
  
  int role;
  MDI_Get_role(&role);
  if (role != MDI_DRIVER) 
    error->all(FLERR,"Must invoke LAMMPS as an MDI driver to use fix mdi/aimd");

  // storage for all atoms

  buf3 = buf3all = nullptr;
  maxbuf = 0;

  // set unit conversion factors

  if (strcmp(update->unit_style, "real") == 0) lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0) lmpunits = METAL;
  else lmpunits = NATIVE;

  unit_conversions();

  // connect to MDI engine

  MDI_Accept_communicator(&mdicomm);
  if (mdicomm <= 0) error->all(FLERR, "Unable to connect to MDI engine");

  nprocs = comm->nprocs;
}

/* ---------------------------------------------------------------------- */

FixMDIAimd::~FixMDIAimd()
{
  // send exit command to engine

  int ierr = MDI_Send_command("EXIT",mdicomm);
  if (ierr) error->all(FLERR,"MDI: EXIT command");

  // clean up

  memory->destroy(buf3);
  memory->destroy(buf3all);
}

/* ---------------------------------------------------------------------- */

int FixMDIAimd::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to request energy
------------------------------------------------------------------------- */

void FixMDIAimd::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::post_force(int vflag)
{
  int ilocal,ierr;
  double cell[9];

  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // if simulation box dynamically changes, send current box to MDI engine

  if (domain->box_change_size || domain->box_change_shape) {
    ierr = MDI_Send_command(">CELL_DISPL",mdicomm);
    if (ierr) error->all(FLERR,"MDI: >CELL_DISPL command");
    cell[0] = domain->boxlo[0] * lmp2mdi_length;
    cell[1] = domain->boxlo[1] * lmp2mdi_length;
    cell[2] = domain->boxlo[2] * lmp2mdi_length;
    ierr = MDI_Send(cell,3,MDI_DOUBLE,mdicomm);
    if (ierr) error->all(FLERR,"MDI: >CELL_DISPL data");

    ierr = MDI_Send_command(">CELL",mdicomm);
    if (ierr) error->all(FLERR,"MDI: >CELL command");
    cell[0] = domain->boxhi[0] - domain->boxlo[0];
    cell[1] = 0.0;
    cell[2] = 0.0;
    cell[3] = domain->xy;
    cell[4] = domain->boxhi[1] - domain->boxlo[1];
    cell[5] = 0.0;
    cell[6] = domain->xz;
    cell[7] = domain->yz;
    cell[8] = domain->boxhi[2] - domain->boxlo[2];
    ierr = MDI_Send(cell,9,MDI_DOUBLE,mdicomm);
    if (ierr) error->all(FLERR,"MDI: >CELL data");
  }

  // gather all coords, ordered by atomID

  reallocate();
  memset(buf3,0,3*atom->natoms*sizeof(double));

  double **x = atom->x;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int> (tag[i]) - 1;
    buf3[3*ilocal+0] = x[i][0] * lmp2mdi_length;
    buf3[3*ilocal+1] = x[i][1] * lmp2mdi_length;
    buf3[3*ilocal+2] = x[i][2] * lmp2mdi_length;
  }

  MPI_Reduce(buf3,buf3all,3*atom->natoms,MPI_DOUBLE,MPI_SUM,0,world);

  // send current coords to MDI engine

  ierr = MDI_Send_command(">COORDS",mdicomm);
  if (ierr) error->all(FLERR,"MDI: >COORDS command");
  ierr = MDI_Send(buf3all,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >COORDS data");

  // trigger engine to evaluate forces,energy,pressure for current system

  ierr = MDI_Send_command("EVAL",mdicomm);
  if (ierr) error->all(FLERR,"MDI: EVAL command");

  // request forces from MDI engine

  ierr = MDI_Send_command("<FORCES",mdicomm);
  if (ierr) error->all(FLERR,"MDI: <FORCES command");
  ierr = MDI_Recv(buf3,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <FORCES data");
  MPI_Bcast(buf3,3*atom->natoms,MPI_DOUBLE,0,world);

  // add forces to owned atoms
  // use atomID to index into ordered buf3

  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int> (tag[i]) - 1;
    f[i][0] += buf3[3*ilocal+0] * mdi2lmp_force;
    f[i][1] += buf3[3*ilocal+1] * mdi2lmp_force;
    f[i][2] += buf3[3*ilocal+2] * mdi2lmp_force;
  }

  // optionally request energy from MDI engine
  // divide by nprocs so each proc stores a portion

  if (eflag_global) {
    ierr = MDI_Send_command("<PE",mdicomm);
    if (ierr) error->all(FLERR,"MDI: <PE command");
    ierr = MDI_Recv(&engine_energy,1,MDI_DOUBLE,mdicomm);
    if (ierr) error->all(FLERR,"MDI: <PE data");
    MPI_Bcast(&engine_energy,1,MPI_DOUBLE,0,world);
    engine_energy *= mdi2lmp_energy / nprocs;
  }

  // optionally request pressure tensor from MDI engine, convert to virial
  // divide by nprocs so each proc stores a portion

  if (vflag_global) {
    double ptensor[6];
    ierr = MDI_Send_command("<PTENSOR",mdicomm);
    if (ierr) error->all(FLERR,"MDI: <PTENSOR command");
    ierr = MDI_Recv(ptensor,6,MDI_DOUBLE,mdicomm);
    if (ierr) error->all(FLERR,"MDI: <PTENSOR data");
    MPI_Bcast(ptensor,6,MPI_DOUBLE,0,world);

    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (int i = 0; i < 6; i++) {
      ptensor[i] *= mdi2lmp_pressure;
      virial[i] = ptensor[i] * volume / force->nktv2p / nprocs;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy from MDI engine
------------------------------------------------------------------------- */

double FixMDIAimd::compute_scalar()
{
  return engine_energy;
}

/* ----------------------------------------------------------------------
   reallocate storage for all atoms if necessary
------------------------------------------------------------------------- */

void FixMDIAimd::reallocate()
{
  if (atom->natoms <= maxbuf) return;

  if (3*atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Natoms too large to use with fix mdi/aimd");

  maxbuf = atom->natoms;

  memory->destroy(buf3);
  memory->destroy(buf3all);

  memory->create(buf3,3*maxbuf,"mdi:buf3");
  memory->create(buf3all,3*maxbuf,"mdi:buf3all");
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void FixMDIAimd::unit_conversions()
{
  double angstrom_to_bohr,kelvin_to_hartree,ev_to_hartree;

  MDI_Conversion_factor("angstrom","bohr",&angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy","hartree",&kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt","hartree",&ev_to_hartree);

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

  // force units

  mdi2lmp_force = 1.0;
  lmp2mdi_force = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_force = (kelvin_to_hartree / force->boltz) / angstrom_to_bohr;
    mdi2lmp_force = 1.0 / lmp2mdi_force;
  } else if (lmpunits == METAL) {
    lmp2mdi_force = ev_to_hartree / angstrom_to_bohr;
    mdi2lmp_force = angstrom_to_bohr / ev_to_hartree;
  }

  // pressure units

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;
}
