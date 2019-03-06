/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstring>
#include "fix_client_md.h"
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"

// CSlib interface

#include "cslib.h"

using namespace LAMMPS_NS;
using namespace CSLIB_NS;
using namespace FixConst;

enum{OTHER,REAL,METAL};
enum{SETUP=1,STEP};
enum{DIM=1,PERIODICITY,ORIGIN,BOX,NATOMS,NTYPES,TYPES,COORDS,UNITS,CHARGE};
enum{FORCES=1,ENERGY,PRESSURE,ERROR};

/* ---------------------------------------------------------------------- */

FixClientMD::FixClientMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->clientserver != 1)
    error->all(FLERR,"Fix client/md requires LAMMPS be running as a client");
  if (!atom->map_style) error->all(FLERR,"Fix client/md requires atom map");

  if (sizeof(tagint) != 4)
    error->all(FLERR,"Fix client/md requires 4-byte atom IDs");

  if (strcmp(update->unit_style,"real") == 0) units = REAL;
  else if (strcmp(update->unit_style,"metal") == 0) units = METAL;
  else units = OTHER;

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  virial_flag = 1;
  thermo_virial = 1;

  inv_nprocs = 1.0 / comm->nprocs;
  if (domain->dimension == 2)
    box[0][2] = box[1][2] = box[2][0] = box[2][1] = box[2][2] = 0.0;

  maxatom = 0;
  xpbc = NULL;

  // unit conversion factors for REAL
  // otherwise not needed
  // message received in METAL units, convert to local REAL units

  fconvert = econvert = pconvert = 1.0;
  if (units == REAL) {
    fconvert = econvert = 23.06035;    // eV -> Kcal/mole
    pconvert = 0.986923;               // bars -> atmospheres
  }
}

/* ---------------------------------------------------------------------- */

FixClientMD::~FixClientMD()
{
  memory->destroy(xpbc);

  CSlib *cs = (CSlib *) lmp->cslib;

  // all-done message to server

  cs->send(-1,0);

  int nfield;
  int *fieldID,*fieldtype,*fieldlen;
  int msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);

  // clean-up

  delete cs;
  lmp->cslib = NULL;
}

/* ---------------------------------------------------------------------- */

int FixClientMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixClientMD::init()
{
  if (3*atom->natoms > INT_MAX)
    error->all(FLERR,"Fix client/md max atoms is 1/3 of 2^31");
}

/* ---------------------------------------------------------------------- */

void FixClientMD::setup(int vflag)
{
  CSlib *cs = (CSlib *) lmp->cslib;

  // SETUP send at beginning of each run
  // required fields: DIM, PERIODICTY, ORIGIN, BOX, NATOMS, NTYPES, TYPES, COORDS
  // optional fields: others in enum above

  int nfields = 8;
  if (units == OTHER) nfields++;
  if (atom->q_flag) nfields++;

  cs->send(SETUP,nfields);

  cs->pack_int(DIM,domain->dimension);
  cs->pack(PERIODICITY,1,3,domain->periodicity);

  pack_box();
  cs->pack(ORIGIN,4,3,domain->boxlo);
  cs->pack(BOX,4,9,&box[0][0]);

  cs->pack_int(NATOMS,atom->natoms);
  cs->pack_int(NTYPES,atom->ntypes);

  cs->pack_parallel(TYPES,1,atom->nlocal,atom->tag,1,atom->type);
  pack_coords();
  cs->pack_parallel(COORDS,4,atom->nlocal,atom->tag,3,xpbc);

  if (units == OTHER) cs->pack_string(UNITS,update->unit_style);

  if (atom->q_flag)
    cs->pack_parallel(CHARGE,4,atom->nlocal,atom->tag,1,atom->q);

  // receive initial forces, energy, virial

  receive_fev(vflag);

  if (server_error) {
    char str[64];
    sprintf(str,"Fix client/md received server error %d\n",server_error);
    error->all(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void FixClientMD::min_setup(int vflag)
{
  setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixClientMD::post_force(int vflag)
{
  int i,j,m;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // STEP send every step
  // required fields: COORDS
  // optional fields: ORIGIN, BOX

  // send coords

  CSlib *cs = (CSlib *) lmp->cslib;

  int nfields = 1;
  if (domain->box_change) nfields += 2;

  cs->send(STEP,nfields);

  pack_coords();
  cs->pack_parallel(COORDS,4,atom->nlocal,atom->tag,3,xpbc);

  if (domain->box_change) {
    pack_box();
    cs->pack(ORIGIN,4,3,domain->boxlo);
    cs->pack(BOX,4,9,&box[0][0]);
  }

  // receive forces, energy, virial

  receive_fev(vflag);

  if (server_error) {
    char str[64];
    sprintf(str,"Fix client/md received server error %d\n",server_error);
    error->all(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void FixClientMD::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy from server
------------------------------------------------------------------------- */

double FixClientMD::compute_scalar()
{
  return eng;
}

/* ----------------------------------------------------------------------
   pack local coords into xpbc, enforcing PBC
------------------------------------------------------------------------- */

void FixClientMD::pack_coords()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (nlocal > maxatom) {
    memory->destroy(xpbc);
    maxatom = atom->nmax;
    memory->create(xpbc,3*maxatom,"message:xpbc");
  }

  memcpy(xpbc,&x[0][0],3*nlocal*sizeof(double));

  int j = 0;
  for (int i = 0; i < nlocal; i++) {
    domain->remap(&xpbc[j]);
    j += 3;
  }
}

/* ----------------------------------------------------------------------
   pack box info into box = 3 edge vectors of simulation box
------------------------------------------------------------------------- */

void FixClientMD::pack_box()
{
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  box[0][0] = boxhi[0] - boxlo[0];
  box[0][1] = 0.0;
  box[0][2] = 0.0;
  box[1][0] = domain->xy;
  box[1][1] = boxhi[1] - boxlo[1];
  box[1][2] = 0.0;
  box[2][0] = domain->xz;
  box[2][1] = domain->yz;
  box[2][2] = boxhi[2] - boxlo[2];
}

/* ----------------------------------------------------------------------
   receive message from server
   required fields: FORCES, ENERGY, PRESSURE
   optional field: ERROR
------------------------------------------------------------------------- */

void FixClientMD::receive_fev(int vflag)
{
  CSlib *cs = (CSlib *) lmp->cslib;

  int nfield;
  int *fieldID,*fieldtype,*fieldlen;

  int msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);

  double *forces = (double *) cs->unpack(FORCES);
  double **f = atom->f;
  int nlocal = atom->nlocal;
  bigint natoms = atom->natoms;
  int m;

  int j = 0;
  for (tagint id = 1; id <= natoms; id++) {
    m = atom->map(id);
    if (m < 0 || m >= nlocal) j += 3;
    else {
      f[m][0] += fconvert * forces[j++];
      f[m][1] += fconvert * forces[j++];
      f[m][2] += fconvert * forces[j++];
    }
  }

  eng = econvert * cs->unpack_double(ENERGY);

  if (vflag) {
    double *v = (double *) cs->unpack(PRESSURE);

    double nktv2p = force->nktv2p;
    double volume = domain->xprd * domain->yprd * domain->zprd;
    double factor = inv_nprocs * pconvert * volume / nktv2p;

    for (int i = 0; i < 6; i++)
      virial[i] = factor * v[i];
  }

  // error return

  server_error = 0;
  if (nfield == 4) server_error = cs->unpack_int(ERROR);
}
