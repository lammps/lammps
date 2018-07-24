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
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

// CSlib interface

#include "cslib.h"

using namespace LAMMPS_NS;
using namespace CSLIB_NS;
using namespace FixConst;

enum{SETUP=1,STEP};
enum{UNITS=1,DIM,NATOMS,NTYPES,BOXLO,BOXHI,BOXTILT,TYPES,COORDS,CHARGE};
enum{FORCES=1,ENERGY,VIRIAL};

/* ---------------------------------------------------------------------- */

FixClientMD::FixClientMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->clientserver != 1) 
    error->all(FLERR,"Fix client/md requires LAMMPS be running as a client");
  if (!atom->map_style) error->all(FLERR,"Fix client/md requires atom map");

  if (sizeof(tagint) != 4) 
    error->all(FLERR,"Fix client/md requires 4-byte atom IDs");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  virial_flag = 1;
  thermo_virial = 1;

  maxatom = 0;
  xpbc = NULL;
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

  // required fields: NATOMS, NTYPES, BOXLO, BOXHI, TYPES, COORDS
  // optional fields: others in enum above

  int nfields = 6;
  if (domain->dimension == 2) nfields++;
  if (domain->triclinic) nfields++;
  if (atom->q_flag) nfields++;

  cs->send(SETUP,nfields);

  cs->pack_int(NATOMS,atom->natoms);
  cs->pack_int(NTYPES,atom->ntypes);
  cs->pack(BOXLO,4,3,domain->boxlo);
  cs->pack(BOXHI,4,3,domain->boxhi);
  cs->pack_parallel(TYPES,1,atom->nlocal,atom->tag,1,atom->type);
  pack_coords();
  cs->pack_parallel(COORDS,4,atom->nlocal,atom->tag,3,xpbc);

  if (domain->dimension == 2) cs->pack_int(DIM,domain->dimension);
  if (domain->triclinic) {
    double boxtilt[3];
    boxtilt[0] = domain->xy;
    if (domain->dimension == 3) {
      boxtilt[1] = domain->xz;
      boxtilt[2] = domain->yz;
    } else boxtilt[1] = boxtilt[2] = 0.0;
    cs->pack(BOXTILT,4,3,boxtilt);
  }
  if (atom->q_flag)
    cs->pack_parallel(CHARGE,4,atom->nlocal,atom->tag,1,atom->q);

  // receive initial forces, energy, virial

  receive_fev(vflag);
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

  // required fields: COORDS
  // optional fields: BOXLO, BOXHI, BOXTILT

  // send coords

  CSlib *cs = (CSlib *) lmp->cslib;

  int nfields = 1;
  if (domain->box_change) nfields += 2;
  if (domain->box_change && domain->triclinic) nfields++;;

  cs->send(STEP,nfields);

  pack_coords();
  cs->pack_parallel(COORDS,4,atom->nlocal,atom->tag,3,xpbc);

  if (domain->box_change) {
    cs->pack(BOXLO,4,3,domain->boxlo);
    cs->pack(BOXHI,4,3,domain->boxhi);
    if (domain->triclinic) {
      double boxtilt[3];
      boxtilt[0] = domain->xy;
      if (domain->dimension == 3) {
        boxtilt[1] = domain->xz;
        boxtilt[2] = domain->yz;
      } else boxtilt[1] = boxtilt[2] = 0.0;
      cs->pack(BOXTILT,4,3,boxtilt);
    }
  }

  // recv forces, energy, virial

  receive_fev(vflag);
}

/* ---------------------------------------------------------------------- */

void FixClientMD::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy from QM code
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
   receive message from server with forces, energy, virial
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
      f[m][0] += forces[j++];
      f[m][1] += forces[j++];
      f[m][2] += forces[j++];
    }
  }

  eng = cs->unpack_double(ENERGY);

  if (vflag) {
    double *v = (double *) cs->unpack(VIRIAL);
    double invnprocs = 1.0 / comm->nprocs;
    for (int i = 0; i < 6; i++)
      virial[i] = invnprocs*v[i];
  }
}
