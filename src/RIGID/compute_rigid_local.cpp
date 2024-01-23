// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_rigid_local.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_rigid_small.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{ID,MOL,MASS,X,Y,Z,XU,YU,ZU,VX,VY,VZ,FX,FY,FZ,IX,IY,IZ,
     TQX,TQY,TQZ,OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     QUATW,QUATI,QUATJ,QUATK,INERTIAX,INERTIAY,INERTIAZ};

/* ---------------------------------------------------------------------- */

ComputeRigidLocal::ComputeRigidLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  rstyle(nullptr), idrigid(nullptr), fixrigid(nullptr), vlocal(nullptr), alocal(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal compute rigid/local command");

  local_flag = 1;
  nvalues = narg - 4;

  idrigid = utils::strdup(arg[3]);

  rstyle = new int[nvalues];

  nvalues = 0;
  for (int iarg = 4; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"id") == 0) rstyle[nvalues++] = ID;
    else if (strcmp(arg[iarg],"mol") == 0) rstyle[nvalues++] = MOL;
    else if (strcmp(arg[iarg],"mass") == 0) rstyle[nvalues++] = MASS;
    else if (strcmp(arg[iarg],"x") == 0) rstyle[nvalues++] = X;
    else if (strcmp(arg[iarg],"y") == 0) rstyle[nvalues++] = Y;
    else if (strcmp(arg[iarg],"z") == 0) rstyle[nvalues++] = Z;
    else if (strcmp(arg[iarg],"xu") == 0) rstyle[nvalues++] = XU;
    else if (strcmp(arg[iarg],"yu") == 0) rstyle[nvalues++] = YU;
    else if (strcmp(arg[iarg],"zu") == 0) rstyle[nvalues++] = ZU;
    else if (strcmp(arg[iarg],"vx") == 0) rstyle[nvalues++] = VX;
    else if (strcmp(arg[iarg],"vy") == 0) rstyle[nvalues++] = VY;
    else if (strcmp(arg[iarg],"vz") == 0) rstyle[nvalues++] = VZ;
    else if (strcmp(arg[iarg],"fx") == 0) rstyle[nvalues++] = FX;
    else if (strcmp(arg[iarg],"fy") == 0) rstyle[nvalues++] = FY;
    else if (strcmp(arg[iarg],"fz") == 0) rstyle[nvalues++] = FZ;
    else if (strcmp(arg[iarg],"ix") == 0) rstyle[nvalues++] = IX;
    else if (strcmp(arg[iarg],"iy") == 0) rstyle[nvalues++] = IY;
    else if (strcmp(arg[iarg],"iz") == 0) rstyle[nvalues++] = IZ;
    else if (strcmp(arg[iarg],"tqx") == 0) rstyle[nvalues++] = TQX;
    else if (strcmp(arg[iarg],"tqy") == 0) rstyle[nvalues++] = TQY;
    else if (strcmp(arg[iarg],"tqz") == 0) rstyle[nvalues++] = TQZ;
    else if (strcmp(arg[iarg],"omegax") == 0) rstyle[nvalues++] = OMEGAX;
    else if (strcmp(arg[iarg],"omegay") == 0) rstyle[nvalues++] = OMEGAY;
    else if (strcmp(arg[iarg],"omegaz") == 0) rstyle[nvalues++] = OMEGAZ;
    else if (strcmp(arg[iarg],"angmomx") == 0) rstyle[nvalues++] = ANGMOMX;
    else if (strcmp(arg[iarg],"angmomy") == 0) rstyle[nvalues++] = ANGMOMY;
    else if (strcmp(arg[iarg],"angmomz") == 0) rstyle[nvalues++] = ANGMOMZ;
    else if (strcmp(arg[iarg],"quatw") == 0) rstyle[nvalues++] = QUATW;
    else if (strcmp(arg[iarg],"quati") == 0) rstyle[nvalues++] = QUATI;
    else if (strcmp(arg[iarg],"quatj") == 0) rstyle[nvalues++] = QUATJ;
    else if (strcmp(arg[iarg],"quatk") == 0) rstyle[nvalues++] = QUATK;
    else if (strcmp(arg[iarg],"inertiax") == 0) rstyle[nvalues++] = INERTIAX;
    else if (strcmp(arg[iarg],"inertiay") == 0) rstyle[nvalues++] = INERTIAY;
    else if (strcmp(arg[iarg],"inertiaz") == 0) rstyle[nvalues++] = INERTIAZ;
    else error->all(FLERR,"Invalid keyword in compute rigid/local command");
  }

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

ncount = nmax = 0;
  vlocal = nullptr;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeRigidLocal::~ComputeRigidLocal()
{
  memory->destroy(vlocal);
  memory->destroy(alocal);
  delete [] idrigid;
  delete [] rstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocal::init()
{
  // set fixrigid

  int ifix = modify->find_fix(idrigid);
  if (ifix < 0)
    error->all(FLERR,"FixRigidSmall ID for compute rigid/local does not exist");
  fixrigid = dynamic_cast<FixRigidSmall *>(modify->fix[ifix]);

  int flag = 0;
  if (strstr(fixrigid->style,"rigid/") == nullptr) flag = 1;
  if (strstr(fixrigid->style,"/small") == nullptr) flag = 1;
  if (flag)
    error->all(FLERR,"Compute rigid/local does not use fix rigid/small fix");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_rigid(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_rigid(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_rigid(1);
}

/* ----------------------------------------------------------------------
   count rigid bodies and compute rigid info on this proc
   if flag is set, compute requested info about rigid body
   owning atom of rigid body must be in group
------------------------------------------------------------------------- */

int ComputeRigidLocal::compute_rigid(int flag)
{
  int i,m,n,ibody;
  double *ptr;
  FixRigidSmall::Body *body;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    ibody = fixrigid->bodyown[i];
    if (ibody < 0) continue;
    body = &fixrigid->body[ibody];

    if (flag) {
      if (nvalues == 1) ptr = &vlocal[m];
      else ptr = alocal[m];

      for (n = 0; n < nvalues; n++) {
        switch (rstyle[n]) {
        case ID:
          ptr[n] = tag[body->ilocal];
          break;
        case MOL:
          ptr[n] = molecule[body->ilocal];
          break;
        case MASS:
          ptr[n] = body->mass;
          break;
        case X:
          ptr[n] = body->xcm[0];
          break;
        case Y:
          ptr[n] = body->xcm[1];
          break;
        case Z:
          ptr[n] = body->xcm[2];
          break;
        case XU:
          ptr[n] = body->xcm[0] +
            ((body->image & IMGMASK) - IMGMAX) * xprd;
          break;
        case YU:
          ptr[n] = body->xcm[1] +
            ((body->image >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
          break;
        case ZU:
          ptr[n] = body->xcm[2] +
            ((body->image >> IMG2BITS) - IMGMAX) * zprd;
          break;
        case VX:
          ptr[n] = body->vcm[0];
          break;
        case VY:
          ptr[n] = body->vcm[1];
          break;
        case VZ:
          ptr[n] = body->vcm[2];
          break;
        case FX:
          ptr[n] = body->fcm[0];
          break;
        case FY:
          ptr[n] = body->fcm[1];
          break;
        case FZ:
          ptr[n] = body->fcm[2];
          break;
        case IX:
          ptr[n] = (body->image & IMGMASK) - IMGMAX;
          break;
        case IY:
          ptr[n] = (body->image >> IMGBITS & IMGMASK) - IMGMAX;
          break;
        case IZ:
          ptr[n] = (body->image >> IMG2BITS) - IMGMAX;
          break;
        case TQX:
          ptr[n] = body->torque[0];
          break;
        case TQY:
          ptr[n] = body->torque[1];
          break;
        case TQZ:
          ptr[n] = body->torque[2];
          break;
        case OMEGAX:
          ptr[n] = body->omega[0];
          break;
        case OMEGAY:
          ptr[n] = body->omega[1];
          break;
        case OMEGAZ:
          ptr[n] = body->omega[2];
          break;
        case ANGMOMX:
          ptr[n] = body->angmom[0];
          break;
        case ANGMOMY:
          ptr[n] = body->angmom[1];
          break;
        case ANGMOMZ:
          ptr[n] = body->angmom[2];
          break;
        case QUATW:
          ptr[n] = body->quat[0];
          break;
        case QUATI:
          ptr[n] = body->quat[1];
          break;
        case QUATJ:
          ptr[n] = body->quat[2];
          break;
        case QUATK:
          ptr[n] = body->quat[3];
          break;
        case INERTIAX:
          ptr[n] = body->inertia[0];
          break;
        case INERTIAY:
          ptr[n] = body->inertia[1];
          break;
        case INERTIAZ:
          ptr[n] = body->inertia[2];
          break;
        }
      }
    }

    m++;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocal::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"rigid/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"rigid/local:array_local");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeRigidLocal::memory_usage()
{
  double bytes = (double)nmax*nvalues * sizeof(double);
  return bytes;
}
