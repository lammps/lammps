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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_neb.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNEB::FixNEB(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix neb command");

  kspring = atof(arg[3]);
  if (kspring <= 0.0) error->all(FLERR,"Illegal fix neb command");

  // nreplica = number of partitions
  // ireplica = which world I am in universe
  // procprev,procnext = root proc in adjacent replicas

  nreplica = universe->nworlds;
  ireplica = universe->iworld;

  if (ireplica > 0) procprev = universe->root_proc[ireplica-1];
  else procprev = -1;
  if (ireplica < nreplica-1) procnext = universe->root_proc[ireplica+1];
  else procnext = -1;
  uworld = universe->uworld;

  // create a new compute pe style
  // id = fix-ID + pe, compute group = all

  int n = strlen(id) + 4;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

  xprev = xnext = tangent = NULL;
}

/* ---------------------------------------------------------------------- */

FixNEB::~FixNEB()
{
  modify->delete_compute(id_pe);
  delete [] id_pe;

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
}

/* ---------------------------------------------------------------------- */

int FixNEB::setmask()
{
  int mask = 0;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNEB::init()
{
  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Potential energy ID for fix neb does not exist");
  pe = modify->compute[icompute];

  // turn off climbing mode, NEB command turns it on after init()

  rclimber = -1;

  // setup xprev and xnext arrays

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
  nebatoms = atom->nlocal;
  memory->create(xprev,nebatoms,3,"neb:xprev");
  memory->create(xnext,nebatoms,3,"neb:xnext");
  memory->create(tangent,nebatoms,3,"neb:tangent");
}

/* ---------------------------------------------------------------------- */

void FixNEB::min_setup(int vflag)
{
  min_post_force(vflag);

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixNEB::min_post_force(int vflag)
{
  double vprev,vnext,vmax,vmin;
  double delx,dely,delz;
  double delta1[3],delta2[3];
  MPI_Status status;
  MPI_Request request;

  // veng = PE of this replica
  // vprev,vnext = PEs of adjacent replicas

  veng = pe->compute_scalar();

  if (ireplica < nreplica-1) MPI_Send(&veng,1,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0) MPI_Recv(&vprev,1,MPI_DOUBLE,procprev,0,uworld,&status);

  if (ireplica > 0) MPI_Send(&veng,1,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1)
    MPI_Recv(&vnext,1,MPI_DOUBLE,procnext,0,uworld,&status);

  // xprev,xnext = atom coords of adjacent replicas
  // assume order of atoms in all replicas is the same
  // check that number of atoms hasn't changed

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (nlocal != nebatoms) error->one(FLERR,"Atom count changed in fix neb");

  if (ireplica > 0)
    MPI_Irecv(xprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
  if (ireplica < nreplica-1)
    MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0) MPI_Wait(&request,&status);

  if (ireplica < nreplica-1)
    MPI_Irecv(xnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
  if (ireplica > 0)
    MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1) MPI_Wait(&request,&status);

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);

  // Compute norm of GradV for log output

  double **f = atom->f;
  double fsq = 0.0;
  for (int i = 0; i < nlocal; i++) {
    fsq += f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2];
  }

  MPI_Allreduce(&fsq,&gradvnorm,1,MPI_DOUBLE,MPI_MAX,world);
  gradvnorm = sqrt(gradvnorm);

  // if this is first or last replica, no change to forces, just return

  if (ireplica == 0 || ireplica == nreplica-1) {
    plen = nlen = 0.0;
    return;
  }

  // tangent = unit tangent vector in 3N space
  // based on delta vectors between atoms and their images in adjacent replicas
  // use one or two delta vecs to compute tangent,
  // depending on relative PEs of 3 replicas
  // see Henkelman & Jonsson 2000 paper, eqs 8-11

  if (vnext > veng && veng > vprev) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tangent[i][0] = xnext[i][0] - x[i][0];
        tangent[i][1] = xnext[i][1] - x[i][1];
        tangent[i][2] = xnext[i][2] - x[i][2];
        domain->minimum_image(tangent[i]);
      }
  } else if (vnext < veng && veng < vprev) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tangent[i][0] = x[i][0] - xprev[i][0];
        tangent[i][1] = x[i][1] - xprev[i][1];
        tangent[i][2] = x[i][2] - xprev[i][2];
        domain->minimum_image(tangent[i]);
      }
  } else {
    vmax = MAX(fabs(vnext-veng),fabs(vprev-veng));
    vmin = MIN(fabs(vnext-veng),fabs(vprev-veng));
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delta1[0] = xnext[i][0] - x[i][0];
        delta1[1] = xnext[i][1] - x[i][1];
        delta1[2] = xnext[i][2] - x[i][2];
        domain->minimum_image(delta1);
        delta2[0] = x[i][0] - xprev[i][0];
        delta2[1] = x[i][1] - xprev[i][1];
        delta2[2] = x[i][2] - xprev[i][2];
        domain->minimum_image(delta2);
        if (vnext > vprev) {
          tangent[i][0] = vmax*delta1[0] + vmin*delta2[0];
          tangent[i][1] = vmax*delta1[1] + vmin*delta2[1];
          tangent[i][2] = vmax*delta1[2] + vmin*delta2[2];
        } else {
          tangent[i][0] = vmin*delta1[0] + vmax*delta2[0];
          tangent[i][1] = vmin*delta1[1] + vmax*delta2[1];
          tangent[i][2] = vmin*delta1[2] + vmax*delta2[2];
        }
      }
  }

  // tlen,plen,nlen = lengths of tangent, prev, next vectors

  double tlen = 0.0;
  plen = 0.0;
  nlen = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      tlen += tangent[i][0]*tangent[i][0] + tangent[i][1]*tangent[i][1] +
        tangent[i][2]*tangent[i][2];

      delx = x[i][0] - xprev[i][0];
      dely = x[i][1] - xprev[i][1];
      delz = x[i][2] - xprev[i][2];
      domain->minimum_image(delx,dely,delz);
      plen += delx*delx + dely*dely + delz*delz;

      delx = xnext[i][0] - x[i][0];
      dely = xnext[i][1] - x[i][1];
      delz = xnext[i][2] - x[i][2];
      domain->minimum_image(delx,dely,delz);
      nlen += delx*delx + dely*dely + delz*delz;
    }

  tlen = sqrt(tlen);
  plen = sqrt(plen);
  nlen = sqrt(nlen);

  // normalize tangent vector

  if (tlen > 0.0) {
    double tleninv = 1.0/tlen;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tangent[i][0] *= tleninv;
        tangent[i][1] *= tleninv;
        tangent[i][2] *= tleninv;
      }
  }

  // reset force on each atom in this replica
  // regular NEB for all replicas except rclimber does hill-climbing NEB
  // currently have F = -Grad(V) = -Grad(V)_perp - Grad(V)_parallel
  // want F = -Grad(V)_perp + Fspring for regular NEB
  // thus Fdelta = Grad(V)_parallel + Fspring for regular NEB
  // want F = -Grad(V) + 2 Grad(V)_parallel for hill-climbing NEB
  // thus Fdelta = 2 Grad(V)_parallel for hill-climbing NEB
  // Grad(V)_parallel = (Grad(V) . utan) * utangent = -(F . utan) * utangent
  // Fspring = k (nlen - plen) * utangent
  // see Henkelman & Jonsson 2000 paper, eqs 3,4,12
  // see Henkelman & Jonsson 2000a paper, eq 5

  double dot = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      dot += f[i][0]*tangent[i][0] + f[i][1]*tangent[i][1] +
        f[i][2]*tangent[i][2];
  }

  double prefactor;
  if (ireplica == rclimber) prefactor = -2.0*dot;
  else prefactor = -dot + kspring*(nlen-plen);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] += prefactor*tangent[i][0];
      f[i][1] += prefactor*tangent[i][1];
      f[i][2] += prefactor*tangent[i][2];
    }
}
