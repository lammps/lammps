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

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_neb.h"
#include "universe.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SINGLE_PROC_DIRECT,SINGLE_PROC_MAP,MULTI_PROC};

/* ---------------------------------------------------------------------- */

FixNEB::FixNEB(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_pe(NULL), pe(NULL), xprev(NULL), xnext(NULL), 
  tangent(NULL), xsend(NULL), xrecv(NULL), tagsend(NULL), tagrecv(NULL), 
  xsendall(NULL), xrecvall(NULL), tagsendall(NULL), tagrecvall(NULL), 
  counts(NULL), displacements(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal fix neb command");

  kspring = force->numeric(FLERR,arg[3]);
  if (kspring <= 0.0) error->all(FLERR,"Illegal fix neb command");

  // nreplica = number of partitions
  // ireplica = which world I am in universe
  // nprocs_universe = # of procs in all replicase
  // procprev,procnext = root proc in adjacent replicas

  me = comm->me;
  nprocs = comm->nprocs;

  nprocs_universe = universe->nprocs;
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

  // initialize local storage

  maxlocal = 0;
  ntotal = 0;

  xprev = xnext = tangent = NULL;
  xsend = xrecv = NULL;
  tagsend = tagrecv = NULL;
  xsendall = xrecvall = NULL;
  tagsendall = tagrecvall = NULL;
  counts = displacements = NULL;
}

/* ---------------------------------------------------------------------- */

FixNEB::~FixNEB()
{
  modify->delete_compute(id_pe);
  delete [] id_pe;

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);

  memory->destroy(xsend);
  memory->destroy(xrecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);

  memory->destroy(xsendall);
  memory->destroy(xrecvall);
  memory->destroy(tagsendall);
  memory->destroy(tagrecvall);

  memory->destroy(counts);
  memory->destroy(displacements);
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

  // nebatoms = # of atoms in fix group = atoms with inter-replica forces

  bigint count = group->count(igroup);
  if (count > MAXSMALLINT) error->all(FLERR,"Too many active NEB atoms");
  nebatoms = count;

  // comm mode for inter-replica exchange of coords

  if (nreplica == nprocs_universe &&
      nebatoms == atom->natoms && atom->sortfreq == 0) 
    cmode = SINGLE_PROC_DIRECT;
  else if (nreplica == nprocs_universe) cmode = SINGLE_PROC_MAP;
  else cmode = MULTI_PROC;

  // ntotal = total # of atoms in system, NEB atoms or not

  if (atom->natoms > MAXSMALLINT) error->all(FLERR,"Too many atoms for NEB");
  ntotal = atom->natoms;

  if (atom->nlocal > maxlocal) reallocate();

  if (MULTI_PROC && counts == NULL) {
    memory->create(xsendall,ntotal,3,"neb:xsendall");
    memory->create(xrecvall,ntotal,3,"neb:xrecvall");
    memory->create(tagsendall,ntotal,"neb:tagsendall");
    memory->create(tagrecvall,ntotal,"neb:tagrecvall");
    memory->create(counts,nprocs,"neb:counts");
    memory->create(displacements,nprocs,"neb:displacements");
  }
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

  // veng = PE of this replica
  // vprev,vnext = PEs of adjacent replicas
  // only proc 0 in each replica communicates

  vprev = vnext = veng = pe->compute_scalar();

  if (ireplica < nreplica-1 && me == 0) 
    MPI_Send(&veng,1,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0 && me == 0) 
    MPI_Recv(&vprev,1,MPI_DOUBLE,procprev,0,uworld,MPI_STATUS_IGNORE);

  if (ireplica > 0 && me == 0)
    MPI_Send(&veng,1,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1 && me == 0)
    MPI_Recv(&vnext,1,MPI_DOUBLE,procnext,0,uworld,MPI_STATUS_IGNORE);

  if (cmode == MULTI_PROC) {
    MPI_Bcast(&vprev,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&vnext,1,MPI_DOUBLE,0,world);
  }

  // communicate atoms to/from adjacent replicas to fill xprev,xnext

  inter_replica_comm();

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);

  // compute norm of GradV for log output

  double **f = atom->f;
  int nlocal = atom->nlocal;

  double fsq = 0.0;
  for (int i = 0; i < nlocal; i++)
    fsq += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];

  MPI_Allreduce(&fsq,&gradvnorm,1,MPI_DOUBLE,MPI_SUM,world);
  gradvnorm = sqrt(gradvnorm);

  // first or last replica has no change to forces, just return

  if (ireplica == 0 || ireplica == nreplica-1) {
    plen = nlen = 0.0;
    return;
  }

  // tangent = unit tangent vector in 3N space
  // based on delta vectors between atoms and their images in adjacent replicas
  // use one or two delta vecs to compute tangent,
  // depending on relative PEs of 3 replicas
  // see Henkelman & Jonsson 2000 paper, eqs 8-11

  double **x = atom->x;
  int *mask = atom->mask;

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

  double lenall;
  MPI_Allreduce(&tlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  tlen = sqrt(lenall);

  MPI_Allreduce(&plen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  plen = sqrt(lenall);

  MPI_Allreduce(&nlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  nlen = sqrt(lenall);

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

  double dotall;
  MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);

  double prefactor;
  if (ireplica == rclimber) prefactor = -2.0*dotall;
  else prefactor = -dotall + kspring*(nlen-plen);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] += prefactor*tangent[i][0];
      f[i][1] += prefactor*tangent[i][1];
      f[i][2] += prefactor*tangent[i][2];
    }
}

/* ----------------------------------------------------------------------
   send/recv NEB atoms to/from adjacent replicas
   received atoms matching my local atoms are stored in xprev,xnext
   replicas 0 and N-1 send but do not receive any atoms
------------------------------------------------------------------------- */

void FixNEB::inter_replica_comm()
{
  int i,m;
  MPI_Request request;
  MPI_Request requests[2];
  MPI_Status statuses[2];

  // reallocate memory if necessary

  if (atom->nlocal > maxlocal) reallocate();

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // -----------------------------------------------------
  // 3 cases: two for single proc per replica
  //          one for multiple procs per replica
  // -----------------------------------------------------

  // single proc per replica
  // all atoms are NEB atoms and no atom sorting
  // direct comm of x -> xprev and x -> xnext

  if (cmode == SINGLE_PROC_DIRECT) {
    if (ireplica > 0)
      MPI_Irecv(xprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
    if (ireplica < nreplica-1)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
    if (ireplica > 0) MPI_Wait(&request,MPI_STATUS_IGNORE);
    
    if (ireplica < nreplica-1)
      MPI_Irecv(xnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
    if (ireplica > 0)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
    if (ireplica < nreplica-1) MPI_Wait(&request,MPI_STATUS_IGNORE);

    return;
  }

  // single proc per replica
  // but only some atoms are NEB atoms or atom sorting is enabled
  // send atom IDs and coords of only NEB atoms to prev/next proc
  // recv procs use atom->map() to match received coords to owned atoms

  if (cmode == SINGLE_PROC_MAP) {
    m = 0;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tagsend[m] = tag[i];
        xsend[m][0] = x[i][0];
        xsend[m][1] = x[i][1];
        xsend[m][2] = x[i][2];
        m++;
      }

    if (ireplica > 0) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,&requests[1]);
    }
    if (ireplica < nreplica-1) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
    }

    if (ireplica > 0) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xprev[m][0] = xrecv[i][0];
        xprev[m][1] = xrecv[i][1];
        xprev[m][2] = xrecv[i][2];
      }
    }
      
    if (ireplica < nreplica-1) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,&requests[1]);
    }
    if (ireplica > 0) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
    }

    if (ireplica < nreplica-1) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xnext[m][0] = xrecv[i][0];
        xnext[m][1] = xrecv[i][1];
        xnext[m][2] = xrecv[i][2];
      }
    }

    return;
  }

  // multiple procs per replica
  // MPI_Gather all coords and atom IDs to root proc of each replica
  // send to root of adjacent replicas
  // bcast within each replica
  // each proc extracts info for atoms it owns via atom->map()

  m = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      tagsend[m] = tag[i];
      xsend[m][0] = x[i][0];
      xsend[m][1] = x[i][1];
      xsend[m][2] = x[i][2];
      m++;
    }

  MPI_Gather(&m,1,MPI_INT,counts,1,MPI_INT,0,world);
  displacements[0] = 0;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  MPI_Gatherv(tagsend,m,MPI_LMP_TAGINT,
              tagsendall,counts,displacements,MPI_LMP_TAGINT,0,world);
  for (i = 0; i < nprocs; i++) counts[i] *= 3;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  if (xsend)
    MPI_Gatherv(xsend[0],3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  else
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);

  if (ireplica > 0 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,
              &requests[1]);
  }
  if (ireplica < nreplica-1 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
  }

  if (ireplica > 0) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xprev[m][0] = xrecvall[i][0];
      xprev[m][1] = xrecvall[i][1];
      xprev[m][2] = xrecvall[i][2];
    }
  }

  if (ireplica < nreplica-1 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,
              &requests[1]);
  }
  if (ireplica > 0 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
  }

  if (ireplica < nreplica-1) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xnext[m][0] = xrecvall[i][0];
      xnext[m][1] = xrecvall[i][1];
      xnext[m][2] = xrecvall[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate xprev,xnext,tangent arrays if necessary
   reallocate communication arrays if necessary
------------------------------------------------------------------------- */

void FixNEB::reallocate()
{
  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);

  if (cmode != SINGLE_PROC_DIRECT) {
    memory->destroy(xsend);
    memory->destroy(xrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
  }

  maxlocal = atom->nmax;

  memory->create(xprev,maxlocal,3,"neb:xprev");
  memory->create(xnext,maxlocal,3,"neb:xnext");
  memory->create(tangent,maxlocal,3,"neb:tangent");

  if (cmode != SINGLE_PROC_DIRECT) {
    memory->create(xsend,maxlocal,3,"neb:xsend");
    memory->create(xrecv,maxlocal,3,"neb:xrecv");
    memory->create(tagsend,maxlocal,"neb:tagsend");
    memory->create(tagrecv,maxlocal,"neb:tagrecv");
  }
}
