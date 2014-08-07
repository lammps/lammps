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

/* ----------------------------------------------------------------------
   Contributing authors: Mike Parks (SNL), Ezwanur Rahman, J.T. Foster (UTSA)
------------------------------------------------------------------------- */

#include "math.h"
#include "fix_peri_neigh.h"
#include "pair_peri_pmb.h"
#include "pair_peri_lps.h"
#include "pair_peri_ves.h"
#include "pair_peri_eps.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPeriNeigh::FixPeriNeigh(LAMMPS *lmp,int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  isPMB = isLPS = isVES = isEPS = 0;
  if (force->pair_match("peri/pmb",1)) isPMB = 1;
  if (force->pair_match("peri/lps",1)) isLPS = 1;
  if (force->pair_match("peri/ves",1)) isVES = 1;
  if (force->pair_match("peri/eps",1)) isEPS = 1;
  
  restart_global = 1;
  restart_peratom = 1;
  first = 1;

  // perform initial allocation of atom-based arrays
  // register with atom class
  // set maxpartner = 1 as placeholder

  maxpartner = 1;
  npartner = NULL;
  partner = NULL;
  deviatorextention = NULL;
  deviatorBackextention = NULL;
  deviatorPlasticextension = NULL;
  lambdaValue = NULL;
  r0 = NULL;
  vinter = NULL;
  wvolume = NULL;

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so atom migration is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;

  // set comm sizes needed by this fix

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixPeriNeigh::~FixPeriNeigh()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(npartner);
  memory->destroy(partner);
  memory->destroy(deviatorextention);
  memory->destroy(deviatorBackextention);
  memory->destroy(deviatorPlasticextension);
  memory->destroy(lambdaValue);
  memory->destroy(r0);
  memory->destroy(vinter);
  memory->destroy(wvolume);
}

/* ---------------------------------------------------------------------- */

int FixPeriNeigh::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPeriNeigh::init()
{
  if (!first) return;

  // need a full neighbor list once

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  // compute PD scale factor, stored in Atom class, used by DumpCFG

  int nlocal = atom->nlocal;
  double vone = 0.0;
  for (int i = 0; i < nlocal; i++) vone += atom->vfrac[i];
  double vave;
  MPI_Allreduce(&vone,&vave,1,MPI_DOUBLE,MPI_SUM,world);
  if (atom->natoms) vave /= atom->natoms;
  if (vave > 0.0) atom->pdscale = 1.44 / pow(vave,1.0/3.0);
  else atom->pdscale = 1.0;
}

/* ---------------------------------------------------------------------- */

void FixPeriNeigh::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   For minimization: setup as with dynamics
------------------------------------------------------------------------- */

void FixPeriNeigh::min_setup(int vflag)
{
  setup(vflag);
}

/* ----------------------------------------------------------------------
   create initial list of neighbor partners via call to neighbor->build()
   must be done in setup (not init) since fix init comes before neigh init
------------------------------------------------------------------------- */

void FixPeriNeigh::setup(int vflag)
{
  int i,j,ii,jj,itype,jtype,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh;
  int **firstneigh;

  double **x = atom->x;
  double *vfrac = atom->vfrac;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

 // only build list of bonds on very first run

  if (!first) return;
  first = 0;

  // build full neighbor list, will copy or build as necessary

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // scan neighbor list to set maxpartner

  Pair *anypair = force->pair_match("peri",0);
  double **cutsq = anypair->cutsq;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq <= cutsq[itype][jtype]) npartner[i]++;
    }
  }

  maxpartner = 0;
  for (i = 0; i < nlocal; i++) maxpartner = MAX(maxpartner,npartner[i]);
  int maxall;
  MPI_Allreduce(&maxpartner,&maxall,1,MPI_INT,MPI_MAX,world);
  maxpartner = maxall;

  // realloc arrays with correct value for maxpartner

  memory->destroy(partner);
  memory->destroy(deviatorextention);
  memory->destroy(deviatorBackextention);
  memory->destroy(deviatorPlasticextension);
  memory->destroy(lambdaValue);  
  memory->destroy(r0);
  memory->destroy(npartner);

  npartner = NULL;
  partner = NULL;
  deviatorextention = NULL;
  deviatorBackextention = NULL;
  deviatorPlasticextension = NULL;
  lambdaValue = NULL;
  r0 = NULL;   
  grow_arrays(atom->nmax);

  // create partner list and r0 values from neighbor list
  // compute vinter for each atom

  for (i = 0; i < nlocal; i++) {
    npartner[i] = 0;
    vinter[i] = 0.0;
    wvolume[i] = 0.0;
    if (isEPS) lambdaValue[i] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq <= cutsq[itype][jtype]) {
        partner[i][npartner[i]] = tag[j];
        if (isVES)
          deviatorextention[i][npartner[i]] = 
            deviatorBackextention[i][npartner[i]] = 0.0;
        if (isEPS)
           deviatorPlasticextension[i][npartner[i]] = 0.0;
        r0[i][npartner[i]] = sqrt(rsq);   
        npartner[i]++;
        vinter[i] += vfrac[j];
      }
    }
  }

  // sanity check: does any atom appear twice in any neigborlist?
  // should only be possible if using pbc and domain < 2*delta

  if (domain->xperiodic || domain->yperiodic || domain->zperiodic) {
    for (i = 0; i < nlocal; i++) {
      jnum = npartner[i];
      for (jj = 0; jj < jnum; jj++) {
        for (int kk = jj+1; kk < jnum; kk++) {
          if (partner[i][jj] == partner[i][kk])
            error->one(FLERR,"Duplicate particle in PeriDynamic bond - "
                       "simulation box is too small");
        }
      }
    }
  }

  // compute wvolume for each atom

  double **x0 = atom->x0;
  double half_lc = 0.5*(domain->lattice->xlattice);
  double vfrac_scale;
  PairPeriLPS *pairlps = static_cast<PairPeriLPS*>(anypair);
  PairPeriVES *pairves = static_cast<PairPeriVES*>(anypair);
  PairPeriEPS *paireps = static_cast<PairPeriEPS*>(anypair);

  for (i = 0; i < nlocal; i++) {
    double xtmp0 = x0[i][0];
    double ytmp0 = x0[i][1];
    double ztmp0 = x0[i][2];
    jnum = npartner[i];
    itype = type[i];

    // loop over partners of particle i

    for (jj = 0; jj < jnum; jj++) {

      // if bond already broken, skip this partner

      if (partner[i][jj] == 0) continue;

      // lookup local index of partner particle

      j = atom->map(partner[i][jj]);

      // skip if particle is "lost"

      if (j < 0) continue;

      double delx0 = xtmp0 - x0[j][0];
      double dely0 = ytmp0 - x0[j][1];
      double delz0 = ztmp0 - x0[j][2];
            
      double rsq0 = delx0*delx0 + dely0*dely0 + delz0*delz0;

      jtype = type[j];
      double delta = sqrt(cutsq[itype][jtype]);

      // scale vfrac[j] if particle j near the horizon

      if ((fabs(r0[i][jj] - delta)) <= half_lc)
        vfrac_scale = (-1.0/(2*half_lc))*(r0[i][jj]) +
          (1.0 + ((delta - half_lc)/(2*half_lc) ) );
      else vfrac_scale = 1.0;

      // for PMB, influence = 1.0, otherwise invoke influence function
      if (isPMB) 
        wvolume[i] += 1.0 * rsq0 * vfrac[j] * vfrac_scale; 
      else if (isLPS)
        wvolume[i] += pairlps->influence_function(delx0,dely0,delz0) *
          rsq0 * vfrac[j] * vfrac_scale;
      else if (isVES)
        wvolume[i] += pairves->influence_function(delx0,dely0,delz0) *
          rsq0 * vfrac[j] * vfrac_scale;
      else if (isEPS)
        wvolume[i] += paireps->influence_function(delx0,dely0,delz0) *
          rsq0 * vfrac[j] * vfrac_scale;    
    }
  }

  // communicate wvolume to ghosts

  comm->forward_comm_fix(this);

  // bond statistics

  int n = 0;
  for (i = 0; i < nlocal; i++) n += npartner[i];
  int nall;
  MPI_Allreduce(&n,&nall,1,MPI_INT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Peridynamic bonds:\n");
      fprintf(screen,"  total # of bonds = %d\n",nall);
      fprintf(screen,"  bonds/atom = %g\n",(double)nall/atom->natoms);
    }
    if (logfile) {
      fprintf(logfile,"Peridynamic bonds:\n");
      fprintf(logfile,"  total # of bonds = %d\n",nall);
      fprintf(logfile,"  bonds/atom = %g\n",(double)nall/atom->natoms);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPeriNeigh::memory_usage()
{ 
  int nmax = atom->nmax;
  int bytes = nmax * sizeof(int);
  bytes += nmax*maxpartner * sizeof(tagint);
  bytes += nmax*maxpartner * sizeof(double);
  if (isVES) {
    bytes += nmax*maxpartner * sizeof(double);
    bytes += nmax*maxpartner * sizeof(double);
  }  
  if (isEPS) {
    bytes += nmax*maxpartner * sizeof(double);
    bytes += nmax * sizeof(double);
  }  
  bytes += nmax * sizeof(double);
  bytes += nmax * sizeof(double);
  return bytes; 
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPeriNeigh::grow_arrays(int nmax)
{
   memory->grow(npartner,nmax,"peri_neigh:npartner");
   memory->grow(partner,nmax,maxpartner,"peri_neigh:partner");
   if (isVES) {
     memory->grow(deviatorextention,nmax,maxpartner,
                  "peri_neigh:deviatorextention");
     memory->grow(deviatorBackextention,nmax,maxpartner,
                  "peri_neigh:deviatorBackextention");
   }
   if (isEPS) memory->grow(deviatorPlasticextension,nmax,maxpartner,
                           "peri_neigh:deviatorPlasticextension");
   memory->grow(r0,nmax,maxpartner,"peri_neigh:r0");
   if (isEPS) memory->grow(lambdaValue,nmax,"peri_neigh:lambdaValue");   
   memory->grow(vinter,nmax,"peri_neigh:vinter");
   memory->grow(wvolume,nmax,"peri_neigh:wvolume");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPeriNeigh::copy_arrays(int i, int j, int delflag)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    if (isVES) {
      deviatorextention[j][m] = deviatorextention[i][m];
      deviatorBackextention[j][m] = deviatorBackextention[i][m];
    }  
    if (isEPS)
      deviatorPlasticextension[j][m] = deviatorPlasticextension[i][m];
    r0[j][m] = r0[i][m];
  }
  if (isEPS) lambdaValue[j] = lambdaValue[i];
  vinter[j] = vinter[i];
  wvolume[j] = wvolume[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPeriNeigh::pack_exchange(int i, double *buf)
{
  // compact list by eliminating partner = 0 entries
  // set buf[0] after compaction

  int m = 1;
  for (int n = 0; n < npartner[i]; n++) {
    if (partner[i][n] == 0) continue;
    buf[m++] = partner[i][n];
    if (isVES) {
      buf[m++] = deviatorextention[i][n];
      buf[m++] = deviatorBackextention[i][n];
    } 
    if (isEPS) buf[m++] = deviatorPlasticextension[i][n];
    buf[m++] = r0[i][n];
  }
  if (isVES) buf[0] = m/4;
  else if (isEPS) buf[0] = m/3;
  else buf[0] = m/2;
  if (isEPS) buf[m++] = lambdaValue[i]; 
  buf[m++] = vinter[i];
  buf[m++] = wvolume[i];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixPeriNeigh::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (buf[m++]);
    if (isVES) {   
      deviatorextention[nlocal][n] = buf[m++];
      deviatorBackextention[nlocal][n] = buf[m++];
    }
    if (isEPS) deviatorPlasticextension[nlocal][n] = buf[m++];
    r0[nlocal][n] = buf[m++];     
  }
  if (isEPS) lambdaValue[nlocal] = buf[m++];
  vinter[nlocal] = buf[m++];
  wvolume[nlocal] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

int FixPeriNeigh::pack_forward_comm(int n, int *list, double *buf,
                                    int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = wvolume[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixPeriNeigh::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    wvolume[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixPeriNeigh::write_restart(FILE *fp)
{
  int n = 0;
  double list[2];
  list[n++] = first;
  list[n++] = maxpartner;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixPeriNeigh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  first = static_cast<int> (list[n++]);
  maxpartner = static_cast<int> (list[n++]);

  // grow 2D arrays now, cannot change size of 2nd array index later

  grow_arrays(atom->nmax);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPeriNeigh::pack_restart(int i, double *buf)
{
  int m = 0;
  if (isVES) buf[m++] = 4*npartner[i] + 4;
  else if (isEPS) buf[m++] = 3*npartner[i] + 5;
  else buf[m++] = 2*npartner[i] + 4;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    if (isVES) { 
      buf[m++] = deviatorextention[i][n];
      buf[m++] = deviatorBackextention[i][n];
    }  
    if (isEPS) buf[m++] = deviatorPlasticextension[i][n];
    buf[m++] = r0[i][n];
  }
  if (isEPS) buf[m++] = lambdaValue[i];
  buf[m++] = vinter[i];
  buf[m++] = wvolume[i];
  return m;  
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPeriNeigh::unpack_restart(int nlocal, int nth)
{

  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (extra[nlocal][m++]);
    if (isVES) { 
      deviatorextention[nlocal][n] = extra[nlocal][m++];
      deviatorBackextention[nlocal][n] = extra[nlocal][m++];
    }  
    if (isEPS) deviatorPlasticextension[nlocal][n] = extra[nlocal][m++];
    r0[nlocal][n] = extra[nlocal][m++];
  }
  if (isEPS) lambdaValue[nlocal] = extra[nlocal][m++];
  vinter[nlocal] = extra[nlocal][m++];
  wvolume[nlocal] = extra[nlocal][m++];  
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPeriNeigh::maxsize_restart()
{
  if (isVES) return 4*maxpartner + 4;
  if (isEPS) return 3*maxpartner + 5;
  return 2*maxpartner + 4;  
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPeriNeigh::size_restart(int nlocal)
{
  if (isVES) return 4*npartner[nlocal] + 4;
  if (isEPS) return 3*npartner[nlocal] + 5;
  return 2*npartner[nlocal] + 4; 
}
