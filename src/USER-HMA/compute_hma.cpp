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
#include <string.h>
#include "compute_hma.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"

#include <vector>


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHMA::ComputeHMA(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_temp(NULL), deltaR(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute hma command");                    
  if (igroup) error->all(FLERR,"Compute hma must use group all");
  if (strcmp(arg[3],"NULL") == 0) {error->all(FLERR,"fix ID specifying the set temperature of canonical simulation is required");}       
  else {
    int n = strlen(arg[3]) + 1;                  
    id_temp = new char[n];   
    strcpy(id_temp,arg[3]);                         
  }
 
  create_attribute = 1;                  
  extscalar = 1;                         
  timeflag = 1;                         

  // (from compute displace/atom) create a new fix STORE style  
  // our new fix's id (id_fix)= compute-ID + COMPUTE_STORE 
  // our new fix's group = same as compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];     
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";              
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);        
  fix = (FixStore *) modify->fix[modify->nfix-1];    
 
  delete [] newarg;                               

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;   
  else {
    double **xoriginal = fix->astore;            
    double **x = atom->x;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      domain->unmap(x[i],image[i],xoriginal[i]);          
  }

  vector_flag = 1;
  extvector = -1;
  comm_forward = 0; // 3 if 2nd derivative needed

  computeU = computeP = computeCv = -1;
  returnAnharmonic = 0;
  std::vector<int> extvec;
  for (int iarg=4; iarg<narg; iarg++) {
    if (!strcasecmp(arg[iarg], "u")) {
      computeU = extvec.size();
      extvec.push_back(1);
    }
    else if (!strcasecmp(arg[iarg], "p")) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix hma command");
      computeP = extvec.size();
      deltaPcap = force->numeric(FLERR, arg[iarg+1]);
      extvec.push_back(0);
      iarg++;
    }
    else if (!strcasecmp(arg[iarg], "cv")) {
      computeCv = extvec.size();
      comm_forward = 3;
      extvec.push_back(1);
    }
    else if (!strcasecmp(arg[iarg], "anharmonic")) {
      // the first time we're called, we'll grab lattice pressure and energy
      returnAnharmonic = -1;
    }
    else {
      error->all(FLERR,"Illegal fix hma command");
    }
  }

  if (extvec.size() == 0) {
    error->all(FLERR,"Illegal fix hma command");
  }
  size_vector = extvec.size();
  memory->create(vector, size_vector, "hma::vector");
  extlist = new int[size_vector];
  for (int i=0; i<size_vector; i++) {
    extlist[i] = extvec[i];
  }

  if (computeU>-1 || computeCv>-1) {
   peflag = 1;                             
  }
  if (computeP>-1) {
    pressflag = 1;
  }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeHMA::~ComputeHMA()
{
  // check nfix in case all fixes have already been deleted
  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] id_temp;
  delete [] extlist;
  memory->destroy(vector);
  memory->destroy(deltaR);
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::init() {
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0; 
  neighbor->requests[irequest]->compute = 1; 
  neighbor->requests[irequest]->occasional = 1; 
}

void ComputeHMA::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

void ComputeHMA::setup()
{
  int dummy=0;
  int ifix = modify->find_fix(id_temp);
  if (ifix < 0) error->all(FLERR,"Could not find compute hma temperature ID");
  double * temperat = (double *) modify->fix[ifix]->extract("t_target",dummy);      
  if (temperat==NULL) error->all(FLERR,"Could not find compute hma temperature ID");
  finaltemp = * temperat;       


  // set fix which stores original atom coords

  int ifix2 = modify->find_fix(id_fix);
  if (ifix2 < 0) error->all(FLERR,"Could not find compute hma ID");
  fix = (FixStore *) modify->fix[ifix2];             
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::compute_vector()
{
  invoked_vector = update->ntimestep;         

  // grow deltaR array if necessary
  if (comm_forward>0 && atom->nmax > nmax) {
    memory->destroy(deltaR);
    nmax = atom->nmax;
    memory->create(deltaR,nmax,3,"hma:deltaR");
  }

  double **xoriginal = fix->astore;
  double fdr = 0.0;
  double **x = atom->x;
  double **f = atom->f;

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;   
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double u = 0.0;
  if (computeU>-1 || computeCv>-1) {
    if (force->pair) u += force->pair->eng_vdwl + force->pair->eng_coul;
    if (force->bond) u += force->bond->energy;
    if (force->angle) u += force->angle->energy;
    if (force->dihedral) u += force->dihedral->energy;
    if (force->improper) u += force->improper->energy;
  }

  int dimension = domain->dimension;
  double p = 0, vol = 0;
  if (computeP>-1) {
    p = virial_compute(dimension);
    vol = xprd * yprd;
    if (dimension == 3) vol *= zprd;
    p *= force->nktv2p / (dimension*vol);
    if (returnAnharmonic == -1) {
      pLat = p;
    }
  }

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++) {
      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      double dx = x[i][0] + xbox*xprd - xoriginal[i][0];
      double dy = x[i][1] + ybox*yprd - xoriginal[i][1];
      double dz = x[i][2] + zbox*zprd - xoriginal[i][2];
      if (comm_forward>0) {
        deltaR[i][0] = dx;
        deltaR[i][1] = dy;
        deltaR[i][2] = dz;
      }
      fdr += dx*f[i][0] + dy*f[i][1] + dz*f[i][2];
    }
  }
  else {
    for (int i = 0; i < nlocal; i++) {
      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      double dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
      double dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
      double dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
      if (comm_forward>0) {
        deltaR[i][0] = dx;
        deltaR[i][1] = dy;
        deltaR[i][2] = dz;
      }
      fdr += ((dx*f[i][0])+(dy*f[i][1])+(dz*f[i][2]));
    }
  }

  double phiSum = 0.0;
  if (computeCv>-1) {
    comm->forward_comm_compute(this);
    int *type = atom->type;
    double** cutsq = force->pair->cutsq;
    if (force->pair) {
      double **x = atom->x;
      double **f = atom->f;
      int *type = atom->type;
      int nlocal = atom->nlocal;
      double *special_lj = force->special_lj;
      double *special_coul = force->special_coul;
      int newton_pair = force->newton_pair;

      if (update->firststep == update->ntimestep) neighbor->build_one(list,1);
      else neighbor->build_one(list);
      int inum = list->inum;
      int *ilist = list->ilist;
      int *numneigh = list->numneigh;
      int **firstneigh = list->firstneigh;

      for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        double fac = (newton_pair || i < nlocal) ? 1.0 : 0.5;
        double* ix = x[i];
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        double *idr = deltaR[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (!newton_pair && j>=nlocal) fac -= 0.5;
          double factor_lj = special_lj[sbmask(j)];
          double factor_coul = special_coul[sbmask(j)];
          j &= NEIGHMASK;
          double* jx = x[j];
          double delr[3];
          delr[0] = ix[0] - jx[0];
          delr[1] = ix[1] - jx[1];
          delr[2] = ix[2] - jx[2];
          double rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];
          int jtype = type[j];
          if (rsq < cutsq[itype][jtype]) {
            double* jdr = deltaR[j];
            double fforce, d2u[6];
            force->pair->single2(i, j, itype, jtype, rsq, delr, factor_coul, factor_lj, fforce, d2u);
            int m = 0;
            for (int k=0; k<3; k++) {
              double a = fac;
              for (int l=k; l<3; l++) {
                phiSum += a*(idr[k]*jdr[l]+jdr[k]*idr[l])*d2u[m];
                phiSum -= a*(idr[k]*idr[l]*d2u[m] + jdr[k]*jdr[l]*d2u[m]);
                m++;
                if (k==l) a *= 2;
              }
            }
          }
        }
      }
    }
  }

  // compute and sum up properties across processors

  double fdrTotal;
  MPI_Allreduce(&fdr,&fdrTotal,1,MPI_DOUBLE,MPI_SUM,world);
  double uTotal;
  if (computeU>-1 || computeCv>-1) {
    MPI_Allreduce(&u,&uTotal,1,MPI_DOUBLE,MPI_SUM,world);
    if (returnAnharmonic == -1) {
      uLat = uTotal;
    }
    if (computeU>-1) {
      if (returnAnharmonic) {
        vector[computeU] = uTotal - uLat + 0.5*fdrTotal;
      }
      else {
        vector[computeU] = uTotal + 0.5*fdrTotal + 0.5*dimension*(atom->natoms - 1)*force->boltz*finaltemp;
      }
    }
  }

  if (computeP>-1) {
    double fv = ((deltaPcap)-(force->boltz*finaltemp*force->nktv2p*atom->natoms/vol))/(force->boltz*finaltemp*dimension*(atom->natoms - 1));
    if (returnAnharmonic) {
      vector[computeP] = p - pLat + (fv*fdrTotal);
    }
    else {
      vector[computeP] = p + (fv*fdrTotal) + deltaPcap;
    }
  }

  if (computeCv>-1) {
    if (computeU==-1) MPI_Allreduce(&u,&uTotal,1,MPI_DOUBLE,MPI_SUM,world);
    double buTot;
    if (returnAnharmonic) {
      buTot = (uTotal - uLat + 0.5*fdrTotal)/finaltemp;
    }
    else {
      buTot = (uTotal + 0.5*fdrTotal)/finaltemp + 0.5*dimension*(atom->natoms - 1)*force->boltz;
    }
    double one = -0.25*(fdr + phiSum)/finaltemp;
    double Cv;
    MPI_Allreduce(&one,&Cv,1,MPI_DOUBLE,MPI_SUM,world);
    vector[computeCv] = Cv + buTot*buTot;
    if (!returnAnharmonic) {
      vector[computeCv] += 0.5*dimension*(atom->natoms-1);
    }
  }
  if (returnAnharmonic == -1) {
    // we have our lattice properties
    returnAnharmonic = 1;
  }
}

double ComputeHMA::virial_compute(int n)
{
  double v = 0;

  // sum contributions to virial from forces and fixes

  if (force->pair) v += sumVirial(n, force->pair->virial);
  if (force->bond) v += sumVirial(n, force->bond->virial);
  if (force->angle) v += sumVirial(n, force->angle->virial);
  if (force->dihedral) v += sumVirial(n, force->dihedral->virial);
  if (force->improper) v += sumVirial(n, force->improper->virial);
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->thermo_virial) v += sumVirial(n, modify->fix[i]->virial);

  // sum virial across procs

  double virial;
  MPI_Allreduce(&v,&virial,1,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (force->kspace)
    for (int i = 0; i < n; i++) virial += force->kspace->virial[i];
  return virial;
}

/* ---------------------------------------------------------------------- */

int ComputeHMA::pack_forward_comm(int n, int *list, double *buf,
                                        int pbc_flag, int *pbc)
{
  double **xoriginal = fix->astore;
  imageint *image = atom->image;
  double **x = atom->x;
  double *h = domain->h;
  double xprd = domain->xprd;   
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int m = 0;
  for (int ii = 0; ii < n; ii++) {
    int i = list[ii];
    buf[m++] = deltaR[i][0];
    buf[m++] = deltaR[i][1];
    buf[m++] = deltaR[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::unpack_forward_comm(int n, int first, double *buf)
{
  double **xoriginal = fix->astore;
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    deltaR[i][0] = buf[m++];
    deltaR[i][1] = buf[m++];
    deltaR[i][2] = buf[m++];
  }
}


/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeHMA::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

double ComputeHMA::memory_usage()
{
  double bytes = nmax * 3 * sizeof(double);
  return bytes;
}
