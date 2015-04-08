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
   Contributing authors: Amit Kumar and Michael Bybee (UIUC)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_brownian.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairBrownian::PairBrownian(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairBrownian::~PairBrownian()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
  }
  delete random;
}

/* ---------------------------------------------------------------------- */

void PairBrownian::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  double randr;
  double prethermostat;
  double xl[3],a_sq,a_sh,a_pu,Fbmag;
  double p1[3],p2[3],p3[3];
  int overlaps = 0;

  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2){ // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)){
        double wallhi[3], walllo[3];
        for (int j = 0; j < 3; j++){
          wallhi[j] = domain->prd[j];
          walllo[j] = 0;
        }
        for (int m = 0; m < wallfix->nwall; m++){
          int dim = wallfix->wallwhich[m] / 2;
          int side = wallfix->wallwhich[m] % 2;
          if (wallfix->xstyle[m] == VARIABLE){
            wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
          }
          else wallcoord = wallfix->coord0[m];
          if (side == 0) walllo[dim] = wallcoord;
          else wallhi[dim] = wallcoord;
        }
        for (int j = 0; j < 3; j++)
          dims[j] = wallhi[j] - walllo[j];
      }
      double vol_T = dims[0]*dims[1]*dims[2];
      double vol_f = vol_P/vol_T;
      if (flaglog == 0) {
        R0  = 6*MY_PI*mu*rad*(1.0 + 2.16*vol_f);
        RT0 = 8*MY_PI*mu*cube(rad);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        R0  = 6*MY_PI*mu*rad*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*MY_PI*mu*cube(rad)*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }

  // scale factor for Brownian moments

  prethermostat = sqrt(24.0*force->boltz*t_target/update->dt);
  prethermostat *= sqrt(force->vxmu2f/force->ftm2v/force->mvv2e);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // FLD contribution to force and torque due to isotropic terms

    if (flagfld) {
      f[i][0] += prethermostat*sqrt(R0)*(random->uniform()-0.5);
      f[i][1] += prethermostat*sqrt(R0)*(random->uniform()-0.5);
      f[i][2] += prethermostat*sqrt(R0)*(random->uniform()-0.5);
      if (flaglog) {
        torque[i][0] += prethermostat*sqrt(RT0)*(random->uniform()-0.5);
        torque[i][1] += prethermostat*sqrt(RT0)*(random->uniform()-0.5);
        torque[i][2] += prethermostat*sqrt(RT0)*(random->uniform()-0.5);
      }
    }

    if (!flagHI) continue;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // scalar resistances a_sq and a_sh

        h_sep = r - 2.0*radi;

        // check for overlaps

        if (h_sep < 0.0) overlaps++;

        // if less than minimum gap, use minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;

        // scale h_sep by radi

        h_sep = h_sep/radi;

        // scalar resistances

        if (flaglog) {
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1.0/h_sep));
          a_sh = 6.0*MY_PI*mu*radi*(1.0/6.0*log(1.0/h_sep));
          a_pu = 8.0*MY_PI*mu*cube(radi)*(3.0/160.0*log(1.0/h_sep));
        } else
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep);

        // generate the Pairwise Brownian Force: a_sq

        Fbmag = prethermostat*sqrt(a_sq);

        // generate a random number

        randr = random->uniform()-0.5;

        // contribution due to Brownian motion

        fx = Fbmag*randr*delx/r;
        fy = Fbmag*randr*dely/r;
        fz = Fbmag*randr*delz/r;

        // add terms due to a_sh

        if (flaglog) {

          // generate two orthogonal vectors to the line of centers

          p1[0] = delx/r; p1[1] = dely/r; p1[2] = delz/r;
          set_3_orthogonal_vectors(p1,p2,p3);

          // magnitude

          Fbmag = prethermostat*sqrt(a_sh);

          // force in each of the two directions

          randr = random->uniform()-0.5;
          fx += Fbmag*randr*p2[0];
          fy += Fbmag*randr*p2[1];
          fz += Fbmag*randr*p2[2];

          randr = random->uniform()-0.5;
          fx += Fbmag*randr*p3[0];
          fy += Fbmag*randr*p3[1];
          fz += Fbmag*randr*p3[2];
        }

        // scale forces to appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;

        // sum to total force

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;

        if (newton_pair || j < nlocal) {
          //randr = random->uniform()-0.5;
          //fx = Fbmag*randr*delx/r;
          //fy = Fbmag*randr*dely/r;
          //fz = Fbmag*randr*delz/r;

          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;
        }

        // torque due to the Brownian Force

        if (flaglog) {

          // location of the point of closest approach on I from its center

          xl[0] = -delx/r*radi;
          xl[1] = -dely/r*radi;
          xl[2] = -delz/r*radi;

          // torque = xl_cross_F

          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;

          // torque is same on both particles

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;

          if (newton_pair || j < nlocal) {
            torque[j][0] -= tx;
            torque[j][1] -= ty;
            torque[j][2] -= tz;
          }

          // torque due to a_pu

          Fbmag = prethermostat*sqrt(a_pu);

          // force in each direction

          randr = random->uniform()-0.5;
          tx = Fbmag*randr*p2[0];
          ty = Fbmag*randr*p2[1];
          tz = Fbmag*randr*p2[2];

          randr = random->uniform()-0.5;
          tx += Fbmag*randr*p3[0];
          ty += Fbmag*randr*p3[1];
          tz += Fbmag*randr*p3[2];

          // torque has opposite sign on two particles

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;

          if (newton_pair || j < nlocal) {
            torque[j][0] += tx;
            torque[j][1] += ty;
            torque[j][2] += tz;
          }
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }

  int print_overlaps = 0;
  if (print_overlaps && overlaps)
    printf("Number of overlaps=%d\n",overlaps);

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBrownian::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBrownian::settings(int narg, char **arg)
{
  if (narg != 7 && narg != 9) error->all(FLERR,"Illegal pair_style command");

  mu = force->numeric(FLERR,arg[0]);
  flaglog = force->inumeric(FLERR,arg[1]);
  flagfld = force->inumeric(FLERR,arg[2]);
  cut_inner_global = force->numeric(FLERR,arg[3]);
  cut_global = force->numeric(FLERR,arg[4]);
  t_target = force->numeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);

  flagHI = flagVF = 1;
  if (narg == 9) {
    flagHI = force->inumeric(FLERR,arg[7]);
    flagVF = force->inumeric(FLERR,arg[8]);
  }

  if (flaglog == 1 && flagHI == 0) {
    error->warning(FLERR,"Cannot include log terms without 1/r terms; "
                   "setting flagHI to 1");
    flagHI = 1;
  }

  // initialize Marsaglia RNG with processor-unique seed

  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBrownian::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;

  if (narg == 4) {
    cut_inner_one = force->numeric(FLERR,arg[2]);
    cut_one = force->numeric(FLERR,arg[3]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBrownian::init_style()
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair brownian requires atom style sphere");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR,
                   "Pair brownian needs newton pair on for "
                   "momentum conservation");

  neighbor->request(this,instance_me);

  // insure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair brownian requires extended particles");

  // require monodisperse system with same radii for all types

  double radtype;
  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->radius_consistency(i,radtype))
      error->all(FLERR,"Pair brownian requires monodisperse particles");
    if (i > 1 && radtype != rad)
      error->all(FLERR,"Pair brownian requires monodisperse particles");
    rad = radtype;
  }

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume
  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.
  //
  // If available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.

  flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0)
      flagdeform = 1;
    else if (strstr(modify->fix[i]->style,"wall") != NULL) {
      if (flagwall) 
        error->all(FLERR,
                   "Cannot use multiple fix wall commands with pair brownian");
      flagwall = 1; // Walls exist
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }
  }

  // set the isotropic constants depending on the volume fraction
  // vol_T = total volumeshearing = flagdeform = flagwall = 0;

  double vol_T, wallcoord;
  if (!flagwall) vol_T = domain->xprd*domain->yprd*domain->zprd;
  else {
    double wallhi[3], walllo[3];
    for (int j = 0; j < 3; j++){
      wallhi[j] = domain->prd[j];
      walllo[j] = 0;
    }
    for (int m = 0; m < wallfix->nwall; m++){
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      if (wallfix->xstyle[m] == VARIABLE){
        wallfix->xindex[m] = input->variable->find(wallfix->xstr[m]);
        // Since fix->wall->init happens after pair->init_style
        wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
      }

      else wallcoord = wallfix->coord0[m];

      if (side == 0) walllo[dim] = wallcoord;
      else wallhi[dim] = wallcoord;
    }
    vol_T = (wallhi[0] - walllo[0]) * (wallhi[1] - walllo[1]) *
      (wallhi[2] - walllo[2]);
  }

  // vol_P = volume of particles, assuming mono-dispersity
  // vol_f = volume fraction

  vol_P = atom->natoms*(4.0/3.0)*MY_PI*cube(rad);

  double vol_f = vol_P/vol_T;

  // set isotropic constants
  if (!flagVF) vol_f = 0;

  if (flaglog == 0) {
    R0  = 6*MY_PI*mu*rad*(1.0 + 2.16*vol_f);
    RT0 = 8*MY_PI*mu*cube(rad);  // not actually needed
  } else {
    R0  = 6*MY_PI*mu*rad*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*cube(rad)*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBrownian::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBrownian::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBrownian::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cut_inner[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBrownian::write_restart_settings(FILE *fp)
{
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&flaglog,sizeof(int),1,fp);
  fwrite(&flagfld,sizeof(int),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&t_target,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&flagHI,sizeof(int),1,fp);
  fwrite(&flagVF,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBrownian::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&mu,sizeof(double),1,fp);
    fread(&flaglog,sizeof(int),1,fp);
    fread(&flagfld,sizeof(int),1,fp);
    fread(&cut_inner_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&t_target, sizeof(double),1,fp);
    fread(&seed, sizeof(int),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&flagHI,sizeof(int),1,fp);
    fread(&flagVF,sizeof(int),1,fp);
  }
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&flaglog,1,MPI_INT,0,world);
  MPI_Bcast(&flagfld,1,MPI_INT,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_target,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&flagHI,1,MPI_INT,0,world);
  MPI_Bcast(&flagVF,1,MPI_INT,0,world);

  // additional setup based on restart parameters

  delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------*/

void PairBrownian::set_3_orthogonal_vectors(double p1[3],
                                            double p2[3], double p3[3])
{
  double norm;
  int ix,iy,iz;

  // find the index of maximum magnitude and store it in iz

  if (fabs(p1[0]) > fabs(p1[1])) {
    iz=0;
    ix=1;
    iy=2;
  } else {
    iz=1;
    ix=2;
    iy=0;
  }

  if (iz==0) {
    if (fabs(p1[0]) < fabs(p1[2])) {
      iz = 2;
      ix = 0;
      iy = 1;
    }
  } else {
    if (fabs(p1[1]) < fabs(p1[2])) {
      iz = 2;
      ix = 0;
      iy = 1;
    }
  }

  // set p2 arbitrarily such that it's orthogonal to p1

  p2[ix]=1.0;
  p2[iy]=1.0;
  p2[iz] = -(p1[ix]*p2[ix] + p1[iy]*p2[iy])/p1[iz];

  // normalize p2

  norm = sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]);

  p2[0] = p2[0]/norm;
  p2[1] = p2[1]/norm;
  p2[2] = p2[2]/norm;

  // Set p3 by taking the cross product p3=p2xp1

  p3[0] = p1[1]*p2[2] - p1[2]*p2[1];
  p3[1] = p1[2]*p2[0] - p1[0]*p2[2];
  p3[2] = p1[0]*p2[1] - p1[1]*p2[0];
}
