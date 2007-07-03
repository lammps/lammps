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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_gran_diag.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{NO_HISTORY,HISTORY,HERTZIAN};

/* ---------------------------------------------------------------------- */

FixGranDiag::FixGranDiag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix gran/diag command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix gran/diag command");
  first = 1;

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all("Fix gran/diag requires atom attributes radius, rmass, omega");

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    char *file = new char[128];

    sprintf(file,"%s.den",arg[4]);
    fpden = fopen(file,"w");
    if (fpden == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix gran/diag file %s",file);
      error->one(str);
    }

    sprintf(file,"%s.vel",arg[4]);
    fpvel = fopen(file,"w");
    if (fpvel == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix gran/diag file %s",file);
      error->one(str);
    }

    sprintf(file,"%s.str",arg[4]);
    fpstr = fopen(file,"w");
    if (fpstr == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix gran/diag file %s",file);
      error->one(str);
    }

    delete [] file;
  }

  step = atof(arg[5]);
  stepinv = 1.0/step;
  PI = 4.0*atan(1.0);

  maxlayers = 0;
}

/* ---------------------------------------------------------------------- */

FixGranDiag::~FixGranDiag()
{
  deallocate();

  if (me == 0) {
    fclose(fpden);
    fclose(fpvel);
    fclose(fpstr);
  }
}

/* ---------------------------------------------------------------------- */

int FixGranDiag::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::init()
{
  // insure use of granular pair_style
  // set local values from Pair values

  Pair *pair = force->pair_match("gran");
  if (pair == NULL)
    error->all("Fix gran/diag is incompatible with Pair style");
  int dampflag;
  double gamman;
  pair->extract_gran(&xkk,&gamman,&xmu,&dampflag);

  // same initialization as in pair_gran_history::init_style()

  xkkt = xkk * 2.0/7.0;
  dt = update->dt;
  double gammas = 0.5*gamman;
  if (dampflag == 0) gammas = 0.0;
  gamman_dl = gamman/dt;
  gammas_dl = gammas/dt;

  // set pairstyle from granular pair style

  if (force->pair_match("gran/no_history")) pairstyle = NO_HISTORY;
  else if (force->pair_match("gran/history")) pairstyle = HISTORY;
  else if (force->pair_match("gran/hertzian")) pairstyle = HERTZIAN;

  // check for freeze Fix and set freeze_group_bit

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::setup()
{
  if (first) end_of_step();
  first = 0;
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::end_of_step()
{
  int i,m;

  // setup for 2d vs 3d

  int dim;
  if (domain->dimension == 3) dim = 2;
  else dim = 1;

  // set bottom of box for binning purposes
  
  boxlo = domain->boxlo[dim];

  // update ghost atom info
  // else ghost x/v is out-of-date at end of timestep

  comm->communicate();

  // insure accumulator arrays are correct length
  // add 10 for buffer

  nlayers = static_cast<int> (domain->prd[dim] / step) + 10;
  if (nlayers > maxlayers) {
    deallocate();
    maxlayers = nlayers;
    allocate();
  }

  // zero all accumulators

  for (i = 0; i < nlayers; i++) {
    numdens[i] = 0;
    dendens[i] = 0.0;
    velx[i] = vely[i] = velz[i] = 0.0;
    velxx[i] = velyy[i] = velzz[i] = velxy[i] = velxz[i] = velyz[i] = 0.0;
    sigxx[i] = sigyy[i] = sigzz[i] = sigxy[i] = sigxz[i] = sigyz[i] = 0.0;
  }

  // density/velocity accumulation by atom
  // assign to layer based on atom distance from bottom of domain

  int overflow = 0;

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      m = static_cast<int> ((x[i][dim]-boxlo) * stepinv);
      if (m >= 0 && m < nlayers) {
	numdens[m]++;
	dendens[m] += rmass[i];
	velx[m] += v[i][0];
	vely[m] += v[i][1];
	velz[m] += v[i][2];
	velxx[m] += v[i][0]*v[i][0];
	velyy[m] += v[i][1]*v[i][1];
	velzz[m] += v[i][2]*v[i][2];
	velxy[m] += v[i][0]*v[i][1];
	velxz[m] += v[i][0]*v[i][2];
	velyz[m] += v[i][1]*v[i][2];
      } else overflow++;
    }

  // m = largest layer # with any counts
  // nmax = # of layers up to m
      
  for (m = nlayers-1; m >= 0; m--) if (numdens[m]) break;
  int nmax = m + 1;

  int tmp = nmax;
  MPI_Allreduce(&tmp,&nmax,1,MPI_INT,MPI_MAX,world);

  // overflow = total # of atoms out-of-bounds of layer arrays
  
  tmp = overflow;
  MPI_Allreduce(&tmp,&overflow,1,MPI_INT,MPI_SUM,world);

  // sum contributions across procs

  int *isum = new int[nmax];
  double *dsum = new double[nmax];
  
  MPI_Allreduce(numdens,isum,nmax,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nmax; i++) numdens[i] = isum[i];

  MPI_Allreduce(dendens,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) dendens[i] = dsum[i];
  
  MPI_Allreduce(velx,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velx[i] = dsum[i];
  MPI_Allreduce(vely,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) vely[i] = dsum[i];
  MPI_Allreduce(velz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velz[i] = dsum[i];
  
  MPI_Allreduce(velxx,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velxx[i] = dsum[i];
  MPI_Allreduce(velyy,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velyy[i] = dsum[i];
  MPI_Allreduce(velzz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velzz[i] = dsum[i];
  MPI_Allreduce(velxy,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velxy[i] = dsum[i];
  MPI_Allreduce(velxz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velxz[i] = dsum[i];
  MPI_Allreduce(velyz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) velyz[i] = dsum[i];

  // compute contribution to stress by every atom pair

  if (pairstyle == NO_HISTORY) stress_no_history();
  else if (pairstyle == HISTORY) stress_history();
  else if (pairstyle == HERTZIAN) stress_hertzian();

  // sum contributions across procs

  MPI_Allreduce(sigxx,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigxx[i] = dsum[i];
  MPI_Allreduce(sigyy,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigyy[i] = dsum[i];
  MPI_Allreduce(sigzz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigzz[i] = dsum[i];

  MPI_Allreduce(sigxy,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigxy[i] = dsum[i];
  MPI_Allreduce(sigxz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigxz[i] = dsum[i];
  MPI_Allreduce(sigyz,dsum,nmax,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmax; i++) sigyz[i] = dsum[i];

  delete [] isum;
  delete [] dsum;

  // density/velocity/stress by layer

  double velxxd1,velyyd1,velzzd1,velxyd1,velxzd1,velyzd1;
  double vollayer = domain->xprd;
  if (domain->dimension == 3) vollayer *= domain->yprd;

  for (m = 0; m < nmax; m++) {
    if (numdens[m] == 0) numdens[m] = 1;
    dendens[m] = stepinv*dendens[m]/vollayer;

    velx11[m] = velx[m]/numdens[m];
    vely11[m] = vely[m]/numdens[m];
    velz11[m] = velz[m]/numdens[m];
    velxx11[m] = velxx[m]/numdens[m];
    velyy11[m] = velyy[m]/numdens[m];
    velzz11[m] = velzz[m]/numdens[m];
    velxy11[m] = velxy[m]/numdens[m];
    velxz11[m] = velxz[m]/numdens[m];
    velyz11[m] = velyz[m]/numdens[m];

    velxxd1 = velxx11[m] - velx11[m]*velx11[m];
    velyyd1 = velyy11[m] - vely11[m]*vely11[m];
    velzzd1 = velzz11[m] - velz11[m]*velz11[m];
    velxyd1 = velxy11[m] - velx11[m]*vely11[m];
    velxzd1 = velxz11[m] - velx11[m]*velz11[m];
    velyzd1 = velyz11[m] - vely11[m]*velz11[m];

    velfxx[m] = velxxd1 * dendens[m];
    velfyy[m] = velyyd1 * dendens[m];
    velfzz[m] = velzzd1 * dendens[m];
    velfxy[m] = velxyd1 * dendens[m];
    velfxz[m] = velxzd1 * dendens[m];
    velfyz[m] = velyzd1 * dendens[m];

    sigx2[m] = sigxx[m]/(2.0*vollayer*step) + velxxd1*dendens[m];
    sigy2[m] = sigyy[m]/(2.0*vollayer*step) + velyyd1*dendens[m];
    sigz2[m] = sigzz[m]/(2.0*vollayer*step) + velzzd1*dendens[m];
    sigxy2[m] = sigxy[m]/(2.0*vollayer*step) + velxyd1*dendens[m];
    sigxz2[m] = sigxz[m]/(2.0*vollayer*step) + velxzd1*dendens[m];
    sigyz2[m] = sigyz[m]/(2.0*vollayer*step) + velyzd1*dendens[m];
  }

  // write out density profile

  if (me == 0) {
    fprintf(fpden,"ITEM: TIMESTEP\n");
    fprintf(fpden,"%d\n",update->ntimestep);
    fprintf(fpden,"ITEM: NUMBER OF LAYERS / OVERFLOWS\n");
    fprintf(fpden,"%d %d\n",nmax,overflow);
    fprintf(fpden,"ITEM: DENSITY BY LAYER\n");
    for (m = 0; m < nmax; m++)
      fprintf(fpden,"%d %g %g\n",m+1,(m+1)*step+boxlo,dendens[m]*PI/6.0);
  }

  // write out velocity profile

  if (me == 0) {
    fprintf(fpvel,"ITEM: TIMESTEP\n");
    fprintf(fpvel,"%d\n",update->ntimestep);
    fprintf(fpvel,"ITEM: NUMBER OF LAYERS / OVERFLOWS\n");
    fprintf(fpvel,"%d %d\n",nmax,overflow);
    fprintf(fpvel,"ITEM: VELOCITY BY LAYER\n");
    for (m = 0; m < nmax; m++)
      fprintf(fpvel,"%d %g %g %g %g %g %g %g %g %g %g\n",
	      m+1,(m+1)*step+boxlo,
	      velx11[m],vely11[m],velz11[m],
	      velxx11[m],velyy11[m],velzz11[m],
	      velxy11[m],velxz11[m],velyz11[m]);
  }

  // write out stress profile

  if (me == 0) {
    fprintf(fpstr,"ITEM: TIMESTEP\n");
    fprintf(fpstr,"%d\n",update->ntimestep);
    fprintf(fpstr,"ITEM: NUMBER OF LAYERS / OVERFLOWS\n");
    fprintf(fpstr,"%d %d\n",nmax,overflow);
    fprintf(fpstr,"ITEM: STRESS BY LAYER\n");
    for (m = 0; m < nmax; m++)
      fprintf(fpstr,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      m+1,(m+1)*step+boxlo,
              sigx2[m],sigy2[m],sigz2[m],
              sigxy2[m],sigxz2[m],sigyz2[m],
              velfxx[m],velfyy[m],velfzz[m],
              velfxy[m],velfxz[m],velfyz[m]);
  }
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::allocate()
{
  numdens = new int[maxlayers];
  dendens = new double[maxlayers];

  velx = new double[maxlayers];
  vely = new double[maxlayers];
  velz = new double[maxlayers];
  velxx = new double[maxlayers];
  velyy = new double[maxlayers];
  velzz = new double[maxlayers];
  velxy = new double[maxlayers];
  velxz = new double[maxlayers];
  velyz = new double[maxlayers];

  velx11 = new double[maxlayers];
  vely11 = new double[maxlayers];
  velz11 = new double[maxlayers];
  velxx11 = new double[maxlayers];
  velyy11 = new double[maxlayers];
  velzz11 = new double[maxlayers];
  velxy11 = new double[maxlayers];
  velxz11 = new double[maxlayers];
  velyz11 = new double[maxlayers];

  sigxx = new double[maxlayers];
  sigyy = new double[maxlayers];
  sigzz = new double[maxlayers];
  sigxy = new double[maxlayers];
  sigxz = new double[maxlayers];
  sigyz = new double[maxlayers];

  sigx2 = new double[maxlayers];
  sigy2 = new double[maxlayers];
  sigz2 = new double[maxlayers];
  sigxy2 = new double[maxlayers];
  sigxz2 = new double[maxlayers];
  sigyz2 = new double[maxlayers];

  velfxx = new double[maxlayers];
  velfyy = new double[maxlayers];
  velfzz = new double[maxlayers];
  velfxy = new double[maxlayers];
  velfxz = new double[maxlayers];
  velfyz = new double[maxlayers];
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::deallocate()
{
  if (maxlayers == 0) return;

  delete [] numdens;  delete [] dendens;

  delete [] velx;  delete [] vely;  delete [] velz;
  delete [] velxx; delete [] velyy; delete [] velzz;
  delete [] velxy; delete [] velxz; delete [] velyz;

  delete [] velx11;   delete [] vely11;   delete [] velz11;
  delete [] velxx11;  delete [] velyy11;  delete [] velzz11;
  delete [] velxy11;  delete [] velxz11;  delete [] velyz11;

  delete [] sigxx;  delete [] sigyy;  delete [] sigzz;
  delete [] sigxy;  delete [] sigxz;  delete [] sigyz;

  delete [] sigx2;   delete [] sigy2;   delete [] sigz2;
  delete [] sigxy2;  delete [] sigxz2;  delete [] sigyz2;

  delete [] velfxx;  delete [] velfyy;  delete [] velfzz;
  delete [] velfxy;  delete [] velfxz;  delete [] velfyz;
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::stress_no_history()
{
  int i,j,k,m,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double xmeff,damp,ccel,ccelx,ccely,ccelz;
  double fn,fs,ft,fs1,fs2,fs3;
  int *neighs;

  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // loop over all neighbors of my atoms
  // store stress for both atoms i and j

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      // skip if neither atom is in fix group

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq <  radsum*radsum) {

	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	vr1 *= dt;
	vr2 *= dt;
	vr3 *= dt;

	//  normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	//  relative rotational velocity

	wr1 = radi*omega[i][0] + radj*omega[j][0];
	wr2 = radi*omega[i][1] + radj*omega[j][1];
	wr3 = radi*omega[i][2] + radj*omega[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// force normalization

	fn = xmu * fabs(ccel*r);
	fs = xmeff*gammas_dl*vrel;
	if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
	else ft = 0.0;

	// shear friction forces

	fs1 = -ft*vtr1;
	fs2 = -ft*vtr2;
	fs3 = -ft*vtr3;

	// forces

	ccelx = delx*ccel + fs1;
	ccely = dely*ccel + fs2;
	ccelz = delz*ccel + fs3;

	// stress contribution of atom pair to z-layers
	// atom i always contributes
	// atom j contributes if newton_pair is on or if owned by this proc

	m = static_cast<int> ((x[i][dim]-boxlo) * stepinv);
	if (m >= 0 && m < nlayers) {
	  sigxx[m] += delx*ccelx;
	  sigyy[m] += dely*ccely;
	  sigzz[m] += delz*ccelz;
	  sigxy[m] += delx*ccely;
	  sigxz[m] += delx*ccelz;
	  sigyz[m] += dely*ccelz;
	}
	
	if (newton_pair || j < nlocal) {
	  m = static_cast<int> ((x[j][dim]-boxlo) * stepinv);
	  if (m >= 0 && m < nlayers) {
	    sigxx[m] += delx*ccelx;
	    sigyy[m] += dely*ccely;
	    sigzz[m] += delz*ccelz;
	    sigxy[m] += delx*ccely;
	    sigxz[m] += delx*ccelz;
	    sigyz[m] += dely*ccelz;
	  }
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::stress_history()
{
  int i,j,k,m,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel,shrx,shry,shrz;
  double xmeff,damp,ccel,ccelx,ccely,ccelz;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *neighs;
  double *firstshear,*shear;

  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // loop over all neighbors of my atoms
  // store stress on both atoms i and j

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    neighs = neighbor->firstneigh[i];
    firstshear = neighbor->firstshear[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      // skip if neither atom is in fix group

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq <  radsum*radsum) {

	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	vr1 *= dt;
	vr2 *= dt;
	vr3 *= dt;

	//  normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	//  relative rotational velocity

	wr1 = radi*omega[i][0] + radj*omega[j][0];
	wr2 = radi*omega[i][1] + radj*omega[j][1];
	wr3 = radi*omega[i][2] + radj*omega[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// shear history effects
	// shrmag = magnitude of shear
	// do not update shear history since not timestepping

	shear = &firstshear[3*k];
	shrx = shear[0] + vtr1;
	shry = shear[1] + vtr2;
	shrz = shear[2] + vtr3;
	shrmag = sqrt(shrx*shrx + shry*shry + shrz*shrz);

	// rotate shear displacements correctly

	rsht = shrx*delx + shry*dely + shrz*delz;
	rsht /= rsq;
        shrx -= rsht*delx;
        shry -= rsht*dely;
        shrz -= rsht*delz;

	// tangential forces

	fs1 = - (xkkt*shrx + xmeff*gammas_dl*vtr1);
	fs2 = - (xkkt*shry + xmeff*gammas_dl*vtr2);
	fs3 = - (xkkt*shrz + xmeff*gammas_dl*vtr3);

	// force normalization
	// rescale frictional displacements and forces if needed

	fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
	fn = xmu * fabs(ccel*r);

	if (fs > fn) {
	  if (shrmag != 0.0) {
	    fs1 *= fn/fs;
	    fs2 *= fn/fs;
	    fs3 *= fn/fs;
	  } else {
	    fs1 = 0.0;
	    fs2 = 0.0;
	    fs3 = 0.0;
	  }
	}

	// forces

	ccelx = delx*ccel + fs1;
	ccely = dely*ccel + fs2;
	ccelz = delz*ccel + fs3;

	// stress contribution of atom pair to z-layers
	// atom i always contributes
	// atom j contributes if newton_pair is on or if owned by this proc

	m = static_cast<int> ((x[i][dim]-boxlo) * stepinv);
	if (m >= 0 && m < nlayers) {
	  sigxx[m] += delx*ccelx;
	  sigyy[m] += dely*ccely;
	  sigzz[m] += delz*ccelz;
	  sigxy[m] += delx*ccely;
	  sigxz[m] += delx*ccelz;
	  sigyz[m] += dely*ccelz;
	}
	
	if (newton_pair || j < nlocal) {
	  m = static_cast<int> ((x[j][dim]-boxlo) * stepinv);
	  if (m >= 0 && m < nlayers) {
	    sigxx[m] += delx*ccelx;
	    sigyy[m] += dely*ccely;
	    sigzz[m] += delz*ccelz;
	    sigxy[m] += delx*ccely;
	    sigxz[m] += delx*ccelz;
	    sigyz[m] += dely*ccelz;
	  }
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGranDiag::stress_hertzian()
{
  int i,j,k,m,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel,shrx,shry,shrz;
  double xmeff,damp,ccel,ccelx,ccely,ccelz;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht,rhertz;
  int *neighs;
  double *firstshear,*shear;

  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // loop over all neighbors of my atoms
  // store stress on both atoms i and j

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    neighs = neighbor->firstneigh[i];
    firstshear = neighbor->firstshear[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      // skip if neither atom is in fix group

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq <  radsum*radsum) {

	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	vr1 *= dt;
	vr2 *= dt;
	vr3 *= dt;

	//  normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	//  relative rotational velocity

	wr1 = radi*omega[i][0] + radj*omega[j][0];
	wr2 = radi*omega[i][1] + radj*omega[j][1];
	wr3 = radi*omega[i][2] + radj*omega[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;
	rhertz = sqrt(radsum - r);
	ccel = rhertz * ccel;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// shear history effects
	// shrmag = magnitude of shear
	// do not update shear history since not timestepping

	shear = &firstshear[3*k];
	shrx = shear[0] + vtr1;
	shry = shear[1] + vtr2;
	shrz = shear[2] + vtr3;
	shrmag = sqrt(shrx*shrx + shry*shry + shrz*shrz);

	// rotate shear displacements correctly

	rsht = shrx*delx + shry*dely + shrz*delz;
	rsht /= rsq;
        shrx -= rsht*delx;
        shry -= rsht*dely;
        shrz -= rsht*delz;

	// tangential forces

        fs1 = -rhertz * (xkkt*shrx + xmeff*gammas_dl*vtr1);
        fs2 = -rhertz * (xkkt*shry + xmeff*gammas_dl*vtr2);
        fs3 = -rhertz * (xkkt*shrz + xmeff*gammas_dl*vtr3);

	// force normalization
	// rescale frictional displacements and forces if needed

	fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
	fn = xmu * fabs(ccel*r);

	if (fs > fn) {
	  if (shrmag != 0.0) {
	    fs1 *= fn/fs;
	    fs2 *= fn/fs;
	    fs3 *= fn/fs;
	  } else {
	    fs1 = 0.0;
	    fs2 = 0.0;
	    fs3 = 0.0;
	  }
	}

	// forces

	ccelx = delx*ccel + fs1;
	ccely = dely*ccel + fs2;
	ccelz = delz*ccel + fs3;

	// stress contribution of atom pair to z-layers
	// atom i always contributes
	// atom j contributes if newton_pair is on or if owned by this proc

	m = static_cast<int> ((x[i][dim]-boxlo) * stepinv);
	if (m >= 0 && m < nlayers) {
	  sigxx[m] += delx*ccelx;
	  sigyy[m] += dely*ccely;
	  sigzz[m] += delz*ccelz;
	  sigxy[m] += delx*ccely;
	  sigxz[m] += delx*ccelz;
	  sigyz[m] += dely*ccelz;
	}
	
	if (newton_pair || j < nlocal) {
	  m = static_cast<int> ((x[j][dim]-boxlo) * stepinv);
	  if (m >= 0 && m < nlayers) {
	    sigxx[m] += delx*ccelx;
	    sigyy[m] += dely*ccely;
	    sigzz[m] += delz*ccelz;
	    sigxy[m] += delx*ccely;
	    sigxz[m] += delx*ccelz;
	    sigyz[m] += dely*ccelz;
	  }
	}
      }
    }
  }
}
