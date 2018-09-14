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
   Contributing authors: Julien Tranchida (SNL)
   			 Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "ewald_dipole.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

EwaldDipole::EwaldDipole(LAMMPS *lmp, int narg, char **arg) : Ewald(lmp, narg, arg)
{
  ewaldflag = dipoleflag = 1;
  group_group_enable = 0;
  muk = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldDipole::~EwaldDipole()
{
  memory->destroy(muk);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void EwaldDipole::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"EwaldDipole initialization ...\n");
    if (logfile) fprintf(logfile,"EwaldDipole initialization ...\n");
  }

  // error check
  
  dipoleflag = atom->mu?1:0;
  qsum_qsq(0); // q[i] might not be declared ?

  if (dipoleflag && q2)
    error->all(FLERR,"Cannot (yet) use charges with Kspace style EwaldDipole");

  triclinic_check();
  
  // no triclinic ewald dipole (yet)
  
  triclinic = domain->triclinic;
  if (triclinic)
    error->all(FLERR,"Cannot (yet) use EwaldDipole with triclinic box");
  
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use EwaldDipole with 2d simulation");

  if (!atom->mu) error->all(FLERR,"Kspace style requires atom attribute mu");
//if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

  if (dipoleflag && strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");
  
  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with EwaldDipole");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab EwaldDipole");
  }

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  if (triclinic)
    error->all(FLERR,"Cannot yet use triclinic cells with EwaldDipole");

  pair_check();

  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  double cutoff = *p_cutoff;

  // kspace TIP4P not yet supported
  // qdist = offset only for TIP4P fictitious charge

  //qdist = 0.0;
  if (tip4pflag)
    error->all(FLERR,"Cannot yet use TIP4P with EwaldDipole");

  // compute musum & musqsum and warn if no dipole

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  musum_musq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup K-space resolution

  bigint natoms = atom->natoms;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab EwaldDipole
  // 3d EwaldDipole just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  if (!gewaldflag) {
    if (accuracy <= 0.0)
      error->all(FLERR,"KSpace accuracy must be > 0");
    if (q2 == 0.0)
      error->all(FLERR,"Must use 'kspace_modify gewald' for uncharged system");
    g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
    if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
    else g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // setup EwaldDipole coefficients so can print stats

  setup();

  // final RMS accuracy

  double lprx = rms(kxmax_orig,xprd,natoms,q2);
  double lpry = rms(kymax_orig,yprd,natoms,q2);
  double lprz = rms(kzmax_orig,zprd_slab,natoms,q2);
  double lpr = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  double q2_over_sqrt = q2 / sqrt(natoms*cutoff*xprd*yprd*zprd_slab);
  double spr = 2.0 *q2_over_sqrt * exp(-g_ewald*g_ewald*cutoff*cutoff);
  double tpr = estimate_table_accuracy(q2_over_sqrt,spr);
  double estimated_accuracy = sqrt(lpr*lpr + spr*spr + tpr*tpr);

  // stats

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  G vector (1/distance) = %g\n",g_ewald);
      fprintf(screen,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(screen,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(screen,"  KSpace vectors: actual max1d max3d = %d %d %d\n",
              kcount,kmax,kmax3d);
      fprintf(screen,"                  kxmax kymax kzmax  = %d %d %d\n",
              kxmax,kymax,kzmax);
    }
    if (logfile) {
      fprintf(logfile,"  G vector (1/distance) = %g\n",g_ewald);
      fprintf(logfile,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(logfile,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(logfile,"  KSpace vectors: actual max1d max3d = %d %d %d\n",
              kcount,kmax,kmax3d);
      fprintf(logfile,"                  kxmax kymax kzmax  = %d %d %d\n",
              kxmax,kymax,kzmax);
    }
  }
}

/* ----------------------------------------------------------------------
   adjust EwaldDipole coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void EwaldDipole::setup()
{
  // volume-dependent factors

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // adjustment of z dimension for 2d slab EwaldDipole
  // 3d EwaldDipole just uses zprd since slab_volfactor = 1.0

  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  unitk[0] = 2.0*MY_PI/xprd;
  unitk[1] = 2.0*MY_PI/yprd;
  unitk[2] = 2.0*MY_PI/zprd_slab;

  int kmax_old = kmax;

  if (kewaldflag == 0) {

    // determine kmax
    // function of current box size, accuracy, G_ewald (short-range cutoff)

    bigint natoms = atom->natoms;
    double err;
    kxmax = 1;
    kymax = 1;
    kzmax = 1;	

    // set kmax in 3 directions to respect accuracy

    err = rms_dipole(kxmax,xprd,natoms);
    while (err > accuracy) {
      kxmax++;
      err = rms_dipole(kxmax,xprd,natoms);
    }

    err = rms_dipole(kxmax,xprd,natoms);
    while (err > accuracy) {
      kymax++;
      err = rms_dipole(kxmax,xprd,natoms);
    }

    err = rms_dipole(kxmax,xprd,natoms);
    while (err > accuracy) {
      kzmax++;
      err = rms_dipole(kxmax,xprd,natoms);
    }

    kmax = MAX(kxmax,kymax);
    kmax = MAX(kmax,kzmax);
    kmax3d = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax;

    double gsqxmx = unitk[0]*unitk[0]*kxmax*kxmax;
    double gsqymx = unitk[1]*unitk[1]*kymax*kymax;
    double gsqzmx = unitk[2]*unitk[2]*kzmax*kzmax;
    gsqmx = MAX(gsqxmx,gsqymx);
    gsqmx = MAX(gsqmx,gsqzmx);

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

  } else {

    kxmax = kx_ewald;
    kymax = ky_ewald;
    kzmax = kz_ewald;

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

    kmax = MAX(kxmax,kymax);
    kmax = MAX(kmax,kzmax);
    kmax3d = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax;

    double gsqxmx = unitk[0]*unitk[0]*kxmax*kxmax;
    double gsqymx = unitk[1]*unitk[1]*kymax*kymax;
    double gsqzmx = unitk[2]*unitk[2]*kzmax*kzmax;
    gsqmx = MAX(gsqxmx,gsqymx);
    gsqmx = MAX(gsqmx,gsqzmx);
  }

  gsqmx *= 1.00001;

  // if size has grown, reallocate k-dependent and nlocal-dependent arrays

  if (kmax > kmax_old) {
    deallocate();
    allocate();
    group_allocate_flag = 0;

    memory->destroy(ek);
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    memory->destroy(muk);
    nmax = atom->nmax;
    memory->create(ek,nmax,3,"ewald:ek");
    memory->create3d_offset(cs,-kmax,kmax,3,nmax,"ewald:cs");
    memory->create3d_offset(sn,-kmax,kmax,3,nmax,"ewald:sn");
    memory->create(muk,kmax3d,nmax,"ewald:muk");
    kmax_created = kmax;
  }

  // pre-compute EwaldDipole coefficients

  coeffs();
}

/* ----------------------------------------------------------------------
   compute dipole RMS accuracy for a dimension
------------------------------------------------------------------------- */

double EwaldDipole::rms_dipole(int km, double prd, bigint natoms)
{
  if (natoms == 0) natoms = 1;   // avoid division by zero

  // error from eq.(46), Wang et al., JCP 115, 6351 (2001)
  
  double value = 8*MY_PI*mu2*g_ewald/volume *
    sqrt(2*MY_PI*km*km*km/(15.0*natoms)) *
    exp(-MY_PI*MY_PI*km*km/(g_ewald*g_ewald*prd*prd));

  return value;
}

/* ----------------------------------------------------------------------
   compute the EwaldDipole long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldDipole::compute(int eflag, int vflag)
{
  int i,j,k;

  // set energy/virial flags

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    musum_musq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ek);
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    memory->destroy(muk);
    nmax = atom->nmax;
    memory->create(ek,nmax,3,"ewald:ek");
    memory->create3d_offset(cs,-kmax,kmax,3,nmax,"ewald:cs");
    memory->create3d_offset(sn,-kmax,kmax,3,nmax,"ewald:sn");
    memory->create(muk,kmax3d,nmax,"ewald:muk");
    kmax_created = kmax;
  }

  // partial structure factors on each processor
  // total structure factor by summing over procs

  eik_dot_r();

  MPI_Allreduce(sfacrl,sfacrl_all,kcount,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sfacim,sfacim_all,kcount,MPI_DOUBLE,MPI_SUM,world);

  // K-space portion of electric field
  // double loop over K-vectors and local atoms
  // perform per-atom calculations if needed

  double **f = atom->f;
  //double *q = atom->q;
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  int kx,ky,kz;
  double cypz,sypz,exprl,expim,partial,partial_peratom;

  for (i = 0; i < nlocal; i++) {
    ek[i][0] = 0.0;
    ek[i][1] = 0.0;
    ek[i][2] = 0.0;
  }

  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];

    for (i = 0; i < nlocal; i++) {

      // calculating  exp(i*k*ri)

      cypz = cs[ky][1][i]*cs[kz][2][i] - sn[ky][1][i]*sn[kz][2][i];
      sypz = sn[ky][1][i]*cs[kz][2][i] + cs[ky][1][i]*sn[kz][2][i];
      exprl = cs[kx][0][i]*cypz - sn[kx][0][i]*sypz;
      expim = sn[kx][0][i]*cypz + cs[kx][0][i]*sypz;

      // taking im-part of struct_fact x exp(i*k*ri) (for force calc.)

      partial = expim*sfacrl_all[k] - exprl*sfacim_all[k];
      //ek[i][0] += partial*eg[k][0];
      //ek[i][1] += partial*eg[k][1];
      //ek[i][2] += partial*eg[k][2];
      ek[i][0] += kx*partial*eg[k][0];
      ek[i][1] += ky*partial*eg[k][1];
      ek[i][2] += kz*partial*eg[k][2];

      if (evflag_atom) {

	// taking re-part of struct_fact x exp(i*k*ri) (for energy calc.)

        partial_peratom = exprl*sfacrl_all[k] + expim*sfacim_all[k];
        //if (eflag_atom) eatom[i] += q[i]*ug[k]*partial_peratom;
        if (eflag_atom) eatom[i] += muk[k][i]*ug[k]*partial_peratom;
        if (vflag_atom)
          for (j = 0; j < 6; j++)
	    // to be done
            vatom[i][j] += ug[k]*vg[k][j]*partial_peratom;
      }
    }
  }

  // convert E-field to force

  //const double qscale = qqrd2e * scale;
  const double muscale = qqrd2e * scale;

  for (i = 0; i < nlocal; i++) {
    //f[i][0] += qscale * q[i]*ek[i][0];
    //f[i][1] += qscale * q[i]*ek[i][1];
    //if (slabflag != 2) f[i][2] += qscale * q[i]*ek[i][2];
    f[i][0] += muscale * mu[i][0] * ek[i][0];
    f[i][1] += muscale * mu[i][1] * ek[i][1];
    if (slabflag != 2) f[i][2] += muscale * mu[i][2] * ek[i][2];
  }

  // sum global energy across Kspace vevs and add in volume-dependent term

  if (eflag_global) {
    for (k = 0; k < kcount; k++) {

      // taking the re-part of struct_fact_i x struct_fact_j

      energy += ug[k] * (sfacrl_all[k]*sfacrl_all[k] +
                         sfacim_all[k]*sfacim_all[k]);
    }

    // substracting self energy and scaling

    energy -= musqsum*2.0*g3/3.0/MY_PIS;
    energy *= muscale;
  }

  // global virial

  if (vflag_global) {
    double uk;
    for (k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
      for (j = 0; j < 6; j++) virial[j] += uk*vg[k][j];
    }
    //for (j = 0; j < 6; j++) virial[j] *= qscale;
    for (j = 0; j < 6; j++) virial[j] *= muscale;
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] -= (mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2])
	  *2.0*g3/3.0/MY_PIS;
        eatom[i] *= muscale;
      }
    }

    if (vflag_atom)
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= q[i]*qscale;
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();
}

/* ----------------------------------------------------------------------
   compute the 
------------------------------------------------------------------------- */

void EwaldDipole::eik_dot_r()
{
  int i,k,l,m,n,ic;
  double cstr1,sstr1,cstr2,sstr2,cstr3,sstr3,cstr4,sstr4;
  double sqk,clpm,slpm;
  double mux, muy, muz;

  double **x = atom->x;
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  n = 0;
  mux = muy = muz = 0.0;

  // loop on different k-directions
  // loop on n kpoints and nlocal atoms
  // store (n x nlocal) tab. of values of (mu_i dot k)
  // store n values of sum_j[ (mu_j dot k) exp(-k dot r_j) ]

  // (k,0,0), (0,l,0), (0,0,m)
  
  // loop 1: k=1, l=1, m=1
  // define first val. of cos and sin

  for (ic = 0; ic < 3; ic++) {
    sqk = (unitk[ic] * unitk[ic]);
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
        sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
        cs[-1][ic][i] = cs[1][0][i];
        sn[-1][ic][i] = -sn[1][0][i];
        muk[n][i] = (mu[i][ic]*unitk[ic]);
        cstr1 += muk[n][i]*cs[1][ic][i];
        sstr1 += muk[n][i]*sn[1][ic][i];
      }
      sfacrl[n] = cstr1;
      sfacim[n++] = sstr1;
    }
  }

  // loop 2: k>1, l>1, m>1

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m*unitk[ic] * m*unitk[ic];
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cs[m][ic][i] = cs[m-1][ic][i]*cs[1][ic][i] -
            sn[m-1][ic][i]*sn[1][ic][i];
          sn[m][ic][i] = sn[m-1][ic][i]*cs[1][ic][i] +
            cs[m-1][ic][i]*sn[1][ic][i];
          cs[-m][ic][i] = cs[m][ic][i];
          sn[-m][ic][i] = -sn[m][ic][i];
	  muk[n][i] = (mu[i][ic]*m*unitk[ic]);
          cstr1 += muk[n][i]*cs[1][ic][i];
          sstr1 += muk[n][i]*sn[1][ic][i];
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
	  mux = mu[i][0];
	  muy = mu[i][1];

	  // dir 1: (k,l,0)
	  muk[n][i] = (mux*k*unitk[0] + muy*l*unitk[1]);
          cstr1 += muk[n][i]*(cs[k][0][i]*cs[l][1][i]-sn[k][0][i]*sn[l][1][i]);
          sstr1 += muk[n][i]*(sn[k][0][i]*cs[l][1][i]+cs[k][0][i]*sn[l][1][i]);
	  
	  // dir 2: (k,-l,0)
	  muk[n+1][i] = (mux*k*unitk[0] - muy*l*unitk[1]);
          cstr2 += muk[n+1][i]*(cs[k][0][i]*cs[l][1][i]+sn[k][0][i]*sn[l][1][i]);
          sstr2 += muk[n+1][i]*(sn[k][0][i]*cs[l][1][i]-cs[k][0][i]*sn[l][1][i]);
	}
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (l*unitk[1] * l*unitk[1]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
	  muy = mu[i][1];
	  muz = mu[i][2];

	  // dir 1: (0,l,m)
	  muk[n][i] = (muy*l*unitk[1] + muz*m*unitk[2]); 
          cstr1 += muk[n][i]*(cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i]);
          sstr1 += muk[n][i]*(sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i]);
	  
	  // dir 2: (0,l,-m)
	  muk[n+1][i] = (muy*l*unitk[1] - muz*m*unitk[2]); 
          cstr2 += muk[n+1][i]*(cs[l][1][i]*cs[m][2][i]+sn[l][1][i]*sn[m][2][i]);
          sstr2 += muk[n+1][i]*(sn[l][1][i]*cs[m][2][i]-cs[l][1][i]*sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }
  
  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
	  mux = mu[i][0];
	  muz = mu[i][2];

	  // dir 1: (k,0,m)
	  muk[n][i] = (mux*k*unitk[0] + muz*m*unitk[2]); 
          cstr1 += muk[n][i]*(cs[k][0][i]*cs[m][2][i]-sn[k][0][i]*sn[m][2][i]);
          sstr1 += muk[n][i]*(sn[k][0][i]*cs[m][2][i]+cs[k][0][i]*sn[m][2][i]);
	  
	  // dir 2: (k,0,-m)
	  muk[n+1][i] = (mux*k*unitk[0] - muz*m*unitk[2]); 
          cstr2 += muk[n+1][i]*(cs[k][0][i]*cs[m][2][i]+sn[k][0][i]*sn[m][2][i]);
          sstr2 += muk[n+1][i]*(sn[k][0][i]*cs[m][2][i]-cs[k][0][i]*sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]) +
          (m*unitk[2] * m*unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          cstr3 = 0.0;
          sstr3 = 0.0;
          cstr4 = 0.0;
          sstr4 = 0.0;
          for (i = 0; i < nlocal; i++) {
	    mux = mu[i][0];
	    muy = mu[i][1];
	    muz = mu[i][2];

	    // dir 1: (k,l,m)
	    muk[n][i] = (mux*k*unitk[0] + muy*l*unitk[1] + muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr1 += muk[n][i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr1 += muk[n][i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 2: (k,-l,m)
	    muk[n+1][i] = (mux*k*unitk[0] - muy*l*unitk[1] + muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr2 += muk[n+1][i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr2 += muk[n+1][i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 3: (k,l,-m)
	    muk[n+2][i] = (mux*k*unitk[0] + muy*l*unitk[1] - muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr3 += muk[n+2][i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr3 += muk[n+2][i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 4: (k,-l,-m)
	    muk[n+3][i] = (mux*k*unitk[0] - muy*l*unitk[1] - muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr4 += muk[n+3][i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr4 += muk[n+3][i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
          }
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
          sfacrl[n] = cstr2;
          sfacim[n++] = sstr2;
          sfacrl[n] = cstr3;
          sfacim[n++] = sstr3;
          sfacrl[n] = cstr4;
          sfacim[n++] = sstr4;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D EwaldDipole if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void EwaldDipole::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd*zprd/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd*zprd/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += ffact * q[i]*(dipole_all - qsum*x[i][2]);
}

/* ----------------------------------------------------------------------
   compute musum,musqsum,mu2
   called initially, when particle count changes, when dipoles are changed
------------------------------------------------------------------------- */

void EwaldDipole::musum_musq()
{
  const int nlocal = atom->nlocal;

  musum = musqsum = mu2 = 0.0;
  if (atom->mu_flag) {
    double** mu = atom->mu;
    double musum_local(0.0), musqsum_local(0.0);

    for (int i = 0; i < nlocal; i++) {
      musum_local += mu[i][0] + mu[i][1] + mu[i][2];
      musqsum_local += mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
    }

    MPI_Allreduce(&musum_local,&musum,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&musqsum_local,&musqsum,1,MPI_DOUBLE,MPI_SUM,world);

    mu2 = musqsum * force->qqrd2e;
  }

  if (mu2 == 0 && comm->me == 0)
    error->all(FLERR,"Using kspace solver PPPMDipole on system with no dipoles");
}
