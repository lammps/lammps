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
   Contributing authors: Julien Tranchida (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "ewald_dipole.h"
#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

EwaldDipole::EwaldDipole(LAMMPS *lmp) : Ewald(lmp),
  tk(NULL), vc(NULL)
{
  ewaldflag = dipoleflag = 1;
  group_group_enable = 0;
  tk = NULL;
  vc = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldDipole::~EwaldDipole()
{
  memory->destroy(tk);
  memory->destroy(vc);
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
    
    // initial guess with old method
    
    g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*mu2);
    if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
    else g_ewald = sqrt(-log(g_ewald)) / cutoff;
    
    // try Newton solver

    double g_ewald_new =
      NewtonSolve(g_ewald,cutoff,natoms,xprd*yprd*zprd,mu2);
    if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
    else error->warning(FLERR,"Ewald/disp Newton solver failed, "
                        "using old method to estimate g_ewald");  
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

    err = rms_dipole(kymax,yprd,natoms);
    while (err > accuracy) {
      kymax++;
      err = rms_dipole(kymax,yprd,natoms);
    }

    err = rms_dipole(kzmax,zprd,natoms);
    while (err > accuracy) {
      kzmax++;
      err = rms_dipole(kzmax,zprd,natoms);
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
    memory->destroy(tk);
    memory->destroy(vc);
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    nmax = atom->nmax;
    memory->create(ek,nmax,3,"ewald_dipole:ek");
    memory->create(tk,nmax,3,"ewald_dipole:tk");
    memory->create(vc,kmax3d,6,"ewald_dipole:tk");
    memory->create3d_offset(cs,-kmax,kmax,3,nmax,"ewald_dipole:cs");
    memory->create3d_offset(sn,-kmax,kmax,3,nmax,"ewald_dipole:sn");
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
  const double g3 = g_ewald*g_ewald*g_ewald;

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

  if (musqsum == 0.0) return;

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ek);
    memory->destroy(tk);
    memory->destroy(vc);
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    nmax = atom->nmax;
    memory->create(ek,nmax,3,"ewald_dipole:ek");
    memory->create(tk,nmax,3,"ewald_dipole:tk");
    memory->create(vc,kmax3d,6,"ewald_dipole:tk");
    memory->create3d_offset(cs,-kmax,kmax,3,nmax,"ewald_dipole:cs");
    memory->create3d_offset(sn,-kmax,kmax,3,nmax,"ewald_dipole:sn");
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
  double **t = atom->torque;
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  int kx,ky,kz;
  double cypz,sypz,exprl,expim;
  double partial,partial2,partial_peratom;
  double vcik[6];
  double mudotk;

  for (i = 0; i < nlocal; i++) {
    ek[i][0] = ek[i][1] = ek[i][2] = 0.0;
    tk[i][0] = tk[i][1] = tk[i][2] = 0.0;
  }

  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];
    for (j = 0; j<6; j++) vc[k][j] = 0.0;

    for (i = 0; i < nlocal; i++) {

      for (j = 0; j<6; j++) vcik[j] = 0.0;
      
      // re-evaluating mu dot k
      mudotk = mu[i][0]*kx*unitk[0] + mu[i][1]*ky*unitk[1] + mu[i][2]*kz*unitk[2];

      // calculating  re and im of exp(i*k*ri)

      cypz = cs[ky][1][i]*cs[kz][2][i] - sn[ky][1][i]*sn[kz][2][i];
      sypz = sn[ky][1][i]*cs[kz][2][i] + cs[ky][1][i]*sn[kz][2][i];
      exprl = cs[kx][0][i]*cypz - sn[kx][0][i]*sypz;
      expim = sn[kx][0][i]*cypz + cs[kx][0][i]*sypz;

      // taking im of struct_fact x exp(i*k*ri) (for force calc.)

      partial = (mudotk)*(expim*sfacrl_all[k] - exprl*sfacim_all[k]);
      ek[i][0] += partial * eg[k][0];
      ek[i][1] += partial * eg[k][1];
      ek[i][2] += partial * eg[k][2];

      // compute field for torque calculation

      partial_peratom = exprl*sfacrl_all[k] + expim*sfacim_all[k];
      tk[i][0] += partial_peratom * eg[k][0];
      tk[i][1] += partial_peratom * eg[k][1];
      tk[i][2] += partial_peratom * eg[k][2];

      // total and per-atom virial correction (dipole only)

      vc[k][0] += vcik[0] = -(partial_peratom * mu[i][0] * eg[k][0]);
      vc[k][1] += vcik[1] = -(partial_peratom * mu[i][1] * eg[k][1]);
      vc[k][2] += vcik[2] = -(partial_peratom * mu[i][2] * eg[k][2]);
      vc[k][3] += vcik[3] = -(partial_peratom * mu[i][0] * eg[k][1]);
      vc[k][4] += vcik[4] = -(partial_peratom * mu[i][0] * eg[k][2]);
      vc[k][5] += vcik[5] = -(partial_peratom * mu[i][1] * eg[k][2]);
      
      // taking re-part of struct_fact x exp(i*k*ri) 
      // (for per-atom energy and virial calc.)

      if (evflag_atom) {
        if (eflag_atom) eatom[i] += mudotk*ug[k]*partial_peratom;
        if (vflag_atom)
          for (j = 0; j < 6; j++)
	    vatom[i][j] += (ug[k]*mudotk*vg[k][j]*partial_peratom - vcik[j]);
      }
    }
  }

  // force and torque calculation

  const double muscale = qqrd2e * scale;

  for (i = 0; i < nlocal; i++) {
    f[i][0] += muscale * ek[i][0];
    f[i][1] += muscale * ek[i][1];
    if (slabflag != 2) f[i][2] += muscale * ek[i][2];
    t[i][0] -= muscale * (mu[i][1]*tk[i][2] - mu[i][2]*tk[i][1]);
    t[i][1] -= muscale * (mu[i][2]*tk[i][0] - mu[i][0]*tk[i][2]);
    if (slabflag != 2) t[i][2] -= muscale * (mu[i][0]*tk[i][1] - mu[i][1]*tk[i][0]);
  }

  // sum global energy across Kspace vevs and add in volume-dependent term
  // taking the re-part of struct_fact_i x struct_fact_j
  // substracting self energy and scaling

  if (eflag_global) {
    for (k = 0; k < kcount; k++) {
      energy += ug[k] * (sfacrl_all[k]*sfacrl_all[k] +
                         sfacim_all[k]*sfacim_all[k]);
    }
    energy -= musqsum*2.0*g3/3.0/MY_PIS;
    energy *= muscale;
  }

  // global virial

  if (vflag_global) {
    double uk, vk;
    for (k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
      for (j = 0; j < 6; j++) virial[j] += uk*vg[k][j] - vc[k][j];
    }
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
        for (j = 0; j < 6; j++) vatom[i][j] *= muscale;
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();
}

/* ----------------------------------------------------------------------
   compute the struc. factors and mu dot k products
------------------------------------------------------------------------- */

void EwaldDipole::eik_dot_r()
{
  int i,k,l,m,n,ic;
  double cstr1,sstr1,cstr2,sstr2,cstr3,sstr3,cstr4,sstr4;
  double sqk,clpm,slpm;
  double mux, muy, muz;
  double mudotk;

  double **x = atom->x;
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  n = 0;
  mux = muy = muz = 0.0;

  // loop on different k-directions
  // loop on n kpoints and nlocal atoms

  // (k,0,0), (0,l,0), (0,0,m)
  
  // loop 1: k=1, l=1, m=1
  // define first val. of cos and sin

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic]*unitk[ic];
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
        sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
        cs[-1][ic][i] = cs[1][ic][i];
        sn[-1][ic][i] = -sn[1][ic][i];
        mudotk = (mu[i][ic]*unitk[ic]);
        cstr1 += mudotk*cs[1][ic][i];
        sstr1 += mudotk*sn[1][ic][i];
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
	  mudotk = (mu[i][ic]*m*unitk[ic]);
          cstr1 += mudotk*cs[m][ic][i];
          sstr1 += mudotk*sn[m][ic][i];
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
	  mudotk = (mux*k*unitk[0] + muy*l*unitk[1]);
          cstr1 += mudotk*(cs[k][0][i]*cs[l][1][i]-sn[k][0][i]*sn[l][1][i]);
          sstr1 += mudotk*(sn[k][0][i]*cs[l][1][i]+cs[k][0][i]*sn[l][1][i]);
	  
	  // dir 2: (k,-l,0)
	  mudotk = (mux*k*unitk[0] - muy*l*unitk[1]);
          cstr2 += mudotk*(cs[k][0][i]*cs[l][1][i]+sn[k][0][i]*sn[l][1][i]);
          sstr2 += mudotk*(sn[k][0][i]*cs[l][1][i]-cs[k][0][i]*sn[l][1][i]);
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
      	  mudotk = (muy*l*unitk[1] + muz*m*unitk[2]); 
          cstr1 += mudotk*(cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i]);
          sstr1 += mudotk*(sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i]);
	  
	  // dir 2: (0,l,-m)
	  mudotk = (muy*l*unitk[1] - muz*m*unitk[2]); 
          cstr2 += mudotk*(cs[l][1][i]*cs[m][2][i]+sn[l][1][i]*sn[m][2][i]);
          sstr2 += mudotk*(sn[l][1][i]*cs[m][2][i]-cs[l][1][i]*sn[m][2][i]);
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
	  mudotk = (mux*k*unitk[0] + muz*m*unitk[2]); 
          cstr1 += mudotk*(cs[k][0][i]*cs[m][2][i]-sn[k][0][i]*sn[m][2][i]);
          sstr1 += mudotk*(sn[k][0][i]*cs[m][2][i]+cs[k][0][i]*sn[m][2][i]);
	  
	  // dir 2: (k,0,-m)
	  mudotk = (mux*k*unitk[0] - muz*m*unitk[2]); 
          cstr2 += mudotk*(cs[k][0][i]*cs[m][2][i]+sn[k][0][i]*sn[m][2][i]);
          sstr2 += mudotk*(sn[k][0][i]*cs[m][2][i]-cs[k][0][i]*sn[m][2][i]);
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
	    mudotk = (mux*k*unitk[0] + muy*l*unitk[1] + muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr1 += mudotk*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr1 += mudotk*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 2: (k,-l,m)
	    mudotk = (mux*k*unitk[0] - muy*l*unitk[1] + muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr2 += mudotk*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr2 += mudotk*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 3: (k,l,-m)
	    mudotk = (mux*k*unitk[0] + muy*l*unitk[1] - muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr3 += mudotk*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr3 += mudotk*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

	    // dir 4: (k,-l,-m)
	    mudotk = (mux*k*unitk[0] - muy*l*unitk[1] - muz*m*unitk[2]); 
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr4 += mudotk*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
            sstr4 += mudotk*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
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

  double **x = atom->x;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  double **mu = atom->mu;
  for (int i = 0; i < nlocal; i++) dipole += mu[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  if (eflag_atom || fabs(qsum) > SMALL) {

    error->all(FLERR,"Cannot (yet) use kspace slab correction with "
      "long-range dipoles and non-neutral systems or per-atom energy");
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume/12.0;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * mu[i][2]*dipole_all;
  }

  // add on torque corrections

  if (atom->torque) {
    double ffact = qscale * (-4.0*MY_PI/volume);
    double **mu = atom->mu;
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++) {
      torque[i][0] += ffact * dipole_all * mu[i][1];
      torque[i][1] += -ffact * dipole_all * mu[i][0];
    }
  }
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
    error->all(FLERR,"Using kspace solver EwaldDipole on system with no dipoles");
}

/* ----------------------------------------------------------------------
   Newton solver used to find g_ewald for LJ systems
------------------------------------------------------------------------- */

double EwaldDipole::NewtonSolve(double x, double Rc,
                              bigint natoms, double vol, double b2)
{
  double dx,tol;
  int maxit;

  maxit = 10000; //Maximum number of iterations
  tol = 0.00001; //Convergence tolerance

  //Begin algorithm

  for (int i = 0; i < maxit; i++) {
    dx = f(x,Rc,natoms,vol,b2) / derivf(x,Rc,natoms,vol,b2);
    x = x - dx; //Update x
    if (fabs(dx) < tol) return x;
    if (x < 0 || x != x) // solver failed
      return -1;
  }
  return -1;
}

/* ----------------------------------------------------------------------
 Calculate f(x)
 ------------------------------------------------------------------------- */

double EwaldDipole::f(double x, double Rc, bigint natoms, double vol, double b2)
{
  double a = Rc*x;
  double f = 0.0;

  double rg2 = a*a;
  double rg4 = rg2*rg2;
  double rg6 = rg4*rg2;
  double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
  double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
  f = (b2/(sqrt(vol*powint(x,4)*powint(Rc,9)*natoms)) *
    sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) *
    exp(-rg2)) - accuracy;

  return f;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x)
 ------------------------------------------------------------------------- */

double EwaldDipole::derivf(double x, double Rc,
                         bigint natoms, double vol, double b2)
{
  double h = 0.000001;  //Derivative step-size
  return (f(x + h,Rc,natoms,vol,b2) - f(x,Rc,natoms,vol,b2)) / h;
}

