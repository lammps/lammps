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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "pair_hybrid.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair_spin_cubic.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpinCubic::PairSpinCubic(LAMMPS *lmp) : PairSpin(lmp),
lockfixnvespin(NULL)
{
  single_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairSpinCubic::~PairSpinCubic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_cubic);
    memory->destroy(K1_mag);
    memory->destroy(K2_mag);
    memory->destroy(K1_mech);
    memory->destroy(K2_mech);
    memory->destroy(cutsq); // to be implemented
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinCubic::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");

  cut_spin_cubic_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_cubic[i][j] = cut_spin_cubic_global;
        }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs
------------------------------------------------------------------------- */

void PairSpinCubic::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // check if args correct

  if (strcmp(arg[2],"cubic") != 0)
    error->all(FLERR,"Incorrect args in pair_style command");
  if (narg != 6)
    error->all(FLERR,"Incorrect args in pair_style command");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // get cubic aniso arguments from input command

  const double rc = force->numeric(FLERR,arg[3]);
  const double k1 = force->numeric(FLERR,arg[4]);
  const double k2 = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_cubic[i][j] = rc;
      K1_mag[i][j] = k1/hbar;
      K2_mag[i][j] = k2/hbar;
      K1_mech[i][j] = k1;
      K2_mech[i][j] = k2;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args in pair_style command");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinCubic::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    ifix++;
  }
  if ((ifix == modify->nfix) && (comm->me == 0))
    error->warning(FLERR,"Using pair/spin style without nve/spin");

  // get the lattice_flag from nve/spin

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nve/spin") == 0) {
      lockfixnvespin = (FixNVESpin *) modify->fix[i];
      lattice_flag = lockfixnvespin->lattice_flag;
    }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinCubic::init_one(int i, int j)
{

   if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  K1_mag[j][i] = K1_mag[i][j];
  K2_mag[j][i] = K2_mag[i][j];
  K1_mech[j][i] = K1_mech[i][j];
  K2_mech[j][i] = K2_mech[i][j];
  cut_spin_cubic[j][i] = cut_spin_cubic[i][j];

  return cut_spin_cubic_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinCubic::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_cubic_global;
  return NULL;
}

/* ---------------------------------------------------------------------- */
// compute cubic aniso

// Warning: procedure does not work if not cubic crystal
// and if orient command was used (has to be in canonical 
// xyz coordinates).
// Also, issue with surface => to be understand

void PairSpinCubic::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl, ecoul;
  double xi[3], eij[3], rij[3];
  double spi[3], spj[3];
  double fi[3], fmi[3];
  double delx,dely,delz;
  double delx2,dely2,delz2;
  double ea1[3], ea2[3], ea3[3];
  double local_cut2;
  double rsq, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *emag = atom->emag;
  double **fm = atom->fm;
  double **sp = atom->sp;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // define dx as an entry param
  double dx2 = 0.2;

  // computation of the cubic interaction
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];

    ea1[0] = ea1[1] = ea1[2] = 0.0;
    ea2[0] = ea2[1] = ea2[2] = 0.0;
    ea3[0] = ea3[1] = ea3[2] = 0.0;

    evdwl = 0.0;


    set_axis(i,ea1,ea2,ea3);


    // first loop over neighbors j 
    // finding the anisotropy axes

    //for (jj = 0; jj < jnum; jj++) {
    //  j = jlist[jj];
    //  j &= NEIGHMASK;
    //  jtype = type[j];

    //  rij[0] = x[j][0] - xi[0];
    //  rij[1] = x[j][1] - xi[1];
    //  rij[2] = x[j][2] - xi[2];
    //  rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    //  inorm = 1.0/sqrt(rsq);
    //  eij[0] = rij[0]*inorm;
    //  eij[1] = rij[1]*inorm;
    //  eij[2] = rij[2]*inorm;
    //  //delx = xi[0] - x[j][0];
    //  //dely = xi[1] - x[j][1];
    //  //delz = xi[2] - x[j][2];

    //  local_cut2 = cut_spin_cubic[itype][jtype]*cut_spin_cubic[itype][jtype];

    //  // finding anisotropy axes

    //  //delx2 = delx*delx;
    //  //dely2 = dely*dely;
    //  //delz2 = delz*delz;
    //  delx2 = rij[0]*rij[0];
    //  dely2 = rij[1]*rij[1];
    //  delz2 = rij[2]*rij[2];

    //  if (delx2 > dx2 && delx2 <= local_cut2) {
    //    ea1[0] += eij[0];
    //    ea1[1] += eij[1];
    //    ea1[2] += eij[2];
    //  } 

    //  if (dely2 > dx2 && dely2 <= local_cut2) {
    //    ea2[0] += eij[0];
    //    ea2[1] += eij[1];
    //    ea2[2] += eij[2];
    //  }

    //  if (delz2 > dx2 && delz2 <= local_cut2) {
    //    ea3[0] += eij[0];
    //    ea3[1] += eij[1];
    //    ea3[2] += eij[2];
    //  }
    //}

    //// normalizing the three aniso axes

    //double inorm1,inorm2,inorm3;
    //inorm1 = 1.0/(ea1[0]*ea1[0]+ea1[1]*ea1[1]+ea1[2]*ea1[2]);
    //ea1[0] *= inorm1;
    //ea1[1] *= inorm1;
    //ea1[2] *= inorm1;
    //inorm2 = 1.0/(ea2[0]*ea2[0]+ea2[1]*ea2[1]+ea2[2]*ea2[2]);
    //ea2[0] *= inorm2;
    //ea2[1] *= inorm2;
    //ea2[2] *= inorm2;
    //inorm3 = 1.0/(ea3[0]*ea3[0]+ea3[1]*ea3[1]+ea3[2]*ea3[2]);
    //ea3[0] *= inorm3;
    //ea3[1] *= inorm3;
    //ea3[2] *= inorm3;

    compute_cubic(i,fmi,spi,ea1,ea2,ea3);
    if (lattice_flag) {
      compute_cubic_mech(i,eij,fi,spi,ea1,ea2,ea3);
    }

    f[i][0] += fi[0];
    f[i][1] += fi[1];
    f[i][2] += fi[2];
    fm[i][0] += fmi[0];
    fm[i][1] += fmi[1];
    fm[i][2] += fmi[2];

    if (newton_pair || j < nlocal) {
      f[j][0] -= fi[0];
      f[j][1] -= fi[1];
      f[j][2] -= fi[2];
    }

    if (eflag) {
      //evdwl -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
      evdwl -= compute_cubic_energy(i,spi,ea1,ea2,ea3);
      evdwl *= hbar;
      emag[i] += evdwl;
    } else evdwl = 0.0;

    //if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
    //    evdwl,ecoul,fi[0],fi[1],fi[2],delx,dely,delz);
  }


  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinCubic::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3];
  double eij[3],rij[3];
  double inorm;

  double delx2,dely2,delz2;
  double ea1[3], ea2[3], ea3[3];
  double spi[3],spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

  double rsq;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // check if interaction applies to type of ii

  itype = type[ii];
  ntypes = atom->ntypes;
  locflag = 0;
  k = 1;
  while (k <= ntypes) {
    if (k <= itype) {
      if (setflag[k][itype] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else if (k > itype) {
      if (setflag[itype][k] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else error->all(FLERR,"Wrong type number");
  }
  
  // define dx as an entry param
  double dx2 = 0.2;

  // if interaction applies to type of ii,
  // locflag = 1 and compute pair interaction

  if (locflag == 1) {

    spi[0] = sp[ii][0];
    spi[1] = sp[ii][1];
    spi[2] = sp[ii][2];
    
    //xi[0] = x[ii][0];
    //xi[1] = x[ii][1];
    //xi[2] = x[ii][2];
    //eij[0] = eij[1] = eij[2] = 0.0;

    //jlist = firstneigh[ii];
    //jnum = numneigh[ii];

    set_axis(ii,ea1,ea2,ea3);

    //for (int jj = 0; jj < jnum; jj++) {

    //  j = jlist[jj];
    //  j &= NEIGHMASK;
    //  jtype = type[j];
    //  local_cut2 = cut_spin_cubic[itype][jtype]*cut_spin_cubic[itype][jtype];

    //  spj[0] = sp[j][0];
    //  spj[1] = sp[j][1];
    //  spj[2] = sp[j][2];

    //  rij[0] = x[j][0] - xi[0];
    //  rij[1] = x[j][1] - xi[1];
    //  rij[2] = x[j][2] - xi[2];
    //  rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    //  inorm = 1.0/sqrt(rsq);
    //  eij[0] = inorm*rij[0];
    //  eij[1] = inorm*rij[1];
    //  eij[2] = inorm*rij[2];

    //  // finding anisotropy axes

    //  //delx2 = delx*delx;
    //  //dely2 = dely*dely;
    //  //delz2 = delz*delz;
    //  delx2 = rij[0]*rij[0];
    //  dely2 = rij[1]*rij[1];
    //  delz2 = rij[2]*rij[2];

    //  if (delx2 > dx2 && delx2 <= local_cut2) {
    //    ea1[0] += eij[0];
    //    ea1[1] += eij[1];
    //    ea1[2] += eij[2];
    //  } 

    //  if (dely2 > dx2 && dely2 <= local_cut2) {
    //    ea2[0] += eij[0];
    //    ea2[1] += eij[1];
    //    ea2[2] += eij[2];
    //  }

    //  if (delz2 > dx2 && delz2 <= local_cut2) {
    //    ea3[0] += eij[0];
    //    ea3[1] += eij[1];
    //    ea3[2] += eij[2];
    //  }
    //}
    //
    //// normalizing the three aniso axes

    //double inorm1,inorm2,inorm3;
    //inorm1 = 1.0/(ea1[0]*ea1[0]+ea1[1]*ea1[1]+ea1[2]*ea1[2]);
    //ea1[0] *= inorm1;
    //ea1[1] *= inorm1;
    //ea1[2] *= inorm1;
    //inorm2 = 1.0/(ea2[0]*ea2[0]+ea2[1]*ea2[1]+ea2[2]*ea2[2]);
    //ea2[0] *= inorm2;
    //ea2[1] *= inorm2;
    //ea2[2] *= inorm2;
    //inorm3 = 1.0/(ea3[0]*ea3[0]+ea3[1]*ea3[1]+ea3[2]*ea3[2]);
    //ea3[0] *= inorm3;
    //ea3[1] *= inorm3;
    //ea3[2] *= inorm3;
    
    compute_cubic(ii,fmi,spi,ea1,ea2,ea3);
    
  }
}

/* ----------------------------------------------------------------------
   compute cubic aniso interaction between spins i and j
------------------------------------------------------------------------- */

void PairSpinCubic::compute_cubic(int i, double fmi[3], double spi[3], double ea1[3], double ea2[3], double ea3[3])
{
  int *type = atom->type;
  int itype, jtype;
  double ke1,ke2;
  double skx,sky,skz;
  double skx2,sky2,skz2;
 
  double four1,four2,four3;
  double fourx,foury,fourz;
  double six1,six2,six3;
  double sixx,sixy,sixz;

  // not dealing with different types for now 
  
  itype = type[i];
  //jtype = type[j];
  jtype = type[i];

  ke1 = K1_mech[itype][jtype];
  ke2 = K2_mech[itype][jtype];

  skx = spi[0]*ea1[0]+spi[1]*ea1[1]+spi[2]*ea1[2];
  sky = spi[0]*ea2[0]+spi[1]*ea2[1]+spi[2]*ea2[2];
  skz = spi[0]*ea3[0]+spi[1]*ea3[1]+spi[2]*ea3[2];

  skx2 = skx*skx;
  sky2 = sky*sky;
  skz2 = skz*skz;

  four1 = 2.0*skx*(sky2+skz2);
  four2 = 2.0*sky*(skx2+skz2);
  four3 = 2.0*skz*(skx2+skz2);

  fourx = ke1*(ea1[0]*four1 + ea2[0]*four2 + ea3[0]*four3);
  foury = ke1*(ea1[1]*four1 + ea2[1]*four2 + ea3[1]*four3);
  fourz = ke1*(ea1[2]*four1 + ea2[2]*four2 + ea3[2]*four3);

  six1 = 2.0*skx*sky2*skz2;
  six2 = 2.0*sky*skx2*skz2;
  six3 = 2.0*skz*skx2*sky2;
  
  sixx = ke2*(ea1[0]*six1 + ea2[0]*six2 + ea3[0]*six3);
  sixy = ke2*(ea1[1]*six1 + ea2[1]*six2 + ea3[1]*six3);
  sixz = ke2*(ea1[2]*six1 + ea2[2]*six2 + ea3[2]*six3);


  //ra = rsq/J3[itype][jtype]/J3[itype][jtype];
  //Jex = 4.0*J1_mag[itype][jtype]*ra;
  //Jex *= (1.0-J2[itype][jtype]*ra);
  //Jex *= exp(-ra);

  //fmi[0] += Jex*spj[0];
  //fmi[1] += Jex*spj[1];
  //fmi[2] += Jex*spj[2];
  
  fmi[0] += fourx + sixx;
  fmi[1] += foury + sixy;
  fmi[2] += fourz + sixz;
}

/* ----------------------------------------------------------------------
   compute cubic aniso interaction between spins i and j
------------------------------------------------------------------------- */

double PairSpinCubic::compute_cubic_energy(int i, double spi[3], double ea1[3], double ea2[3], double ea3[3])
{
  int *type = atom->type;
  int itype, jtype;
  double ke1,ke2,energy;
  double skx,sky,skz;

  // not dealing with different types for now

  energy = 0.0;
  itype = type[i];
  //jtype = type[j];
  jtype = type[i];

  ke1 = K1_mech[itype][jtype];
  ke2 = K2_mech[itype][jtype];

  skx = spi[0]*ea1[0]+spi[1]*ea1[1]+spi[2]*ea1[2];
  sky = spi[0]*ea2[0]+spi[1]*ea2[1]+spi[2]*ea2[2];
  skz = spi[0]*ea3[0]+spi[1]*ea3[1]+spi[2]*ea3[2];

  energy = ke1*(skx*skx*sky*sky + sky*sky*skz*skz + skx*skx*skz*skz);
  energy += ke2*skx*skx*sky*sky*skz*skz;

  //printf("test energy: %g\n",energy);
  //printf("test k1,k2: %g %g\n",ke1,ke2);
  //printf("test skx,sky,sky: %g %g %g\n",skx,sky,skz);
  //printf("test ea1,ea2,ea3: %g %g %g\n",ea1[0],ea2[0],ea3[0]);

  return energy;
}

/* ----------------------------------------------------------------------
   compute the mechanical force due to the cubic aniso between atom i and atom j
------------------------------------------------------------------------- */

void PairSpinCubic::compute_cubic_mech(int i, double eij[3], double fi[3],  
    double spi[3], double ea1[3], double ea2[3], double ea3[3])
{
  int *type = atom->type;
  int itype, jtype;
  double ke1,ke2;
  
  // not dealing with different types for now
  
  itype = type[i];
  //jtype = type[j];
  jtype = type[i];

  ke1 = K1_mech[itype][jtype];
  ke2 = K2_mech[itype][jtype];

  //Jex = J1_mech[itype][jtype];
  //iJ3 = 1.0/(J3[itype][jtype]*J3[itype][jtype]);

  //ra = rsq*iJ3;
  //rr = sqrt(rsq)*iJ3;

  //Jex_mech = 1.0-ra-J2[itype][jtype]*ra*(2.0-ra);
  //Jex_mech *= 8.0*Jex*rr*exp(-ra);
  //Jex_mech *= (spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2]);

  //fi[0] -= Jex_mech*eij[0];
  //fi[1] -= Jex_mech*eij[1];
  //fi[2] -= Jex_mech*eij[2];
}

/* ----------------------------------------------------------------------
   set three cubic axis
------------------------------------------------------------------------- */

void PairSpinCubic::set_axis(int ii, double ea1[3], double ea2[3], double ea3[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3];
  double eij[3],rij[3];
  double inorm,rsq;

  double delx2,dely2,delz2;
  double spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;


  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // check if interaction applies to type of ii

  itype = type[ii];
  ntypes = atom->ntypes;

  jlist = firstneigh[ii];
  jnum = numneigh[ii];
  
  for (int jj = 0; jj < jnum; jj++) {

    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = type[j];
    local_cut2 = cut_spin_cubic[itype][jtype]*cut_spin_cubic[itype][jtype];

    spj[0] = sp[j][0];
    spj[1] = sp[j][1];
    spj[2] = sp[j][2];

    rij[0] = x[j][0] - xi[0];
    rij[1] = x[j][1] - xi[1];
    rij[2] = x[j][2] - xi[2];
    rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    inorm = 1.0/sqrt(rsq);
    eij[0] = inorm*rij[0];
    eij[1] = inorm*rij[1];
    eij[2] = inorm*rij[2];

    // finding anisotropy axes

    //delx2 = delx*delx;
    //dely2 = dely*dely;
    //delz2 = delz*delz;
    delx2 = rij[0]*rij[0];
    dely2 = rij[1]*rij[1];
    delz2 = rij[2]*rij[2];
  
    // define dx as an entry param
    double dx2 = 0.2;

    if (delx2 > dx2 && delx2 <= local_cut2) {
      ea1[0] += eij[0];
      ea1[1] += eij[1];
      ea1[2] += eij[2];
    } 

    if (dely2 > dx2 && dely2 <= local_cut2) {
      ea2[0] += eij[0];
      ea2[1] += eij[1];
      ea2[2] += eij[2];
    }

    if (delz2 > dx2 && delz2 <= local_cut2) {
      ea3[0] += eij[0];
      ea3[1] += eij[1];
      ea3[2] += eij[2];
    }
  }
  
  // normalizing the three aniso axes

  double inorm1,inorm2,inorm3;
  inorm1 = 1.0/(ea1[0]*ea1[0]+ea1[1]*ea1[1]+ea1[2]*ea1[2]);
  ea1[0] *= inorm1;
  ea1[1] *= inorm1;
  ea1[2] *= inorm1;
  inorm2 = 1.0/(ea2[0]*ea2[0]+ea2[1]*ea2[1]+ea2[2]*ea2[2]);
  ea2[0] *= inorm2;
  ea2[1] *= inorm2;
  ea2[2] *= inorm2;
  inorm3 = 1.0/(ea3[0]*ea3[0]+ea3[1]*ea3[1]+ea3[2]*ea3[2]);
  ea3[0] *= inorm3;
  ea3[1] *= inorm3;
  ea3[2] *= inorm3;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinCubic::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_cubic,n+1,n+1,"pair/spin/cubic:cut_spin_cubic");
  memory->create(K1_mag,n+1,n+1,"pair/spin/cubic:J1_mag");
  memory->create(K2_mag,n+1,n+1,"pair/spin/cubic:J1_mag");
  memory->create(K1_mech,n+1,n+1,"pair/spin/cubic:J1_mech");
  memory->create(K2_mech,n+1,n+1,"pair/spin/cubic:J1_mech");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}



/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinCubic::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&K1_mag[i][j],sizeof(double),1,fp);
        fwrite(&K2_mag[i][j],sizeof(double),1,fp);
        fwrite(&K1_mech[i][j],sizeof(double),1,fp);
        fwrite(&K2_mech[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_cubic[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinCubic::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&K1_mag[i][j],sizeof(double),1,fp);
          fread(&K2_mag[i][j],sizeof(double),1,fp);
          fread(&K1_mech[i][j],sizeof(double),1,fp);
          fread(&K2_mech[i][j],sizeof(double),1,fp);
          fread(&cut_spin_cubic[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&K1_mag[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K2_mag[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K2_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_cubic[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinCubic::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_cubic_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinCubic::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_spin_cubic_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_cubic_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}


