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
   Contributing author: Tod A Pascal (Caltech)
------------------------------------------------------------------------- */

#include "pair_hbond_dreiding_lj.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "domain.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.001
#define CHUNK 8

/* ---------------------------------------------------------------------- */

PairHbondDreidingLJ::PairHbondDreidingLJ(LAMMPS *lmp) : Pair(lmp)
{
  // hbond cannot compute virial as F dot r
  // due to using map() to find bonded H atoms which are not near donor atom

  no_virial_fdotr_compute = 1;
  restartinfo = 0;

  nparams = maxparam = 0;
  params = NULL;

  nextra = 2;
  pvector = new double[2];
}

/* ---------------------------------------------------------------------- */

PairHbondDreidingLJ::~PairHbondDreidingLJ()
{
  memory->sfree(params);
  delete [] pvector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    delete [] donor;
    delete [] acceptor;
    memory->destroy(type2param);
  }
}

/* ---------------------------------------------------------------------- */

void PairHbondDreidingLJ::compute(int eflag, int vflag)
{
  int i,j,k,m,ii,jj,kk,inum,jnum,knum,itype,jtype,ktype,iatom,imol;
  tagint tagprev;
  double delx,dely,delz,rsq,rsq1,rsq2,r1,r2;
  double factor_hb,force_angle,force_kernel,evdwl,eng_lj,ehbond,force_switch;
  double c,s,a,b,ac,a11,a12,a22,vx1,vx2,vy1,vy2,vz1,vz2,d;
  double fi[3],fj[3],delr1[3],delr2[3];
  double r2inv,r10inv;
  double switch1,switch2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  tagint *klist;

  evdwl = ehbond = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int *type = atom->type;
  double *special_lj = force->special_lj;
  int molecular = atom->molecular;
  Molecule **onemols = atom->avec->onemols;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // ii = loop over donors
  // jj = loop over acceptors
  // kk = loop over hydrogens bonded to donor

  int hbcount = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (!donor[itype]) continue;
    if (molecular == 1) {
      klist = special[i];
      knum = nspecial[i][0];
    } else {
      if (molindex[i] < 0) continue;
      imol = molindex[i];
      iatom = molatom[i];
      klist = onemols[imol]->special[iatom];
      knum = onemols[imol]->nspecial[iatom][0];
      tagprev = tag[i] - iatom - 1;
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_hb = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      jtype = type[j];
      if (!acceptor[jtype]) continue;

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      for (kk = 0; kk < knum; kk++) {
        if (molecular == 1) k = atom->map(klist[kk]);
        else k = atom->map(klist[kk]+tagprev);
        if (k < 0) continue;
        ktype = type[k];
        m = type2param[itype][jtype][ktype];
        if (m < 0) continue;
        const Param &pm = params[m];

        if (rsq < pm.cut_outersq) {
          delr1[0] = x[i][0] - x[k][0];
          delr1[1] = x[i][1] - x[k][1];
          delr1[2] = x[i][2] - x[k][2];
          domain->minimum_image(delr1);
          rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
          r1 = sqrt(rsq1);

          delr2[0] = x[j][0] - x[k][0];
          delr2[1] = x[j][1] - x[k][1];
          delr2[2] = x[j][2] - x[k][2];
          domain->minimum_image(delr2);
          rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
          r2 = sqrt(rsq2);

          // angle (cos and sin)

          c = delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2];
          c /= r1*r2;
          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
          ac = acos(c);

          if (ac > pm.cut_angle && ac < (2.0*MY_PI - pm.cut_angle)) {
            s = sqrt(1.0 - c*c);
            if (s < SMALL) s = SMALL;

            // LJ-specific kernel

            r2inv = 1.0/rsq;
            r10inv = r2inv*r2inv*r2inv*r2inv*r2inv;
            force_kernel = r10inv*(pm.lj1*r2inv - pm.lj2)*r2inv *
              powint(c,pm.ap);
            force_angle = pm.ap * r10inv*(pm.lj3*r2inv - pm.lj4) *
              powint(c,pm.ap-1)*s;

            eng_lj = r10inv*(pm.lj3*r2inv - pm.lj4);

            force_switch=0.0;

            if (rsq > pm.cut_innersq) {
              switch1 = (pm.cut_outersq-rsq) * (pm.cut_outersq-rsq) *
                        (pm.cut_outersq + 2.0*rsq - 3.0*pm.cut_innersq) /
                        pm.denom_vdw;
              switch2 = 12.0*rsq * (pm.cut_outersq-rsq) *
                        (rsq-pm.cut_innersq) / pm.denom_vdw;

              force_kernel *= switch1;
              force_angle  *= switch1;
              force_switch  = eng_lj*switch2/rsq;
              eng_lj       *= switch1;
            }

            if (eflag) {
              evdwl = eng_lj * powint(c,pm.ap);
              evdwl *= factor_hb;
              ehbond += evdwl;
            }

            a = factor_hb*force_angle/s;
            b = factor_hb*force_kernel;
            d = factor_hb*force_switch;

            a11 = a*c / rsq1;
            a12 = -a / (r1*r2);
            a22 = a*c / rsq2;

            vx1 = a11*delr1[0] + a12*delr2[0];
            vx2 = a22*delr2[0] + a12*delr1[0];
            vy1 = a11*delr1[1] + a12*delr2[1];
            vy2 = a22*delr2[1] + a12*delr1[1];
            vz1 = a11*delr1[2] + a12*delr2[2];
            vz2 = a22*delr2[2] + a12*delr1[2];

            fi[0] = vx1 + b*delx + d*delx;
            fi[1] = vy1 + b*dely + d*dely;
            fi[2] = vz1 + b*delz + d*delz;
            fj[0] = vx2 - b*delx - d*delx;
            fj[1] = vy2 - b*dely - d*dely;
            fj[2] = vz2 - b*delz - d*delz;

            f[i][0] += fi[0];
            f[i][1] += fi[1];
            f[i][2] += fi[2];

            f[j][0] += fj[0];
            f[j][1] += fj[1];
            f[j][2] += fj[2];

            f[k][0] -= vx1 + vx2;
            f[k][1] -= vy1 + vy2;
            f[k][2] -= vz1 + vz2;

            // KIJ instead of IJK b/c delr1/delr2 are both with respect to k

            if (evflag) ev_tally3(k,i,j,evdwl,0.0,fi,fj,delr1,delr2);

            hbcount++;
          }
        }
      }
    }
  }

  if (eflag_global) {
    pvector[0] = hbcount;
    pvector[1] = ehbond;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairHbondDreidingLJ::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  // mark all setflag as set, since don't require pair_coeff of all I,J

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  donor = new int[n+1];
  acceptor = new int[n+1];
  memory->create(type2param,n+1,n+1,n+1,"pair:type2param");

  int i,j,k;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      for (k = 1; k <= n; k++)
        type2param[i][j][k] = -1;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairHbondDreidingLJ::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  ap_global = force->inumeric(FLERR,arg[0]);
  cut_inner_global = force->numeric(FLERR,arg[1]);
  cut_outer_global = force->numeric(FLERR,arg[2]);
  cut_angle_global = force->numeric(FLERR,arg[3]) * MY_PI/180.0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHbondDreidingLJ::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 10)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,klo,khi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
  force->bounds(FLERR,arg[2],atom->ntypes,klo,khi);

  int donor_flag;
  if (strcmp(arg[3],"i") == 0) donor_flag = 0;
  else if (strcmp(arg[3],"j") == 0) donor_flag = 1;
  else error->all(FLERR,"Incorrect args for pair coefficients");

  double epsilon_one = force->numeric(FLERR,arg[4]);
  double sigma_one = force->numeric(FLERR,arg[5]);

  int ap_one = ap_global;
  if (narg > 6) ap_one = force->inumeric(FLERR,arg[6]);
  double cut_inner_one = cut_inner_global;
  double cut_outer_one = cut_outer_global;
  if (narg > 8) {
    cut_inner_one = force->numeric(FLERR,arg[7]);
    cut_outer_one = force->numeric(FLERR,arg[8]);
  }
  if (cut_inner_one>cut_outer_one)
    error->all(FLERR,"Pair inner cutoff >= Pair outer cutoff");
  double cut_angle_one = cut_angle_global;
  if (narg == 10) cut_angle_one = force->numeric(FLERR,arg[9]) * MY_PI/180.0;
  // grow params array if necessary

  if (nparams == maxparam) {
    maxparam += CHUNK;
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                        "pair:params");
  }

  params[nparams].epsilon = epsilon_one;
  params[nparams].sigma = sigma_one;
  params[nparams].ap = ap_one;
  params[nparams].cut_inner = cut_inner_one;
  params[nparams].cut_outer = cut_outer_one;
  params[nparams].cut_innersq = cut_inner_one*cut_inner_one;
  params[nparams].cut_outersq = cut_outer_one*cut_outer_one;
  params[nparams].cut_angle = cut_angle_one;
  params[nparams].denom_vdw =
    (params[nparams].cut_outersq-params[nparams].cut_innersq) *
    (params[nparams].cut_outersq-params[nparams].cut_innersq) *
    (params[nparams].cut_outersq-params[nparams].cut_innersq);

  // flag type2param with either i,j = D,A or j,i = D,A

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo,i); j <= jhi; j++)
      for (int k = klo; k <= khi; k++) {
        if (donor_flag == 0) type2param[i][j][k] = nparams;
        else type2param[j][i][k] = nparams;
        count++;
      }
  nparams++;

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHbondDreidingLJ::init_style()
{
  // molecular system required to use special list to find H atoms
  // tags required to use special list
  // pair newton on required since are looping over D atoms
  //   and computing forces on A,H which may be on different procs

  if (atom->molecular == 0)
    error->all(FLERR,"Pair style hbond/dreiding requires molecular system");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style hbond/dreiding requires atom IDs");
  if (atom->map_style == 0)
    error->all(FLERR,"Pair style hbond/dreiding requires an atom map, "
               "see atom_modify");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style hbond/dreiding requires newton pair on");

  // set donor[M]/acceptor[M] if any atom of type M is a donor/acceptor

  int anyflag = 0;
  int n = atom->ntypes;
  for (int m = 1; m <= n; m++) donor[m] = acceptor[m] = 0;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      for (int k = 1; k <= n; k++)
        if (type2param[i][j][k] >= 0) {
          anyflag = 1;
          donor[i] = 1;
          acceptor[j] = 1;
        }

  if (!anyflag) error->all(FLERR,"No pair hbond/dreiding coefficients set");

  // set additional param values
  // offset is for LJ only, angle term is not included

  for (int m = 0; m < nparams; m++) {
    params[m].lj1 = 60.0*params[m].epsilon*pow(params[m].sigma,12.0);
    params[m].lj2 = 60.0*params[m].epsilon*pow(params[m].sigma,10.0);
    params[m].lj3 = 5.0*params[m].epsilon*pow(params[m].sigma,12.0);
    params[m].lj4 = 6.0*params[m].epsilon*pow(params[m].sigma,10.0);

    /*
    if (offset_flag) {
      double ratio = params[m].sigma / params[m].cut_outer;
      params[m].offset = params[m].epsilon *
        ((2.0*pow(ratio,9.0)) - (3.0*pow(ratio,6.0)));
    } else params[m].offset = 0.0;
    */
  }

  // full neighbor list request

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHbondDreidingLJ::init_one(int i, int j)
{
  int m;

  // return maximum cutoff for any K with I,J = D,A or J,I = D,A
  // donor/acceptor is not symmetric, IJ interaction != JI interaction

  double cut = 0.0;
  for (int k = 1; k <= atom->ntypes; k++) {
    m = type2param[i][j][k];
    if (m >= 0) cut = MAX(cut,params[m].cut_outer);
    m = type2param[j][i][k];
    if (m >= 0) cut = MAX(cut,params[m].cut_outer);
  }
  return cut;
}

/* ---------------------------------------------------------------------- */

double PairHbondDreidingLJ::single(int i, int j, int itype, int jtype,
                                   double rsq,
                                   double /*factor_coul*/, double /*factor_lj*/,
                                   double &fforce)
{
  int k,kk,ktype,knum,m;
  tagint tagprev;
  double eng,eng_lj,force_kernel,force_angle;
  double rsq1,rsq2,r1,r2,c,s,ac,r2inv,r10inv,factor_hb;
  double switch1,switch2;
  double delr1[3],delr2[3];
  tagint *klist;

  double **x = atom->x;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  eng = 0.0;
  fforce = 0;

  // sanity check

  if (!donor[itype]) return 0.0;
  if (!acceptor[jtype]) return 0.0;

  int molecular = atom->molecular;
  if (molecular == 1) {
    klist = atom->special[i];
    knum = atom->nspecial[i][0];
  } else {
    if (atom->molindex[i] < 0) return 0.0;
    int imol = atom->molindex[i];
    int iatom = atom->molatom[i];
    Molecule **onemols = atom->avec->onemols;
    klist = onemols[imol]->special[iatom];
    knum = onemols[imol]->nspecial[iatom][0];
    tagprev = atom->tag[i] - iatom - 1;
  }

  factor_hb = special_lj[sbmask(j)];

  for (kk = 0; kk < knum; kk++) {
    if (molecular == 1) k = atom->map(klist[kk]);
    else k = atom->map(klist[kk]+tagprev);

    if (k < 0) continue;
    ktype = type[k];
    m = type2param[itype][jtype][ktype];
    if (m < 0) continue;
    const Param &pm = params[m];

    delr1[0] = x[i][0] - x[k][0];
    delr1[1] = x[i][1] - x[k][1];
    delr1[2] = x[i][2] - x[k][2];
    domain->minimum_image(delr1);
    rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
    r1 = sqrt(rsq1);

    delr2[0] = x[j][0] - x[k][0];
    delr2[1] = x[j][1] - x[k][1];
    delr2[2] = x[j][2] - x[k][2];
    domain->minimum_image(delr2);
    rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2];
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    ac = acos(c);

    if (ac < pm.cut_angle || ac > (2.0*MY_PI - pm.cut_angle)) return 0.0;
    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;

    // LJ-specific kernel

    r2inv = 1.0/rsq;
    r10inv = r2inv*r2inv*r2inv*r2inv*r2inv;
    force_kernel = r10inv*(pm.lj1*r2inv - pm.lj2)*r2inv * powint(c,pm.ap);
    force_angle = pm.ap * r10inv*(pm.lj3*r2inv - pm.lj4) *
      powint(c,pm.ap-1)*s;

    // only lj part for now

    eng_lj = r10inv*(pm.lj3*r2inv - pm.lj4);
    if (rsq > pm.cut_innersq) {
      switch1 = (pm.cut_outersq-rsq) * (pm.cut_outersq-rsq) *
                (pm.cut_outersq + 2.0*rsq - 3.0*pm.cut_innersq) / pm.denom_vdw;
      switch2 = 12.0*rsq * (pm.cut_outersq-rsq) *
                (rsq-pm.cut_innersq) / pm.denom_vdw;
      force_kernel = force_kernel*switch1 + eng_lj*switch2;
      eng_lj *= switch1;
    }

    fforce += force_kernel*powint(c,pm.ap) + eng_lj*force_angle;
    eng += eng_lj * powint(c,pm.ap) * factor_hb;
  }

  return eng;
}
