// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tod A Pascal (Caltech)
------------------------------------------------------------------------- */

#include "pair_hbond_dreiding_morse.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.001
#define CHUNK 8

/* ---------------------------------------------------------------------- */

PairHbondDreidingMorse::PairHbondDreidingMorse(LAMMPS *lmp) :
  PairHbondDreidingLJ(lmp) {}

/* ---------------------------------------------------------------------- */

void PairHbondDreidingMorse::compute(int eflag, int vflag)
{
  int i,j,k,m,ii,jj,kk,inum,jnum,knum,itype,jtype,ktype,imol,iatom;
  tagint tagprev;
  double delx,dely,delz,rsq,rsq1,rsq2,r1,r2;
  double factor_hb,force_angle,force_kernel,force_switch,evdwl,ehbond;
  double c,s,a,b,d,ac,a11,a12,a22,vx1,vx2,vy1,vy2,vz1,vz2;
  double fi[3],fj[3],delr1[3],delr2[3];
  double r,dr,dexp,eng_morse,switch1,switch2;
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
    if (molecular == Atom::MOLECULAR) {
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
        if (molecular == Atom::MOLECULAR) k = atom->map(klist[kk]);
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

            // Morse-specific kernel

            r = sqrt(rsq);
            dr = r - pm.r0;
            dexp = exp(-pm.alpha * dr);
            eng_morse = pm.d0 * (dexp*dexp - 2.0*dexp);
            force_kernel = pm.morse1*(dexp*dexp - dexp)/r * powint(c,pm.ap);
            force_angle = pm.ap * eng_morse * powint(c,pm.ap-1)*s;
            force_switch = 0.0;

            if (rsq > pm.cut_innersq) {
              switch1 = (pm.cut_outersq-rsq) * (pm.cut_outersq-rsq) *
                        (pm.cut_outersq + 2.0*rsq - 3.0*pm.cut_innersq) /
                        pm.denom_vdw;
              switch2 = 12.0*rsq * (pm.cut_outersq-rsq) *
                        (rsq-pm.cut_innersq) / pm.denom_vdw;

              force_kernel *= switch1;
              force_angle  *= switch1;
              force_switch  = eng_morse*switch2/rsq;
              eng_morse    *= switch1;
            }

            if (eflag) {
              evdwl = eng_morse * powint(c,pm.ap);
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

            fi[0] = vx1 + (b+d)*delx;
            fi[1] = vy1 + (b+d)*dely;
            fi[2] = vz1 + (b+d)*delz;
            fj[0] = vx2 - (b+d)*delx;
            fj[1] = vy2 - (b+d)*dely;
            fj[2] = vz2 - (b+d)*delz;

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
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHbondDreidingMorse::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 11)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,klo,khi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
  utils::bounds(FLERR,arg[2],1,atom->ntypes,klo,khi,error);

  int donor_flag;
  if (strcmp(arg[3],"i") == 0) donor_flag = 0;
  else if (strcmp(arg[3],"j") == 0) donor_flag = 1;
  else error->all(FLERR,"Incorrect args for pair coefficients");

  double d0_one = utils::numeric(FLERR,arg[4],false,lmp);
  double alpha_one = utils::numeric(FLERR,arg[5],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[6],false,lmp);

  int ap_one = ap_global;
  if (narg > 7) ap_one = utils::inumeric(FLERR,arg[7],false,lmp);
  double cut_inner_one = cut_inner_global;
  double cut_outer_one = cut_outer_global;
  if (narg > 9) {
    cut_inner_one = utils::numeric(FLERR,arg[8],false,lmp);
    cut_outer_one = utils::numeric(FLERR,arg[9],false,lmp);
  }
  if (cut_inner_one>cut_outer_one)
    error->all(FLERR,"Pair inner cutoff >= Pair outer cutoff");
  double cut_angle_one = cut_angle_global;
  if (narg > 10) cut_angle_one = utils::numeric(FLERR,arg[10],false,lmp) * MY_PI/180.0;

  // grow params array if necessary

  if (nparams == maxparam) {
    maxparam += CHUNK;
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                        "pair:params");

    // make certain all addional allocated storage is initialized
    // to avoid false positives when checking with valgrind

    memset(params + nparams, 0, CHUNK*sizeof(Param));
  }

  params[nparams].d0 = d0_one;
  params[nparams].alpha = alpha_one;
  params[nparams].r0 = r0_one;
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

void PairHbondDreidingMorse::init_style()
{
  // molecular system required to use special list to find H atoms
  // tags required to use special list
  // pair newton on required since are looping over D atoms
  //   and computing forces on A,H which may be on different procs

  if (atom->molecular == Atom::ATOMIC)
    error->all(FLERR,"Pair style hbond/dreiding requires molecular system");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style hbond/dreiding requires atom IDs");
  if (atom->map_style == Atom::MAP_NONE)
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
  // offset is for Morse only, angle term is not included

  for (int m = 0; m < nparams; m++) {
    params[m].morse1 = 2.0*params[m].d0*params[m].alpha;

    /*
    if (offset_flag) {
      double alpha_dr = -params[m].alpha * (params[m].cut - params[m].r0);
      params[m].offset = params[m].d0 *
        ((exp(2.0*alpha_dr)) - (2.0*exp(alpha_dr)));
    } else params[m].offset = 0.0;
    */
  }

  // full neighbor list request

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairHbondDreidingMorse::single(int i, int j, int itype, int jtype,
                                     double rsq,
                                     double /*factor_coul*/, double /*factor_lj*/,
                                     double &fforce)
{
  int k,kk,ktype,knum,m;
  tagint tagprev;
  double eng,eng_morse,force_kernel,force_angle;
  double rsq1,rsq2,r1,r2,c,s,ac,r,dr,dexp,factor_hb;
  double switch1,switch2;
  double delr1[3],delr2[3];
  tagint *klist;

  double **x = atom->x;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  eng = 0.0;
  fforce = 0;

  //sanity check

  if (!donor[itype]) return 0.0;
  if (!acceptor[jtype]) return 0.0;

  int molecular = atom->molecular;
  if (molecular == Atom::MOLECULAR) {
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
    if (molecular == Atom::MOLECULAR) k = atom->map(klist[kk]);
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

    // Morse-specific kernel

    r = sqrt(rsq);
    dr = r - pm.r0;
    dexp = exp(-pm.alpha * dr);
    eng_morse = pm.d0 * (dexp*dexp - 2.0*dexp);  //<-- BUGFIX 2012-11-14
    force_kernel = pm.morse1*(dexp*dexp - dexp)/r * powint(c,pm.ap);
    force_angle = pm.ap * eng_morse * powint(c,pm.ap-1)*s;

    if (rsq > pm.cut_innersq) {
      switch1 = (pm.cut_outersq-rsq) * (pm.cut_outersq-rsq) *
                (pm.cut_outersq + 2.0*rsq - 3.0*pm.cut_innersq) /
                pm.denom_vdw;
      switch2 = 12.0*rsq * (pm.cut_outersq-rsq) *
                (rsq-pm.cut_innersq) / pm.denom_vdw;
      force_kernel = force_kernel*switch1 + eng_morse*switch2;
      eng_morse *= switch1;
    }

    eng += eng_morse * powint(c,pm.ap)* factor_hb;
    fforce += force_kernel*powint(c,pm.ap) + eng_morse*force_angle;
  }

  return eng;
}
