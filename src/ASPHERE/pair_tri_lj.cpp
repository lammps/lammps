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

#include "pair_tri_lj.h"
#include <cmath>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_tri.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

#define DELTA 20

/* ---------------------------------------------------------------------- */

PairTriLJ::PairTriLJ(LAMMPS *lmp) : Pair(lmp)
{
  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  single_enable = 0;
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairTriLJ::~PairTriLJ()
{
  memory->sfree(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairTriLJ::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,term1,term2,sig,sig3,forcelj;
  double dxi,dxj,dyi,dyj,dzi,dzj;
  double xi[3],xj[3],fi[3],fj[3],ti[3],tj[3],p[3][3];
  double dc1[3],dc2[3],dc3[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  AtomVecTri::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *tri = atom->tri;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // grow discrete list if necessary and initialize

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->create(dnum,nall,"pair:dnum");
    memory->create(dfirst,nall,"pair:dfirst");
  }
  for (i = 0; i < nall; i++) dnum[i] = 0;
  ndiscrete = 0;

  // loop over neighbors of my atoms

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

      if (rsq >= cutsq[itype][jtype]) continue;

      // tri/tri interactions = NxN particles
      // c1,c2,c3 = corner pts of triangle I or J

      evdwl = 0.0;
      if (tri[i] >= 0 && tri[j] >= 0) {
        if (dnum[i] == 0) {
          MathExtra::quat_to_mat(bonus[tri[i]].quat,p);
          MathExtra::matvec(p,bonus[tri[i]].c1,dc1);
          MathExtra::matvec(p,bonus[tri[i]].c2,dc2);
          MathExtra::matvec(p,bonus[tri[i]].c3,dc3);
          dfirst[i] = ndiscrete;
          discretize(i,sigma[itype][itype],dc1,dc2,dc3);
          dnum[i] = ndiscrete - dfirst[i];
        }
        npi = dnum[i];
        ifirst = dfirst[i];

        if (dnum[j] == 0) {
          MathExtra::quat_to_mat(bonus[tri[j]].quat,p);
          MathExtra::matvec(p,bonus[tri[j]].c1,dc1);
          MathExtra::matvec(p,bonus[tri[j]].c2,dc2);
          MathExtra::matvec(p,bonus[tri[j]].c3,dc3);
          dfirst[j] = ndiscrete;
          discretize(j,sigma[jtype][jtype],dc1,dc2,dc3);
          dnum[j] = ndiscrete - dfirst[j];
        }
        npj = dnum[j];
        jfirst = dfirst[j];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;
          dzi = discrete[ifirst+ni].dz;

          for (nj = 0; nj < npj; nj++) {
            dxj = discrete[jfirst+nj].dx;
            dyj = discrete[jfirst+nj].dy;
            dzj = discrete[jfirst+nj].dz;

            xi[0] = x[i][0] + dxi;
            xi[1] = x[i][1] + dyi;
            xi[2] = x[i][2] + dzi;
            xj[0] = x[j][0] + dxj;
            xj[1] = x[j][1] + dyj;
            xj[2] = x[j][2] + dzj;

            delx = xi[0] - xj[0];
            dely = xi[1] - xj[1];
            delz = xi[2] - xj[2];
            rsq = delx*delx + dely*dely + delz*delz;

            sig = 0.5 * (discrete[ifirst+ni].sigma+discrete[jfirst+nj].sigma);
            sig3 = sig*sig*sig;
            term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
            term1 = 2.0 * term2 * sig3*sig3;
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (term1*r6inv - term2);
            fpair = forcelj*r2inv;

            if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

            fi[0] = delx*fpair;
            fi[1] = dely*fpair;
            fi[2] = delz*fpair;
            f[i][0] += fi[0];
            f[i][1] += fi[1];
            f[i][2] += fi[2];
            ti[0] = dyi*fi[2] - dzi*fi[1];
            ti[1] = dzi*fi[0] - dxi*fi[2];
            ti[2] = dxi*fi[1] - dyi*fi[0];
            torque[i][0] += ti[0];
            torque[i][1] += ti[1];
            torque[i][2] += ti[2];

            if (newton_pair || j < nlocal) {
              fj[0] = -delx*fpair;
              fj[1] = -dely*fpair;
              fj[2] = -delz*fpair;
              f[j][0] += fj[0];
              f[j][1] += fj[1];
              f[j][2] += fj[2];
              tj[0] = dyj*fj[2] - dzj*fj[1];
              tj[1] = dzj*fj[0] - dxj*fj[2];
              tj[2] = dxj*fj[1] - dyj*fj[0];
              torque[j][0] += tj[0];
              torque[j][1] += tj[1];
              torque[j][2] += tj[2];
            }
          }
        }

      // tri/particle interaction = Nx1 particles
      // c1,c2,c3 = corner pts of triangle I

      } else if (tri[i] >= 0) {

        if (dnum[i] == 0) {
          MathExtra::quat_to_mat(bonus[tri[i]].quat,p);
          MathExtra::matvec(p,bonus[tri[i]].c1,dc1);
          MathExtra::matvec(p,bonus[tri[i]].c2,dc2);
          MathExtra::matvec(p,bonus[tri[i]].c3,dc3);
          dfirst[i] = ndiscrete;
          discretize(i,sigma[itype][itype],dc1,dc2,dc3);
          dnum[i] = ndiscrete - dfirst[i];
        }
        npi = dnum[i];
        ifirst = dfirst[i];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;
          dzi = discrete[ifirst+ni].dz;

          xi[0] = x[i][0] + dxi;
          xi[1] = x[i][1] + dyi;
          xi[2] = x[i][2] + dzi;
          xj[0] = x[j][0];
          xj[1] = x[j][1];
          xj[2] = x[j][2];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          sig = 0.5 * (discrete[ifirst+ni].sigma+sigma[jtype][jtype]);
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] = delx*fpair;
          fi[1] = dely*fpair;
          fi[2] = delz*fpair;
          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          ti[0] = dyi*fi[2] - dzi*fi[1];
          ti[1] = dzi*fi[0] - dxi*fi[2];
          ti[2] = dxi*fi[1] - dyi*fi[0];
          torque[i][2] += ti[0];
          torque[i][1] += ti[1];
          torque[i][2] += ti[2];

          if (newton_pair || j < nlocal) {
            fj[0] = -delx*fpair;
            fj[1] = -dely*fpair;
            fj[2] = -delz*fpair;
            f[j][0] += fj[0];
            f[j][1] += fj[1];
            f[j][2] += fj[2];
          }
        }

      // particle/tri interaction = Nx1 particles
      // c1,c2,c3 = corner pts of triangle J

      } else if (tri[j] >= 0) {
        if (dnum[j] == 0) {
          MathExtra::quat_to_mat(bonus[tri[j]].quat,p);
          MathExtra::matvec(p,bonus[tri[j]].c1,dc1);
          MathExtra::matvec(p,bonus[tri[j]].c2,dc2);
          MathExtra::matvec(p,bonus[tri[j]].c3,dc3);
          dfirst[j] = ndiscrete;
          discretize(j,sigma[jtype][jtype],dc1,dc2,dc3);
          dnum[j] = ndiscrete - dfirst[j];
        }
        npj = dnum[j];
        jfirst = dfirst[j];

        for (nj = 0; nj < npj; nj++) {
          dxj = discrete[jfirst+nj].dx;
          dyj = discrete[jfirst+nj].dy;
          dzj = discrete[jfirst+nj].dz;

          xi[0] = x[i][0];
          xi[1] = x[i][1];
          xi[2] = x[i][2];
          xj[0] = x[j][0] + dxj;
          xj[1] = x[j][1] + dyj;
          xj[2] = x[j][2] + dzj;

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          sig = 0.5 * (sigma[itype][itype]+discrete[jfirst+nj].sigma);
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] = delx*fpair;
          fi[1] = dely*fpair;
          fi[2] = delz*fpair;
          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];

          if (newton_pair || j < nlocal) {
            fj[0] = -delx*fpair;
            fj[1] = -dely*fpair;
            fj[2] = -delz*fpair;
            f[j][0] += fj[0];
            f[j][1] += fj[1];
            f[j][2] += fj[2];
            tj[0] = dyj*fj[2] - dzj*fj[1];
            tj[1] = dzj*fj[0] - dxj*fj[2];
            tj[2] = dxj*fj[1] - dyj*fj[0];
            torque[j][0] += tj[0];
            torque[j][1] += tj[1];
            torque[j][2] += tj[2];
          }
        }

      // particle/particle interaction = 1x1 particles

      } else {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = forcelj*r2inv;

        if (eflag)
          evdwl += r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTriLJ::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTriLJ::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTriLJ::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTriLJ::init_style()
{
  avec = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  if (!avec) error->all(FLERR,"Pair tri/lj requires atom style tri");

  neighbor->add_request(this,NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTriLJ::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   recursively discretize triangle I with displaced corners c1,c2,c3
   into N sub-tris no more than sigma in size
   recurse by making 2 tris via bisecting longest side
   store new discrete particles in Discrete list
------------------------------------------------------------------------- */

void PairTriLJ::discretize(int i, double sigma,
                           double *c1, double *c2, double *c3)
{
  double centroid[3],dc1[3],dc2[3],dc3[3];

  centroid[0] = (c1[0] + c2[0] + c3[0]) / 3.0;
  centroid[1] = (c1[1] + c2[1] + c3[1]) / 3.0;
  centroid[2] = (c1[2] + c2[2] + c3[2]) / 3.0;

  MathExtra::sub3(c1,centroid,dc1);
  MathExtra::sub3(c2,centroid,dc2);
  MathExtra::sub3(c3,centroid,dc3);

  double sigmasq = 0.25 * sigma*sigma;
  double len1sq = MathExtra::lensq3(dc1);
  double len2sq = MathExtra::lensq3(dc2);
  double len3sq = MathExtra::lensq3(dc3);

  // if sigma sphere overlaps all corner points, add particle at centroid

  if ((len1sq <= sigmasq) && (len2sq <= sigmasq) && (len3sq <= sigmasq)) {
    if (ndiscrete == dmax) {
      dmax += DELTA;
      discrete = (Discrete *)
        memory->srealloc(discrete,dmax*sizeof(Discrete),"pair:discrete");
    }
    discrete[ndiscrete].dx = centroid[0];
    discrete[ndiscrete].dy = centroid[1];
    discrete[ndiscrete].dz = centroid[2];
    sigmasq = MAX(len1sq,len2sq);
    sigmasq = MAX(sigmasq,len3sq);
    discrete[ndiscrete].sigma = 2.0 * sqrt(sigmasq);
    ndiscrete++;
    return;
  }

  // else break triangle into 2 sub-triangles and recurse

  double c12[3],c23[3],c13[3],mid[3];

  MathExtra::sub3(c2,c3,c23);
  len1sq = MathExtra::lensq3(c23);
  MathExtra::sub3(c1,c3,c13);
  len2sq = MathExtra::lensq3(c13);
  MathExtra::sub3(c1,c2,c12);
  len3sq = MathExtra::lensq3(c12);

  double maxsq = MAX(len1sq,len2sq);
  maxsq = MAX(maxsq,len3sq);

  if (len1sq == maxsq) {
    MathExtra::add3(c2,c3,mid);
    MathExtra::scale3(0.5,mid);
    discretize(i,sigma,c1,c2,mid);
    discretize(i,sigma,c1,c3,mid);
  } else if (len2sq == maxsq) {
    MathExtra::add3(c1,c3,mid);
    MathExtra::scale3(0.5,mid);
    discretize(i,sigma,c2,c1,mid);
    discretize(i,sigma,c2,c3,mid);
  } else {
    MathExtra::add3(c1,c2,mid);
    MathExtra::scale3(0.5,mid);
    discretize(i,sigma,c3,c1,mid);
    discretize(i,sigma,c3,c2,mid);
  }
}

/* ----------------------------------------------------------------------
   recursively discretize triangle I with displaced corners c1,c2,c3
   into N sub-tris no more than sigma in size
   recurse by making 6 tris via centroid
   store new discrete particles in Discrete list
------------------------------------------------------------------------- */

/*
void PairTriLJ::discretize(int i, double sigma,
                          double *c1, double *c2, double *c3)
{
  double centroid[3],dc1[3],dc2[3],dc3[3];

  centroid[0] = (c1[0] + c2[0] + c3[0]) / 3.0;
  centroid[1] = (c1[1] + c2[1] + c3[1]) / 3.0;
  centroid[2] = (c1[2] + c2[2] + c3[2]) / 3.0;

  MathExtra::sub3(c1,centroid,dc1);
  MathExtra::sub3(c2,centroid,dc2);
  MathExtra::sub3(c3,centroid,dc3);

  double sigmasq = 0.25 * sigma*sigma;
  double len1sq = MathExtra::lensq3(dc1);
  double len2sq = MathExtra::lensq3(dc2);
  double len3sq = MathExtra::lensq3(dc3);

  // if sigma sphere overlaps all corner points, add particle at centroid

  if (len1sq <= sigmasq && len2sq <= sigmasq & len3sq <= sigmasq) {
    if (ndiscrete == dmax) {
      dmax += DELTA;
      discrete = (Discrete *)
        memory->srealloc(discrete,dmax*sizeof(Discrete),"pair:discrete");
    }
    discrete[ndiscrete].dx = centroid[0];
    discrete[ndiscrete].dy = centroid[1];
    discrete[ndiscrete].dz = centroid[2];
    sigmasq = MAX(len1sq,len2sq);
    sigmasq = MAX(sigmasq,len3sq);
    discrete[ndiscrete].sigma = 2.0 * sqrt(sigmasq);
    ndiscrete++;
    return;
  }

  // else break triangle into 6 sub-triangles and recurse

  double c1c2mid[3],c2c3mid[3],c1c3mid[3];

  MathExtra::add3(c1,c2,c1c2mid);
  MathExtra::scale3(0.5,c1c2mid);
  MathExtra::add3(c2,c3,c2c3mid);
  MathExtra::scale3(0.5,c2c3mid);
  MathExtra::add3(c1,c3,c1c3mid);
  MathExtra::scale3(0.5,c1c3mid);

  discretize(i,sigma,c1,c1c2mid,centroid);
  discretize(i,sigma,c1,c1c3mid,centroid);
  discretize(i,sigma,c2,c2c3mid,centroid);
  discretize(i,sigma,c2,c1c2mid,centroid);
  discretize(i,sigma,c3,c1c3mid,centroid);
  discretize(i,sigma,c3,c2c3mid,centroid);
}

*/
