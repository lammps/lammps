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

#include "pair_exp6_rx.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "modify.h"
#include "fix.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

#define MAXLINE 1024
#define DELTA 4

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define oneFluidApproxParameter (-1)
#define isOneFluidApprox(_site) ( (_site) == oneFluidApproxParameter )

#define exp6PotentialType (1)
#define isExp6PotentialType(_type) ( (_type) == exp6PotentialType )

// Create a structure to hold the parameter data for all
// local and neighbor particles. Pack inside this struct
// to avoid any name clashes.
struct PairExp6ParamDataType
{
   int n;
   double *epsilon1, *alpha1, *rm1, *mixWtSite1,
          *epsilon2, *alpha2, *rm2, *mixWtSite2,
          *epsilonOld1, *alphaOld1, *rmOld1, *mixWtSite1old,
          *epsilonOld2, *alphaOld2, *rmOld2, *mixWtSite2old;

   // Default constructor -- nullify everything.
   PairExp6ParamDataType(void)
      : n(0), epsilon1(NULL), alpha1(NULL), rm1(NULL), mixWtSite1(NULL),
              epsilon2(NULL), alpha2(NULL), rm2(NULL), mixWtSite2(NULL),
              epsilonOld1(NULL), alphaOld1(NULL), rmOld1(NULL), mixWtSite1old(NULL),
              epsilonOld2(NULL), alphaOld2(NULL), rmOld2(NULL), mixWtSite2old(NULL)
   {}
};

/* ---------------------------------------------------------------------- */

PairExp6rx::PairExp6rx(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;

  nspecies = 0;
  nparams = maxparam = 0;
  params = NULL;
  mol2param = NULL;
  fractionalWeighting = true;
}

/* ---------------------------------------------------------------------- */

PairExp6rx::~PairExp6rx()
{
  if (copymode) return;

  if (params != NULL) {
    for (int i=0; i < nparams; ++i) {
      delete[] params[i].name;
      delete[] params[i].potential;
    }
    memory->destroy(params);
  }
  memory->destroy(mol2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  if(scalingFlag == POLYNOMIAL){
    memory->destroy(coeffAlpha);
    memory->destroy(coeffEps);
    memory->destroy(coeffRm);
  }
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq,r2inv,r6inv,forceExp6,factor_lj;
  double rCut,rCutInv,rCut2inv,rCut6inv,rCutExp,urc,durc;
  double rm2ij,rm6ij;
  double r,rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwlOld = 0.0;
  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  double alphaOld12_ij, rmOld12_ij, epsilonOld12_ij;
  double alphaOld21_ij, rmOld21_ij, epsilonOld21_ij;
  double alpha12_ij, rm12_ij, epsilon12_ij;
  double alpha21_ij, rm21_ij, epsilon21_ij;
  double rminv, buck1, buck2;
  double epsilonOld1_i,alphaOld1_i,rmOld1_i;
  double epsilonOld1_j,alphaOld1_j,rmOld1_j;
  double epsilonOld2_i,alphaOld2_i,rmOld2_i;
  double epsilonOld2_j,alphaOld2_j,rmOld2_j;
  double epsilon1_i,alpha1_i,rm1_i;
  double epsilon1_j,alpha1_j,rm1_j;
  double epsilon2_i,alpha2_i,rm2_i;
  double epsilon2_j,alpha2_j,rm2_j;
  double evdwlOldEXP6_12, evdwlOldEXP6_21, fpairOldEXP6_12, fpairOldEXP6_21;
  double evdwlEXP6_12, evdwlEXP6_21;
  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  const int nRep = 12;
  const double shift = 1.05;
  double rin1, aRep, uin1, win1, uin1rep, rin1exp, rin6, rin6inv;

  // Initialize the Exp6 parameter data for both the local
  // and ghost atoms. Make the parameter data persistent
  // and exchange like any other atom property later.

  PairExp6ParamDataType PairExp6ParamData;

  {
     const int np_total = nlocal + atom->nghost;

     memory->create( PairExp6ParamData.epsilon1     , np_total, "PairExp6ParamData.epsilon1");
     memory->create( PairExp6ParamData.alpha1       , np_total, "PairExp6ParamData.alpha1");
     memory->create( PairExp6ParamData.rm1          , np_total, "PairExp6ParamData.rm1");
     memory->create( PairExp6ParamData.mixWtSite1    , np_total, "PairExp6ParamData.mixWtSite1");
     memory->create( PairExp6ParamData.epsilon2     , np_total, "PairExp6ParamData.epsilon2");
     memory->create( PairExp6ParamData.alpha2       , np_total, "PairExp6ParamData.alpha2");
     memory->create( PairExp6ParamData.rm2          , np_total, "PairExp6ParamData.rm2");
     memory->create( PairExp6ParamData.mixWtSite2    , np_total, "PairExp6ParamData.mixWtSite2");
     memory->create( PairExp6ParamData.epsilonOld1  , np_total, "PairExp6ParamData.epsilonOld1");
     memory->create( PairExp6ParamData.alphaOld1    , np_total, "PairExp6ParamData.alphaOld1");
     memory->create( PairExp6ParamData.rmOld1       , np_total, "PairExp6ParamData.rmOld1");
     memory->create( PairExp6ParamData.mixWtSite1old , np_total, "PairExp6ParamData.mixWtSite1old");
     memory->create( PairExp6ParamData.epsilonOld2  , np_total, "PairExp6ParamData.epsilonOld2");
     memory->create( PairExp6ParamData.alphaOld2    , np_total, "PairExp6ParamData.alphaOld2");
     memory->create( PairExp6ParamData.rmOld2       , np_total, "PairExp6ParamData.rmOld2");
     memory->create( PairExp6ParamData.mixWtSite2old , np_total, "PairExp6ParamData.mixWtSite2old");

     for (i = 0; i < np_total; ++i)
     {
        getMixingWeights (i, PairExp6ParamData.epsilon1[i],
                          PairExp6ParamData.alpha1[i],
                          PairExp6ParamData.rm1[i],
                          PairExp6ParamData.mixWtSite1[i],
                          PairExp6ParamData.epsilon2[i],
                          PairExp6ParamData.alpha2[i],
                          PairExp6ParamData.rm2[i],
                          PairExp6ParamData.mixWtSite2[i],
                          PairExp6ParamData.epsilonOld1[i],
                          PairExp6ParamData.alphaOld1[i],
                          PairExp6ParamData.rmOld1[i],
                          PairExp6ParamData.mixWtSite1old[i],
                          PairExp6ParamData.epsilonOld2[i],
                          PairExp6ParamData.alphaOld2[i],
                          PairExp6ParamData.rmOld2[i],
                          PairExp6ParamData.mixWtSite2old[i]);
     }
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    {
       epsilon1_i     = PairExp6ParamData.epsilon1[i];
       alpha1_i       = PairExp6ParamData.alpha1[i];
       rm1_i          = PairExp6ParamData.rm1[i];
       mixWtSite1_i    = PairExp6ParamData.mixWtSite1[i];
       epsilon2_i     = PairExp6ParamData.epsilon2[i];
       alpha2_i       = PairExp6ParamData.alpha2[i];
       rm2_i          = PairExp6ParamData.rm2[i];
       mixWtSite2_i    = PairExp6ParamData.mixWtSite2[i];
       epsilonOld1_i  = PairExp6ParamData.epsilonOld1[i];
       alphaOld1_i    = PairExp6ParamData.alphaOld1[i];
       rmOld1_i       = PairExp6ParamData.rmOld1[i];
       mixWtSite1old_i = PairExp6ParamData.mixWtSite1old[i];
       epsilonOld2_i  = PairExp6ParamData.epsilonOld2[i];
       alphaOld2_i    = PairExp6ParamData.alphaOld2[i];
       rmOld2_i       = PairExp6ParamData.rmOld2[i];
       mixWtSite2old_i = PairExp6ParamData.mixWtSite2old[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;

        r = sqrt(rsq);
        rCut2inv = 1.0/cutsq[itype][jtype];
        rCut6inv = rCut2inv*rCut2inv*rCut2inv;
        rCut = sqrt(cutsq[itype][jtype]);
        rCutInv = 1.0/rCut;

        //
        // A. Compute the exp-6 potential
        //

        // A1.  Get alpha, epsilon and rm for particle j

        {
           epsilon1_j     = PairExp6ParamData.epsilon1[j];
           alpha1_j       = PairExp6ParamData.alpha1[j];
           rm1_j          = PairExp6ParamData.rm1[j];
           mixWtSite1_j    = PairExp6ParamData.mixWtSite1[j];
           epsilon2_j     = PairExp6ParamData.epsilon2[j];
           alpha2_j       = PairExp6ParamData.alpha2[j];
           rm2_j          = PairExp6ParamData.rm2[j];
           mixWtSite2_j    = PairExp6ParamData.mixWtSite2[j];
           epsilonOld1_j  = PairExp6ParamData.epsilonOld1[j];
           alphaOld1_j    = PairExp6ParamData.alphaOld1[j];
           rmOld1_j       = PairExp6ParamData.rmOld1[j];
           mixWtSite1old_j = PairExp6ParamData.mixWtSite1old[j];
           epsilonOld2_j  = PairExp6ParamData.epsilonOld2[j];
           alphaOld2_j    = PairExp6ParamData.alphaOld2[j];
           rmOld2_j       = PairExp6ParamData.rmOld2[j];
           mixWtSite2old_j = PairExp6ParamData.mixWtSite2old[j];
        }

        // A2.  Apply Lorentz-Berthelot mixing rules for the i-j pair
        alphaOld12_ij = sqrt(alphaOld1_i*alphaOld2_j);
        rmOld12_ij = 0.5*(rmOld1_i + rmOld2_j);
        epsilonOld12_ij = sqrt(epsilonOld1_i*epsilonOld2_j);
        alphaOld21_ij = sqrt(alphaOld2_i*alphaOld1_j);
        rmOld21_ij = 0.5*(rmOld2_i + rmOld1_j);
        epsilonOld21_ij = sqrt(epsilonOld2_i*epsilonOld1_j);

        alpha12_ij = sqrt(alpha1_i*alpha2_j);
        rm12_ij = 0.5*(rm1_i + rm2_j);
        epsilon12_ij = sqrt(epsilon1_i*epsilon2_j);
        alpha21_ij = sqrt(alpha2_i*alpha1_j);
        rm21_ij = 0.5*(rm2_i + rm1_j);
        epsilon21_ij = sqrt(epsilon2_i*epsilon1_j);

        evdwlOldEXP6_12 = 0.0;
        evdwlOldEXP6_21 = 0.0;
        evdwlEXP6_12 = 0.0;
        evdwlEXP6_21 = 0.0;
        fpairOldEXP6_12 = 0.0;
        fpairOldEXP6_21 = 0.0;

        if(rmOld12_ij!=0.0 && rmOld21_ij!=0.0){
          if(alphaOld21_ij == 6.0 || alphaOld12_ij == 6.0)
            error->all(FLERR,"alpha_ij is 6.0 in pair exp6");

          // A3.  Compute some convenient quantities for evaluating the force
          rminv = 1.0/rmOld12_ij;
          buck1 = epsilonOld12_ij / (alphaOld12_ij - 6.0);
          rexp = expValue(alphaOld12_ij*(1.0-r*rminv));
          rm2ij = rmOld12_ij*rmOld12_ij;
          rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          rCutExp = expValue(alphaOld12_ij*(1.0-rCut*rminv));
          buck2 = 6.0*alphaOld12_ij;
          urc = buck1*(6.0*rCutExp - alphaOld12_ij*rm6ij*rCut6inv);
          durc = -buck1*buck2*(rCutExp* rminv - rCutInv*rm6ij*rCut6inv);
          rin1 = shift*rmOld12_ij*func_rin(alphaOld12_ij);
          if(r < rin1){
            rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            rin6inv = 1.0/rin6;

            rin1exp = expValue(alphaOld12_ij*(1.0-rin1*rminv));

            uin1 = buck1*(6.0*rin1exp - alphaOld12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            aRep = win1*powint(rin1,nRep)/nRep;

            uin1rep = aRep/powint(rin1,nRep);

            forceExp6 = double(nRep)*aRep/powint(r,nRep);
            fpairOldEXP6_12 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_12 = uin1 - uin1rep + aRep/powint(r,nRep);
          } else {
            forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
            fpairOldEXP6_12 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_12 = buck1*(6.0*rexp - alphaOld12_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }

          // A3.  Compute some convenient quantities for evaluating the force
          rminv = 1.0/rmOld21_ij;
          buck1 = epsilonOld21_ij / (alphaOld21_ij - 6.0);
          buck2 = 6.0*alphaOld21_ij;
          rexp = expValue(alphaOld21_ij*(1.0-r*rminv));
          rm2ij = rmOld21_ij*rmOld21_ij;
          rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          rCutExp = expValue(alphaOld21_ij*(1.0-rCut*rminv));
          buck2 = 6.0*alphaOld21_ij;
          urc = buck1*(6.0*rCutExp - alphaOld21_ij*rm6ij*rCut6inv);
          durc = -buck1*buck2*(rCutExp* rminv - rCutInv*rm6ij*rCut6inv);
          rin1 = shift*rmOld21_ij*func_rin(alphaOld21_ij);

          if(r < rin1){
            rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            rin6inv = 1.0/rin6;

            rin1exp = expValue(alphaOld21_ij*(1.0-rin1*rminv));

            uin1 = buck1*(6.0*rin1exp - alphaOld21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            aRep = win1*powint(rin1,nRep)/nRep;

            uin1rep = aRep/powint(rin1,nRep);

            forceExp6 = double(nRep)*aRep/powint(r,nRep);
            fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_21 = uin1 - uin1rep + aRep/powint(r,nRep);
          } else {
            forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
            fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_21 = buck1*(6.0*rexp - alphaOld21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }

          if (isite1 == isite2)
            evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12;
          else
            evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*evdwlOldEXP6_21;

          evdwlOld *= factor_lj;

          uCG[i] += 0.5*evdwlOld;
          if (newton_pair || j < nlocal)
            uCG[j] += 0.5*evdwlOld;
        }

        if(rm12_ij!=0.0 && rm21_ij!=0.0){
          if(alpha21_ij == 6.0 || alpha12_ij == 6.0)
            error->all(FLERR,"alpha_ij is 6.0 in pair exp6");

          // A3.  Compute some convenient quantities for evaluating the force
          rminv = 1.0/rm12_ij;
          buck1 = epsilon12_ij / (alpha12_ij - 6.0);
          buck2 = 6.0*alpha12_ij;
          rexp = expValue(alpha12_ij*(1.0-r*rminv));
          rm2ij = rm12_ij*rm12_ij;
          rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          rCutExp = expValue(alpha12_ij*(1.0-rCut*rminv));
          urc = buck1*(6.0*rCutExp - alpha12_ij*rm6ij*rCut6inv);
          durc = -buck1*buck2*(rCutExp*rminv - rCutInv*rm6ij*rCut6inv);
          rin1 = shift*rm12_ij*func_rin(alpha12_ij);

          if(r < rin1){
            rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            rin6inv = 1.0/rin6;

            rin1exp = expValue(alpha12_ij*(1.0-rin1*rminv));

            uin1 = buck1*(6.0*rin1exp - alpha12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            aRep = win1*powint(rin1,nRep)/nRep;

            uin1rep = aRep/powint(rin1,nRep);

            evdwlEXP6_12 = uin1 - uin1rep + aRep/powint(r,nRep);
          } else {
            evdwlEXP6_12 = buck1*(6.0*rexp - alpha12_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }

          rminv = 1.0/rm21_ij;
          buck1 = epsilon21_ij / (alpha21_ij - 6.0);
          buck2 = 6.0*alpha21_ij;
          rexp = expValue(alpha21_ij*(1.0-r*rminv));
          rm2ij = rm21_ij*rm21_ij;
          rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          rCutExp = expValue(alpha21_ij*(1.0-rCut*rminv));
          urc = buck1*(6.0*rCutExp - alpha21_ij*rm6ij*rCut6inv);
          durc = -buck1*buck2*(rCutExp*rminv - rCutInv*rm6ij*rCut6inv);
          rin1 = shift*rm21_ij*func_rin(alpha21_ij);

          if(r < rin1){
            rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            rin6inv = 1.0/rin6;

            rin1exp = expValue(alpha21_ij*(1.0-rin1*rminv));

            uin1 = buck1*(6.0*rin1exp - alpha21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            aRep = win1*powint(rin1,nRep)/nRep;

            uin1rep = aRep/powint(rin1,nRep);

            evdwlEXP6_21 = uin1 - uin1rep + aRep/powint(r,nRep);
          } else {
            evdwlEXP6_21 = buck1*(6.0*rexp - alpha21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }
        }

        //
        // Apply Mixing Rule to get the overall force for the CG pair
        //
        if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12;
        else fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*fpairOldEXP6_21;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (isite1 == isite2) evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12;
        else evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12 + sqrt(mixWtSite2_i*mixWtSite1_j)*evdwlEXP6_21;
        evdwl *= factor_lj;

        uCGnew[i]   += 0.5*evdwl;
        if (newton_pair || j < nlocal)
          uCGnew[j] += 0.5*evdwl;
        evdwl = evdwlOld;
        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();

  // Release the local parameter data.
  {
     if (PairExp6ParamData.epsilon1    ) memory->destroy(PairExp6ParamData.epsilon1);
     if (PairExp6ParamData.alpha1      ) memory->destroy(PairExp6ParamData.alpha1);
     if (PairExp6ParamData.rm1         ) memory->destroy(PairExp6ParamData.rm1);
     if (PairExp6ParamData.mixWtSite1   ) memory->destroy(PairExp6ParamData.mixWtSite1);
     if (PairExp6ParamData.epsilon2    ) memory->destroy(PairExp6ParamData.epsilon2);
     if (PairExp6ParamData.alpha2      ) memory->destroy(PairExp6ParamData.alpha2);
     if (PairExp6ParamData.rm2         ) memory->destroy(PairExp6ParamData.rm2);
     if (PairExp6ParamData.mixWtSite2   ) memory->destroy(PairExp6ParamData.mixWtSite2);
     if (PairExp6ParamData.epsilonOld1 ) memory->destroy(PairExp6ParamData.epsilonOld1);
     if (PairExp6ParamData.alphaOld1   ) memory->destroy(PairExp6ParamData.alphaOld1);
     if (PairExp6ParamData.rmOld1      ) memory->destroy(PairExp6ParamData.rmOld1);
     if (PairExp6ParamData.mixWtSite1old) memory->destroy(PairExp6ParamData.mixWtSite1old);
     if (PairExp6ParamData.epsilonOld2 ) memory->destroy(PairExp6ParamData.epsilonOld2);
     if (PairExp6ParamData.alphaOld2   ) memory->destroy(PairExp6ParamData.alphaOld2);
     if (PairExp6ParamData.rmOld2      ) memory->destroy(PairExp6ParamData.rmOld2);
     if (PairExp6ParamData.mixWtSite2old) memory->destroy(PairExp6ParamData.mixWtSite2old);
  }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairExp6rx::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut_lj");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairExp6rx::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // optional keywords

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fractional") == 0) fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0) fractionalWeighting = false;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  allocated = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExp6rx::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 9) error->all(FLERR,"Incorrect args for pair coefficients");

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairExp6rx requires a fix rx command.");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  int n;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  nspecies = atom->nspecies_dpd;
  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");
  read_file(arg[2]);

  n = strlen(arg[3]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[3]);

  int ispecies;
  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site1,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site1,"1fluid") != 0)
    error->all(FLERR,"Site1 name not recognized in pair coefficients");

  n = strlen(arg[4]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[4]);

  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site2,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site2,"1fluid") != 0)
    error->all(FLERR,"Site2 name not recognized in pair coefficients");

  {
    // Set isite1 and isite2 parameters based on site1 and site2 strings.

    if (strcmp(site1,"1fluid") == 0)
      isite1 = oneFluidApproxParameter;
    else
      {
        int isp;
        for (isp = 0; isp < nspecies; isp++)
          if (strcmp(site1, &atom->dname[isp][0]) == 0) break;

        if (isp == nspecies)
          error->all(FLERR,"Site1 name not recognized in pair coefficients");
        else
          isite1 = isp;
      }

    if (strcmp(site2,"1fluid") == 0)
      isite2 = oneFluidApproxParameter;
    else
      {
        int isp;
        for (isp = 0; isp < nspecies; isp++)
        if (strcmp(site2, &atom->dname[isp][0]) == 0) break;

        if (isp == nspecies)
          error->all(FLERR,"Site2 name not recognized in pair coefficients");
        else
          isite2 = isp;
      }

    // Set the interaction potential type to the enumerated type.
    for (int iparam = 0; iparam < nparams; ++iparam)
      {
        if (strcmp( params[iparam].potential, "exp6") == 0)
          params[iparam].potentialType = exp6PotentialType;
        else
          error->all(FLERR,"params[].potential type unknown");
      }
  }
  delete[] site1;
  delete[] site2;
  site1 = site2 = NULL;

  setup();

  double cut_one = cut_global;
  if (strcmp(arg[5],"exponent") == 0){
    scalingFlag = EXPONENT;
    exponentR = force->numeric(FLERR,arg[6]);
    exponentEpsilon = force->numeric(FLERR,arg[7]);
    if (narg > 9) error->all(FLERR,"Incorrect args for pair coefficients");
    if (narg == 9) cut_one = force->numeric(FLERR,arg[8]);
  } else if (strcmp(arg[5],"polynomial") == 0){
    scalingFlag = POLYNOMIAL;
    memory->create(coeffAlpha,6,"pair:coeffAlpha");
    memory->create(coeffEps,6,"pair:coeffEps");
    memory->create(coeffRm,6,"pair:coeffRm");
    read_file2(arg[6]);
    if (narg > 8) error->all(FLERR,"Incorrect args for pair coefficients");
    if (narg == 8) cut_one = force->numeric(FLERR,arg[7]);
  } else if (strcmp(arg[5],"none") == 0){
    scalingFlag = NONE;
    if (narg > 7) error->all(FLERR,"Incorrect args for pair coefficients");
    if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);
  } else {
    error->all(FLERR,"Incorrect args for pair coefficients");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExp6rx::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::read_file(char *file)
{
  int params_per_line = 5;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open exp6/rx potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines

  int n,nwords,ispecies;
  char line[MAXLINE],*ptr;
  int eof = 0;
  char *r_token;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in exp6/rx potential file");

    // words = ptrs to all words in line

    r_token = line;
    nwords = 0;
    words[nwords++] = utils::strtok_r(r_token," \t\n\r\f",&r_token);
    while ((words[nwords++] = utils::strtok_r(NULL," \t\n\r\f",&r_token))) continue;

    for (ispecies = 0; ispecies < nspecies; ispecies++)
      if (strcmp(words[0],&atom->dname[ispecies][0]) == 0) break;
    if (ispecies == nspecies) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ispecies = ispecies;

    n = strlen(&atom->dname[ispecies][0]) + 1;
    params[nparams].name = new char[n];
    strcpy(params[nparams].name,&atom->dname[ispecies][0]);

    n = strlen(words[1]) + 1;
    params[nparams].potential = new char[n];
    strcpy(params[nparams].potential,words[1]);
    if (strcmp(params[nparams].potential,"exp6") == 0){
      params[nparams].alpha = atof(words[2]);
      params[nparams].epsilon = atof(words[3]);
      params[nparams].rm = atof(words[4]);
      if (params[nparams].epsilon <= 0.0 || params[nparams].rm <= 0.0 ||
          params[nparams].alpha < 0.0)
        error->all(FLERR,"Illegal exp6/rx parameters.  Rm and Epsilon must be greater than zero.  Alpha cannot be negative.");
    } else {
      error->all(FLERR,"Illegal exp6/rx parameters.  Interaction potential does not exist.");
    }
    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::read_file2(char *file)
{
  int params_per_line = 7;
  char **words = new char*[params_per_line+1];

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open polynomial file %s",file);
      error->one(FLERR,str);
    }
  }

  // one set of params can span multiple lines
  int n,nwords;
  char line[MAXLINE],*ptr;
  int eof = 0;
  char *r_token;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in polynomial file");

    // words = ptrs to all words in line

    r_token = line;
    nwords = 0;
    words[nwords++] = utils::strtok_r(r_token," \t\n\r\f",&r_token);
    while ((words[nwords++] = utils::strtok_r(NULL," \t\n\r\f",&r_token))) continue;

    if (strcmp(words[0],"alpha") == 0){
      for (int ii=1; ii<params_per_line; ii++)
        coeffAlpha[ii-1] = atof(words[ii]);
    }
    if (strcmp(words[0],"epsilon") == 0){
      for (int ii=1; ii<params_per_line; ii++)
        coeffEps[ii-1] = atof(words[ii]);
    }
    if (strcmp(words[0],"rm") == 0){
      for (int ii=1; ii<params_per_line; ii++)
        coeffRm[ii-1] = atof(words[ii]);
    }
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::setup()
{
  int i,j,n;

  // set mol2param for all combinations
  // must be a single exact match to lines read from file

  memory->destroy(mol2param);
  memory->create(mol2param,nspecies,"pair:mol2param");

  for (i = 0; i < nspecies; i++) {
    n = -1;
    for (j = 0; j < nparams; j++) {
      if (i == params[j].ispecies) {
        if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
        n = j;
      }
    }
    mol2param[i] = n;
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExp6rx::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExp6rx::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExp6rx::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExp6rx::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::getMixingWeights(int id,double &epsilon1,double &alpha1,double &rm1, double &mixWtSite1,double &epsilon2,double &alpha2,double &rm2,double &mixWtSite2,double &epsilon1_old,double &alpha1_old,double &rm1_old, double &mixWtSite1old,double &epsilon2_old,double &alpha2_old,double &rm2_old,double &mixWtSite2old) const
{
  int iparam, jparam;
  double rmi, rmj, rmij, rm3ij;
  double epsiloni, epsilonj, epsilonij;
  double alphai, alphaj, alphaij;
  double epsilon_old, rm3_old, alpha_old;
  double epsilon, rm3, alpha;
  double xMolei, xMolej, xMolei_old, xMolej_old;

  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalOld;

  rm3 = 0.0;
  epsilon = 0.0;
  alpha = 0.0;
  epsilon_old = 0.0;
  rm3_old = 0.0;
  alpha_old = 0.0;
  fractionOFA = 0.0;
  fractionOFAold = 0.0;
  nMoleculesOFA = 0.0;
  nMoleculesOFAold = 0.0;
  nTotal = 0.0;
  nTotalOld = 0.0;

  // Compute the total number of molecules in the old and new CG particle as well as the total number of molecules in the fluid portion of the old and new CG particle
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    nTotal += atom->dvector[ispecies][id];
    nTotalOld += atom->dvector[ispecies+nspecies][id];

    iparam = mol2param[ispecies];

    if (iparam < 0 || params[iparam].potentialType != exp6PotentialType ) continue;
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == params[iparam].ispecies || isite2 == params[iparam].ispecies) continue;
      nMoleculesOFAold += atom->dvector[ispecies+nspecies][id];
      nMoleculesOFA += atom->dvector[ispecies][id];
    }
  }
  if(nTotal < MY_EPSILON || nTotalOld < MY_EPSILON)
    error->all(FLERR,"The number of molecules in CG particle is less than 10*DBL_EPSILON.");

  // Compute the mole fraction of molecules within the fluid portion of the particle (One Fluid Approximation)
  fractionOFAold = nMoleculesOFAold / nTotalOld;
  fractionOFA = nMoleculesOFA / nTotal;

  for (int ispecies = 0; ispecies < nspecies; ispecies++) {
    iparam = mol2param[ispecies];
    if (iparam < 0 || params[iparam].potentialType != exp6PotentialType ) continue;

    // If Site1 matches a pure species, then grab the parameters
    if (isite1 == params[iparam].ispecies){
      rm1_old = params[iparam].rm;
      rm1 = params[iparam].rm;
      epsilon1_old = params[iparam].epsilon;
      epsilon1 = params[iparam].epsilon;
      alpha1_old = params[iparam].alpha;
      alpha1 = params[iparam].alpha;

      // Compute the mole fraction of Site1
      nMoleculesOld1 = atom->dvector[ispecies+nspecies][id];
      nMolecules1 = atom->dvector[ispecies][id];
      fractionOld1 = nMoleculesOld1/nTotalOld;
      fraction1 = nMolecules1/nTotal;
    }

    // If Site2 matches a pure species, then grab the parameters
    if (isite2 == params[iparam].ispecies){
      rm2_old = params[iparam].rm;
      rm2 = params[iparam].rm;
      epsilon2_old = params[iparam].epsilon;
      epsilon2 = params[iparam].epsilon;
      alpha2_old = params[iparam].alpha;
      alpha2 = params[iparam].alpha;

      // Compute the mole fraction of Site2
      nMoleculesOld2 = atom->dvector[ispecies+nspecies][id];
      nMolecules2 = atom->dvector[ispecies][id];
      fractionOld2 = atom->dvector[ispecies+nspecies][id]/nTotalOld;
      fraction2 = atom->dvector[ispecies][id]/nTotal;
    }

    // If Site1 or Site2 matches is a fluid, then compute the parameters
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == params[iparam].ispecies || isite2 == params[iparam].ispecies) continue;
      rmi = params[iparam].rm;
      epsiloni = params[iparam].epsilon;
      alphai = params[iparam].alpha;
      if(nMoleculesOFA<MY_EPSILON) xMolei = 0.0;
      else xMolei = atom->dvector[ispecies][id]/nMoleculesOFA;
      if(nMoleculesOFAold<MY_EPSILON) xMolei_old = 0.0;
      else xMolei_old = atom->dvector[ispecies+nspecies][id]/nMoleculesOFAold;

      for (int jspecies = 0; jspecies < nspecies; jspecies++) {
        jparam = mol2param[jspecies];
        if (jparam < 0 || params[jparam].potentialType != exp6PotentialType ) continue;
        if (isite1 == params[jparam].ispecies || isite2 == params[jparam].ispecies) continue;
        rmj = params[jparam].rm;
        epsilonj = params[jparam].epsilon;
        alphaj = params[jparam].alpha;
        if(nMoleculesOFA<MY_EPSILON) xMolej = 0.0;
        else xMolej = atom->dvector[jspecies][id]/nMoleculesOFA;
        if(nMoleculesOFAold<MY_EPSILON) xMolej_old = 0.0;
        else xMolej_old = atom->dvector[jspecies+nspecies][id]/nMoleculesOFAold;

        rmij = (rmi+rmj)/2.0;
        rm3ij = rmij*rmij*rmij;
        epsilonij = sqrt(epsiloni*epsilonj);
        alphaij = sqrt(alphai*alphaj);

        if(fractionOFAold > 0.0){
          rm3_old += xMolei_old*xMolej_old*rm3ij;
          epsilon_old += xMolei_old*xMolej_old*rm3ij*epsilonij;
          alpha_old += xMolei_old*xMolej_old*rm3ij*epsilonij*alphaij;
        }
        if(fractionOFA > 0.0){
          rm3 += xMolei*xMolej*rm3ij;
          epsilon += xMolei*xMolej*rm3ij*epsilonij;
          alpha += xMolei*xMolej*rm3ij*epsilonij*alphaij;
        }
      }
    }
  }

  if (isOneFluidApprox(isite1)){
    rm1 = cbrt(rm3);
    if(rm1 < MY_EPSILON) {
      rm1 = 0.0;
      epsilon1 = 0.0;
      alpha1 = 0.0;
    } else {
      epsilon1 = epsilon / rm3;
      alpha1 = alpha / epsilon1 / rm3;
    }
    nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
    fraction1 = fractionOFA;

    rm1_old = cbrt(rm3_old);
    if(rm1_old < MY_EPSILON) {
      rm1_old = 0.0;
      epsilon1_old = 0.0;
      alpha1_old = 0.0;
    } else {
      epsilon1_old = epsilon_old / rm3_old;
      alpha1_old = alpha_old / epsilon1_old / rm3_old;
    }
    nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
    fractionOld1 = fractionOFAold;

    if(scalingFlag == EXPONENT){
      exponentScaling(nMoleculesOFA,epsilon1,rm1);
      exponentScaling(nMoleculesOFAold,epsilon1_old,rm1_old);
    } else if(scalingFlag == POLYNOMIAL){
      polynomialScaling(nMoleculesOFA,alpha1,epsilon1,rm1);
      polynomialScaling(nMoleculesOFAold,alpha1_old,epsilon1_old,rm1_old);
    }
  }

  if (isOneFluidApprox(isite2)){
    rm2 = cbrt(rm3);
    if(rm2 < MY_EPSILON) {
      rm2 = 0.0;
      epsilon2 = 0.0;
      alpha2 = 0.0;
    } else {
      epsilon2 = epsilon / rm3;
      alpha2 = alpha / epsilon2 / rm3;
    }
    nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
    fraction2 = fractionOFA;

    rm2_old = cbrt(rm3_old);
    if(rm2_old < MY_EPSILON) {
      rm2_old = 0.0;
      epsilon2_old = 0.0;
      alpha2_old = 0.0;
    } else {
      epsilon2_old = epsilon_old / rm3_old;
      alpha2_old = alpha_old / epsilon2_old / rm3_old;
    }
    nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
    fractionOld2 = fractionOFAold;

    if(scalingFlag == EXPONENT){
      exponentScaling(nMoleculesOFA,epsilon2,rm2);
      exponentScaling(nMoleculesOFAold,epsilon2_old,rm2_old);
    } else if(scalingFlag == POLYNOMIAL){
      polynomialScaling(nMoleculesOFA,alpha2,epsilon2,rm2);
      polynomialScaling(nMoleculesOFAold,alpha2_old,epsilon2_old,rm2_old);
    }
  }

  // Check that no fractions are less than zero
  if(fraction1 < 0.0 || nMolecules1 < 0.0){
    if(fraction1 < -MY_EPSILON || nMolecules1 < -MY_EPSILON){
      error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");
    }
    nMolecules1 = 0.0;
    fraction1 = 0.0;
  }
  if(fraction2 < 0.0 || nMolecules2 < 0.0){
    if(fraction2 < -MY_EPSILON || nMolecules2 < -MY_EPSILON){
      error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");
    }
    nMolecules2 = 0.0;
    fraction2 = 0.0;
  }
  if(fractionOld1 < 0.0 || nMoleculesOld1 < 0.0){
    if(fractionOld1 < -MY_EPSILON || nMoleculesOld1 < -MY_EPSILON){
      error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");
    }
    nMoleculesOld1 = 0.0;
    fractionOld1 = 0.0;
  }
  if(fractionOld2 < 0.0 || nMoleculesOld2 < 0.0){
    if(fractionOld2 < -MY_EPSILON || nMoleculesOld2 < -MY_EPSILON){
      error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");
    }
    nMoleculesOld2 = 0.0;
    fractionOld2 = 0.0;
  }

  if(fractionalWeighting){
    mixWtSite1old = fractionOld1;
    mixWtSite1 = fraction1;
    mixWtSite2old = fractionOld2;
    mixWtSite2 = fraction2;
  } else {
    mixWtSite1old = nMoleculesOld1;
    mixWtSite1 = nMolecules1;
    mixWtSite2old = nMoleculesOld2;
    mixWtSite2 = nMolecules2;
  }
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::exponentScaling(double phi, double &epsilon, double &rm) const
{
  double powfuch;

  if(exponentEpsilon < 0.0){
    powfuch = pow(phi,-exponentEpsilon);
    if(powfuch<MY_EPSILON) epsilon = 0.0;
    else epsilon *= 1.0/powfuch;
  } else {
    epsilon *= pow(phi,exponentEpsilon);
  }

  if(exponentR < 0.0){
    powfuch = pow(phi,-exponentR);
    if(powfuch<MY_EPSILON) rm = 0.0;
    else rm *= 1.0/powfuch;
  } else {
    rm *= pow(phi,exponentR);
  }
}

/* ---------------------------------------------------------------------- */

void PairExp6rx::polynomialScaling(double phi, double &alpha, double &epsilon, double &rm) const
{
    double phi2 = phi*phi;
    double phi3 = phi2*phi;
    double phi4 = phi2*phi2;
    double phi5 = phi2*phi3;

    alpha = (coeffAlpha[0]*phi5 + coeffAlpha[1]*phi4 + coeffAlpha[2]*phi3 + coeffAlpha[3]*phi2 + coeffAlpha[4]*phi + coeffAlpha[5]);
    epsilon *= (coeffEps[0]*phi5 + coeffEps[1]*phi4 + coeffEps[2]*phi3 + coeffEps[3]*phi2 + coeffEps[4]*phi + coeffEps[5]);
    rm *= (coeffRm[0]*phi5 + coeffRm[1]*phi4 + coeffRm[2]*phi3 + coeffRm[3]*phi2 + coeffRm[4]*phi + coeffRm[5]);
}

/* ---------------------------------------------------------------------- */

inline double PairExp6rx::func_rin(const double &alpha) const
{
  double function;

  const double a = 3.7682065;
  const double b = -1.4308614;

  function = a+b*sqrt(alpha);
  function = expValue(function);

  return function;
}

/* ---------------------------------------------------------------------- */

inline double PairExp6rx::expValue(double value) const
{
  double returnValue;
  if(value < DBL_MIN_EXP) returnValue = 0.0;
  else returnValue = exp(value);

  return returnValue;
}
