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
   Contributing authors: Amit Kumar and Michael Bybee (UIUC)
                         Dave Heine (Corning), polydispersity
------------------------------------------------------------------------- */

#include "pair_brownian_poly.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "math_special.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstdint>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

/* ----------------------------------------------------------------------
   deterministic, traversal- and MPI-rank-independent uniform random number
   in [0,1), keyed on the *unordered* atom-tag pair, the timestep, the style
   seed, and a stream index k.  Both atoms of a pair (and either MPI rank that
   owns them under "newton off") draw the identical value, so the pairwise
   Brownian force can be applied equal-and-opposite and obeys Newton's 3rd law
   (conserves momentum).  See GitHub issue #2933.  A per-rank RNG sequence
   cannot be used here because the two halves of a ghosted pair are evaluated
   on different ranks.  Mixing uses a boost-style hash_combine followed by the
   splitmix64 finalizer for good avalanche.
------------------------------------------------------------------------- */

double PairBrownianPoly::pair_uniform(tagint ti, tagint tj, bigint step, int seed, int k)
{
  uint64_t lo = (uint64_t) (ti < tj ? ti : tj);
  uint64_t hi = (uint64_t) (ti < tj ? tj : ti);
  uint64_t h = (uint64_t) seed * 0x9E3779B97F4A7C15ULL;
  h ^= lo + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  h ^= hi + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  h ^= (uint64_t) step + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  h ^= (uint64_t) k + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  h += 0x9E3779B97F4A7C15ULL;
  h = (h ^ (h >> 30)) * 0xBF58476D1CE4E5B9ULL;
  h = (h ^ (h >> 27)) * 0x94D049BB133111EBULL;
  h = h ^ (h >> 31);
  return (double) (h >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------------------------------------------------------------------- */

PairBrownianPoly::PairBrownianPoly(LAMMPS *lmp) : PairBrownian(lmp)
{
  no_virial_fdotr_compute = 1;
  rad = 0.0; // set to a default value
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBrownianPoly::settings(int narg, char **arg)
{
  PairBrownian::settings(narg, arg);
  // NOTE: the code for volume fraction correction was copied from pair style brownian,
  // which requires a uniform radius (stored in the variable rad). For a polydisperse
  // system that is not correct and the variable rad unset. Thus we stop here with an error.
  if (flagVF)
    error->all(FLERR, "Pair style brownian/poly does not support volume fraction corrections");
}

/* ---------------------------------------------------------------------- */

void PairBrownianPoly::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  bigint step = update->ntimestep;

  double vxmu2f = force->vxmu2f;
  double prethermostat;
  double a_sq,a_sh,a_pu,Fbmag;
  double p1[3],p2[3],p3[3];

  // this section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2) { // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)) {
        double wallhi[3], walllo[3];
        for (j = 0; j < 3; j++) {
          wallhi[j] = domain->prd[j];
          walllo[j] = 0;
        }
        for (int m = 0; m < wallfix->nwall; m++) {
          int dim = wallfix->wallwhich[m] / 2;
          int side = wallfix->wallwhich[m] % 2;
          if (wallfix->xstyle[m] == FixWall::VARIABLE) {
            wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
          }
          else wallcoord = wallfix->coord0[m];
          if (side == 0) walllo[dim] = wallcoord;
          else wallhi[dim] = wallcoord;
        }
        for (j = 0; j < 3; j++)
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
      f[i][0] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      f[i][1] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      f[i][2] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      if (flaglog) {
        const double radi3 = radi*radi*radi;
        torque[i][0] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
        torque[i][1] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
        torque[i][2] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
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
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // canonical (tag-ordered) description of the pair so that the random
        // Brownian force is computed identically no matter which atom -- or,
        // under "newton off", which MPI rank -- evaluates the pair.  All
        // geometry, resistances and random draws below use the lower-tag atom
        // as the reference particle.  csgn = +1 if the local atom i is that
        // reference, -1 otherwise; (ex,ey,ez) is the unit line-of-centers
        // pointing from the higher-tag toward the lower-tag atom.  This makes
        // the pairwise force exactly equal and opposite (Newton's 3rd law),
        // fixing the momentum/energy injection of GitHub issue #2933.

        double csgn = (tag[i] < tag[j]) ? 1.0 : -1.0;
        double rref = (csgn > 0.0) ? radi : radj;   // lower-tag (reference) radius
        double roth = (csgn > 0.0) ? radj : radi;   // higher-tag radius
        double ex = csgn*delx/r;
        double ey = csgn*dely/r;
        double ez = csgn*delz/r;

        // scalar resistances a_sq, a_sh, a_pu

        h_sep = r - radi-radj;

        // if less than minimum gap, use minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // scale h_sep by the reference radius

        h_sep = h_sep/rref;
        beta0 = roth/rref;
        beta1 = 1.0 + beta0;

        // scalar resistances.  By the fluctuation-dissipation theorem the
        // Brownian variance equals the lubrication friction, so these use the
        // same symmetric-gap Jeffrey & Onishi expansion as pair lubricate/poly
        // (xi = 2*gap/(radi+radj) = 2*h_sep/beta1) after the GitHub issue
        // #1933 fix; using bare h_sep / dropping the beta0 prefactor would
        // make the magnitude depend on which particle is the reference.

        if (flaglog) {
          double xi = 2.0*h_sep/beta1;
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            beta0*(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/xi);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *
                   pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/xi);
          a_sq *= 6.0*MY_PI*mu*rref;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/xi);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/xi);
          a_sh *= 6.0*MY_PI*mu*rref;
          a_pu = 2.0*beta0/5.0/beta1*log(1.0/xi);
          a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/xi);
          a_pu *= 8.0*MY_PI*mu*pow(rref,3.0);
        } else a_sq = 6.0*MY_PI*mu*rref*(beta0*beta0/beta1/beta1/h_sep);

        // build the random Brownian force on the *reference* (lower-tag) atom.
        // random scalars are drawn from the order-independent pair RNG, one
        // independent stream (index k) per component, so both halves of the
        // pair agree on the value.

        // squeeze term (along line of centers): a_sq

        Fbmag = prethermostat*sqrt(a_sq);
        double s0 = pair_uniform(tag[i],tag[j],step,seed,0)-0.5;
        fx = Fbmag*s0*ex;
        fy = Fbmag*s0*ey;
        fz = Fbmag*s0*ez;

        // shear terms: a_sh

        if (flaglog) {

          // two orthogonal vectors to the (canonical) line of centers

          p1[0] = ex; p1[1] = ey; p1[2] = ez;
          set_3_orthogonal_vectors(p1,p2,p3);

          Fbmag = prethermostat*sqrt(a_sh);

          double s2 = pair_uniform(tag[i],tag[j],step,seed,1)-0.5;
          fx += Fbmag*s2*p2[0];
          fy += Fbmag*s2*p2[1];
          fz += Fbmag*s2*p2[2];

          double s3 = pair_uniform(tag[i],tag[j],step,seed,2)-0.5;
          fx += Fbmag*s3*p3[0];
          fy += Fbmag*s3*p3[1];
          fz += Fbmag*s3*p3[2];
        }

        // scale to force units

        fx *= vxmu2f;
        fy *= vxmu2f;
        fz *= vxmu2f;

        // (fx,fy,fz) is the negated Brownian force on the reference atom.
        // apply equal and opposite to the local atom i via csgn; the other
        // atom accumulates its (opposite) share when it is processed as a
        // local atom, so f[j] is intentionally not touched (newton off).

        f[i][0] -= csgn*fx;
        f[i][1] -= csgn*fy;
        f[i][2] -= csgn*fz;

        if (flaglog) {

          // torque on the local atom from the Brownian force acting at the
          // point of closest approach: tau = (radius_i * e) x F_i, which here
          // reduces to radi*(e x (fx,fy,fz)) for either pair member (only the
          // shear part contributes; the squeeze part is parallel to e).

          tx = radi*(ey*fz - ez*fy);
          ty = radi*(ez*fx - ex*fz);
          tz = radi*(ex*fy - ey*fx);

          torque[i][0] += tx;
          torque[i][1] += ty;
          torque[i][2] += tz;

          // pumping (rotational) Brownian torque: a_pu.  This is an
          // antisymmetric pair torque (opposite sign on the two particles),
          // so it carries csgn.  Note: as in the original code it is not
          // scaled by vxmu2f.

          Fbmag = prethermostat*sqrt(a_pu);

          double t2 = pair_uniform(tag[i],tag[j],step,seed,3)-0.5;
          tx = Fbmag*t2*p2[0];
          ty = Fbmag*t2*p2[1];
          tz = Fbmag*t2*p2[2];

          double t3 = pair_uniform(tag[i],tag[j],step,seed,4)-0.5;
          tx += Fbmag*t3*p3[0];
          ty += Fbmag*t3*p3[1];
          tz += Fbmag*t3*p3[2];

          torque[i][0] -= csgn*tx;
          torque[i][1] -= csgn*ty;
          torque[i][2] -= csgn*tz;
        }

        // tally only the local atom's contribution (j = nlocal, newton 0)

        if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,
                                 0.0,0.0,-csgn*fx,-csgn*fy,-csgn*fz,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBrownianPoly::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair brownian/poly requires newton pair off");
  if (!atom->radius_flag)
    error->all(FLERR,"Pair brownian/poly requires atom attribute radius");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair brownian/poly requires atom IDs");

  // ensure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair brownian/poly requires extended particles");

  neighbor->add_request(this, NeighConst::REQ_FULL);

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
  wallfix = nullptr;

  if (modify->get_fix_by_style("^deform").size() > 0) flagdeform = 1;
  auto fixes = modify->get_fix_by_style("^wall");
  if (fixes.size() > 1)
    error->all(FLERR, "Cannot use multiple fix wall commands with pair brownian/poly");
  else if (fixes.size() == 1) {
    wallfix = dynamic_cast<FixWall *>(fixes[0]);
    if (!wallfix)
      error->all(FLERR, "Fix {} is not compatible with pair brownian/poly", fixes[0]->style);
    flagwall = 1;
    if (wallfix->xflag) flagwall = 2; // Moving walls exist
  }

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume

  double vol_T, wallcoord;
  if (!flagwall) vol_T = domain->xprd*domain->yprd*domain->zprd;
  else {
    double wallhi[3], walllo[3];
    for (int j = 0; j < 3; j++) {
      wallhi[j] = domain->prd[j];
      walllo[j] = 0;
    }
    for (int m = 0; m < wallfix->nwall; m++) {
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      if (wallfix->xstyle[m] == FixWall::VARIABLE) {
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

  double volP = 0.0;
  for (int i = 0; i < nlocal; i++)
    volP += (4.0/3.0)*MY_PI*pow(atom->radius[i],3.0);
  MPI_Allreduce(&volP,&vol_P,1,MPI_DOUBLE,MPI_SUM,world);

  double vol_f = vol_P/vol_T;

  if (!flagVF) vol_f = 0;
  // set isotropic constants

  if (flaglog == 0) {
    R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
    RT0 = 8*MY_PI*mu;
  } else {
    R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBrownianPoly::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];
  return cut[i][j];
}
