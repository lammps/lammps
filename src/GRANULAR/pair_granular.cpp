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
   Contributing authors:
   Dan Bolintineanu (SNL), Ishan Srivastava (SNL), Jeremy Lechman(SNL)
   Leo Silbert (SNL), Gary Grest (SNL)
----------------------------------------------------------------------- */

#include "pair_granular.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS (4.0/3.0)                 // 4/3
#define THREEQUARTERS 0.75                 // 3/4

#define EPSILON 1e-10

enum {HOOKE, HERTZ, HERTZ_MATERIAL, DMT, JKR};
enum {VELOCITY, MASS_VELOCITY, VISCOELASTIC, TSUJI};
enum {TANGENTIAL_NOHISTORY, TANGENTIAL_HISTORY,
      TANGENTIAL_MINDLIN, TANGENTIAL_MINDLIN_RESCALE,
      TANGENTIAL_MINDLIN_FORCE, TANGENTIAL_MINDLIN_RESCALE_FORCE};
enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
enum {ROLL_NONE, ROLL_SDS};

/* ---------------------------------------------------------------------- */

PairGranular::PairGranular(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  finitecutflag = 1;

  single_extra = 12;
  svector = new double[single_extra];

  neighprev = 0;

  nmax = 0;
  mass_rigid = nullptr;

  onerad_dynamic = nullptr;
  onerad_frozen = nullptr;
  maxrad_dynamic = nullptr;
  maxrad_frozen = nullptr;

  limit_damping = nullptr;

  history_transfer_factors = nullptr;

  dt = update->dt;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  use_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;
  tangential_history_index = 0;
  roll_history_index = twist_history_index = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>(modify->add_fix("NEIGH_HISTORY_GRANULAR_DUMMY all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairGranular::~PairGranular()
{
  delete[] svector;
  delete[] history_transfer_factors;

  if (!fix_history) modify->delete_fix("NEIGH_HISTORY_GRANULAR_DUMMY");
  else modify->delete_fix("NEIGH_HISTORY_GRANULAR");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutoff_type);

    memory->destroy(normal_coeffs);
    memory->destroy(tangential_coeffs);
    memory->destroy(roll_coeffs);
    memory->destroy(twist_coeffs);

    memory->destroy(Emod);
    memory->destroy(poiss);

    memory->destroy(normal_model);
    memory->destroy(damping_model);
    memory->destroy(tangential_model);
    memory->destroy(roll_model);
    memory->destroy(twist_model);
    memory->destroy(limit_damping);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }

  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

void PairGranular::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,nx,ny,nz;
  double radi,radj,radsum,rsq,r,rinv,factor_lj;
  double Reff, delta, dR, dR2, dist_to_contact;

  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;

  double knfac, damp_normal=0.0, damp_normal_prefactor;
  double k_tangential, damp_tangential;
  double Fne, Ft, Fdamp, Fntot, Fncrit, Fscrit, Frcrit;
  double fs, fs1, fs2, fs3, tor1, tor2, tor3;

  double mi,mj,meff;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3;

  // for JKR
  double R2, coh, F_pulloff, delta_pulloff, dist_pulloff, a, a2, E;
  double t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  // rolling
  double k_roll, damp_roll;
  double torroll1, torroll2, torroll3;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;

  // twisting
  double k_twist, damp_twist, mu_twist;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double tortwist1, tortwist2, tortwist3;

  double shrmag,rsht,prjmag;
  bool frameupdate;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;
  const bool historyupdate = update->setupflag == 0;

  ev_init(eflag,vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      E = normal_coeffs[itype][jtype][0];
      Reff = radi*radj/radsum;
      touchflag = false;

      if (normal_model[itype][jtype] == JKR) {
        E *= THREEQUARTERS;
        if (touch[jj]) {
          R2 = Reff*Reff;
          coh = normal_coeffs[itype][jtype][3];
          a = cbrt(9.0*MY_PI*coh*R2/(4*E));
          delta_pulloff = a*a/Reff - 2*sqrt(MY_PI*coh*a/E);
          dist_pulloff = radsum-delta_pulloff;
          touchflag = (rsq < dist_pulloff*dist_pulloff);
        } else {
          touchflag = (rsq < radsum*radsum);
        }
      } else {
        touchflag = (rsq < radsum*radsum);
      }

      if (!touchflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history*jj];
          for (int k = 0; k < size_history; k++) history[k] = 0.0;
        }
      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;

        nx = delx*rinv;
        ny = dely*rinv;
        nz = delz*rinv;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*nx + vr2*ny + vr3*nz; //v_R . n
        vn1 = nx*vnnr;
        vn2 = ny*vnnr;
        vn3 = nz*vnnr;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi*mj / (mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        delta = radsum - r;
        dR = delta*Reff;

        if (normal_model[itype][jtype] == JKR) {
          touch[jj] = 1;
          R2=Reff*Reff;
          coh = normal_coeffs[itype][jtype][3];
          dR2 = dR*dR;
          t0 = coh*coh*R2*R2*E;
          t1 = PI27SQ*t0;
          t2 = 8*dR*dR2*E*E*E;
          t3 = 4*dR2*E;
          // in case sqrt(0) < 0 due to precision issues
          sqrt1 = MAX(0, t0*(t1+2*t2));
          t4 = cbrt(t1+t2+THREEROOT3*MY_PI*sqrt(sqrt1));
          t5 = t3/t4 + t4/E;
          sqrt2 = MAX(0, 2*dR + t5);
          t6 = sqrt(sqrt2);
          sqrt3 = MAX(0, 4*dR - t5 + SIXROOT6*coh*MY_PI*R2/(E*t6));
          a = INVROOT6*(t6 + sqrt(sqrt3));
          a2 = a*a;
          knfac = normal_coeffs[itype][jtype][0]*a;
          Fne = knfac*a2/Reff - MY_2PI*a2*sqrt(4*coh*E/(MY_PI*a));
        } else {
          knfac = E; // Hooke
          Fne = knfac*delta;
          a = sqrt(dR);
          if (normal_model[itype][jtype] != HOOKE) {
            Fne *= a;
            knfac *= a;
          }
          if (normal_model[itype][jtype] == DMT)
            Fne -= 4*MY_PI*normal_coeffs[itype][jtype][3]*Reff;
        }

        // NOTE: consider restricting Hooke to only have
        // 'velocity' as an option for damping?

        if (damping_model[itype][jtype] == VELOCITY) {
          damp_normal = 1;
        } else if (damping_model[itype][jtype] == MASS_VELOCITY) {
          damp_normal = meff;
        } else if (damping_model[itype][jtype] == VISCOELASTIC) {
          damp_normal = a*meff;
        } else if (damping_model[itype][jtype] == TSUJI) {
          damp_normal = sqrt(meff*knfac);
        } else damp_normal = 0.0;

        damp_normal_prefactor = normal_coeffs[itype][jtype][1]*damp_normal;
        Fdamp = -damp_normal_prefactor*vnnr;

        Fntot = Fne + Fdamp;
        if (limit_damping[itype][jtype] && (Fntot < 0.0)) Fntot = 0.0;

        //****************************************
        // tangential force, including history effects
        //****************************************

        // For linear, mindlin, mindlin_rescale:
        // history = cumulative tangential displacement
        //
        // For mindlin/force, mindlin_rescale/force:
        // history = cumulative tangential elastic force

        // tangential component
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity
        wr1 = (radi*omega[i][0] + radj*omega[j][0]);
        wr2 = (radi*omega[i][1] + radj*omega[j][1]);
        wr3 = (radi*omega[i][2] + radj*omega[j][2]);

        // relative tangential velocities
        vtr1 = vt1 - (nz*wr2-ny*wr3);
        vtr2 = vt2 - (nx*wr3-nz*wr1);
        vtr3 = vt3 - (ny*wr1-nx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // if any history is needed
        if (use_history) {
          touch[jj] = 1;
          history = &allhistory[size_history*jj];
        }

        if (normal_model[itype][jtype] == JKR) {
          F_pulloff = 3*MY_PI*coh*Reff;
          Fncrit = fabs(Fne + 2*F_pulloff);
        } else if (normal_model[itype][jtype] == DMT) {
          F_pulloff = 4*MY_PI*coh*Reff;
          Fncrit = fabs(Fne + 2*F_pulloff);
        } else {
          Fncrit = fabs(Fntot);
        }
        Fscrit = tangential_coeffs[itype][jtype][2] * Fncrit;

        //------------------------------
        // tangential forces
        //------------------------------
        k_tangential = tangential_coeffs[itype][jtype][0];
        damp_tangential = tangential_coeffs[itype][jtype][1] *
          damp_normal_prefactor;

        if (tangential_history) {
          if (tangential_model[itype][jtype] == TANGENTIAL_MINDLIN ||
              tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_FORCE) {
            k_tangential *= a;
          } else if (tangential_model[itype][jtype] ==
                     TANGENTIAL_MINDLIN_RESCALE ||
                     tangential_model[itype][jtype] ==
                     TANGENTIAL_MINDLIN_RESCALE_FORCE) {
            k_tangential *= a;
            // on unloading, rescale the shear displacements/force
            if (a < history[3]) {
              double factor = a/history[3];
              history[0] *= factor;
              history[1] *= factor;
              history[2] *= factor;
            }
          }
          // rotate and update displacements / force.
          // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
          if (historyupdate) {
            rsht = history[0]*nx + history[1]*ny + history[2]*nz;
            if (tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_FORCE ||
                tangential_model[itype][jtype] ==
                TANGENTIAL_MINDLIN_RESCALE_FORCE)
              frameupdate = fabs(rsht) > EPSILON*Fscrit;
            else
              frameupdate = fabs(rsht)*k_tangential > EPSILON*Fscrit;
            if (frameupdate) {
              shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
                                               history[2]*history[2]);
              // projection
              history[0] -= rsht*nx;
              history[1] -= rsht*ny;
              history[2] -= rsht*nz;

              // also rescale to preserve magnitude
              prjmag = sqrt(history[0]*history[0] + history[1]*history[1] +
                                               history[2]*history[2]);
              if (prjmag > 0) scalefac = shrmag/prjmag;
              else scalefac = 0;
              history[0] *= scalefac;
              history[1] *= scalefac;
              history[2] *= scalefac;
            }
            // update history
            if (tangential_model[itype][jtype] == TANGENTIAL_HISTORY ||
                tangential_model[itype][jtype] == TANGENTIAL_MINDLIN ||
                tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_RESCALE) {
              // tangential displacement
              history[0] += vtr1*dt;
              history[1] += vtr2*dt;
              history[2] += vtr3*dt;
            } else {
              // tangential force
              // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
              history[0] -= k_tangential*vtr1*dt;
              history[1] -= k_tangential*vtr2*dt;
              history[2] -= k_tangential*vtr3*dt;
            }
            if (tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_RESCALE ||
                tangential_model[itype][jtype] ==
                TANGENTIAL_MINDLIN_RESCALE_FORCE)
              history[3] = a;
          }

          // tangential forces = history + tangential velocity damping
          if (tangential_model[itype][jtype] == TANGENTIAL_HISTORY ||
              tangential_model[itype][jtype] == TANGENTIAL_MINDLIN ||
              tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_RESCALE) {
            fs1 = -k_tangential*history[0] - damp_tangential*vtr1;
            fs2 = -k_tangential*history[1] - damp_tangential*vtr2;
            fs3 = -k_tangential*history[2] - damp_tangential*vtr3;
          } else {
            fs1 = history[0] - damp_tangential*vtr1;
            fs2 = history[1] - damp_tangential*vtr2;
            fs3 = history[2] - damp_tangential*vtr3;
          }

          // rescale frictional displacements and forces if needed
          fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
          if (fs > Fscrit) {
            shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
                                    history[2]*history[2]);
            if (shrmag != 0.0) {
              if (tangential_model[itype][jtype] == TANGENTIAL_HISTORY ||
                  tangential_model[itype][jtype] == TANGENTIAL_MINDLIN ||
                  tangential_model[itype][jtype] ==
                  TANGENTIAL_MINDLIN_RESCALE) {
                history[0] = -1.0/k_tangential*(Fscrit*fs1/fs +
                                                damp_tangential*vtr1);
                history[1] = -1.0/k_tangential*(Fscrit*fs2/fs +
                                                damp_tangential*vtr2);
                history[2] = -1.0/k_tangential*(Fscrit*fs3/fs +
                                                damp_tangential*vtr3);
              } else {
                history[0] = Fscrit*fs1/fs + damp_tangential*vtr1;
                history[1] = Fscrit*fs2/fs + damp_tangential*vtr2;
                history[2] = Fscrit*fs3/fs + damp_tangential*vtr3;
              }
              fs1 *= Fscrit/fs;
              fs2 *= Fscrit/fs;
              fs3 *= Fscrit/fs;
            } else fs1 = fs2 = fs3 = 0.0;
          }
        } else { // classic pair gran/hooke (no history)
          fs = damp_tangential*vrel;
          if (vrel != 0.0) Ft = MIN(Fscrit,fs) / vrel;
          else Ft = 0.0;
          fs1 = -Ft*vtr1;
          fs2 = -Ft*vtr2;
          fs3 = -Ft*vtr3;
        }

        if (roll_model[itype][jtype] != ROLL_NONE ||
            twist_model[itype][jtype] != TWIST_NONE) {
          relrot1 = omega[i][0] - omega[j][0];
          relrot2 = omega[i][1] - omega[j][1];
          relrot3 = omega[i][2] - omega[j][2];
        }
        //****************************************
        // rolling resistance
        //****************************************

        if (roll_model[itype][jtype] != ROLL_NONE) {
          // rolling velocity,
          // see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
          // this is different from the Marshall papers,
          // which use the Bagi/Kuhn formulation
          // for rolling velocity (see Wang et al for why the latter is wrong)
          vrl1 = Reff*(relrot2*nz - relrot3*ny);
          vrl2 = Reff*(relrot3*nx - relrot1*nz);
          vrl3 = Reff*(relrot1*ny - relrot2*nx);

          int rhist0 = roll_history_index;
          int rhist1 = rhist0 + 1;
          int rhist2 = rhist1 + 1;

          k_roll = roll_coeffs[itype][jtype][0];
          damp_roll = roll_coeffs[itype][jtype][1];
          Frcrit = roll_coeffs[itype][jtype][2] * Fncrit;

          if (historyupdate) {
            rolldotn = history[rhist0]*nx + history[rhist1]*ny + history[rhist2]*nz;
            frameupdate = fabs(rolldotn)*k_roll > EPSILON*Frcrit;
            if (frameupdate) { // rotate into tangential plane
              rollmag = sqrt(history[rhist0]*history[rhist0] +
                            history[rhist1]*history[rhist1] +
                            history[rhist2]*history[rhist2]);
              // projection
              history[rhist0] -= rolldotn*nx;
              history[rhist1] -= rolldotn*ny;
              history[rhist2] -= rolldotn*nz;
              // also rescale to preserve magnitude
              prjmag = sqrt(history[rhist0]*history[rhist0] +
                            history[rhist1]*history[rhist1] +
                            history[rhist2]*history[rhist2]);
              if (prjmag > 0) scalefac = rollmag/prjmag;
              else scalefac = 0;
              history[rhist0] *= scalefac;
              history[rhist1] *= scalefac;
              history[rhist2] *= scalefac;
            }
            history[rhist0] += vrl1*dt;
            history[rhist1] += vrl2*dt;
            history[rhist2] += vrl3*dt;
          }

          fr1 = -k_roll*history[rhist0] - damp_roll*vrl1;
          fr2 = -k_roll*history[rhist1] - damp_roll*vrl2;
          fr3 = -k_roll*history[rhist2] - damp_roll*vrl3;

          // rescale frictional displacements and forces if needed

          fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
          if (fr > Frcrit) {
            rollmag = sqrt(history[rhist0]*history[rhist0] +
                                    history[rhist1]*history[rhist1] +
                                    history[rhist2]*history[rhist2]);
            if (rollmag != 0.0) {
              history[rhist0] = -1.0/k_roll*(Frcrit*fr1/fr + damp_roll*vrl1);
              history[rhist1] = -1.0/k_roll*(Frcrit*fr2/fr + damp_roll*vrl2);
              history[rhist2] = -1.0/k_roll*(Frcrit*fr3/fr + damp_roll*vrl3);
              fr1 *= Frcrit/fr;
              fr2 *= Frcrit/fr;
              fr3 *= Frcrit/fr;
            } else fr1 = fr2 = fr3 = 0.0;
          }
        }

        //****************************************
        // twisting torque, including history effects
        //****************************************

        if (twist_model[itype][jtype] != TWIST_NONE) {
          // omega_T (eq 29 of Marshall)
          magtwist = relrot1*nx + relrot2*ny + relrot3*nz;
          if (twist_model[itype][jtype] == TWIST_MARSHALL) {
            k_twist = 0.5*k_tangential*a*a;; // eq 32 of Marshall paper
            damp_twist = 0.5*damp_tangential*a*a;
            mu_twist = TWOTHIRDS*a*tangential_coeffs[itype][jtype][2];
          } else {
            k_twist = twist_coeffs[itype][jtype][0];
            damp_twist = twist_coeffs[itype][jtype][1];
            mu_twist = twist_coeffs[itype][jtype][2];
          }
          if (historyupdate) {
            history[twist_history_index] += magtwist*dt;
          }
          magtortwist = -k_twist*history[twist_history_index] -
            damp_twist*magtwist; // M_t torque (eq 30)
          signtwist = (magtwist > 0) - (magtwist < 0);
          Mtcrit = mu_twist*Fncrit; // critical torque (eq 44)
          if (fabs(magtortwist) > Mtcrit) {
            history[twist_history_index] = 1.0/k_twist*(Mtcrit*signtwist -
                                                        damp_twist*magtwist);
            magtortwist = -Mtcrit * signtwist; // eq 34
          }
        }

        // apply forces & torques

        fx = nx*Fntot + fs1;
        fy = ny*Fntot + fs2;
        fz = nz*Fntot + fs3;
        fx *= factor_lj;
        fy *= factor_lj;
        fz *= factor_lj;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = ny*fs3 - nz*fs2;
        tor2 = nz*fs1 - nx*fs3;
        tor3 = nx*fs2 - ny*fs1;
        tor1 *= factor_lj;
        tor2 *= factor_lj;
        tor3 *= factor_lj;

        dist_to_contact = radi-0.5*delta;
        torque[i][0] -= dist_to_contact*tor1;
        torque[i][1] -= dist_to_contact*tor2;
        torque[i][2] -= dist_to_contact*tor3;

        if (twist_model[itype][jtype] != TWIST_NONE) {
          tortwist1 = magtortwist * nx * factor_lj;
          tortwist2 = magtortwist * ny * factor_lj;
          tortwist3 = magtortwist * nz * factor_lj;

          torque[i][0] += tortwist1;
          torque[i][1] += tortwist2;
          torque[i][2] += tortwist3;
        }

        if (roll_model[itype][jtype] != ROLL_NONE) {
          torroll1 = Reff*(ny*fr3 - nz*fr2) * factor_lj; // n cross fr
          torroll2 = Reff*(nz*fr1 - nx*fr3) * factor_lj;
          torroll3 = Reff*(nx*fr2 - ny*fr1) * factor_lj;

          torque[i][0] += torroll1;
          torque[i][1] += torroll2;
          torque[i][2] += torroll3;
        }

        if (force->newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

          dist_to_contact = radj-0.5*delta;
          torque[j][0] -= dist_to_contact*tor1;
          torque[j][1] -= dist_to_contact*tor2;
          torque[j][2] -= dist_to_contact*tor3;

          if (twist_model[itype][jtype] != TWIST_NONE) {
            torque[j][0] -= tortwist1;
            torque[j][1] -= tortwist2;
            torque[j][2] -= tortwist3;
          }
          if (roll_model[itype][jtype] != ROLL_NONE) {
            torque[j][0] -= torroll1;
            torque[j][1] -= torroll2;
            torque[j][2] -= torroll3;
          }
        }
        if (evflag) ev_tally_xyz(i,j,nlocal,force->newton_pair,
            0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranular::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutoff_type,n+1,n+1,"pair:cutoff_type");
  memory->create(normal_coeffs,n+1,n+1,4,"pair:normal_coeffs");
  memory->create(tangential_coeffs,n+1,n+1,3,"pair:tangential_coeffs");
  memory->create(roll_coeffs,n+1,n+1,3,"pair:roll_coeffs");
  memory->create(twist_coeffs,n+1,n+1,3,"pair:twist_coeffs");

  memory->create(Emod,n+1,n+1,"pair:Emod");
  memory->create(poiss,n+1,n+1,"pair:poiss");

  memory->create(normal_model,n+1,n+1,"pair:normal_model");
  memory->create(damping_model,n+1,n+1,"pair:damping_model");
  memory->create(tangential_model,n+1,n+1,"pair:tangential_model");
  memory->create(roll_model,n+1,n+1,"pair:roll_model");
  memory->create(twist_model,n+1,n+1,"pair:twist_model");
  memory->create(limit_damping,n+1,n+1,"pair:limit_damping");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranular::settings(int narg, char **arg)
{
  if (narg == 1) {
    cutoff_global = utils::numeric(FLERR,arg[0],false,lmp);
  } else {
    cutoff_global = -1; // will be set based on particle sizes, model choice
  }

  normal_history = tangential_history = 0;
  roll_history = twist_history = 0;
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranular::coeff(int narg, char **arg)
{
  int normal_model_one, damping_model_one;
  int tangential_model_one, roll_model_one, twist_model_one;

  double normal_coeffs_one[4];
  double tangential_coeffs_one[4];
  double roll_coeffs_one[4];
  double twist_coeffs_one[4];

  double cutoff_one = -1;

  if (narg < 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  //Defaults
  normal_model_one = tangential_model_one = -1;
  roll_model_one = ROLL_NONE;
  twist_model_one = TWIST_NONE;
  damping_model_one = VISCOELASTIC;
  int ld_flag = 0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "hooke") == 0) {
      if (iarg + 2 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hooke option");
      normal_model_one = HOOKE;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // kn
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      iarg += 3;
    } else if (strcmp(arg[iarg], "hertz") == 0) {
      if (iarg + 2 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz option");
      normal_model_one = HERTZ;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // kn
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      iarg += 3;
    } else if (strcmp(arg[iarg], "hertz/material") == 0) {
      if (iarg + 3 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz/material option");
      normal_model_one = HERTZ_MATERIAL;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      iarg += 4;
    } else if (strcmp(arg[iarg], "dmt") == 0) {
      if (iarg + 4 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz option");
      normal_model_one = DMT;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      normal_coeffs_one[3] = utils::numeric(FLERR,arg[iarg+4],false,lmp); // cohesion
      iarg += 5;
    } else if (strcmp(arg[iarg], "jkr") == 0) {
      if (iarg + 4 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for JKR option");
      beyond_contact = 1;
      normal_model_one = JKR;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      normal_coeffs_one[3] = utils::numeric(FLERR,arg[iarg+4],false,lmp); // cohesion
      iarg += 5;
    } else if (strcmp(arg[iarg], "damping") == 0) {
      if (iarg+1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, "
                   "not enough parameters provided for damping model");
      if (strcmp(arg[iarg+1], "velocity") == 0) {
        damping_model_one = VELOCITY;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "mass_velocity") == 0) {
        damping_model_one = MASS_VELOCITY;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "viscoelastic") == 0) {
        damping_model_one = VISCOELASTIC;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "tsuji") == 0) {
        damping_model_one = TSUJI;
        iarg += 1;
      } else error->all(FLERR, "Illegal pair_coeff command, "
                        "unrecognized damping model");
      iarg += 1;
    } else if (strcmp(arg[iarg], "tangential") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, must specify "
                   "tangential model after tangential keyword");
      if (strcmp(arg[iarg+1], "linear_nohistory") == 0) {
        if (iarg + 3 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for tangential model");
        tangential_model_one = TANGENTIAL_NOHISTORY;
        tangential_coeffs_one[0] = 0;
        // gammat and friction coeff
        tangential_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        tangential_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        iarg += 4;
      } else if ((strcmp(arg[iarg+1], "linear_history") == 0) ||
               (strcmp(arg[iarg+1], "mindlin") == 0) ||
               (strcmp(arg[iarg+1], "mindlin_rescale") == 0) ||
               (strcmp(arg[iarg+1], "mindlin/force") == 0) ||
               (strcmp(arg[iarg+1], "mindlin_rescale/force") == 0)) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for tangential model");
        if (strcmp(arg[iarg+1], "linear_history") == 0)
          tangential_model_one = TANGENTIAL_HISTORY;
        else if (strcmp(arg[iarg+1], "mindlin") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN;
        else if (strcmp(arg[iarg+1], "mindlin_rescale") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_RESCALE;
        else if (strcmp(arg[iarg+1], "mindlin/force") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_FORCE;
        else if (strcmp(arg[iarg+1], "mindlin_rescale/force") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_RESCALE_FORCE;
        tangential_history = 1;
        if ((tangential_model_one == TANGENTIAL_MINDLIN ||
             tangential_model_one == TANGENTIAL_MINDLIN_RESCALE ||
             tangential_model_one == TANGENTIAL_MINDLIN_FORCE ||
             tangential_model_one == TANGENTIAL_MINDLIN_RESCALE_FORCE) &&
            (strcmp(arg[iarg+2], "NULL") == 0)) {
          if (normal_model_one == HERTZ || normal_model_one == HOOKE) {
            error->all(FLERR, "NULL setting for Mindlin tangential "
                       "stiffness requires a normal contact model that "
                       "specifies material properties");
          }
          tangential_coeffs_one[0] = -1;
        } else {
          tangential_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // kt
        }
        // gammat and friction coeff
        tangential_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        tangential_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "tangential model not recognized");
      }
    } else if (strcmp(arg[iarg], "rolling") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0) {
        roll_model_one = ROLL_NONE;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "sds") == 0) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for rolling model");
        roll_model_one = ROLL_SDS;
        roll_history = 1;
        // kR and gammaR and rolling friction coeff
        roll_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        roll_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        roll_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "rolling friction model not recognized");
      }
    } else if (strcmp(arg[iarg], "twisting") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0) {
        twist_model_one = TWIST_NONE;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "marshall") == 0) {
        twist_model_one = TWIST_MARSHALL;
        twist_history = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "sds") == 0) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for twist model");
        twist_model_one = TWIST_SDS;
        twist_history = 1;
        // kt and gammat and friction coeff
        twist_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        twist_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        twist_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "twisting friction model not recognized");
      }
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      cutoff_one = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      ld_flag = 1;
      iarg += 1;
    } else error->all(FLERR, "Illegal pair_coeff command");
  }

  // error not to specify normal or tangential model
  if ((normal_model_one < 0) || (tangential_model_one < 0))
    error->all(FLERR, "Illegal pair_coeff command, "
               "must specify normal or tangential contact model");

  int count = 0;
  double damp;
  if (damping_model_one == TSUJI) {
    double cor;
    cor = normal_coeffs_one[1];
    damp = 1.2728-4.2783*cor+11.087*square(cor)-22.348*cube(cor)+
        27.467*powint(cor,4)-18.022*powint(cor,5)+4.8218*powint(cor,6);
  } else damp = normal_coeffs_one[1];

  if (ld_flag && normal_model_one == JKR)
    error->all(FLERR,"Illegal pair_coeff command, "
        "Cannot limit damping with JKR model");

  if (ld_flag && normal_model_one == DMT)
    error->all(FLERR,"Illegal pair_coeff command, "
        "Cannot limit damping with DMT model");

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      normal_model[i][j] = normal_model[j][i] = normal_model_one;
      normal_coeffs[i][j][1] = normal_coeffs[j][i][1] = damp;
      if (normal_model_one != HERTZ && normal_model_one != HOOKE) {
        Emod[i][j] = Emod[j][i] = normal_coeffs_one[0];
        poiss[i][j] = poiss[j][i] = normal_coeffs_one[2];
        normal_coeffs[i][j][0] = normal_coeffs[j][i][0] =
          FOURTHIRDS*mix_stiffnessE(Emod[i][j],Emod[i][j],
                                    poiss[i][j],poiss[i][j]);
      } else {
        normal_coeffs[i][j][0] = normal_coeffs[j][i][0] = normal_coeffs_one[0];
      }
      if ((normal_model_one == JKR) || (normal_model_one == DMT))
        normal_coeffs[i][j][3] = normal_coeffs[j][i][3] = normal_coeffs_one[3];

      damping_model[i][j] = damping_model[j][i] = damping_model_one;

      tangential_model[i][j] = tangential_model[j][i] = tangential_model_one;
      if (tangential_coeffs_one[0] == -1) {
        tangential_coeffs[i][j][0] = tangential_coeffs[j][i][0] =
          8*mix_stiffnessG(Emod[i][j],Emod[i][j],poiss[i][j],poiss[i][j]);
      } else {
        tangential_coeffs[i][j][0] = tangential_coeffs[j][i][0] =
          tangential_coeffs_one[0];
      }
      for (int k = 1; k < 3; k++)
        tangential_coeffs[i][j][k] = tangential_coeffs[j][i][k] =
          tangential_coeffs_one[k];

      roll_model[i][j] = roll_model[j][i] = roll_model_one;
      if (roll_model_one != ROLL_NONE)
        for (int k = 0; k < 3; k++)
          roll_coeffs[i][j][k] = roll_coeffs[j][i][k] = roll_coeffs_one[k];

      twist_model[i][j] = twist_model[j][i] = twist_model_one;
      if (twist_model_one != TWIST_NONE && twist_model_one != TWIST_MARSHALL)
        for (int k = 0; k < 3; k++)
          twist_coeffs[i][j][k] = twist_coeffs[j][i][k] = twist_coeffs_one[k];

      cutoff_type[i][j] = cutoff_type[j][i] = cutoff_one;

      limit_damping[i][j] = limit_damping[j][i] = ld_flag;

      setflag[i][j] = 1;
      count++;
    }
  }


  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
  init specific to this pair style
------------------------------------------------------------------------- */

void PairGranular::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, rmass");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  // determine whether we need a granular neigh list, how large it needs to be

  use_history = normal_history || tangential_history ||
    roll_history || twist_history;

  // for JKR, will need fix/neigh/history to keep track of touch arrays

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (normal_model[i][j] == JKR) use_history = 1;

  size_history = 3*tangential_history + 3*roll_history + twist_history;

  // determine location of tangential/roll/twist histories in array

  if (roll_history) {
    if (tangential_history) roll_history_index = 3;
    else roll_history_index = 0;
  }
  if (twist_history) {
    if (tangential_history) {
      if (roll_history) twist_history_index = 6;
      else twist_history_index = 3;
    } else {
      if (roll_history) twist_history_index = 3;
      else twist_history_index = 0;
    }
  }
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (tangential_model[i][j] == TANGENTIAL_MINDLIN_RESCALE ||
          tangential_model[i][j] == TANGENTIAL_MINDLIN_RESCALE_FORCE) {
        size_history += 1;
        roll_history_index += 1;
        twist_history_index += 1;
        nondefault_history_transfer = 1;
        history_transfer_factors = new int[size_history];
        for (int ii = 0; ii < size_history; ++ii)
          history_transfer_factors[ii] = -1;
        history_transfer_factors[3] = 1;
        break;
      }

  if (use_history) neighbor->add_request(this, NeighConst::REQ_SIZE|NeighConst::REQ_HISTORY);
  else neighbor->add_request(this, NeighConst::REQ_SIZE);

  dt = update->dt;

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix("NEIGH_HISTORY_GRANULAR_DUMMY",
                                                          "NEIGH_HISTORY_GRANULAR"
                                                          " all NEIGH_HISTORY "
                                                          + std::to_string(size_history),1));
    fix_history->pair = this;
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);

  // set fix which stores history info

  if (size_history > 0) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranular::init_one(int i, int j)
{
  double cutoff=0.0;

  if (setflag[i][j] == 0) {
    if ((normal_model[i][i] != normal_model[j][j]) ||
        (damping_model[i][i] != damping_model[j][j]) ||
        (tangential_model[i][i] != tangential_model[j][j]) ||
        (roll_model[i][i] != roll_model[j][j]) ||
        (twist_model[i][i] != twist_model[j][j])) {
      error->all(FLERR,"Granular pair style functional forms are different, "
                 "cannot mix coefficients for types {} and {}. \n"
                 "This combination must be set explicitly via a "
                 "pair_coeff command",i,j);
    }

    if (normal_model[i][j] == HERTZ || normal_model[i][j] == HOOKE)
      normal_coeffs[i][j][0] = normal_coeffs[j][i][0] =
        mix_geom(normal_coeffs[i][i][0], normal_coeffs[j][j][0]);
    else
      normal_coeffs[i][j][0] = normal_coeffs[j][i][0] =
        mix_stiffnessE(Emod[i][i], Emod[j][j], poiss[i][i], poiss[j][j]);

    normal_coeffs[i][j][1] = normal_coeffs[j][i][1] =
      mix_geom(normal_coeffs[i][i][1], normal_coeffs[j][j][1]);
    if ((normal_model[i][j] == JKR) || (normal_model[i][j] == DMT))
      normal_coeffs[i][j][3] = normal_coeffs[j][i][3] =
        mix_geom(normal_coeffs[i][i][3], normal_coeffs[j][j][3]);

    for (int k = 0; k < 3; k++)
      tangential_coeffs[i][j][k] = tangential_coeffs[j][i][k] =
        mix_geom(tangential_coeffs[i][i][k], tangential_coeffs[j][j][k]);

    if (roll_model[i][j] != ROLL_NONE) {
      for (int k = 0; k < 3; k++)
        roll_coeffs[i][j][k] = roll_coeffs[j][i][k] =
          mix_geom(roll_coeffs[i][i][k], roll_coeffs[j][j][k]);
    }

    if (twist_model[i][j] != TWIST_NONE && twist_model[i][j] != TWIST_MARSHALL) {
      for (int k = 0; k < 3; k++)
        twist_coeffs[i][j][k] = twist_coeffs[j][i][k] =
          mix_geom(twist_coeffs[i][i][k], twist_coeffs[j][j][k]);
    }
  }

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  double pulloff;

  if (cutoff_type[i][j] < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
      pulloff = 0.0;
      if (normal_model[i][j] == JKR) {
        pulloff = pulloff_distance(maxrad_dynamic[i], maxrad_dynamic[j], i, j);
        cutoff += pulloff;
      }

      if (normal_model[i][j] == JKR)
        pulloff = pulloff_distance(maxrad_frozen[i], maxrad_dynamic[j], i, j);
      cutoff = MAX(cutoff, maxrad_frozen[i]+maxrad_dynamic[j]+pulloff);

      if (normal_model[i][j] == JKR)
        pulloff = pulloff_distance(maxrad_dynamic[i], maxrad_frozen[j], i, j);
      cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]+pulloff);

    } else {

      // radius info about either i or j does not exist
      // (i.e. not present and not about to get poured;
      // set to largest value to not interfere with neighbor list)

      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax,2.0*maxrad_dynamic[k]);
        cutmax = MAX(cutmax,2.0*maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  } else if (cutoff_type[i][j] > 0) {
    cutoff = cutoff_type[i][j];
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

  return cutoff;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranular::write_restart(FILE *fp)
{
  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&normal_model[i][j],sizeof(int),1,fp);
        fwrite(&damping_model[i][j],sizeof(int),1,fp);
        fwrite(&tangential_model[i][j],sizeof(int),1,fp);
        fwrite(&roll_model[i][j],sizeof(int),1,fp);
        fwrite(&twist_model[i][j],sizeof(int),1,fp);
        fwrite(&limit_damping[i][j],sizeof(int),1,fp);
        fwrite(normal_coeffs[i][j],sizeof(double),4,fp);
        fwrite(tangential_coeffs[i][j],sizeof(double),3,fp);
        fwrite(roll_coeffs[i][j],sizeof(double),3,fp);
        fwrite(twist_coeffs[i][j],sizeof(double),3,fp);
        fwrite(&cutoff_type[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranular::read_restart(FILE *fp)
{
  allocate();
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&normal_model[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&damping_model[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&tangential_model[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&roll_model[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&twist_model[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&limit_damping[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,normal_coeffs[i][j],sizeof(double),4,fp,nullptr,error);
          utils::sfread(FLERR,tangential_coeffs[i][j],sizeof(double),3,fp,nullptr,error);
          utils::sfread(FLERR,roll_coeffs[i][j],sizeof(double),3,fp,nullptr,error);
          utils::sfread(FLERR,twist_coeffs[i][j],sizeof(double),3,fp,nullptr,error);
          utils::sfread(FLERR,&cutoff_type[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&normal_model[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&damping_model[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&tangential_model[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&roll_model[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&twist_model[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&limit_damping[i][j],1,MPI_INT,0,world);
        MPI_Bcast(normal_coeffs[i][j],4,MPI_DOUBLE,0,world);
        MPI_Bcast(tangential_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(roll_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(twist_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutoff_type[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
                            double rsq, double /* factor_coul */,
                            double /* factor_lj */, double &fforce)
{
  double radi,radj,radsum;
  double r,rinv,delx,dely,delz, nx, ny, nz, Reff;
  double dR, dR2;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3;

  double knfac, damp_normal, damp_normal_prefactor;
  double k_tangential, damp_tangential;
  double Fne, Ft, Fdamp, Fntot, Fncrit, Fscrit, Frcrit;
  double fs, fs1, fs2, fs3;

  // for JKR
  double R2, coh, F_pulloff, delta_pulloff, dist_pulloff, a, a2, E;
  double delta, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  // rolling
  double k_roll, damp_roll;
  double rollmag;
  double fr, fr1, fr2, fr3;

  // twisting
  double k_twist, damp_twist, mu_twist;
  double signtwist, magtwist, magtortwist, Mtcrit;

  double shrmag;
  int jnum;
  int *jlist;
  double *history,*allhistory;

  int nall = atom->nlocal + atom->nghost;
  if ((i >= nall) || (j >= nall))
    error->all(FLERR,"Not enough atoms for pair granular single function");

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;
  Reff = (radsum > 0.0) ? radi*radj/radsum : 0.0;

  bool touchflag;
  E = normal_coeffs[itype][jtype][0];
  if (normal_model[itype][jtype] == JKR) {
    E *= THREEQUARTERS;
    R2 = Reff*Reff;
    coh = normal_coeffs[itype][jtype][3];
    a = cbrt(9.0*MY_PI*coh*R2/(4*E));
    delta_pulloff = a*a/Reff - 2*sqrt(MY_PI*coh*a/E);
    dist_pulloff = radsum+delta_pulloff;
    touchflag = (rsq <= dist_pulloff*dist_pulloff);
  } else touchflag = (rsq <= radsum*radsum);

  if (!touchflag) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  r = sqrt(rsq);
  rinv = 1.0/r;

  nx = delx*rinv;
  ny = dely*rinv;
  nz = delz*rinv;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  vnnr = vr1*nx + vr2*ny + vr3*nz;
  vn1 = nx*vnnr;
  vn2 = ny*vnnr;
  vn3 = nz*vnnr;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi*omega[i][0] + radj*omega[j][0]);
  wr2 = (radi*omega[i][1] + radj*omega[j][1]);
  wr3 = (radi*omega[i][2] + radj*omega[j][2]);

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: ensure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  delta = radsum - r;
  dR = delta*Reff;
  if (normal_model[itype][jtype] == JKR) {
    dR2 = dR*dR;
    t0 = coh*coh*R2*R2*E;
    t1 = PI27SQ*t0;
    t2 = 8*dR*dR2*E*E*E;
    t3 = 4*dR2*E;
    // in case sqrt(0) < 0 due to precision issues
    sqrt1 = MAX(0, t0*(t1+2*t2));
    t4 = cbrt(t1+t2+THREEROOT3*MY_PI*sqrt(sqrt1));
    t5 = t3/t4 + t4/E;
    sqrt2 = MAX(0, 2*dR + t5);
    t6 = sqrt(sqrt2);
    sqrt3 = MAX(0, 4*dR - t5 + SIXROOT6*coh*MY_PI*R2/(E*t6));
    a = INVROOT6*(t6 + sqrt(sqrt3));
    a2 = a*a;
    knfac = normal_coeffs[itype][jtype][0]*a;
    Fne = knfac*a2/Reff - MY_2PI*a2*sqrt(4*coh*E/(MY_PI*a));
  } else {
    knfac = E;
    Fne = knfac*delta;
    a = sqrt(dR);
    if (normal_model[itype][jtype] != HOOKE) {
      Fne *= a;
      knfac *= a;
    }
    if (normal_model[itype][jtype] == DMT)
      Fne -= 4*MY_PI*normal_coeffs[itype][jtype][3]*Reff;
  }

  if (damping_model[itype][jtype] == VELOCITY) {
    damp_normal = 1;
  } else if (damping_model[itype][jtype] == MASS_VELOCITY) {
    damp_normal = meff;
  } else if (damping_model[itype][jtype] == VISCOELASTIC) {
    damp_normal = a*meff;
  } else if (damping_model[itype][jtype] == TSUJI) {
    damp_normal = sqrt(meff*knfac);
  } else damp_normal = 0.0;

  damp_normal_prefactor = normal_coeffs[itype][jtype][1]*damp_normal;
  Fdamp = -damp_normal_prefactor*vnnr;

  Fntot = Fne + Fdamp;
  if (limit_damping[itype][jtype] && (Fntot < 0.0)) Fntot = 0.0;

  jnum = list->numneigh[i];
  jlist = list->firstneigh[i];

  if (use_history) {
    if ((fix_history == nullptr) || (fix_history->firstvalue == nullptr))
      error->one(FLERR,"Pair granular single computation needs history");
    allhistory = fix_history->firstvalue[i];
    for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
    }
    history = &allhistory[size_history*neighprev];
  }

  //****************************************
  // tangential force, including history effects
  //****************************************

  // For linear, mindlin, mindlin_rescale:
  // history = cumulative tangential displacement
  //
  // For mindlin/force, mindlin_rescale/force:
  // history = cumulative tangential elastic force

  // tangential component
  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity
  wr1 = (radi*omega[i][0] + radj*omega[j][0]);
  wr2 = (radi*omega[i][1] + radj*omega[j][1]);
  wr3 = (radi*omega[i][2] + radj*omega[j][2]);

  // relative tangential velocities
  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  if (normal_model[itype][jtype] == JKR) {
    F_pulloff = 3*MY_PI*coh*Reff;
    Fncrit = fabs(Fne + 2*F_pulloff);
  } else if (normal_model[itype][jtype] == DMT) {
    F_pulloff = 4*MY_PI*coh*Reff;
    Fncrit = fabs(Fne + 2*F_pulloff);
  } else {
    Fncrit = fabs(Fntot);
  }
  Fscrit = tangential_coeffs[itype][jtype][2] * Fncrit;

  //------------------------------
  // tangential forces
  //------------------------------
  k_tangential = tangential_coeffs[itype][jtype][0];
  damp_tangential = tangential_coeffs[itype][jtype][1]*damp_normal_prefactor;

  if (tangential_history) {
    if (tangential_model[itype][jtype] != TANGENTIAL_HISTORY) {
      k_tangential *= a;
    }

    shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
        history[2]*history[2]);

    // tangential forces = history + tangential velocity damping
    if (tangential_model[itype][jtype] == TANGENTIAL_HISTORY ||
        tangential_model[itype][jtype] == TANGENTIAL_MINDLIN ||
        tangential_model[itype][jtype] == TANGENTIAL_MINDLIN_RESCALE) {
      fs1 = -k_tangential*history[0] - damp_tangential*vtr1;
      fs2 = -k_tangential*history[1] - damp_tangential*vtr2;
      fs3 = -k_tangential*history[2] - damp_tangential*vtr3;
    } else {
      fs1 = history[0] - damp_tangential*vtr1;
      fs2 = history[1] - damp_tangential*vtr2;
      fs3 = history[2] - damp_tangential*vtr3;
    }

    // rescale frictional forces if needed
    fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    if (fs > Fscrit) {
      if (shrmag != 0.0) {
        fs1 *= Fscrit/fs;
        fs2 *= Fscrit/fs;
        fs3 *= Fscrit/fs;
        fs *= Fscrit/fs;
      } else fs1 = fs2 = fs3 = fs = 0.0;
    }

  // classic pair gran/hooke (no history)
  } else {
    fs = damp_tangential*vrel;
    if (vrel != 0.0) Ft = MIN(Fscrit,fs) / vrel;
    else Ft = 0.0;
    fs1 = -Ft*vtr1;
    fs2 = -Ft*vtr2;
    fs3 = -Ft*vtr3;
    fs = Ft*vrel;
  }

  //****************************************
  // rolling resistance
  //****************************************

  if ((roll_model[itype][jtype] != ROLL_NONE)
      || (twist_model[itype][jtype] != TWIST_NONE)) {
    relrot1 = omega[i][0] - omega[j][0];
    relrot2 = omega[i][1] - omega[j][1];
    relrot3 = omega[i][2] - omega[j][2];

    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // this is different from the Marshall papers,
    // which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)

    vrl1 = Reff*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
    vrl2 = Reff*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
    vrl3 = Reff*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;

    int rhist0 = roll_history_index;
    int rhist1 = rhist0 + 1;
    int rhist2 = rhist1 + 1;

    // rolling displacement
    rollmag = sqrt(history[rhist0]*history[rhist0] +
        history[rhist1]*history[rhist1] +
        history[rhist2]*history[rhist2]);

    k_roll = roll_coeffs[itype][jtype][0];
    damp_roll = roll_coeffs[itype][jtype][1];
    fr1 = -k_roll*history[rhist0] - damp_roll*vrl1;
    fr2 = -k_roll*history[rhist1] - damp_roll*vrl2;
    fr3 = -k_roll*history[rhist2] - damp_roll*vrl3;

    // rescale frictional displacements and forces if needed
    Frcrit = roll_coeffs[itype][jtype][2] * Fncrit;

    fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
    if (fr > Frcrit) {
      if (rollmag != 0.0) {
        fr1 *= Frcrit/fr;
        fr2 *= Frcrit/fr;
        fr3 *= Frcrit/fr;
                fr *= Frcrit/fr;
      } else fr1 = fr2 = fr3 = fr = 0.0;
    }
  } else fr1 = fr2 = fr3 = fr = 0.0;

  //****************************************
  // twisting torque, including history effects
  //****************************************

  if (twist_model[itype][jtype] != TWIST_NONE) {
    // omega_T (eq 29 of Marshall)
    magtwist = relrot1*nx + relrot2*ny + relrot3*nz;
    if (twist_model[itype][jtype] == TWIST_MARSHALL) {
      k_twist = 0.5*k_tangential*a*a;; //eq 32
      damp_twist = 0.5*damp_tangential*a*a;
      mu_twist = TWOTHIRDS*a*tangential_coeffs[itype][jtype][2];;
    } else {
      k_twist = twist_coeffs[itype][jtype][0];
      damp_twist = twist_coeffs[itype][jtype][1];
      mu_twist = twist_coeffs[itype][jtype][2];
    }
    // M_t torque (eq 30)
    magtortwist = -k_twist*history[twist_history_index] - damp_twist*magtwist;
    signtwist = (magtwist > 0) - (magtwist < 0);
    Mtcrit = mu_twist*Fncrit; // critical torque (eq 44)
    if (fabs(magtortwist) > Mtcrit) {
      magtortwist = -Mtcrit * signtwist; // eq 34
    } else magtortwist = 0.0;
  } else magtortwist = 0.0;

  // set force and return no energy

  fforce = Fntot*rinv;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = fr1;
  svector[5] = fr2;
  svector[6] = fr3;
  svector[7] = fr;
  svector[8] = magtortwist;
  svector[9] = delx;
  svector[10] = dely;
  svector[11] = delz;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranular::pack_forward_comm(int n, int *list, double *buf,
                                    int /* pbc_flag */, int * /* pbc */)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranular::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGranular::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessE(double Eii, double Ejj,
                                    double poisii, double poisjj)
{
  return 1/((1-poisii*poisii)/Eii+(1-poisjj*poisjj)/Ejj);
}

/* ----------------------------------------------------------------------
   mixing of shear modulus (G)
------------------------------------------------------------------------ */

double PairGranular::mix_stiffnessG(double Eii, double Ejj,
                                    double poisii, double poisjj)
{
  return 1/((2*(2-poisii)*(1+poisii)/Eii) + (2*(2-poisjj)*(1+poisjj)/Ejj));
}

/* ----------------------------------------------------------------------
   mixing of everything else
------------------------------------------------------------------------- */

double PairGranular::mix_geom(double valii, double valjj)
{
  return sqrt(valii*valjj);
}


/* ----------------------------------------------------------------------
   compute pull-off distance (beyond contact) for a given radius and atom type
------------------------------------------------------------------------- */

double PairGranular::pulloff_distance(double radi, double radj,
                                      int itype, int jtype)
{
  double E, coh, a, Reff;
  Reff = radi*radj/(radi+radj);
  if (Reff <= 0) return 0;
  coh = normal_coeffs[itype][jtype][3];
  E = normal_coeffs[itype][jtype][0]*THREEQUARTERS;
  a = cbrt(9*MY_PI*coh*Reff*Reff/(4*E));
  return a*a/Reff - 2*sqrt(MY_PI*coh*a/E);
}

/* ----------------------------------------------------------------------
   transfer history during fix/neigh/history exchange
   only needed if any history entries i-j are not just negative of j-i entries
------------------------------------------------------------------------- */

void PairGranular::transfer_history(double* source, double* target)
{
  for (int i = 0; i < size_history; i++)
    target[i] = history_transfer_factors[i]*source[i];
}

/* ----------------------------------------------------------------------
   self-interaction range of particle
------------------------------------------------------------------------- */

double PairGranular::atom2cut(int i)
{
  double cut;

  cut = atom->radius[i]*2;
  if(beyond_contact) {
    int itype = atom->type[i];
    if(normal_model[itype][itype] == JKR) {
      cut += pulloff_distance(cut, cut, itype, itype);
    }
  }

  return cut;
}

/* ----------------------------------------------------------------------
   maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairGranular::radii2cut(double r1, double r2)
{
  double cut = 0.0;

  if(beyond_contact) {
    int n = atom->ntypes;
    double temp;

    // Check all combinations of i and j to find theoretical maximum pull off distance
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        if(normal_model[i][j] == JKR) {
          temp = pulloff_distance(r1, r2, i, j);
          if(temp > cut) cut = temp;
        }
      }
    }
  }

  cut += r1 + r2;

  return cut;
}
