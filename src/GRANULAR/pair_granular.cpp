/* ----------------------------------------------------------------------
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing authors: Leo Silbert (SNL), Gary Grest (SNL),
		      Jeremy Lechman (SNL), Dan Bolintineanu (SNL), Ishan Srivastava (SNL)
----------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_granular.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define ONETHIRD 0.33333333333333333
#define TWOTHIRDS 0.66666666666666666
#define POW6ONE 0.550321208149104 //6^(-1/3)
#define POW6TWO 0.30285343213869  //6^(-2/3)

#define EPSILON 1e-10

enum {STIFFNESS, MATERIAL};
enum {VELOCITY, VISCOELASTIC, TSUJI};
enum {HOOKE, HERTZ, DMT, JKR};
enum {TANGENTIAL_MINDLIN, TANGENTIAL_NOHISTORY};
enum {NONE, TWISTING_NOHISTORY, TWISTING_SDS, TWISTING_MARSHALL};
enum {NONE, ROLLING_NOHISTORY, ROLLING_SDS};

/* ---------------------------------------------------------------------- */

PairGranular::PairGranular(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  use_history = 1;
  fix_history = NULL;

  single_extra = 10;
  svector = new double[10];

  neighprev = 0;

  nmax = 0;
  mass_rigid = NULL;

  onerad_dynamic = NULL;
  onerad_frozen = NULL;
  maxrad_dynamic = NULL;
  maxrad_frozen = NULL;

  dt = update->dt;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  beyond_contact = 0;
  rolling_history_index = twisting_history_index = 0;
  tangential_history_index = -1;
}

/* ---------------------------------------------------------------------- */
PairGranular::~PairGranular()
{
  delete [] svector;
  if (fix_history) modify->delete_fix("NEIGH_HISTORY");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(E);
    memory->destroy(G);
    memory->destroy(normaldamp);
    memory->destroy(rollingdamp);
    memory->destroy(alpha);
    memory->destroy(gamman);
    memory->destroy(muS);
    memory->destroy(Ecoh);
    memory->destroy(kR);
    memory->destroy(muR);
    memory->destroy(etaR);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }
  memory->destroy(mass_rigid);
}

void PairGranular::compute(int eflag, int vflag){
  /*
#ifdef TEMPLATED_PAIR_GRANULAR
    if (normal == 0){
      if (damping == 0){
        if (tangential == 0){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,0,0,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,0,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,0,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,0,0,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,0,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,0,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,0,0,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,0,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,0,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 1){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,0,1,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,1,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,1,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,0,1,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,1,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,1,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,0,1,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,1,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,1,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 2){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,0,2,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,2,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,2,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,0,2,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,2,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,2,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,0,2,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,0,2,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,0,2,2,2>(eflag, vflag);
          }
        }
      }
      else if (damping == 1){
        if (tangential == 0){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,1,0,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,0,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,0,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,1,0,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,0,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,0,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,1,0,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,0,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,0,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 1){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,1,1,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,1,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,1,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,1,1,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,1,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,1,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,1,1,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,1,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,1,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 2){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,1,2,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,2,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,2,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,1,2,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,2,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,2,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,1,2,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,1,2,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,1,2,2,2>(eflag, vflag);
          }
        }
      }
      else if (damping == 2){
        if (tangential == 0){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,2,0,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,0,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,0,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,2,0,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,0,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,0,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,2,0,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,0,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,0,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 1){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,2,1,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,1,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,1,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,2,1,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,1,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,1,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,2,1,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,1,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,1,2,2>(eflag, vflag);
          }
        }
        else if (tangential == 2){
          if (rolling == 0){
            if (twisting == 0)      compute_templated<0,0,2,2,0,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,2,0,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,2,0,2>(eflag, vflag);
          }
          else if (rolling == 1){
            if (twisting == 0)      compute_templated<0,0,2,2,1,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,2,1,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,2,1,2>(eflag, vflag);
          }
          else if (rolling == 2){
            if (twisting == 0)      compute_templated<0,0,2,2,2,0>(eflag, vflag);
            else if (twisting == 1) compute_templated<0,0,2,2,2,1>(eflag, vflag);
            else if (twisting == 2) compute_templated<0,0,2,2,2,2>(eflag, vflag);
          }
        }
      }
    }


  }
#else
#endif
*/
  compute_untemplated(Tp_coeff_types, Tp_normal, Tp_damping, Tp_tangential,
      Tp_rolling, Tp_twisting, eflag, vflag);
}
/* ---------------------------------------------------------------------- */
/*#ifdef TEMPLATED_PAIR_GRANULAR
template < int Tp_coeff_types,
           int Tp_normal, int Tp_damping, int Tp_tangential,
           int Tp_rolling, int Tp_twisting >
void PairGranular::compute_templated(int eflag, int vflag)
#else
*/
void PairGranular::compute_untemplated
  (int Tp_coeff_types,
   int Tp_normal, int Tp_damping, int Tp_tangential,
   int Tp_rolling, int Tp_twisting,
   int eflag, int vflag)
//#endif
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,nx,ny,nz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,R,a;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double kn, kt, k_Q, k_R, eta_N, eta_T, eta_Q, eta_R;
  double Fne, Fdamp, Fntot, Fscrit, Frcrit;

  //For JKR
  double R, R2, coh, delta_pulloff, dist_pulloff, a, E;
  double overlap, olapsq, olapcubed, sqrtterm, tmp, a0;
  double keyterm, keyterm2, keyterm3, aovera0, foverFc;

  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double fs,fs1,fs2,fs3,roll1,roll2,roll3,torroll1,torroll2,torroll3;
  double tortwist1, tortwist2, tortwist3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool untouchflag;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0){
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
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
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  double coh;
  double **cohesion;
  double **stiffness;
  double **damping;
  if (Tp_normal == JKR){
    cohesion = coeffs[normal_coeff_inds[4]];
  }


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++){
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;
      untouchflag = (rsq >= radsum*radsum);

      if (normal[itype][jtype] == JKR){
        R = radi*radj/(radi+radj);
        R2 = R*R;
        coh = cohesion[itype][jtype];
        a = cbrt(9.0*M_PI*coh*R2/(4*E));
        delta_pulloff = a*a/R - 2*sqrt(M_PI*coh*a/E);
        dist_pulloff = radsum+delta_pulloff;
        untouchflag = (rsq >= (dist_pulloff)*(dist_pulloff));
      }

      if (untouchflag){
        // unset non-touching neighbors
        touch[jj] = 0;
        history = &allhistory[size_history*jj];
        for (int k = 0; k < size_history; k++) history[k] = 0.0;
      }
      else{
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

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

        //****************************************
        //Normal force = JKR-adjusted Hertzian contact + damping
        //****************************************
        if (normal[itype][jtype] == JKR){

        }
        //Damping
        kn = 4.0/3.0*E[itype][jtype]*a;
        if (normaldamp[itype][jtype] == BRILLIANTOV) eta_N = a*meff*gamman[itype][jtype];
        else if (normaldamp[itype][jtype] == TSUJI) eta_N=alpha[itype][jtype]*sqrt(meff*kn);

        Fdamp = -eta_N*vnnr; //F_nd eq 23 and Zhao eq 19

        Fntot = Fne + Fdamp;
        //if (screen) fprintf(screen,"%d %d %16.16g %16.16g  \n",itype,jtype,Ecoh[itype][jtype],E[itype][jtype]);
        //if (logfile) fprintf(logfile,"%d %d %16.16g %16.16g \n",itype,jtype,Ecoh[itype][jtype],E[itype][jtype]);

        //****************************************
        //Tangential force, including history effects
        //****************************************

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

        // history effects
        touch[jj] = 1;
        history = &allhistory[size_history*jj];
        shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
            history[2]*history[2]);

        // Rotate and update displacements.
        // See e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
        if (historyupdate) {
          rsht = history[0]*nx + history[1]*ny + history[2]*nz;
          if (fabs(rsht) < EPSILON) rsht = 0;
          if (rsht > 0){
            scalefac = shrmag/(shrmag - rsht); //if rhst == shrmag, contacting pair has rotated 90 deg. in one step, in which case you deserve a crash!
            history[0] -= rsht*nx;
            history[1] -= rsht*ny;
            history[2] -= rsht*nz;
            //Also rescale to preserve magnitude
            history[0] *= scalefac;
            history[1] *= scalefac;
            history[2] *= scalefac;
          }
          //Update history
          history[0] += vtr1*dt;
          history[1] += vtr2*dt;
          history[2] += vtr3*dt;
        }

        // tangential forces = history + tangential velocity damping
        // following Zhao and Marshall Phys Fluids v20, p043302 (2008)
        kt=8.0*G[itype][jtype]*a;

        eta_T = eta_N; //Based on discussion in Marshall; eta_T can also be an independent parameter
        fs1 = -kt*history[0] - eta_T*vtr1; //eq 26
        fs2 = -kt*history[1] - eta_T*vtr2;
        fs3 = -kt*history[2] - eta_T*vtr3;

        // rescale frictional displacements and forces if needed
        Fscrit = muS[itype][jtype] * fabs(Fne + 2*F_C);
        // For JKR, use eq 43 of Marshall. For DMT, use Fne instead

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        if (fs > Fscrit) {
          if (shrmag != 0.0) {
            //history[0] = (Fcrit/fs) * (history[0] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
            //history[1] = (Fcrit/fs) * (history[1] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
            //history[2] = (Fcrit/fs) * (history[2] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
            history[0] = -1.0/kt*(Fscrit*fs1/fs + eta_T*vtr1); //Same as above, but simpler (check!)
            history[1] = -1.0/kt*(Fscrit*fs2/fs + eta_T*vtr2);
            history[2] = -1.0/kt*(Fscrit*fs3/fs + eta_T*vtr3);
            fs1 *= Fscrit/fs;
            fs2 *= Fscrit/fs;
            fs3 *= Fscrit/fs;
          } else fs1 = fs2 = fs3 = 0.0;
        }

        //****************************************
        // Rolling force, including history history effects
        //****************************************

        relrot1 = omega[i][0] - omega[j][0];
        relrot2 = omega[i][1] - omega[j][1];
        relrot3 = omega[i][2] - omega[j][2];

        // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
        // This is different from the Marshall papers, which use the Bagi/Kuhn formulation
        // for rolling velocity (see Wang et al for why the latter is wrong)
        vrl1 = R*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
        vrl2 = R*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
        vrl3 = R*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;
        vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
        if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
        else vrlmaginv = 0.0;

        // Rolling displacement
        rollmag = sqrt(history[3]*history[3] + history[4]*history[4] + history[5]*history[5]);
        rolldotn = history[3]*nx + history[4]*ny + history[5]*nz;

        if (historyupdate) {
          if (fabs(rolldotn) < EPSILON) rolldotn = 0;
          if (rolldotn > 0){ //Rotate into tangential plane
            scalefac = rollmag/(rollmag - rolldotn);
            history[3] -= rolldotn*nx;
            history[4] -= rolldotn*ny;
            history[5] -= rolldotn*nz;
            //Also rescale to preserve magnitude
            history[3] *= scalefac;
            history[4] *= scalefac;
            history[5] *= scalefac;
          }
          history[3] += vrl1*dt;
          history[4] += vrl2*dt;
          history[5] += vrl3*dt;
        }

        k_R = kR[itype][jtype]*4.0*F_C*pow(aovera0,1.5);
        if (rollingdamp[itype][jtype] == INDEP) eta_R = etaR[itype][jtype];
        else if (rollingdamp[itype][jtype] == BRILLROLL) eta_R = muR[itype][jtype]*fabs(Fne);
        fr1 = -k_R*history[3] - eta_R*vrl1;
        fr2 = -k_R*history[4] - eta_R*vrl2;
        fr3 = -k_R*history[5] - eta_R*vrl3;

        // rescale frictional displacements and forces if needed
        Frcrit = muR[itype][jtype] * fabs(Fne + 2*F_C);

        fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
        if (fr > Frcrit) {
          if (rollmag != 0.0) {
            history[3] = -1.0/k_R*(Frcrit*fr1/fr + eta_R*vrl1);
            history[4] = -1.0/k_R*(Frcrit*fr2/fr + eta_R*vrl2);
            history[5] = -1.0/k_R*(Frcrit*fr3/fr + eta_R*vrl3);
            fr1 *= Frcrit/fr;
            fr2 *= Frcrit/fr;
            fr3 *= Frcrit/fr;
          } else fr1 = fr2 = fr3 = 0.0;
        }


        //****************************************
        // Twisting torque, including history effects
        //****************************************
        magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
        history[6] += magtwist*dt;
        k_Q = 0.5*kt*a*a;; //eq 32
        eta_Q = 0.5*eta_T*a*a;
        magtortwist = -k_Q*history[6] - eta_Q*magtwist;//M_t torque (eq 30)

        signtwist = (magtwist > 0) - (magtwist < 0);
        Mtcrit=TWOTHIRDS*a*Fscrit;//critical torque (eq 44)
        if (fabs(magtortwist) > Mtcrit) {
          //history[6] = Mtcrit/k_Q*magtwist/fabs(magtwist);
          history[6] = 1.0/k_Q*(Mtcrit*signtwist - eta_Q*magtwist);
          magtortwist = -Mtcrit * signtwist; //eq 34
        }

        // Apply forces & torques

        fx = nx*Fntot + fs1;
        fy = ny*Fntot + fs2;
        fz = nz*Fntot + fs3;

        //if (screen) fprintf(screen,"%16.16g %16.16g %16.16g %16.16g %16.16g %16.16g %16.16g \n",fs1,fs2,fs3,Fntot,nx,ny,nz);
        //if (logfile) fprintf(logfile,"%16.16g %16.16g %16.16g %16.16g %16.16g %16.16g %16.16g \n",fs1,fs2,fs3,Fntot,nx,ny,nz);

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = ny*fs3 - nz*fs2;
        tor2 = nz*fs1 - nx*fs3;
        tor3 = nx*fs2 - ny*fs1;

        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        tortwist1 = magtortwist * nx;
        tortwist2 = magtortwist * ny;
        tortwist3 = magtortwist * nz;

        torque[i][0] += tortwist1;
        torque[i][1] += tortwist2;
        torque[i][2] += tortwist3;

        torroll1 = R*(ny*fr3 - nz*fr2); //n cross fr
        torroll2 = R*(nz*fr1 - nx*fr3);
        torroll3 = R*(nx*fr2 - ny*fr1);

        torque[i][0] += torroll1;
        torque[i][1] += torroll2;
        torque[i][2] += torroll3;

        if (force->newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

          torque[j][0] -= radj*tor1;
          torque[j][1] -= radj*tor2;
          torque[j][2] -= radj*tor3;

          torque[j][0] -= tortwist1;
          torque[j][1] -= tortwist2;
          torque[j][2] -= tortwist3;

          torque[j][0] -= torroll1;
          torque[j][1] -= torroll2;
          torque[j][2] -= torroll3;
        }
        if (evflag) ev_tally_xyz(i,j,nlocal,0,
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
  memory->create(normal_coeffs,n+1,n+1,4,"pair:normal_coeffs");
  memory->create(tangential_coeffs,n+1,n+1,3,"pair:tangential_coeffs");
  memory->create(rolling_coeffs,n+1,n+1,3,"pair:rolling_coeffs");
  memory->create(twisting_coeffs,n+1,n+1,3,"pair:twisting_coeffs");

  memory->create(normal,n+1,n+1,"pair:normal");
  memory->create(tangential,n+1,n+1,"pair:tangential");
  memory->create(rolling,n+1,n+1,"pair:rolling");
  memory->create(twisting,n+1,n+1,"pair:twisting");

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
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  int iarg = 0;

  //Some defaults
  damping_global = VISCOELASTIC;
  coeff_types = STIFFNESS;
  tangential_global = -1; //Needs to be explicitly set, since it requires parameters
  rolling_global = NONE;
  twisting_global = NONE;
  tangential_history = rolling_history = twisting_history = 0;

  if (strcmp(arg[iarg], "material") == 0){
    coeff_types = MATERIAL;
    iarg += 1;
  }
  while (iarg < narg){
    if (strcmp(arg[iarg], "hooke") == 0){
      if (coeff_types == MATERIAL) error->all(FLERR,"Illegal pair_coeff command, 'stiffness' coefficients required for Hooke");
      if (iarg + 2 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for Hooke option");
      normal_global = HOOKE;
      memory->create(normal_coeffs_global, 2, "pair:normal_coeffs_global");
      normal_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //kn
      normal_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      iarg += 3;
    }
    else if (strcmp(arg[iarg], "hertz") == 0){
      int num_coeffs = 2;
      if (coeff_types == MATERIAL) num_coeffs += 1;
      if (iarg + offset >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for Hertz option");
      normal_global = HERTZ;
      memory->create(normal_coeffs_global, num_coeffs, "pair:normal_coeffs_global");
      normal_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //kn or E
      normal_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      if (coeff_types == MATERIAL) normal_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //G (if 'material')
      iarg += num_coeffs+1;
    }
    else if (strcmp(arg[iarg], "dmt") == 0){
      if (coeff_types == STIFFNESS) error->all(FLERR,"Illegal pair_style command, 'material' coefficients required for DMT");
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for Hertz option");
      normal_global = DMT;
      memory->create(normal_coeffs_global, 4, "pair:normal_coeffs_global");
      normal_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //E
      normal_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_global[3] = force->numeric(FLERR,arg[iarg+3]); //cohesion
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "jkr") == 0){
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for JKR option");
      if (coeff_types == STIFFNESS) error->all(FLERR,"Illegal pair_style command, 'material' coefficients required for JKR");
      beyond_contact = 1;
      normal_global = JKR;
      memory->create(normal_coeffs_global, 4, "pair:normal_coeffs_global");
      normal_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //E
      normal_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_global[3] = force->numeric(FLERR,arg[iarg+4]); //cohesion
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "damp_velocity") == 0){
      damping_global = VELOCITY;
      iarg += 1;
    }
    else if (strcmp(arg[iarg], "damp_viscoelastic") == 0){
      damping_global = VISCOELASTIC;
      iarg += 1;
    }
    else if (strcmp(arg[iarg], "damp_tsuji") == 0){
      damping_global = TSUJI;
      iarg += 1;
    }
    else if (strstr(arg[iarg], "tangential") != NULL){
      if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for tangential model");
      if (strstr(arg[iarg], "nohistory") != NULL){
        tangential_global = TANGENTIAL_NOHISTORY;
      }
      else{
        tangential_global = TANGENTIAL_MINDLIN;
        tangential_history = 1;
      }
      memory->create(tangential_coeffs_global, 3, "pair:tangential_coeffs_global");
      tangential_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //kt
      tangential_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
      tangential_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
      iarg += 4;
    }
    else if (strstr(arg[iarg], "rolling") != NULL){
      if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for rolling model");
      if (strstr(arg[iarg], "nohistory") != NULL){
        rolling_global = ROLLING_NOHISTORY;
      }
      else{
        rolling_global = ROLLING_SDS;
        rolling_history = 1;
      }
      memory->create(rolling_coeffs_global, 3, "pair:rolling_coeffs_global");
      rolling_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //kt
      rolling_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
      rolling_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
      iarg += 4;
    }
    else if (strstr(arg[iarg], "twisting") != NULL){
      if (strstr(arg[iarg], "marshall") != NULL){
        twisting_global = TWISTING_MARSHALL;
        twisting_history = 1;
        memory->create(twisting_coeffs_global, 3, "pair:twisting_coeffs_global"); //To be filled later
      }
      else{
        if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_style command, not enough parameters provided for twisting model");
        if (strstr(arg[iarg], "nohistory") != NULL){
          twisting_global = TWISTING_NOHISTORY;
        }
        else{
          twisting_global = TWISTING_SDS;
          twisting_history = 1;
        }
        memory->create(twisting_coeffs_global, 3, "pair:twisting_coeffs_global");
        twisting_coeffs_global[0] = force->numeric(FLERR,arg[iarg+1]); //kt
        twisting_coeffs_global[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
        twisting_coeffs_global[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
        iarg += 4;
      }
    }
  }

  //Set all i-i entrie, which may be replaced by pair coeff commands
  //It may also make sense to consider removing all of the above, and only
  // having the option for pair_coeff to set the parameters, similar to most LAMMPS pair styles
  // The reason for the current setup is to keep true to existing pair gran/hooke etc. syntax,
  // where coeffs are set in the pair_style command, and a pair_coeff * * command is issued.
  allocate();
  double damp;
  for (int i = 1; i <= atom->ntypes; i++){
    normal[i][i] = normal_global;
    damping[i][i] = damping_global;
    tangential[i][i] = tangential_global;
    rolling[i][i] = rolling_global;
    twisting[i][i] = twisting_global;

    if (damping_global == TSUJI){
      double cor = normal_coeffs_global[1];
      damp = 1.2728-4.2783*cor+11.087*pow(cor,2)-22.348*pow(cor,3)+
          27.467*pow(cor,4)-18.022*pow(cor,5)+
          4.8218*pow(cor,6);
    }
    else damp = normal_coeffs_global[1];
    normal_coeffs[i][i][0] = normal_coeffs_global[0];
    normal_coeffs[i][i][1] = damp;
    if (coeff_types == MATERIAL) normal_coeffs[i][i][2] = normal_coeffs_global[2];
    if ((normal_global == JKR) || (normal_global == DMT))
      normal_coeffs[i][i][3] = normal_coeffs_global[3];

    tangential[i][i] = tangential_global;
    if (tangential_global != NONE)
      for (int k = 0; k < 3; k++)
        tangential_coeffs[i][i][k] = tangential_coeffs_global[k];

    rolling[i][i] = rolling_global;
    if (rolling_global != NONE)
      for (int k = 0; k < 3; k++)
        rolling_coeffs[i][i][k] = rolling_coeffs_global[k];

    twisting[i][i] = twisting_global;
    if (twisting_global != NONE)
      for (int k = 0; k < 3; k++)
        twisting_coeffs[i][i][k] = twisting_coeffs_local[k];

    setflag[i][i] = 1;
  }

  //Additional checks
  if (tangential_global == -1){
    error->all(FLERR, "Illegal pair_style command: must specify tangential model");
  }
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranular::coeff(int narg, char **arg)
{
  int normal_local, damping_local, tangential_local, rolling_local, twisting_local;
  double *normal_coeffs_local;
  double *tangential_coeffs_local;
  double *rolling_coeffs_local;
  double *twisting_coeffs_local;

  normal_coeffs_local = new double[4];
  tangential_coeffs_local = new double[4];
  rolling_coeffs_local = new double[4];
  twisting_coeffs_local = new double[4];

  if (narg < 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int iarg = 2;
  while (iarg < narg){
    if (strcmp(arg[iarg], "hooke") == 0){
      if (coeff_types == MATERIAL) error->all(FLERR,"Illegal pair_coeff command, 'stiffness' coefficients required for Hooke");
      if (iarg + 2 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hooke option");
      normal_local = HOOKE;
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kn
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      iarg += 3;
    }
    else if (strcmp(arg[iarg], "hertz") == 0){
      int num_coeffs = 2;
      if (coeff_types == MATERIAL) num_coeffs += 1;
      if (iarg + offset >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hertz option");
      normal_local = HERTZ;
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kn or E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      if (coeff_types == MATERIAL) normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G (if 'material')
      iarg += num_coeffs+1;
    }
    else if (strcmp(arg[iarg], "dmt") == 0){
      if (coeff_types == STIFFNESS) error->all(FLERR,"Illegal pair_coeff command, 'material' coefficients required for DMT");
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hertz option");
      normal_local = DMT;
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_local[3] = force->numeric(FLERR,arg[iarg+3]); //cohesion
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "jkr") == 0){
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for JKR option");
      if (coeff_types == STIFFNESS) error->all(FLERR,"Illegal pair_coeff command, 'material' coefficients required for JKR");
      beyond_contact = 1;
      normal_local = JKR;
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_local[3] = force->numeric(FLERR,arg[iarg+4]); //cohesion
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "damp_velocity") == 0){
      damping_local = VELOCITY;
      iarg += 1;
    }
    else if (strcmp(arg[iarg], "damp_viscoelastic") == 0){
      damping_local = VISCOELASTIC;
      iarg += 1;
    }
    else if (strcmp(arg[iarg], "damp_tsuji") == 0){
      damping_local = TSUJI;
      iarg += 1;
    }
    else if (strstr(arg[iarg], "tangential") != NULL){
      if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for tangential model");
      if (strstr(arg[iarg], "nohistory") != NULL){
        tangential_local = TANGENTIAL_NOHISTORY;
      }
      else{
        tangential_local = TANGENTIAL_MINDLIN;
        tangential_history = 1;
      }
      tangential_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kt
      tangential_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
      tangential_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
      iarg += 4;
    }
    else if (strstr(arg[iarg], "rolling") != NULL){
      if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for rolling model");
      if (strstr(arg[iarg], "nohistory") != NULL){
        rolling_local = ROLLING_NOHISTORY;
      }
      else{
        rolling_local = ROLLING_SDS;
        rolling_history = 1;
      }
      rolling_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kt
      rolling_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
      rolling_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
      iarg += 4;
    }
    else if (strstr(arg[iarg], "twisting") != NULL){
      if (strstr(arg[iarg], "marshall") != NULL){
        twisting_local = TWISTING_MARSHALL;
        twisting_history = 1;
      }
      else{
        if (iarg + 3 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for twisting model");
        if (strstr(arg[iarg], "nohistory") != NULL){
          twisting_local = TWISTING_NOHISTORY;
        }
        else{
          twisting_local = TWISTING_SDS;
          twisting_history = 1;
          size_history += 1;
        }
        twisting_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kt
        twisting_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //gammat
        twisting_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //friction coeff.
        iarg += 4;
      }
    }
  }

  int count = 0;
  double damp;
  if (damping_local == TSUJI){
    double cor = normal_coeffs_local[1];
    damp = 1.2728-4.2783*cor+11.087*pow(cor,2)-22.348*pow(cor,3)+
        27.467*pow(cor,4)-18.022*pow(cor,5)+
        4.8218*pow(cor,6);
  }
  else damp = normal_coeffs_local[1];

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      normal[i][j] = normal_local;
      normal_coeffs[i][j][0] = normal_coeffs_local[0];
      normal_coeffs[i][j][1] = damp;
      if (coeff_types == MATERIAL) normal_coeffs[i][j][2] = normal_coeffs_local[2];
      if ((normal_local == JKR) || (normal_local == DMT))
          normal_coeffs[i][j][3] = normal_coeffs_local[3];

      tangential[i][j] = tangential_local;
      if (tangential_local != NONE)
        for (int k = 0; k < 3; k++)
          tangential_coeffs[i][j][k] = tangential_coeffs_local[k];

      rolling[i][j] = rolling_local;
      if (rolling_local != NONE)
        for (int k = 0; k < 3; k++)
          rolling_coeffs[i][j][k] = rolling_coeffs_local[k];

      twisting[i][j] = twisting_local;
      if (twisting_local != NONE)
        for (int k = 0; k < 3; k++)
          twisting_coeffs[i][j][k] = twisting_coeffs_local[k];

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

  // need a granular neigh list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->size = 1;
  if (history) neighbor->requests[irequest]->history = 1;

  dt = update->dt;

  // if history is stored:
  // if first init, create Fix needed for storing history

  if (history && fix_history == NULL) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",size_history);
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "NEIGH_HISTORY";
    fixarg[3] = dnumstr;
    modify->add_fix(4,fixarg,1);
    delete [] fixarg;
    fix_history = (FixNeighHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
  }

  // check for FixFreeze and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = NULL;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];

  // check for FixPour and FixDeposit so can extract particle radii

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      onerad_dynamic[i] =
          *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      onerad_dynamic[i] =
          *((double *) modify->fix[idep]->extract("radius",itype));
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) 
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
      MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
      MPI_DOUBLE,MPI_MAX,world);

  // set fix which stores history info

  if (size_history > 0) {
    int ifix = modify->find_fix("NEIGH_HISTORY");
    if (ifix < 0) error->all(FLERR,"Could not find pair fix neigh history ID");
    fix_history = (FixNeighHistory *) modify->fix[ifix];
  }
}

/* ----------------------------------------------------------------------
	 init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranular::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    if ((normal[i] != normal[j]) ||
        (damping[i] != damping[j]) ||
        (tangential[i] != tangential[j]) ||
        (rolling[i] != rolling[j]) ||
        (twisting[i] != twisting[j])){

      char str[512];
      sprintf(str,"Granular pair style functional forms are different, cannot mix coefficients for types %d and %d. \nThis combination must be set explicitly via pair_coeff command.",i,j);
      error->one(FLERR,str);
    }

    if (coeff_types == MATERIAL){
      normal_coeffs[i][j][0] = mix_stiffnessE(normal_coeffs[i][i][0], normal_coeffs[j][j][0],
          normal_coeffs[i][i][2], normal_coeffs[j][j][2]);
      normal_coeffs[i][j][2] = mix_stiffnessG(normal_coeffs[i][i][0], normal_coeffs[j][j][0],
          normal_coeffs[i][i][2], normal_coeffs[j][j][2]);
    }
    else{
      normal_coeffs[i][j][0] = mix_geom(normal_coeffs[i][i][0], normal_coeffs[j][j][0];)
    }

    normal_coeffs[i][j][1] = mix_geom(normal_coeffs[i][i][1], normal_coeffs[j][j][1];)
    if ((normal[i][i] == JKR) || (normal[i][i] == DMT))
      normal_coeffs[i][j][3] = mix_geom(normal_coeffs[i][i][3], normal_coeffs[j][j][3]);

    if (tangential[i][i] != NONE){
      for (int k = 0; k < 3; k++)
        tangential_coeffs[i][j][k] = mix_geom(tangential_coeffs[i][i][k], tangential_coeffs[j][j][k]);
    }

    if (rolling[i][i] != NONE){
      for (int k = 0; k < 3; k++)
        rolling_coeffs[i][j][k] = mix_geom(rolling_coeffs[i][i][k], rolling_coeffs[j][j][k]);
    }

    if (twisting[i][i] != NONE){
      for (int k = 0; k < 3; k++)
        twisting_coeffs[i][j][k] = mix_geom(twisting_coeffs[i][i][k], twisting_coeffs[j][j][k]);
    }

  double cutoff = cut[i][j];

  // It is likely that cut[i][j] at this point is still 0.0. This can happen when 
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue,for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = min(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  if (cut[i][j] < 0.0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) || ((maxrad_dynamic[i] > 0.0) && (maxrad_frozen[j] > 0.0)) ||
        ((maxrad_frozen[i] > 0.0) && (maxrad_dynamic[j] > 0.0))) { // radius info about both i and j exist
      cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
      cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
      cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);
    }
    else { // radius info about either i or j does not exist (i.e. not present and not about to get poured; set to largest value to not interfere with neighbor list)
      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax,2.0*maxrad_dynamic[k]);
        cutmax = MAX(cutmax,2.0*maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  }
  return cutoff;
}


/* ----------------------------------------------------------------------
	 proc 0 writes to restart file
	 ------------------------------------------------------------------------- */

void PairGranular::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&E[i][j],sizeof(double),1,fp);
        fwrite(&G[i][j],sizeof(double),1,fp);
        fwrite(&normaldamp[i][j],sizeof(int),1,fp);
        fwrite(&rollingdamp[i][j],sizeof(int),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&gamman[i][j],sizeof(double),1,fp);
        fwrite(&muS[i][j],sizeof(double),1,fp);
        fwrite(&Ecoh[i][j],sizeof(double),1,fp);
        fwrite(&kR[i][j],sizeof(double),1,fp);
        fwrite(&muR[i][j],sizeof(double),1,fp);
        fwrite(&etaR[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
	 proc 0 reads from restart file, bcasts
	 ------------------------------------------------------------------------- */

void PairGranular::read_restart(FILE *fp)
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
          fread(&E[i][j],sizeof(double),1,fp);
          fread(&G[i][j],sizeof(double),1,fp);
          fread(&normaldamp[i][j],sizeof(int),1,fp);
          fread(&rollingdamp[i][j],sizeof(int),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&gamman[i][j],sizeof(double),1,fp);
          fread(&muS[i][j],sizeof(double),1,fp);
          fread(&Ecoh[i][j],sizeof(double),1,fp);
          fread(&kR[i][j],sizeof(double),1,fp);
          fread(&muR[i][j],sizeof(double),1,fp);
          fread(&etaR[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&E[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&G[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&normaldamp[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&rollingdamp[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamman[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&muS[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&Ecoh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kR[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&muR[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&etaR[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
	 proc 0 writes to restart file
	 ------------------------------------------------------------------------- */

void PairGranular::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
	 proc 0 reads from restart file, bcasts
	 ------------------------------------------------------------------------- */

void PairGranular::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce)
{
  //	feenableexcept(FE_INVALID | FE_OVERFLOW);
  double radi,radj,radsum;
  double r,rinv,rsqinv,delx,dely,delz, nx, ny, nz, R;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double overlap, a;
  double mi,mj,meff,damp,kn,kt;
  double Fdamp,Fne,Fntot,Fscrit;
  double eta_N,eta_T;
  double vtr1,vtr2,vtr3,vrel;
  double fs1,fs2,fs3,fs;
  double shrmag;
  double F_C, delta_C, olapsq, olapcubed, sqrtterm, tmp, a0;
  double keyterm, keyterm2, keyterm3, aovera0, foverFc;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;
  R = radi*radj/(radi+radj);
  a0 = pow(9.0*M_PI*Ecoh[itype][jtype]*R*R/E[itype][jtype],ONETHIRD);
  delta_C = 0.5*a0*a0*POW6ONE/R;

  int *touch = fix_history->firstflag[i];
  if ((rsq >= (radsum+delta_C)*(radsum+delta_C) )||
      (rsq >= radsum*radsum && touch[j])){
    fforce = 0.0;
    svector[0] = svector[1] = svector[2] = svector[3] = 0.0;
    return 0.0;
  }

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  nx = delx*rinv;
  ny = dely*rinv;
  nz = delz*rinv;


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
  int *type = atom->type;
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


  // normal force = JKR
  F_C = 3.0*R*M_PI*Ecoh[itype][jtype];
  overlap = radsum - r;
  olapsq = overlap*overlap;
  olapcubed = olapsq*olapsq;
  sqrtterm = sqrt(1.0 + olapcubed);
  tmp = 2.0 + olapcubed + 2.0*sqrtterm;
  keyterm = pow(tmp,ONETHIRD);
  keyterm2 = olapsq/keyterm;
  keyterm3 = sqrt(overlap + keyterm2 + keyterm);
  aovera0 = POW6TWO * (keyterm3 +
      sqrt(2.0*overlap - keyterm2 - keyterm + 4.0/keyterm3));// eq 41
  a = aovera0*a0;
  foverFc = 4.0*((aovera0*aovera0*aovera0) - pow(aovera0,1.5));//F_ne/F_C (eq 40)

  Fne = F_C*foverFc;

  //Damping
  kn = 4.0/3.0*E[itype][jtype]*a;
  if (normaldamp[itype][jtype] == BRILLIANTOV) eta_N = a*meff*gamman[itype][jtype];
  else if (normaldamp[itype][jtype] == TSUJI) eta_N=alpha[itype][jtype]*sqrt(meff*kn);

  Fdamp = -eta_N*vnnr; //F_nd eq 23 and Zhao eq 19

  Fntot = Fne + Fdamp;

  // relative velocities

  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  double *allhistory = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  double *history = &allhistory[3*neighprev];
  shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
      history[2]*history[2]);

  // tangential forces = history + tangential velocity damping
  kt=8.0*G[itype][jtype]*a;

  eta_T = eta_N; 
  fs1 = -kt*history[0] - eta_T*vtr1;
  fs2 = -kt*history[1] - eta_T*vtr2;
  fs3 = -kt*history[2] - eta_T*vtr3;

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  Fscrit= muS[itype][jtype] * fabs(Fne + 2*F_C);

  if (fs > Fscrit) {
    if (shrmag != 0.0) {
      fs1 *= Fscrit/fs;
      fs2 *= Fscrit/fs;
      fs3 *= Fscrit/fs;
      fs *= Fscrit/fs;
    } else fs1 = fs2 = fs3 = fs = 0.0;
  }

  // set all forces and return no energy

  fforce = Fntot;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranular::pack_forward_comm(int n, int *list, double *buf,
    int pbc_flag, int *pbc)
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
  double bytes = nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
 mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessE(double Eii, double Ejj, double Gii, double Gjj)
{
  double poisii = Eii/(2.0*Gii) - 1.0;
  double poisjj = Ejj/(2.0*Gjj) - 1.0;
  return 1/((1-poisii*poisjj)/Eii+(1-poisjj*poisjj)/Ejj);
}

/* ----------------------------------------------------------------------
	 mixing of shear modulus (G)
	 ------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessG(double Eii, double Ejj, double Gii, double Gjj)
{
  double poisii = Eii/(2.0*Gii) - 1.0;
  double poisjj = Ejj/(2.0*Gjj) - 1.0;
  return 1/((2.0 -poisjj)/Gii+(2.0-poisjj)/Gjj);
}

/* ----------------------------------------------------------------------
	 mixing of everything else 
	 ------------------------------------------------------------------------- */

double PairGranular::mix_geom(double valii, double valjj)
{
  return sqrt(valii*valjj); 
}
