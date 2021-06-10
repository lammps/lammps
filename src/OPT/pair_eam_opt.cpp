// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
     James Fischer, High Performance Technologies, Inc.
     Charles Cornwell, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natoli, Stone Ridge Technology
------------------------------------------------------------------------- */

#include "pair_eam_opt.h"
#include <cmath>

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMOpt::PairEAMOpt(LAMMPS *lmp) : PairEAM(lmp) {}

/* ---------------------------------------------------------------------- */

void PairEAMOpt::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_pair) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

/* ---------------------------------------------------------------------- */

template < int EVFLAG, int EFLAG, int NEWTON_PAIR >
void PairEAMOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;

  typedef struct {
    double rhor0i,rhor1i,rhor2i,rhor3i;
    double rhor0j,rhor1j,rhor2j,rhor3j;
  } fast_alpha_t;

  typedef struct {
    double rhor4i,rhor5i,rhor6i;
    double rhor4j,rhor5j,rhor6j;
    double z2r0,z2r1,z2r2,z2r3,z2r4,z2r5,z2r6;
    double _pad[3];
  } fast_gamma_t;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl = 0.0;
  double* _noalias coeff;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(numforce);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    memory->create(numforce,nmax,"pair:numforce");
  }

  double** _noalias x = atom->x;
  double** _noalias f = atom->f;
  int* _noalias type = atom->type;
  int nlocal = atom->nlocal;

  vec3_t* _noalias xx = (vec3_t*)x[0];
  vec3_t* _noalias ff = (vec3_t*)f[0];

  double tmp_cutforcesq = cutforcesq;
  double tmp_rdr = rdr;
  int nr2 = nr-2;
  int nr1 = nr-1;

  inum = list->inum;
  int* _noalias ilist = list->ilist;
  int** _noalias firstneigh = list->firstneigh;
  int* _noalias numneigh = list->numneigh;

  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;

  fast_alpha_t* _noalias fast_alpha =
    (fast_alpha_t*) malloc((size_t)ntypes2*(nr+1)*sizeof(fast_alpha_t));
  for (i = 0; i < ntypes; i++) for (j = 0; j < ntypes; j++) {
    fast_alpha_t* _noalias tab = &fast_alpha[i*ntypes*nr+j*nr];
    if (type2rhor[i+1][j+1] >= 0) {
      for (int m = 1; m <= nr; m++) {
        tab[m].rhor0i =  rhor_spline[type2rhor[i+1][j+1]][m][6];
        tab[m].rhor1i =  rhor_spline[type2rhor[i+1][j+1]][m][5];
        tab[m].rhor2i =  rhor_spline[type2rhor[i+1][j+1]][m][4];
        tab[m].rhor3i =  rhor_spline[type2rhor[i+1][j+1]][m][3];
      }
    }
    if (type2rhor[j+1][i+1] >= 0) {
      for (int m = 1; m <= nr; m++) {
        tab[m].rhor0j =  rhor_spline[type2rhor[j+1][i+1]][m][6];
        tab[m].rhor1j =  rhor_spline[type2rhor[j+1][i+1]][m][5];
        tab[m].rhor2j =  rhor_spline[type2rhor[j+1][i+1]][m][4];
        tab[m].rhor3j =  rhor_spline[type2rhor[j+1][i+1]][m][3];
      }
    }
  }
  fast_alpha_t* _noalias tabeight = fast_alpha;

  fast_gamma_t* _noalias fast_gamma =
    (fast_gamma_t*) malloc((size_t)ntypes2*(nr+1)*sizeof(fast_gamma_t));
  for (i = 0; i < ntypes; i++) for (j = 0; j < ntypes; j++) {
    fast_gamma_t* _noalias tab = &fast_gamma[i*ntypes*nr+j*nr];
    if (type2rhor[i+1][j+1] >= 0) {
      for (int m = 1; m <= nr; m++) {
        tab[m].rhor4i =  rhor_spline[type2rhor[i+1][j+1]][m][2];
        tab[m].rhor5i =  rhor_spline[type2rhor[i+1][j+1]][m][1];
        tab[m].rhor6i =  rhor_spline[type2rhor[i+1][j+1]][m][0];
      }
    }
    if (type2rhor[j+1][i+1] >= 0) {
      for (int m = 1; m <= nr; m++) {
        tab[m].rhor4j =  rhor_spline[type2rhor[j+1][i+1]][m][2];
        tab[m].rhor5j =  rhor_spline[type2rhor[j+1][i+1]][m][1];
        tab[m].rhor6j =  rhor_spline[type2rhor[j+1][i+1]][m][0];
        tab[m].z2r6 =  z2r_spline[type2z2r[i+1][j+1]][m][0];
      }
    }
    if (type2z2r[i+1][j+1] >= 0) {
      for (int m = 1; m <= nr; m++) {
        tab[m].z2r0 =  z2r_spline[type2z2r[i+1][j+1]][m][6];
        tab[m].z2r1 =  z2r_spline[type2z2r[i+1][j+1]][m][5];
        tab[m].z2r2 =  z2r_spline[type2z2r[i+1][j+1]][m][4];
        tab[m].z2r3 =  z2r_spline[type2z2r[i+1][j+1]][m][3];
        tab[m].z2r4 =  z2r_spline[type2z2r[i+1][j+1]][m][2];
        tab[m].z2r5 =  z2r_spline[type2z2r[i+1][j+1]][m][1];
        tab[m].z2r6 =  z2r_spline[type2z2r[i+1][j+1]][m][0];
      }
    }
  }
  fast_gamma_t* _noalias tabss = fast_gamma;

  // zero out density

  if (NEWTON_PAIR) {
    int m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    itype = type[i] - 1;
    int* _noalias jlist = firstneigh[i];
    jnum = numneigh[i];

    double tmprho = rho[i];
    fast_alpha_t* _noalias tabeighti = &tabeight[itype*ntypes*nr];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      double delx = xtmp - xx[j].x;
      double dely = ytmp - xx[j].y;
      double delz = ztmp - xx[j].z;
      double rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < tmp_cutforcesq) {
        jtype = type[j] - 1;

        double p = sqrt(rsq)*tmp_rdr;
        if ((int)p <= nr2) {
          int m = (int)p + 1;
          p -= (double)((int)p);
          fast_alpha_t& a = tabeighti[jtype*nr+m];
          tmprho += ((a.rhor3j*p+a.rhor2j)*p+a.rhor1j)*p+a.rhor0j;
          if (NEWTON_PAIR || j < nlocal) {
            rho[j] += ((a.rhor3i*p+a.rhor2i)*p+a.rhor1i)*p+a.rhor0i;
          }
        } else {
          fast_alpha_t& a = tabeighti[jtype*nr+nr1];
          tmprho += a.rhor3j+a.rhor2j+a.rhor1j+a.rhor0j;
          if (NEWTON_PAIR || j < nlocal) {
            rho[j] += a.rhor3i+a.rhor2i+a.rhor1i+a.rhor0i;
          }
        }
      }
    }
    rho[i] = tmprho;
  }

  // communicate and sum densities

  if (NEWTON_PAIR) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double p = rho[i]*rdrho + 1.0;
    int m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (EFLAG) {
      double phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
      phi *= scale[type[i]][type[i]];
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);
  embedstep = update->ntimestep;

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    int itype1 = type[i] - 1;
    int* _noalias jlist = firstneigh[i];
    jnum = numneigh[i];

    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;

    fast_gamma_t* _noalias tabssi = &tabss[itype1*ntypes*nr];
    double* _noalias scale_i = scale[itype1+1]+1;
    numforce[i] = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      double delx = xtmp - xx[j].x;
      double dely = ytmp - xx[j].y;
      double delz = ztmp - xx[j].z;
      double rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < tmp_cutforcesq) {
        ++numforce[i];
        jtype = type[j] - 1;
        double r = sqrt(rsq);
        double rhoip,rhojp,z2,z2p;
        double p = r*tmp_rdr;
        if ((int)p <= nr2) {
          int m = (int) p + 1;
          m = MIN(m,nr-1);
          p -= (double)((int) p);
          p = MIN(p,1.0);

          fast_gamma_t& a = tabssi[jtype*nr+m];
          rhoip = (a.rhor6i*p + a.rhor5i)*p + a.rhor4i;
          rhojp = (a.rhor6j*p + a.rhor5j)*p + a.rhor4j;
          z2 = ((a.z2r3*p + a.z2r2)*p + a.z2r1)*p + a.z2r0;
          z2p = (a.z2r6*p + a.z2r5)*p + a.z2r4;

        } else {

          fast_gamma_t& a = tabssi[jtype*nr+nr1];
          rhoip = a.rhor6i + a.rhor5i + a.rhor4i;
          rhojp = a.rhor6j + a.rhor5j + a.rhor4j;
          z2 = a.z2r3 + a.z2r2 + a.z2r1 + a.z2r0;
          z2p = a.z2r6 + a.z2r5 + a.z2r4;
        }

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        double recip = 1.0/r;
        double phi = z2*recip;
        double phip = z2p*recip - phi*recip;
        double psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        double fpair = -scale_i[jtype]*psip*recip;

        tmpfx += delx*fpair;
        tmpfy += dely*fpair;
        tmpfz += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          ff[j].x -= delx*fpair;
          ff[j].y -= dely*fpair;
          ff[j].z -= delz*fpair;
        }

        if (EFLAG) evdwl = scale_i[jtype]*phi;

        if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }

    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }

  free(fast_alpha); fast_alpha = 0;
  free(fast_gamma); fast_gamma = 0;

  if (vflag_fdotr) virial_fdotr_compute();
}
