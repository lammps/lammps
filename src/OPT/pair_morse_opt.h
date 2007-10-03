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
   Contributing authors:  
     James Fischer, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natol, Stone Ridge Technology
------------------------------------------------------------------------- */

#ifndef PAIR_MORSE_OPT_H
#define PAIR_MORSE_OPT_H

#include "math.h"
#include "stdlib.h"
#include "pair_morse.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

namespace LAMMPS_NS {

class PairMorseOpt : public PairMorse {
 public:
  PairMorseOpt(class LAMMPS *);
  void compute(int, int);

 private:
  template < int EFLAG, int VFLAG, int NEWTON_PAIR > void eval();
};

template < int EFLAG, int VFLAG, int NEWTON_PAIR >
void PairMorseOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;
  
  typedef struct {
    double cutsq,r0,alpha,morse1,d0,offset;
    double _pad[2];
  } fast_alpha_t;
  
  int i,j,ii,jj,inum,jnum,itype,jtype;
  
  eng_vdwl = 0.0;
  if (VFLAG) for (i = 0; i < 6; i++) virial[i] = 0.0;
  
  double** __restrict__ x = atom->x;
  double** __restrict__ f = atom->f;
  int* __restrict__ type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double* __restrict__ special_lj = force->special_lj;
  
  inum = list->inum;
  int* __restrict__ ilist = list->ilist;
  int** __restrict__ firstneigh = list->firstneigh;
  int* __restrict__ numneigh = list->numneigh;
  
  vec3_t* __restrict__ xx = (vec3_t*)x[0];
  vec3_t* __restrict__ ff = (vec3_t*)f[0];
  
  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;
  
  fast_alpha_t* __restrict__ fast_alpha = 
    (fast_alpha_t*) malloc(ntypes2*sizeof(fast_alpha_t));
  for (i = 0; i < ntypes; i++) for (j = 0; j < ntypes; j++) {
    fast_alpha_t& a = fast_alpha[i*ntypes+j];
    a.cutsq = cutsq[i+1][j+1];
    a.r0 = r0[i+1][j+1];
    a.alpha = alpha[i+1][j+1];
    a.morse1 = morse1[i+1][j+1];
    a.d0 = d0[i+1][j+1];
    a.offset = offset[i+1][j+1];
  }
  fast_alpha_t* __restrict__ tabsix = fast_alpha;
  
  // loop over neighbors of my atoms
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    itype = type[i] - 1;
    int* __restrict__ jlist = firstneigh[i];
    jnum = numneigh[i];
    
    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;
    
    fast_alpha_t* __restrict__ tabsixi = (fast_alpha_t*)&tabsix[itype*ntypes];
   
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      double factor_lj;

      if (j < nall) {
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	double rsq = delx*delx + dely*dely + delz*delz;
	
	jtype = type[j] - 1;
	
	fast_alpha_t& a = tabsixi[jtype];
	if (rsq < a.cutsq) {
	  double r = sqrt(rsq);
	  double dr = r - a.r0;
	  double dexp = exp(-a.alpha * dr);
	  double fforce = a.morse1 * (dexp*dexp - dexp) / r;
	  if (EFLAG) {
	    double phi = a.d0 * (dexp*dexp - 2.0*dexp) - a.offset;
	    if (NEWTON_PAIR || j < nlocal) {
	      eng_vdwl += phi;
	    } else {
	      eng_vdwl += 0.5*phi;
	    }
	  }
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
	  }
	  
	  if (VFLAG == 1) {
	    if (NEWTON_PAIR == 0 && j >= nlocal) fforce *= 0.5;
	    virial[0] += delx*delx*fforce;
	    virial[1] += dely*dely*fforce;
	    virial[2] += delz*delz*fforce;
	    virial[3] += delx*dely*fforce;
	    virial[4] += delx*delz*fforce;
	    virial[5] += dely*delz*fforce;
	  }
	}
      } else {
	factor_lj = special_lj[j/nall];
	j = j % nall;
	
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	double rsq = delx*delx + dely*dely + delz*delz;
	
	jtype = type[j] - 1;
	
	fast_alpha_t& a = tabsixi[jtype];
	if (rsq < a.cutsq) {
	  double r = sqrt(rsq);
	  double dr = r - a.r0;
	  double dexp = exp(-a.alpha * dr);
	  double fforce = factor_lj * a.morse1 * (dexp*dexp - dexp) / r;
	  if (EFLAG) {
	    double phi = a.d0 * (dexp*dexp - 2.0*dexp) - a.offset;
	    if (NEWTON_PAIR || j < nlocal) {
	      eng_vdwl += factor_lj*phi;
	    } else {
	      eng_vdwl += 0.5*factor_lj*phi;
	    }
	  }
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
	  }
	  
	  if (VFLAG == 1) {
	    if (NEWTON_PAIR == 0 && j >= nlocal) fforce *= 0.5;
	    virial[0] += delx*delx*fforce;
	    virial[1] += dely*dely*fforce;
	    virial[2] += delz*delz*fforce;
	    virial[3] += delx*dely*fforce;
	    virial[4] += delx*delz*fforce;
	    virial[5] += dely*delz*fforce;
	  }
	}
      }
    }
    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }
  
  if (VFLAG == 2) virial_compute();
  
  free(fast_alpha); fast_alpha = 0;
  
}

}

#endif
