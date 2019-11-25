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
   Contributing author: Philipp Kloza (University of Cambridge)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "pair_mesocnt.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

#define MAXLINE 1024
#define SELF_CUTOFF 5
#define SMALL 1.0e-6
#define SWITCH 1.0e-6
#define DELTA1 1.0
#define DELTA2 2.0
#define QUADRATURE 100
#define UINF_POINTS 1001
#define GAMMA_POINTS 26
#define PHI_POINTS 1001
#define USEMI_POINTS 1001

/* ---------------------------------------------------------------------- */

PairMesoCNT::PairMesoCNT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  respa_enable = 0;
  one_coeff = 1;
  manybody_flag = 1;
  no_virial_fdotr_compute = 0;
  writedata = 0;
  ghostneigh = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMesoCNT::~PairMesoCNT()
{
  if (allocated) {
    memory->destroy(cutsq);
    memory->destroy(setflag);

    memory->destroy(uinf_coeff);
    memory->destroy(gamma_coeff);
    memory->destroy(phi_coeff);
    memory->destroy(usemi_coeff);

    memory->destroy(redlist);
    memory->destroy(chain);
    memory->destroy(nchain);
    memory->destroy(end);

    memory->destroy(p1);
    memory->destroy(p2);

    memory->destroy(param);

    memory->destroy(flocal);
    memory->destroy(basis);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{
  double *r1,*r2,*q1,*q2,*qe;

  double evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,evflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int *tag = atom->tag;
  int *mol = atom->molecule;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // iterate over all bonds
  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    i1 &= NEIGHMASK;
    i2 &= NEIGHMASK;

    r1 = x[i1];
    r2 = x[i2];

    // reduce neighbors to common list

    int numneigh1,numneigh2;
    int *neighlist1,*neighlist2;
    if (i1 > nlocal-1) numneigh1 = 0;
    else {
      neighlist1 = firstneigh[i1];
      numneigh1 = numneigh[i1];
    }
    if (i2 > nlocal-1) numneigh2 = 0;
    else {
      neighlist2 = firstneigh[i2];
      numneigh2 = numneigh[i2];
    }

    int numneigh_max = numneigh1 + numneigh2;
    if (numneigh_max < 2) continue;
    if (numneigh_max > redlist_size) {
      redlist_size = 2 * numneigh_max;
      memory->grow(redlist,redlist_size,"pair:redlist");
    }

    int numred = 0;
    for (int j = 0; j < numneigh1; j++) {
      int ind = neighlist1[j];
      if (mol[ind] == mol[i1] && abs(tag[ind] - tag[i1]) < SELF_CUTOFF)
	      continue;
      redlist[numred++] = ind;
    }
    int inflag = 0;
    for (int j2 = 0; j2 < numneigh2; j2++) {
      for (int j1 = 0; j1 < numneigh1; j1++) {
        if (neighlist1[j1] == neighlist2[j2]) {
          inflag = 1;
	  break;
	}
      }
      if (inflag) {
        inflag = 0;
	continue;
      }
      int ind = neighlist2[j2];
      if (mol[ind] == mol[i2] && abs(tag[ind] - tag[i2]) < SELF_CUTOFF)
	      continue;
      redlist[numred++] = ind;
    }

    if (numred < 2) continue;

    // sort list according to atom-id

    sort(redlist,numred);

    // split neighbor list into connected chains

    int cid = 0;
    int clen = 0;

    if (numred > chain_size) {
      chain_size = 2 * numred;
      memory->destroy(chain);
      memory->create(chain,chain_size,chain_size,"pair:chain");
      memory->grow(nchain,chain_size,"pair:nchain");
    }

    for (int j = 0; j < numred-1; j++) {
      int j1 = redlist[j];
      int j2 = redlist[j+1];
      chain[cid][clen++] = j1;
      if (tag[j2] - tag[j1] != 1 || mol[j1] != mol[j2]) {
        nchain[cid++] = clen;
	clen = 0;
      }
    }
    chain[cid][clen++] = redlist[numred-1];
    nchain[cid++] = clen;

    // check for chain ends

    if (cid > end_size) {
      end_size = 2 * cid;
      memory->grow(end,end_size,"pair:end");
    }

    for (int j = 0; j < cid; j++) {
      int clen = nchain[j];
      int cstart = chain[j][0];
      int cend = chain[j][clen-1];
      int tagstart = tag[cstart];
      int tagend = tag[cend];
      end[j] = 0;
      if (tagstart == 1) end[j] = 1;
      else {
        int idprev = atom->map(tagstart-1);
	if (idprev == -1 || mol[tagstart] != mol[idprev]) end[j] = 1;
      }
      if (tagend == atom->natoms) end[j] = 2;
      else {
        int idnext = atom->map(tagend+1);
	if (idnext == -1 || mol[tagend] != mol[idnext]) end[j] = 2;
      }
    }

    // iterate over all neighbouring chains

    for (int j = 0; j < cid; j++) {
      if (nchain[j] < 2) continue;
      
      zero3(p1);
      zero3(p2);
      
      double sumw = 0;
      int clen = nchain[j];

      // assign end position

      if (end[j] == 1) qe = x[chain[j][0]];
      else if (end[j] == 2) qe = x[chain[j][clen-1]];

      // compute substitute straight (semi-)infinite CNT

      for (int k = 0; k < clen-1; k++) {
        q1 = x[chain[j][k]];
	q2 = x[chain[j][k+1]];

	double w = weight(r1,r2,q1,q2);

	if (w == 0) {
          if (end[j] == 1 && k == 0) end[j] = 0;
	  else if (end[j] == 2 && k == clen-2) end[j] = 0;
	  continue;
	}
	sumw += w;

	scaleadd3(w,q1,p1,p1);
	scaleadd3(w,q2,p2,p2);
      }

      if (sumw == 0) continue;

      double sumwrec = 1.0 / sumw;
      scale3(sumwrec,p1);
      scale3(sumwrec,p2);

      // compute geometry and forces

      // infinite CNT case

      if (end[j] == 0) {
        geominf(r1,r2,p1,p2,param,basis);
	if (param[0] > cutoff) continue;
	finf(param,evdwl,flocal);
      }

      // semi-infinite CNT case with end at start of chain

      else if (end[j] == 1) {
        geomsemi(r1,r2,p1,p2,p1,param,basis);
	if (param[0] > cutoff) continue;
	fsemi(param,evdwl,flocal);
      }

      // semi-infinite CNT case with end at end of chain

      else {
        geomsemi(r1,r2,p1,p2,p2,param,basis);
	if (param[0] > cutoff) continue;
	fsemi(param,evdwl,flocal);
      }

      // convert forces to global coordinate system

      f[i1][0] += flocal[0][0]*basis[0][0]
	      + flocal[0][1]*basis[1][0]
	      + flocal[0][2]*basis[2][0];
      f[i1][1] += flocal[0][0]*basis[0][1]
	      + flocal[0][1]*basis[1][1]
	      + flocal[0][2]*basis[2][1];
      f[i1][2] += flocal[0][0]*basis[0][2]
	      + flocal[0][1]*basis[1][2]
	      + flocal[0][2]*basis[2][2];
      f[i2][0] += flocal[1][0]*basis[0][0]
	      + flocal[1][1]*basis[1][0]
	      + flocal[1][2]*basis[2][0];
      f[i2][1] += flocal[1][0]*basis[0][1]
	      + flocal[1][1]*basis[1][1]
	      + flocal[1][2]*basis[2][1];
      f[i2][2] += flocal[1][0]*basis[0][2]
	      + flocal[1][1]*basis[1][2]
	      + flocal[1][2]*basis[2][2];

      // compute energy

      if (eflag_either) {
	if (eflag_global) eng_vdwl += 0.5 * evdwl;
	if (eflag_atom) {
          eatom[i1] += 0.25 * evdwl;
	  eatom[i2] += 0.25 * evdwl;
	}
      }

      // compute
      if (vflag_atom) {
        vatom[i1][0] += f[i1][0]*x[i1][0];
        vatom[i1][1] += f[i1][1]*x[i1][1];
        vatom[i1][2] += f[i1][2]*x[i1][2];
        vatom[i1][3] += f[i1][1]*x[i1][0];
        vatom[i1][4] += f[i1][2]*x[i1][0];
        vatom[i1][5] += f[i1][2]*x[i1][1];
        vatom[i2][0] += f[i2][0]*x[i2][0];
        vatom[i2][1] += f[i2][1]*x[i2][1];
        vatom[i2][2] += f[i2][2]*x[i2][2];
        vatom[i2][3] += f[i2][1]*x[i2][0];
        vatom[i2][4] += f[i2][2]*x[i2][0];
        vatom[i2][5] += f[i2][2]*x[i2][1];
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  int ntypes = atom->ntypes;
  redlist_size = 64;
  chain_size = 64;
  end_size = 64;
  
  memory->create(cutsq,ntypes+1,ntypes+1,"pair:cutsq");
  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;

  memory->create(uinf_coeff,uinf_points,4,"pair:uinf_coeff");
  memory->create(gamma_coeff,gamma_points,4,"pair:gamma_coeff");
  memory->create(phi_coeff,phi_points,phi_points,4,4,"pair:phi_coeff");
  memory->create(usemi_coeff,usemi_points,usemi_points,4,4,"pair:usemi_coeff");

  memory->create(redlist,redlist_size,"pair:redlist");
  memory->create(chain,chain_size,chain_size,"pair:chain");
  memory->create(nchain,chain_size,"pair:nchain");
  memory->create(end,end_size,"pair:end");

  memory->create(p1,3,"pair:p1");
  memory->create(p2,3,"pair:p2");

  memory->create(param,7,"pair:param");

  memory->create(flocal,2,3,"pair:flocal");
  memory->create(basis,3,3,"pair:basis");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg == 0) {
    uinf_points = UINF_POINTS;
    gamma_points = GAMMA_POINTS;
    phi_points = PHI_POINTS;
    usemi_points = USEMI_POINTS;
  }
  else if (narg == 2) {
    uinf_points = force->inumeric(FLERR,arg[0]);
    gamma_points = force->inumeric(FLERR,arg[1]);
    phi_points = force->inumeric(FLERR,arg[0]);
    usemi_points = force->inumeric(FLERR,arg[0]);
  }
  else if (narg == 4) {
    uinf_points = force->inumeric(FLERR,arg[0]);
    gamma_points = force->inumeric(FLERR,arg[1]);
    phi_points = force->inumeric(FLERR,arg[2]);
    usemi_points = force->inumeric(FLERR,arg[3]);
  }
  else error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  // CNT constants
  n = force->inumeric(FLERR,arg[2]);
  sig = force->numeric(FLERR,arg[3]);
   
  // file names
  uinf_file = arg[4];
  gamma_file = arg[5];
  phi_file = arg[6];
  usemi_file = arg[7];

  // units
  ang = force->angstrom;
  angrec = 1.0 / ang;
  e = force->qelectron;
  erec = 1.0 / e;
  funit = e * angrec;

  // potential variables
  r = 1.421 * 3 * n / MY_2PI * ang;
  rsq = r * r;
  d = 2 * r;
  d_ang = d * angrec;
  rc = 3.0 * sig;
  cutoff = rc + d;
  cutoffsq = cutoff * cutoff;
  cutoff_ang = cutoff * angrec;
  cutoffsq_ang = cutoff_ang * cutoff_ang;
  comega = 0.275 * (1.0 - 1.0/(1.0 + 0.59*r*angrec));
  ctheta = 0.35 + 0.0226*(r*angrec - 6.785);

  // parse and bcast data
  int me;
  double *uinf_data,*gamma_data,**phi_data,**usemi_data;
  memory->create(uinf_data,uinf_points,"pair:uinf_data");
  memory->create(gamma_data,gamma_points,"pair:gamma_data");
  memory->create(phi_data,phi_points,phi_points,"pair:phi_data");
  memory->create(usemi_data,usemi_points,phi_points,"pair:usemi_data");

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    read_file(uinf_file,uinf_data,hstart_uinf,delh_uinf,uinf_points);
    read_file(gamma_file,gamma_data,hstart_gamma,delh_gamma,gamma_points);
    read_file(phi_file,phi_data,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_points);
    read_file(usemi_file,usemi_data,hstart_usemi,xistart_usemi,
		    delh_usemi,delxi_usemi,usemi_points);
  }

  MPI_Bcast(&hstart_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&hstart_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&hstart_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&psistart_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&hstart_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&xistart_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delpsi_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delxi_usemi,1,MPI_DOUBLE,0,world);

  MPI_Bcast(uinf_data,uinf_points,MPI_DOUBLE,0,world);
  MPI_Bcast(gamma_data,gamma_points,MPI_DOUBLE,0,world);
  for (int i = 0; i < phi_points; i++)
    MPI_Bcast(phi_data[i],phi_points,MPI_DOUBLE,0,world);
  for (int i = 0; i < usemi_points; i++)
    MPI_Bcast(usemi_data[i],usemi_points,MPI_DOUBLE,0,world);

  // compute spline coefficients
  spline_coeff(uinf_data,uinf_coeff,delh_uinf,uinf_points);
  spline_coeff(gamma_data,gamma_coeff,delh_gamma,gamma_points);
  spline_coeff(phi_data,phi_coeff,delh_phi,delpsi_phi,phi_points);
  spline_coeff(usemi_data,usemi_coeff,delh_usemi,delxi_usemi,usemi_points);

  memory->destroy(uinf_data);
  memory->destroy(gamma_data);
  memory->destroy(phi_data);
  memory->destroy(usemi_data);

  int ntypes = atom->ntypes; 
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 1;

  std::ofstream uinffile("uinf.dat");
  std::ofstream finffile("finf.dat");

  double deltah = 2.45e-10;
  double h = d + deltah;
  double alphastart = 0;
  double alphaend = MY_2PI;
  int points = 1001;
  double dalpha = (alphaend - alphastart) / (points - 1);
  double del = 1.0e-15;
  for (int i = 0; i < points; i++) {
    double alpha = alphastart + i*dalpha;
    double p1[3] = {-1.0e-9, 0, 0};
    double p2[3] = {1.0e-9, 0, 0};
    double r1[3] = {-1.0e-9*cos(alpha), -1.0e-9*sin(alpha), h};
    double r2[3] = {1.0e-9*cos(alpha), 1.0e-9*sin(alpha), h};
    
    double rx1_hi[3],ry1_hi[3],rz1_hi[3];
    double rx1_lo[3],ry1_lo[3],rz1_lo[3];
    double rx2_hi[3],ry2_hi[3],rz2_hi[3];
    double rx2_lo[3],ry2_lo[3],rz2_lo[3];
    copy3(r1,rx1_lo);
    copy3(r1,ry1_lo);
    copy3(r1,rz1_lo);
    copy3(r1,rx1_hi);
    copy3(r1,ry1_hi);
    copy3(r1,rz1_hi);
    copy3(r2,rx2_lo);
    copy3(r2,ry2_lo);
    copy3(r2,rz2_lo);
    copy3(r2,rx2_hi);
    copy3(r2,ry2_hi);
    copy3(r2,rz2_hi);
    rx1_lo[0] -= del;
    ry1_lo[1] -= del;
    rz1_lo[2] -= del;
    rx1_hi[0] += del;
    ry1_hi[1] += del;
    rz1_hi[2] += del;
    rx2_lo[0] -= del;
    ry2_lo[1] -= del;
    rz2_lo[2] -= del;
    rx2_hi[0] += del;
    ry2_hi[1] += del;
    rz2_hi[2] += del;
    
    double ux1_lo,uy1_lo,uz1_lo;
    double ux1_hi,uy1_hi,uz1_hi;
    double ux2_lo,uy2_lo,uz2_lo;
    double ux2_hi,uy2_hi,uz2_hi;

    geominf(rx1_lo,r2,p1,p2,param,basis);
    finf(param,ux1_lo,flocal);
    geominf(ry1_lo,r2,p1,p2,param,basis);
    finf(param,uy1_lo,flocal);
    geominf(rz1_lo,r2,p1,p2,param,basis);
    finf(param,uz1_lo,flocal);
    geominf(rx1_hi,r2,p1,p2,param,basis);
    finf(param,ux1_hi,flocal);
    geominf(ry1_hi,r2,p1,p2,param,basis);
    finf(param,uy1_hi,flocal);
    geominf(rz1_hi,r2,p1,p2,param,basis);
    finf(param,uz1_hi,flocal);

    geominf(r1,rx2_lo,p1,p2,param,basis);
    finf(param,ux2_lo,flocal);
    geominf(r1,ry2_lo,p1,p2,param,basis);
    finf(param,uy2_lo,flocal);
    geominf(r1,rz2_lo,p1,p2,param,basis);
    finf(param,uz2_lo,flocal);
    geominf(r1,rx2_hi,p1,p2,param,basis);
    finf(param,ux2_hi,flocal);
    geominf(r1,ry2_hi,p1,p2,param,basis);
    finf(param,uy2_hi,flocal);
    geominf(r1,rz2_hi,p1,p2,param,basis);
    finf(param,uz2_hi,flocal);

    double f1x = -0.5 * (ux1_hi - ux1_lo) / del;
    double f1y = -0.5 * (uy1_hi - uy1_lo) / del;
    double f1z = -0.5 * (uz1_hi - uz1_lo) / del;
    double f2x = -0.5 * (ux2_hi - ux2_lo) / del;
    double f2y = -0.5 * (uy2_hi - uy2_lo) / del;
    double f2z = -0.5 * (uz2_hi - uz2_lo) / del;
    geominf(r1,r2,p1,p2,param,basis);
    double evdwl;
    finf(param,evdwl,flocal);
    double f1x_an = flocal[0][0]*basis[0][0]
	      + flocal[0][1]*basis[1][0]
	      + flocal[0][2]*basis[2][0];
    double f1y_an = flocal[0][0]*basis[0][1]
	      + flocal[0][1]*basis[1][1]
	      + flocal[0][2]*basis[2][1];
    double f1z_an = flocal[0][0]*basis[0][2]
	      + flocal[0][1]*basis[1][2]
	      + flocal[0][2]*basis[2][2];
    double f2x_an = flocal[1][0]*basis[0][0]
	      + flocal[1][1]*basis[1][0]
	      + flocal[1][2]*basis[2][0];
    double f2y_an = flocal[1][0]*basis[0][1]
	      + flocal[1][1]*basis[1][1]
	      + flocal[1][2]*basis[2][1];
    double f2z_an = flocal[1][0]*basis[0][2]
	      + flocal[1][1]*basis[1][2]
	      + flocal[1][2]*basis[2][2];

    uinffile << alpha * 180.0 / MY_PI << " " << evdwl << std::endl;
    finffile << alpha * 180.0 / MY_PI
	    << " " << f1x << " " << f1y << " " << f1z
	    << " " << f2x << " " << f2y << " " << f2z 
	    << " " << f1x_an << " " << f1y_an << " " << f1z_an
            << " " << f2x_an << " " << f2y_an << " " << f2z_an
	    << std::endl;
  }

  uinffile.close();
  finffile.close();

  std::ofstream usemifile("usemi.dat");
  std::ofstream fsemifile("fsemi.dat");

  deltah = 3.15e-10;
  h = d + deltah;
  double xi1 = -1.0e-9;
  double xi2 = 1.0e-9;
  double alpha1 = 0.0 / 180.0 * MY_PI;
  double alpha2 = 10.0 / 180.0 * MY_PI;
  double alpha3 = 30.0 / 180.0 * MY_PI;
  double alpha4 = 90.0 / 180.0 * MY_PI;
  points = 1001;
  double etaestart = -2.0e-9;
  double etaeend = 2.0e-9;
  double detae = (etaeend - etaestart) / (points - 1);
  for (int i = 0; i < points; i++) {
    double etae = etaestart + i*detae;
    double p1[3] = {etae, 0, 0};
    double p2[3] = {5.0e-9, 0, 0};
    double r1[3] = {xi1*cos(alpha1), xi1*sin(alpha1), h};
    double r2[3] = {xi2*cos(alpha1), xi2*sin(alpha1), h};

    double rx1_hi[3],ry1_hi[3],rz1_hi[3];
    double rx1_lo[3],ry1_lo[3],rz1_lo[3];
    double rx2_hi[3],ry2_hi[3],rz2_hi[3];
    double rx2_lo[3],ry2_lo[3],rz2_lo[3];
    copy3(r1,rx1_lo);
    copy3(r1,ry1_lo);
    copy3(r1,rz1_lo);
    copy3(r1,rx1_hi);
    copy3(r1,ry1_hi);
    copy3(r1,rz1_hi);
    copy3(r2,rx2_lo);
    copy3(r2,ry2_lo);
    copy3(r2,rz2_lo);
    copy3(r2,rx2_hi);
    copy3(r2,ry2_hi);
    copy3(r2,rz2_hi);
    rx1_lo[0] -= del;
    ry1_lo[1] -= del;
    rz1_lo[2] -= del;
    rx1_hi[0] += del;
    ry1_hi[1] += del;
    rz1_hi[2] += del;
    rx2_lo[0] -= del;
    ry2_lo[1] -= del;
    rz2_lo[2] -= del;
    rx2_hi[0] += del;
    ry2_hi[1] += del;
    rz2_hi[2] += del;
    
    double ux1_lo,uy1_lo,uz1_lo;
    double ux1_hi,uy1_hi,uz1_hi;
    double ux2_lo,uy2_lo,uz2_lo;
    double ux2_hi,uy2_hi,uz2_hi;

    geomsemi(rx1_lo,r2,p1,p2,p1,param,basis);
    fsemi(param,ux1_lo,flocal);
    geomsemi(ry1_lo,r2,p1,p2,p1,param,basis);
    fsemi(param,uy1_lo,flocal);
    geomsemi(rz1_lo,r2,p1,p2,p1,param,basis);
    fsemi(param,uz1_lo,flocal);
    geomsemi(rx1_hi,r2,p1,p2,p1,param,basis);
    fsemi(param,ux1_hi,flocal);
    geomsemi(ry1_hi,r2,p1,p2,p1,param,basis);
    fsemi(param,uy1_hi,flocal);
    geomsemi(rz1_hi,r2,p1,p2,p1,param,basis);
    fsemi(param,uz1_hi,flocal);

    geomsemi(r1,rx2_lo,p1,p2,p1,param,basis);
    fsemi(param,ux2_lo,flocal);
    geomsemi(r1,ry2_lo,p1,p2,p1,param,basis);
    fsemi(param,uy2_lo,flocal);
    geomsemi(r1,rz2_lo,p1,p2,p1,param,basis);
    fsemi(param,uz2_lo,flocal);
    geomsemi(r1,rx2_hi,p1,p2,p1,param,basis);
    fsemi(param,ux2_hi,flocal);
    geomsemi(r1,ry2_hi,p1,p2,p1,param,basis);
    fsemi(param,uy2_hi,flocal);
    geomsemi(r1,rz2_hi,p1,p2,p1,param,basis);
    fsemi(param,uz2_hi,flocal);

    double f1x = -0.5 * (ux1_hi - ux1_lo) / del;
    double f1y = -0.5 * (uy1_hi - uy1_lo) / del;
    double f1z = -0.5 * (uz1_hi - uz1_lo) / del;
    double f2x = -0.5 * (ux2_hi - ux2_lo) / del;
    double f2y = -0.5 * (uy2_hi - uy2_lo) / del;
    double f2z = -0.5 * (uz2_hi - uz2_lo) / del;
 
    double f1x_an = flocal[0][0]*basis[0][0]
	      + flocal[0][1]*basis[1][0]
	      + flocal[0][2]*basis[2][0];
    double f1y_an = flocal[0][0]*basis[0][1]
	      + flocal[0][1]*basis[1][1]
	      + flocal[0][2]*basis[2][1];
    double f1z_an = flocal[0][0]*basis[0][2]
	      + flocal[0][1]*basis[1][2]
	      + flocal[0][2]*basis[2][2];
    double f2x_an = flocal[1][0]*basis[0][0]
	      + flocal[1][1]*basis[1][0]
	      + flocal[1][2]*basis[2][0];
    double f2y_an = flocal[1][0]*basis[0][1]
	      + flocal[1][1]*basis[1][1]
	      + flocal[1][2]*basis[2][1];
    double f2z_an = flocal[1][0]*basis[0][2]
	      + flocal[1][1]*basis[1][2]
	      + flocal[1][2]*basis[2][2];
    geomsemi(r1,r2,p1,p2,p1,param,basis);

    double evdwl;
    fsemi(param,evdwl,flocal);
    usemifile << etae*angrec << " " << evdwl;
    fsemifile << etae*angrec << " " 
	    << f1x << " " << f1y << " " << f1z << " "
	    << f2x << " " << f2y << " " << f2z << " " 
	    << f1x_an << " " << f1y_an << " " << f1z_an << " "
	    << f2x_an << " " << f2y_an << " " << f2z_an << " "
	    << std::endl;
    
    double r12[3] = {xi1*cos(alpha2), xi1*sin(alpha2), h};
    double r22[3] = {xi2*cos(alpha2), xi2*sin(alpha2), h};
    geomsemi(r12,r22,p1,p2,p1,param,basis);
    fsemi(param,evdwl,flocal);
    usemifile << " " << evdwl;
    
    double r13[3] = {xi1*cos(alpha3), xi1*sin(alpha3), h};
    double r23[3] = {xi2*cos(alpha3), xi2*sin(alpha3), h};
    geomsemi(r13,r23,p1,p2,p1,param,basis);
    fsemi(param,evdwl,flocal);
    usemifile << " " << evdwl;
    
    double r14[3] = {xi1*cos(alpha4), xi1*sin(alpha4), h};
    double r24[3] = {xi2*cos(alpha4), xi2*sin(alpha4), h};
    geomsemi(r14,r24,p1,p2,p1,param,basis);
    fsemi(param,evdwl,flocal);
    usemifile << " " << evdwl << std::endl;
  }

  usemifile.close();
  fsemifile.close();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMesoCNT::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style mesoCNT requires atom IDs");
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair style mesoCNT requires newton pair off");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMesoCNT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutoff;
}

/* ----------------------------------------------------------------------
   insertion sort list according to corresponding atom ID
------------------------------------------------------------------------- */

void PairMesoCNT::sort(int *list, int size)
{
  int i,j,temp1,temp2;
  int *tag = atom->tag;
  for (int i = 1; i < size; i++) {
    j = i;
    temp1 = list[j-1];
    temp2 = list[j];
    while (j > 0 && tag[temp1] > tag[temp2]) {
      list[j] = temp1;
      list[j-1] = temp2;
      j--;
      temp1 = list[j-1];
      temp2 = list[j];
    }
  }
}

/* ----------------------------------------------------------------------
   read 1D data file
------------------------------------------------------------------------- */

void PairMesoCNT::read_file(const char *file, double *data, 
		double &xstart, double &dx, int ninput)
{
  char line[MAXLINE];

  // open file
  
  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    std::string str("Cannot open file ");
    str += file;
    error->one(FLERR,str.c_str());
  }

  // read values from file

  int cerror = 0;
  int serror = 0;
  double x,xtemp,dxtemp;

  for (int i = 0; i < ninput; i++) {
    if (NULL == fgets(line,MAXLINE,fp)) {
      std::string str("Premature end of file in pair table ");
      str += file;
      error->one(FLERR,str.c_str());
    }
    if (i > 0) xtemp = x;
    if (2 != sscanf(line,"%lg %lg",&x,&data[i])) cerror++;
    if (i == 0) xstart = x;
    else {
      dxtemp = x - xtemp;
      if (i == 1) dx = dxtemp;
      if (fabs(dxtemp - dx)/dx > SMALL) serror++;
    }
  }

  // warn if data was read incompletely, e.g. colums were missing

  if (cerror) { 
    char str[128];
    sprintf(str,"%d of %d lines were incomplete\n"
		    "  or could not be parsed completely\n" 
		    "  in pair table ",cerror,ninput);
    std::string errstr = str;
    errstr += file;
    error->warning(FLERR,errstr.c_str());
  }

  // warn if spacing between data points is not constant
  
  if (serror) {
    char str[128];
    sprintf(str,"%d spacings in first column were different\n"
		    "  from first spacing in pair table ",serror);
    std::string errstr = str;
    errstr += file;
    error->warning(FLERR,errstr.c_str());
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   read 2D data file
------------------------------------------------------------------------- */

void PairMesoCNT::read_file(const char *file, double **data, 
		double &xstart, double &ystart, 
		double &dx, double &dy, int ninput)
{
  char line[MAXLINE];

  // open file
  
  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    std::string str("Cannot open file ");
    str += file;
    error->one(FLERR,str.c_str());
  }

  // read values from file

  int cerror = 0;
  int sxerror = 0;
  int syerror = 0;
  double x,y,xtemp,ytemp,dxtemp,dytemp;

  for (int i = 0; i < ninput; i++) {
    if (i > 0) xtemp = x;
    for (int j = 0; j < ninput; j++) {
      if (NULL == fgets(line,MAXLINE,fp)) {
        std::string str("Premature end of file in pair table ");
        str += file;
        error->one(FLERR,str.c_str());
      }
      if (j > 0) ytemp = y;
      if (3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[i][j])) cerror++;
      if (i == 0 && j == 0) ystart = y;
      if (j > 0) {
	dytemp = y - ytemp;
    	if (j == 1) dy = dytemp;
	if (fabs(dytemp - dy)/dy > SMALL) syerror++;
      }
    }
    if (i == 0) xstart = x;
    else {
      dxtemp = x - xtemp;
      if (i == 1) dx = dxtemp;
      if (fabs(dxtemp - dx)/dx > SMALL) sxerror++;
    }
  }

  // warn if data was read incompletely, e.g. colums were missing

  if (cerror) { 
    char str[128];
    sprintf(str,"%d of %d lines were incomplete\n"
		    "  or could not be parsed completely\n" 
		    "  in pair table ",cerror,ninput*ninput);
    std::string errstr = str;
    errstr += file;
    error->warning(FLERR,errstr.c_str());
  }

  // warn if spacing between data points is not constant
  
  if (sxerror) {
    char str[128];
    sprintf(str,"%d spacings in first column were different\n"
		    "  from first spacing in pair table ",sxerror);
    std::string errstr = str;
    errstr += file;
    error->warning(FLERR,errstr.c_str());
  }
  if (syerror) {
    char str[128];
    sprintf(str,"%d spacings in first column were different\n"
		    "  from first spacing in pair table ",syerror);
    std::string errstr = str;
    errstr += file;
    error->warning(FLERR,errstr.c_str());
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   compute cubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double *data, double **coeff,
		double dx, int size)
{
  double *u = data;
  double **g = coeff;
  int n = size;

  double d,*p,*bprime,*dprime,**b;
  memory->create(p,n,"pair:p");
  memory->create(b,n,n,"pair:b");
  memory->create(bprime,n,"pair:bprime");
  memory->create(dprime,n,"pair:dprime");

  double dxrec = 1.0 / dx;
  double dxsqrec = dxrec * dxrec;
  double dxcbrec = dxrec * dxsqrec;

  double ax[4][4] =
  {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {-3*dxsqrec, -2*dxrec, 3*dxsqrec, -dxrec},
    {2*dxcbrec, dxsqrec, -2*dxcbrec, dxsqrec}
  };

  // compute finite difference derivatives at boundaries
  
  p[0] = (u[1] - u[0]) * dxrec;
  p[n-1] = (u[n-1] - u[n-2]) * dxrec;

  // compute derivatives inside domain
  
  for (int i = 1; i < n-1; i++) {
    if (i > 1) b[i][i-1] = dx;
    b[i][i] = 4 * dx;
    if (i < n-2) b[i][i+1] = dx;
  }
  bprime[1] = b[1][1];
  for(int i = 2; i < n-1; i++)
    bprime[i] = b[i][i] - b[i][i-1]*b[i-1][i]/bprime[i-1];

  for (int i = 1; i < n-1; i++) {
    d = 3 * (u[i+1] - u[i-1]);
    if (i == 1) d -= dx * p[i-1];
    if (i == n-2) d -= dx * p[i+1];
    dprime[i] = d;
    if (i != 1) dprime[i] -= b[i][i-1] * dprime[i-1] / bprime[i-1];
  }

  p[n-2] = dprime[n-2] / bprime[n-2];
  for (int i = n-3; i > 0; i--)
    p[i] = (dprime[i] - b[i][i+1]*p[i+1]) / bprime[i];

  // compute spline coefficients

  for (int i = 1; i < n; i++) {
    for (int j = 0; j < 4; j++)
      g[i][j] = 0;

    double k[4] = {u[i-1], p[i-1], u[i], p[i]};

    for (int j = 0; j < 4; j++)
      for (int l = 0; l < 4; l++)
        g[i][j] += ax[j][l] * k[l];
  }

  memory->destroy(p);
  memory->destroy(b);
  memory->destroy(bprime);
  memory->destroy(dprime);
}

/* ----------------------------------------------------------------------
   compute bicubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ****coeff,
		double dx, double dy, int size)
{
  double **u = data;
  double ****g = coeff;
  int n = size;

  double d,*bprime,*dprime,**p,**q,**s,**b;
  memory->create(p,n,n,"pair:p");
  memory->create(q,n,n,"pair:q");
  memory->create(s,n,n,"pair:s");
  memory->create(b,n,n,"pair:b");
  memory->create(bprime,n,"pair:bprime");
  memory->create(dprime,n,"pair:dprime");

  double dxrec = 1.0 / dx;
  double dyrec = 1.0 / dy;
  double dxsqrec = dxrec * dxrec;
  double dysqrec = dyrec * dyrec;
  double dxcbrec = dxrec * dxsqrec;
  double dycbrec = dyrec * dysqrec;

  double ax[4][4] =
  {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {-3*dxsqrec, -2*dxrec, 3*dxsqrec, -dxrec},
    {2*dxcbrec, dxsqrec, -2*dxcbrec, dxsqrec}
  };
  double ay[4][4] =
  {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {-3*dysqrec, -2*dyrec, 3*dysqrec, -dyrec},
    {2*dycbrec, dysqrec, -2*dycbrec, dysqrec}
  };

  // compute finite difference derivatives at boundaries

  for (int i = 0; i < n; i++) {
    p[0][i] = (u[1][i] - u[0][i]) * dxrec;
    p[n-1][i] = (u[n-1][i] - u[n-2][i]) * dxrec;
  }

  for (int i = 0; i < n; i++) {
    q[i][0] = (u[i][1] - u[i][0]) * dyrec;
    q[i][n-1] = (u[i][n-1] - u[i][n-2]) * dyrec;
  }

  s[0][0] = (p[0][1] - p[0][0]) * dyrec;
  s[0][n-1] = (p[0][n-1] - p[0][n-2]) * dyrec;
  s[n-1][0] = (p[n-1][1] - p[n-1][0]) * dyrec;
  s[n-1][n-1] = (p[n-1][n-1] - p[n-1][n-2]) * dyrec;

  // compute derivatives inside domain

  // sweep in x

  for (int i = 1; i < n-1; i++) {
    if (i > 1) b[i][i-1] = dx;
    b[i][i] = 4 * dx;
    if (i < n-2) b[i][i+1] = dx;
  }
  bprime[1] = b[1][1];
  for(int i = 2; i < n-1; i++)
    bprime[i] = b[i][i] - b[i][i-1]*b[i-1][i]/bprime[i-1];

  // compute p
  
  for (int j = 0; j < n; j++) {
    for (int i = 1; i < n-1; i++) {
      d = 3 * (u[i+1][j] - u[i-1][j]);
      if (i == 1) d -= dx * p[i-1][j];
      if (i == n-2) d -= dx * p[i+1][j];
      dprime[i] = d;
      if (i != 1) dprime[i] -= b[i][i-1] * dprime[i-1] / bprime[i-1];
    }

    p[n-2][j] = dprime[n-2] / bprime[n-2];
    for (int i = n-3; i > 0; i--)
      p[i][j] = (dprime[i] - b[i][i+1]*p[i+1][j]) / bprime[i];
  }

  // compute s

  for (int j = 0; j < n; j += n-1) {
    for (int i = 1; i < n-1; i++) {
      d = 3 * (q[i+1][j] - q[i-1][j]);
      if (i == 1) d -= dx * s[i-1][j];
      if (i == n-2) d -= dx * s[i+1][j];
      dprime[i] = d;
      if (i != 1) dprime[i] -= b[i][i-1] * dprime[i-1] / bprime[i-1];
    }

    s[n-2][j] = dprime[n-2] / bprime[n-2];
    for (int i = n-3; i > 0; i--)
      s[i][j] = (dprime[i] - b[i][i+1]*s[i+1][j]) / bprime[i];
  }

  // sweep in y
  
  for (int i = 1; i < n-1; i++) {
    if (i > 1) b[i][i-1] = dy;
    b[i][i] = 4 * dy;
    if (i < n-2) b[i][i+1] = dy;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < n-1; i++)
    bprime[i] = b[i][i] - b[i][i-1]*b[i-1][i]/bprime[i-1];

  // compute q

  for (int i = 0; i < n; i++) {
    for (int j = 1; j < n-1; j++) {
      d = 3 * (u[i][j+1] - u[i][j-1]);
      if(j == 1) d -= dy * q[i][j-1];
      if(j == n-2) d -= dy * q[i][j+1];
      dprime[j] = d;
      if(j != 1) dprime[j] -= b[j][j-1] * dprime[j-1] / bprime[j-1];
    }

    q[i][n-2] = dprime[n-2] / bprime[n-2];
    for (int j = n-3; j > 0; j--)
      q[i][j] = (dprime[j] - b[j][j+1]*q[i][j+1]) / bprime[j];
  }

  // compute s

  for (int i = 0; i < n; i++) {
    for (int j = 1; j < n-1; j++) {
      d = 3 * (p[i][j+1] - p[i][j-1]);
      if(j == 1) d -= dy * s[i][j-1];
      if(j == n-2) d -= dy * s[i][j+1];
      dprime[j] = d;
      if(j != 1) dprime[j] -= b[j][j-1] * dprime[j-1] / bprime[j-1];
    }

    s[i][n-2] = dprime[n-2] / bprime[n-2];
    for (int j = n-3; j > 0; j--)
      s[i][j] = (dprime[j] - b[j][j+1]*s[i][j+1]) / bprime[j];
  }

  for (int i = 1; i < n; i++)
    for (int j = 1; j < n; j++) {
      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++)
	  g[i][j][l][m] = 0;
      
      double k[4][4] =
      {
        {u[i-1][j-1], q[i-1][j-1], u[i-1][j], q[i-1][j]},
        {p[i-1][j-1], s[i-1][j-1], p[i-1][j], s[i-1][j]},
        {u[i][j-1], q[i][j-1], u[i][j], q[i][j]},
        {p[i][j-1], s[i][j-1], p[i][j], s[i][j]}
      };
      
      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
	    for (int o = 0; o < 4; o++)
	      g[i][j][l][m] += ax[l][n] * k[n][o] * ay[m][o];
    }

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(s);
  memory->destroy(b);
  memory->destroy(bprime);
  memory->destroy(dprime);
}


/* ----------------------------------------------------------------------
   cubic spline evaluation
------------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = ceil((x - xstart)/dx);

  // linear extrapolation

  if (i < 1) return coeff[1][0] + coeff[1][1]*(x - xstart);
  
  // constant extrapolation

  else if (i > coeff_size-1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }

  // cubic interpolation

  double xlo = xstart + (i-1)*dx;
  double xbar = x - xlo;

  return coeff[i][0]
	  + xbar*(coeff[i][1] + xbar*(coeff[i][2] + xbar*coeff[i][3]));
}

/* ----------------------------------------------------------------------
   cubic spline derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dspline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = ceil((x - xstart)/dx);

  // constant extrapolation

  if (i < 1) return coeff[1][1];

  // constant extrapolation
  
  else if (i > coeff_size-1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }

  // cubic interpolation

  double xlo = xstart + (i-1)*dx;
  double xbar = x - xlo;

  return coeff[i][1] + xbar*(2*coeff[i][2] + 3*xbar*coeff[i][3]);
}

/* ----------------------------------------------------------------------
   bicubic spline evaluation
------------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart)/dx);
  int j = ceil((y - ystart)/dy);
  
  // constant extrapolation
  
  if (i < 1) {
    i = 1;
    x = xstart;
  }
  else if (i > coeff_size-1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  }
  else if (j > coeff_size-1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  // cubic interpolation

  double xlo = xstart + (i-1)*dx;
  double ylo = ystart + (j-1)*dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][0]
	  + ybar*(coeff[i][j][0][1]
	  + ybar*(coeff[i][j][0][2]
	  + ybar*(coeff[i][j][0][3])));
  double y1 = coeff[i][j][1][0]
	  + ybar*(coeff[i][j][1][1]
	  + ybar*(coeff[i][j][1][2]
	  + ybar*(coeff[i][j][1][3])));
  double y2 = coeff[i][j][2][0]
	  + ybar*(coeff[i][j][2][1]
	  + ybar*(coeff[i][j][2][2]
	  + ybar*(coeff[i][j][2][3])));
  double y3 = coeff[i][j][3][0]
	  + ybar*(coeff[i][j][3][1]
	  + ybar*(coeff[i][j][3][2]
	  + ybar*(coeff[i][j][3][3])));

  return y0 + xbar*(y1 + xbar*(y2 + xbar*y3));
}

/* ----------------------------------------------------------------------
   bicubic spline partial x derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dxspline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart)/dx);
  int j = ceil((y - ystart)/dy);
  
  // constant extrapolation
  
  if (i < 1) {
    i = 1;
    x = xstart;
  }
  else if (i > coeff_size-1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  }
  else if (j > coeff_size-1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  // cubic interpolation

  double xlo = xstart + (i-1)*dx;
  double ylo = ystart + (j-1)*dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][0]
	  + ybar*(coeff[i][j][0][1]
	  + ybar*(coeff[i][j][0][2]
	  + ybar*(coeff[i][j][0][3])));
  double y1 = coeff[i][j][1][0]
	  + ybar*(coeff[i][j][1][1]
	  + ybar*(coeff[i][j][1][2]
	  + ybar*(coeff[i][j][1][3])));
  double y2 = coeff[i][j][2][0]
	  + ybar*(coeff[i][j][2][1]
	  + ybar*(coeff[i][j][2][2]
	  + ybar*(coeff[i][j][2][3])));
  double y3 = coeff[i][j][3][0]
	  + ybar*(coeff[i][j][3][1]
	  + ybar*(coeff[i][j][3][2]
	  + ybar*(coeff[i][j][3][3])));

  return y1 + xbar*(2*y2 + 3*xbar*y3);
}

/* ----------------------------------------------------------------------
   bicubic spline partial y derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dyspline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart)/dx);
  int j = ceil((y - ystart)/dy);
  
  // constant extrapolation
  
  if (i < 1) {
    i = 1;
    x = xstart;
  }
  else if (i > coeff_size-1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  }
  else if (j > coeff_size-1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  // cubic interpolation

  double xlo = xstart + (i-1)*dx;
  double ylo = ystart + (j-1)*dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][1]
	  + ybar*(2*coeff[i][j][0][2]
	  + 3*ybar*coeff[i][j][0][3]);
  double y1 = coeff[i][j][1][1]
	  + ybar*(2*coeff[i][j][1][2]
	  + 3*ybar*coeff[i][j][1][3]);
  double y2 = coeff[i][j][2][1]
	  + ybar*(2*coeff[i][j][2][2]
	  + 3*ybar*coeff[i][j][2][3]);
  double y3 = coeff[i][j][3][1]
	  + ybar*(2*coeff[i][j][3][2]
	  + 3*ybar*coeff[i][j][3][3]);

  return y0 + xbar*(y1 + xbar*(y2 + xbar*y3));
}


/* ----------------------------------------------------------------------
   geometric parameters for infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::geominf(const double *r1, const double *r2, 
		const double *p1, const double *p2, 
		double *param, double **basis)
{
  double r[3],p[3],delr[3],m[3],l[3],rbar[3],pbar[3],delrbar[3];
  double psil[3],psim[3],dell_psim[3],delpsil_m[3];
  double delr1[3],delr2[3],delp1[3],delp2[3];

  double *ex = basis[0];
  double *ey = basis[1];
  double *ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);

  sub3(p,r,delr);

  sub3(r2,r1,l);
  norm3(l);
  sub3(p2,p1,m);
  norm3(m);

  double psi = dot3(l,m);
  if (psi > 1.0) psi = 1.0;
  else if (psi < -1.0) psi = -1.0;
  double denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  double taur,taup;
  
  // parallel case
  
  if (denom < SWITCH) {
    taur = dot3(delr,l);
    taup = 0;
  }

  // non-parallel case

  else {
    double frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);

  double h = len3(delrbar);

  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1.0/h,ex);
  cross3(ez,ex,ey);

  double alpha;
  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  double xi1 = dot3(delr1,l);
  double xi2 = dot3(delr2,l);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
}

/* ----------------------------------------------------------------------
   geometric parameters for semi-infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::geomsemi(const double *r1, const double *r2, 
		const double *p1, const double *p2, const double *qe,
		double *param, double **basis)
{
  double r[3],p[3],delr[3],m[3],l[3],rbar[3],pbar[3],delrbar[3];
  double psil[3],psim[3],dell_psim[3],delpsil_m[3];
  double delr1[3],delr2[3],delp1[3],delp2[3],delpqe[3];

  double *ex = basis[0];
  double *ey = basis[1];
  double *ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);

  sub3(p,r,delr);

  sub3(r2,r1,l);
  norm3(l);
  sub3(p2,p1,m);
  norm3(m);

  double psi = dot3(l,m);
  if (psi > 1.0) psi = 1.0;
  else if (psi < -1.0) psi = -1.0;
  double denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  sub3(p,qe,delpqe);
  double rhoe = dot3(delpqe,m);
  double taur,taup;
  double etae;

  // parallel case
  
  if (denom < SWITCH) {
    taur = dot3(delr,l) - rhoe*psi;
    taup = -rhoe;
    etae = 0;
  }

  // non-parallel case

  else {
    double frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
    etae = -rhoe - taup;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);

  double h = len3(delrbar);

  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1.0/h,ex);
  cross3(ez,ex,ey);

  double alpha;
  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  double xi1 = dot3(delr1,l);
  double xi2 = dot3(delr2,l);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[6] = etae; 
}


/* ----------------------------------------------------------------------
   weight for substitute CNT chain
------------------------------------------------------------------------- */

double PairMesoCNT::weight(const double *r1, const double *r2,
		const double *p1, const double *p2)
{
  double r[3],p[3],delr[3];

  add3(r1,r2,r);
  add3(p1,p2,p);
  scale3(0.5,r);
  scale3(0.5,p);

  double temp = sqrt(0.25*distsq3(r1,r2) + rsq);
  double rhoc = temp + sqrt(0.25*distsq3(p1,p2) + rsq) + rc;
  double rhomin = 0.72 * rhoc;
  double rho = sqrt(distsq3(r,p));

  return s((rho-rhomin)/(rhoc-rhomin));
}

/* ----------------------------------------------------------------------
   forces for infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::finf(const double *param, double &evdwl, double **f)
{
  double h = param[0] * angrec;
  double alpha = param[1];
  double xi1 = param[2] * angrec;
  double xi2 = param[3] * angrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha * sin_alpha;

  // parallel case 

  if(sin_alphasq < SWITCH) {
    double ubar = spline(h,hstart_uinf,delh_uinf,uinf_coeff,uinf_points);
    double delxi = xi2 - xi1;
    f[0][0] = 0.5 * delxi 
	    * dspline(h,hstart_uinf,delh_uinf,uinf_coeff,uinf_points)
	    * funit;
    f[1][0] = f[0][0];
    f[0][1] = 0;
    f[1][1] = 0;
    f[0][2] = ubar * funit;
    f[1][2] = -f[0][2];
    evdwl = ubar * delxi * e;
  }
 
  // non-parallel case
  
  else {
    double sin_alpharec = 1.0 / sin_alpha;
    double sin_alphasqrec = sin_alpharec * sin_alpharec;
    double cos_alpha = cos(alpha);
    double cot_alpha = cos_alpha * sin_alpharec;

    double omega = 1.0 / (1.0 - comega*sin_alphasq);
    double c1 = omega * sin_alpha;
    double c1rec = 1.0 / c1;
    double domega = 2 * comega * cos_alpha * c1 * omega;

    double gamma_orth = 
	    spline(h,hstart_gamma,delh_gamma,gamma_coeff,gamma_points);
    double dgamma_orth = 
	    dspline(h,hstart_gamma,delh_gamma,gamma_coeff,gamma_points);
    double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;
    double gammarec = 1.0 / gamma;
    double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
    double dh_gamma = dgamma_orth * sin_alphasq;

    double zeta1 = xi1 * c1;
    double zeta2 = xi2 * c1;  
    
    double smooth = s5((h-d_ang-DELTA1)/(DELTA2-DELTA1));
    double dsmooth = ds5((h-d_ang-DELTA1)/(DELTA2-DELTA1));
    double g = d_ang + DELTA2;
    double hsq = h * h;

    double zetaminbar;
    if (h >= g) zetaminbar = 0;
    else zetaminbar = sqrt(g*g - hsq);
    double zetamin = smooth * zetaminbar;
    double zetamax = sqrt(cutoffsq_ang - hsq);

    double dzetaminbar;
    if (h >= g) dzetaminbar = 0;
    else dzetaminbar = -h / zetaminbar;
    double dzetamin = 
	    dzetaminbar*smooth + zetaminbar*dsmooth/(DELTA2-DELTA1);
    double dzetamax = -h / zetamax;

    double zeta_rangerec = 1.0 / (zetamax - zetamin);
    double delzeta1 = fabs(zeta1) - zetamin;
    double delzeta2 = fabs(zeta2) - zetamin;

    double psi1 = delzeta1 * zeta_rangerec;
    double psi2 = delzeta2 * zeta_rangerec;

    double phi1 = spline(h,psi1,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);
    double dh_phibar1 = dxspline(h,psi1,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);
    double dpsi_phibar1 = dyspline(h,psi1,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);
    double phi2 = spline(h,psi2,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);
    double dh_phibar2 = dxspline(h,psi2,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);
    double dpsi_phibar2 = dyspline(h,psi2,hstart_phi,psistart_phi,
		    delh_phi,delpsi_phi,phi_coeff,phi_points);

    double dzeta_range = dzetamax - dzetamin;
    double dh_psi1 = -zeta_rangerec * (dzetamin + dzeta_range*psi1);
    double dh_psi2 = -zeta_rangerec * (dzetamin + dzeta_range*psi2);
    double dh_phi1 = dh_phibar1 + dpsi_phibar1*dh_psi1;
    double dh_phi2 = dh_phibar2 + dpsi_phibar2*dh_psi2;

    double dzeta_phi1 = dpsi_phibar1 * zeta_rangerec;
    double dzeta_phi2 = dpsi_phibar2 * zeta_rangerec;

    if (zeta1 < 0) {
      phi1 *= -1;
      dh_phi1 *= -1;
    }
    if (zeta2 < 0) {
      phi2 *= -1;
      dh_phi2 *= -1;
    }

    double deldzeta_phi = dzeta_phi2 - dzeta_phi1;
    
    double c2 = gamma * c1rec;
    double u = c2 * (phi2 - phi1);
    double c3 = u * gammarec;

    double dh_u = dh_gamma*c3 + c2*(dh_phi2 - dh_phi1);
    double dalpha_u = dalpha_gamma*c3 
	    + c1rec*(domega*sin_alpha + omega*cos_alpha)
	    * (gamma*(xi2*dzeta_phi2 - xi1*dzeta_phi1) - u);

    double lrec = 1.0 / (xi2 - xi1);
    double cx = h * gamma * sin_alphasqrec * deldzeta_phi;
    double cy = gamma * cot_alpha * deldzeta_phi;

    f[0][0] = lrec * (xi2*dh_u - cx) * funit;
    f[1][0] = lrec * (-xi1*dh_u + cx) * funit;
    f[0][1] = lrec * (dalpha_u - xi2*cy) * funit;
    f[1][1] = lrec * (-dalpha_u + xi1*cy) * funit;
    f[0][2] = gamma * dzeta_phi1 * funit;
    f[1][2] = -gamma * dzeta_phi2 * funit;
    evdwl = u * e;
  }
}

/* ----------------------------------------------------------------------
   forces for semi-infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::fsemi(const double *param, double &evdwl, double **f)
{
  double h = param[0] * angrec;
  double alpha = param[1];
  double xi1 = param[2] * angrec;
  double xi2 = param[3] * angrec;
  double etae = param[6] * angrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha * sin_alpha;
  double cos_alpha = cos(alpha);

  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double omegasq = omega * omega;
  double domega = 2 * comega * sin_alpha * cos_alpha * omegasq;

  double theta = 1.0 - ctheta*sin_alphasq;
  double dtheta = -2 * ctheta * sin_alpha * cos_alpha;

  double c1 = omega * sin_alpha;
  double c1sq = c1 * c1;
  double c2 = theta * etae;

  double gamma_orth = spline(h,hstart_gamma,delh_gamma,
		  gamma_coeff,gamma_points);
  double dgamma_orth = dspline(h,hstart_gamma,delh_gamma,
		  gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1)*sin_alphasq;
  double gammarec = 1.0 / gamma;
  double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
  double dh_gamma = dgamma_orth * sin_alphasq;

  double delxi = (xi2 - xi1) / (QUADRATURE - 1);
  double c3 = delxi * gamma;

  double jh = 0;
  double jh1 = 0;
  double jh2 = 0;
  double jxi = 0;
  double jxi1 = 0;
  double ubar = 0;

  for (int i = 0; i < QUADRATURE; i++) {
    double xibar = xi1 + i*delxi;
    double g = xibar * c1;
    double hbar = sqrt(h*h + g*g);
    double thetabar = xibar*cos_alpha - c2;

    double c = 1.0;
    if (i == 0 || i == QUADRATURE-1) c = 0.5;

    double u = c * spline(hbar,thetabar,hstart_usemi,xistart_usemi,
		    delh_usemi,delxi_usemi,usemi_coeff,usemi_points);
    double uh;
    if (hbar == 0) uh = 0;
    else uh = c / hbar * dxspline(hbar,thetabar,hstart_usemi,xistart_usemi,
		    delh_usemi,delxi_usemi,usemi_coeff,usemi_points);
    double uxi = c * dyspline(hbar,thetabar,hstart_usemi,xistart_usemi,
		    delh_usemi,delxi_usemi,usemi_coeff,usemi_points);

    double uh1 = xibar * uh;
    jh += uh;
    jh1 += uh1;
    jh2 += xibar * uh1;
    jxi += uxi;
    jxi1 += xibar * uxi;
    ubar += u;
  }

  jh *= c3;
  jh1 *= c3;
  jh2 *= c3;
  jxi *= c3;
  jxi1 *= c3;
  ubar *= c3;

  double c4 = gammarec * ubar;
  double dh_ubar = dh_gamma*c4 + h*jh;
  double dalpha_ubar = dalpha_gamma*c4
	  + c1*(domega*sin_alpha + omega*cos_alpha)*jh2
	  - sin_alpha*jxi1 - dtheta*etae*jxi;

  double cx = h * (omegasq*jh1 + cos_alpha*ctheta*jxi);
  double cy = sin_alpha * (cos_alpha*omegasq*jh1 + (ctheta-1)*jxi);
  double cz1 = c1sq*jh1 + cos_alpha*jxi;
  double cz2 = c1sq*jh2 + cos_alpha*jxi1;

  double lrec = 1.0 / (xi2 - xi1);
  f[0][0] = lrec * (xi2*dh_ubar - cx) * funit;
  f[1][0] = lrec * (cx - xi1*dh_ubar) * funit;
  f[0][1] = lrec * (dalpha_ubar - xi2*cy) * funit;
  f[1][1] = lrec * (xi1*cy - dalpha_ubar) * funit;
  f[0][2] = lrec * (cz2 + ubar - xi2*cz1) * funit;
  f[1][2] = lrec * (xi1*cz1 - cz2 - ubar) * funit;
  evdwl = ubar * e;
}
