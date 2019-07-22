#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include "pair_mesocnt.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define SMALL 1.0e-6
#define LOW 1.0e-6
#define HIGH 1.0e-6
#define DELTA1 1.0
#define DELTA2 2.0
#define DELTAREC 1.0

/* ---------------------------------------------------------------------- */

PairMesoCNT::PairMesoCNT(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  writedata = 1;

  
}

/* ---------------------------------------------------------------------- */

PairMesoCNT::~PairMesoCNT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(gamma_data);
    memory->destroy(uinf_data);
    memory->destroy(starth_usemi);
    memory->destroy(startzeta_phi);
    memory->destroy(delh_usemi);
    memory->destroy(delzeta_phi);

    memory->destroy(usemi_data);
    memory->destroy(phi_data);
    memory->destroy(gamma_coeff);
    memory->destroy(uinf_coeff);
    
    memory->destroy(usemi_coeff);
    memory->destroy(phi_coeff);
 
    memory->destroy(redlist);
    memory->destroy(chain);
    memory->destroy(nchain);
    memory->destroy(end);

    memory->destroy(p1);
    memory->destroy(p2);
    memory->destroy(m);
    memory->destroy(dw1);
    memory->destroy(dw2);
    memory->destroy(sumdw1);
    memory->destroy(sumdw2);
    memory->destroy(fchain11);
    memory->destroy(fchain12);
    memory->destroy(fchain21);
    memory->destroy(fchain22);
    memory->destroy(fchain1);
    memory->destroy(fchain2);

    memory->destroy(param);
    memory->destroy(wvector);

    memory->destroy(dw1q1);
    memory->destroy(dw1q2);
    memory->destroy(dw2q1);
    memory->destroy(dw2q2);
    memory->destroy(sumdw1q1);
    memory->destroy(sumdw1q2);
    memory->destroy(sumdw2q1);
    memory->destroy(sumdw2q2);
    memory->destroy(sumdw1p1);
    memory->destroy(sumdw1p2);
    memory->destroy(sumdw2p1);
    memory->destroy(sumdw2p2);
    memory->destroy(dp1_dr1);
    memory->destroy(dp2_dr1);
    memory->destroy(dp1_dr2);
    memory->destroy(dp2_dr2);

    memory->destroy(flocal);
    memory->destroy(fglobal);
    memory->destroy(basis);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{
  int inum,i,i1,i2,j,j1,j2,jj,jj1,jj2,j1num,j2num,listmax,numred,inflag,n;
  int ind,loc1,loc2,cid,cnum,cn,tag1,tag2,idprev,idnext;
  double w,sumw,sumwreq,evdwl;
  double *r1,*r2,*q1,*q2,*qe;

  int *ilist,*j1list,*j2list,*numneigh,**firstneigh;
 
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int *tag = atom->tag;
  int *mol = atom->molecule;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
 
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for(n = 0; n < nbondlist; n++){
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    i1 &= NEIGHMASK;
    i2 &= NEIGHMASK;
    
    r1 = x[i1];
    r2 = x[i2];
    
    // reduce neighbors to common list
    
    if(i1 > nlocal-1){
      j1num = 0;
    }
    else{
      j1list = firstneigh[i1];
      j1num = numneigh[i1];
    }
    if(i2 > nlocal-1){
      j2num = 0;
    }
    else{
      j2list = firstneigh[i2];
      j2num = numneigh[i2];
    }

    listmax = j1num + j2num;
    if(listmax < 2) continue;
    if(listmax > redlist_size){
      redlist_size = 2 * listmax;
      memory->grow(redlist,redlist_size,"pair:redlist");
    }

    numred = 0;
    for(jj1 = 0; jj1 < j1num; jj1++){
      ind = j1list[jj1];
      if ((mol[ind] == mol[i1]) && (abs(tag[ind] - tag[i1]) < 5)) continue;
      redlist[numred++] = ind;
    }
    for(jj2 = 0; jj2 < j2num; jj2++){
      for(jj1 = 0; jj1 < j1num; jj1++){
        if(j1list[jj1] == j2list[jj2]){
	  inflag = 1;
	  break;
	}
      }
      if(inflag){
        inflag = 0;
	continue;
      }
      ind = j2list[jj2];
      if ((mol[ind] == mol[i2]) && (abs(tag[ind] - tag[i2]) < 5)) continue; 
      redlist[numred++] = ind;
    }
    
    if (numred < 2) continue;
    
    // insertion sort according to atom-id
    
    sort(redlist,numred);
 
    // split into connected chains
    
    cid = 0;
    cnum = 0;
    
    if(numred > chain_size){
      chain_size = 2 * numred;
      memory->destroy(chain);
      memory->create(chain,chain_size,chain_size,"pair:chain");
      memory->grow(nchain,chain_size,"pair:nchain");
    }
    
    for(jj = 0; jj < numred-1; jj++){
      j1 = redlist[jj];
      j2 = redlist[jj+1];
      chain[cid][cnum++] = j1;
      if(abs(tag[j1] - tag[j2]) != 1 || mol[j1] != mol[j2]){
        nchain[cid++] = cnum;
	cnum = 0;
      }
    }
    chain[cid][cnum++] = redlist[numred-1];
    nchain[cid++] = cnum;

    // check for ends
 
    if(cid > end_size){
      end_size = 2 * cid;
      memory->grow(end,end_size,"pair:end_size");
    }
    
    for(i = 0; i < cid; i++){
      cn = nchain[i];
      loc1 = chain[i][0];
      loc2 = chain[i][cn-1];
      tag1 = tag[loc1];
      tag2 = tag[loc2];
      end[i] = 0;
      if(tag1 == 1) end[i] = 1;
      else{
        idprev = atom->map(tag1-1);
	if(idprev == -1 || mol[loc1] != mol[idprev]) end[i] = 1;
      }
      if(tag2 == atom->natoms) end[i] = 2;
      else{
        idnext = atom->map(tag2+1);
	if(idnext == -1 || mol[loc2] != mol[idnext]) end[i] = 2;
      }
    }

    // compute subsitute chains, forces and energies

    using namespace MathExtra;
    
    double chainEnd1,chainEnd2;
    for(i = 0; i < cid; i++){
      if(nchain[i] < 2) continue;

      zero3(p1);
      zero3(p2);
      
      sumw = 0;
      zero3(sumdw1);
      zero3(sumdw2);
      zero(sumdw1q1);
      zero(sumdw1q2);
      zero(sumdw2q1);
      zero(sumdw2q2);

      chainEnd1 = 0;
      chainEnd2 = 0;

      loc1 = chain[i][0];
      loc2 = chain[i][nchain[i]-1];
      if(end[i] == 1) qe = x[loc1];
      if(end[i] == 2) qe = x[loc2];

      for(j = 0; j < nchain[i]-1; j++){
	q1 = x[chain[i][j]];
	q2 = x[chain[i][j+1]];

	weight(r1,r2,q1,q2,wvector);
	w = wvector[0];
	dw1[0] = wvector[1];
	dw1[1] = wvector[2];
	dw1[2] = wvector[3];
	dw2[0] = wvector[4];
	dw2[1] = wvector[5];
	dw2[2] = wvector[6];

	if(w == 0){ 
          if(end[i] == 1 && j == 0) end[i] = 0;
	  else if(end[i] == 2 && j == nchain[i] - 2) end[i] = 0;
	  continue;
	}
	sumw += w;
	add3(dw1,sumdw1,sumdw1);
	add3(dw2,sumdw2,sumdw2);

	chainEnd1 += w * chain[i][j];
	chainEnd2 += w * chain[i][j+1];

	outer(dw1,q1,dw1q1);
	outer(dw1,q2,dw1q2);
	outer(dw2,q1,dw2q1);
	outer(dw2,q2,dw2q2);

	plus(dw1q1,sumdw1q1,sumdw1q1);
	plus(dw1q2,sumdw1q2,sumdw1q2);
	plus(dw2q1,sumdw2q1,sumdw2q1);
	plus(dw2q2,sumdw2q2,sumdw2q2);

	scaleadd3(w,q1,p1,p1);
	scaleadd3(w,q2,p2,p2);
      }
      if (sumw == 0) continue;

      sumwreq = 1.0 / sumw;

      chainEnd1 *= sumwreq;
      chainEnd2 *= sumwreq;

      scale3(sumwreq,p1);
      scale3(sumwreq,p2);
     
      if(end[i] == 0){
        geominf(r1,r2,p1,p2,param,m,basis);
        if(param[0] > cutoff) continue;
        finf(param,r1,r2,p1,p2,m,flocal);
      }
      else if(end[i] == 1){
        geomsemi(r1,r2,p1,p2,qe,param,m,basis);
	if(param[0] > cutoff) continue;
	fsemi(param,r1,r2,p1,p2,m,flocal);
      }
      else{
	geomsemi(r1,r2,p2,p1,qe,param,m,basis);
	if(param[0] > cutoff) continue;
	fsemi(param,r1,r2,p1,p2,m,flocal);
      }
      
      scale3(-sumwreq,sumdw1);
      scale3(-sumwreq,sumdw2);
      
      outer(sumdw1,p1,sumdw1p1);
      outer(sumdw1,p2,sumdw1p2);
      outer(sumdw2,p1,sumdw2p1);
      outer(sumdw2,p2,sumdw2p2);

      scalar(sumwreq,sumdw1q1);
      scalar(sumwreq,sumdw1q2);
      scalar(sumwreq,sumdw2q1);
      scalar(sumwreq,sumdw2q2);

      plus(sumdw1q1,sumdw1p1,dp1_dr1);
      plus(sumdw1q2,sumdw1p2,dp2_dr1);
      plus(sumdw2q1,sumdw2p1,dp1_dr2);
      plus(sumdw2q2,sumdw2p2,dp2_dr2);

      trans_matrix_vector(basis,flocal[0],fglobal[0]);
      trans_matrix_vector(basis,flocal[1],fglobal[1]);
      trans_matrix_vector(basis,flocal[2],fglobal[2]);
      trans_matrix_vector(basis,flocal[3],fglobal[3]);

      matrix_vector(dp1_dr1,fglobal[2],fchain11);
      matrix_vector(dp2_dr1,fglobal[3],fchain12);
      matrix_vector(dp1_dr2,fglobal[2],fchain21);
      matrix_vector(dp2_dr2,fglobal[3],fchain22);

      add3(fchain11,fchain12,fchain1);
      add3(fchain21,fchain22,fchain2);

      f[i1][0] += fglobal[0][0] + fchain1[0];
      f[i1][1] += fglobal[0][1] + fchain1[1];
      f[i1][2] += fglobal[0][2] + fchain1[2];
      f[i2][0] += fglobal[1][0] + fchain2[0];
      f[i2][1] += fglobal[1][1] + fchain2[1];
      f[i2][2] += fglobal[1][2] + fchain2[2];
    
      if(eflag_either){
	if(eflag_global){
          if(end[i] == 0) evdwl = uinf(param);
          else evdwl = usemi(param);
	  eng_vdwl += 0.5 * evdwl;
	}
	if(eflag_atom){
	  eatom[i1] += 0.25 * evdwl;
	  eatom[i2] += 0.25 * evdwl;
	}
      }

      if(vflag_atom){
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

void PairMesoCNT::init_style()
{
  int irequest;
  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::init_one(int i, int j)
{
  return cutoff;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  redlist_size = 64;
  chain_size = 64;
  end_size = 64;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
 
  memory->create(gamma_data,gamma_points,"pair:gamma_data");
  memory->create(uinf_data,uinf_points,"pair:uinf_data");
  memory->create(starth_usemi,usemi_points,"pair:starth_usemi");
  memory->create(startzeta_phi,phi_points,"pair:startzeta_phi");
  memory->create(delh_usemi,usemi_points,"pair:delh_usemi");
  memory->create(delzeta_phi,phi_points,"pair:delzeta_phi");
 
  memory->create(usemi_data,usemi_points,usemi_points,"pair:usemi_data");
  memory->create(phi_data,phi_points,phi_points,"pair:phi_data");
  memory->create(gamma_coeff,gamma_points,4,"pair:gamma_coeff");
  memory->create(uinf_coeff,uinf_points,4,"pair:uinf_coeff");
  
  memory->create(usemi_coeff,usemi_points,usemi_points,4,"pair:usemi_coeff");
  memory->create(phi_coeff,phi_points,phi_points,4,4,"pair:phi_coeff");

  memory->create(redlist,redlist_size,"pair:redlist");
  memory->create(chain,chain_size,chain_size,"pair:chain");
  memory->create(nchain,chain_size,"pair:nchain");
  memory->create(end,end_size,"pair:end");

  memory->create(p1,3,"pair:p1");
  memory->create(p2,3,"pair:p2");
  memory->create(m,3,"pair:m");
  memory->create(dw1,3,"pair:dw1");
  memory->create(dw2,3,"pair:dw2");
  memory->create(sumdw1,3,"pair:sumdw1");
  memory->create(sumdw2,3,"pair:sumdw2");
  memory->create(fchain11,3,"pair:fchain11");
  memory->create(fchain12,3,"pair:fchain12");
  memory->create(fchain21,3,"pair:fchain21");
  memory->create(fchain22,3,"pair:fchain22"); 
  memory->create(fchain1,3,"pair:fchain1");
  memory->create(fchain2,3,"pair:fchain2"); 

  memory->create(param,7,"pair:param");
  memory->create(wvector,7,"pair:wvector");

  memory->create(dw1q1,3,3,"pair:dw1q1");
  memory->create(dw1q2,3,3,"pair:dw1q2");
  memory->create(dw2q1,3,3,"pair:dw2q1");
  memory->create(dw2q2,3,3,"pair:dw2q2");
  memory->create(sumdw1q1,3,3,"pair:sumdw1q1");
  memory->create(sumdw1q2,3,3,"pair:sumdw1q2");
  memory->create(sumdw2q1,3,3,"pair:sumdw2q1");
  memory->create(sumdw2q2,3,3,"pair:sumdw2q2");
  memory->create(sumdw1p1,3,3,"pair:sumdw1p1");
  memory->create(sumdw1p2,3,3,"pair:sumdw1p2");
  memory->create(sumdw2p1,3,3,"pair:sumdw2p1");
  memory->create(sumdw2p2,3,3,"pair:sumdw2p2");
  memory->create(dp1_dr1,3,3,"pair:dp1_dr1");
  memory->create(dp2_dr1,3,3,"pair:dp2_dr1");
  memory->create(dp1_dr2,3,3,"pair:dp1_dr2");
  memory->create(dp2_dr2,3,3,"pair:dp2_dr2");

  memory->create(flocal,4,3,"pair:flocal");
  memory->create(fglobal,4,3,"pair:fglobal");
  memory->create(basis,3,3,"pair:basis");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  gamma_points = force->inumeric(FLERR,arg[0]);
  pot_points = force->inumeric(FLERR,arg[1]);

  uinf_points = 2001;
  usemi_points = 1001;
  phi_points = 2001;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  n = force->inumeric(FLERR,arg[2]);
  sigma = force->numeric(FLERR,arg[3]);
  epsilon = force->numeric(FLERR,arg[4]);
  n_sigma = force->numeric(FLERR,arg[5]);

  gamma_file = arg[6];
  uinf_file = arg[7];
  usemi_file = arg[8];
  phi_file = arg[9];
 
  angstrom = force->angstrom;
  angstromrec = 1 / angstrom;
  qelectron = force->qelectron;
  qelectronrec = 1 / qelectron;
  forceunit = qelectron * angstromrec;

  radius = 1.421*3*n / MY_2PI * angstrom;
  radiussq = radius * radius;
  diameter = 2 * radius;
  rc = 3 * sigma;
  comega = 0.275*(1.0 - 1.0/(1.0 + 0.59*radius*angstromrec));
  ctheta = 0.35 + 0.0226*(radius*angstromrec - 6.785);

  //Parse and bcast data
  MPI_Comm_rank(world,&me);
  if (me == 0) { 
    read_file(gamma_file,gamma_data,&start_gamma,&del_gamma,gamma_points);
    read_file(uinf_file,uinf_data,&start_uinf,&del_uinf,uinf_points);
    read_file(usemi_file,usemi_data,starth_usemi,&startxi_usemi,
		    delh_usemi,&delxi_usemi,usemi_points);
    read_file(phi_file,phi_data,startzeta_phi,&starth_phi,
		    delzeta_phi,&delh_phi,phi_points);
  }
  
  MPI_Bcast(&start_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&del_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&start_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&del_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&startxi_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delxi_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&starth_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_phi,1,MPI_DOUBLE,0,world);

  MPI_Bcast(gamma_data,gamma_points,MPI_DOUBLE,0,world);
  MPI_Bcast(uinf_data,uinf_points,MPI_DOUBLE,0,world);
  MPI_Bcast(starth_usemi,usemi_points,MPI_DOUBLE,0,world);
  MPI_Bcast(startzeta_phi,phi_points,MPI_DOUBLE,0,world); 
  MPI_Bcast(delh_usemi,usemi_points,MPI_DOUBLE,0,world);
  MPI_Bcast(delzeta_phi,phi_points,MPI_DOUBLE,0,world);
  for(int i = 0; i < usemi_points; i++){
    MPI_Bcast(usemi_data[i],usemi_points,MPI_DOUBLE,0,world);
  }
  for(int i = 0; i < phi_points; i++){
    MPI_Bcast(phi_data[i],phi_points,MPI_DOUBLE,0,world);
  }

  spline_coeff(gamma_data,gamma_coeff,gamma_points);
  spline_coeff(uinf_data,uinf_coeff,uinf_points);
  spline_coeff(usemi_data,usemi_coeff,usemi_points);
  spline_coeff(phi_data,phi_coeff,delzeta_phi[0],delh_phi,phi_points);
  printf("Coefficients computed\n");

  int n = atom->ntypes;

  cutoff = rc + diameter;
  cutoffsq = cutoff * cutoff;

  diameter_angstrom = diameter * angstromrec;
  cutoff_angstrom = cutoff * angstromrec;
  cutoffsq_angstrom = cutoff_angstrom * cutoff_angstrom;
  
  for (int i = 1; i <= n; i++){
    for (int j = i; j <= n; j++){
      setflag[i][j] = 1;
      cutsq[i][j] = cutoffsq;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = floor((x - xstart)/dx); 
  if(i < 0){
    i = 0;
    x = xstart;
  }
  else if(i > coeff_size-2){ 
    i = coeff_size - 2;
    x = xstart + (coeff_size-1)*dx;
  }
 
  double xlo = xstart + i*dx;
  double xbar = (x - xlo)/dx;

  return coeff[i][0] + xbar*(coeff[i][1] 
		  + xbar*(coeff[i][2] + xbar*coeff[i][3]));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double y, double *xstart, double ystart,
                double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  else if(i > coeff_size-3){
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }
  
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = spline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = spline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    
    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return a0 + ybar*(a1 + ybar*(a2 + a3*ybar));
}

/* ---------------------------------------------------------------------- */


double PairMesoCNT::spline(double x, double y, double xstart, double ystart,
		double dx, double dy, double ****coeff, int coeff_size)
{
  int j = ceil((y - ystart)/dy);
  if(j < 1){
    j = 1;
    y = ystart;
  }
  else if(j > coeff_size-1){
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  int i = ceil((x - xstart)/dx);
  if(i < 1){
    i = 1;
    x = xstart;
  }
  else if(i > coeff_size-1){
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }
  
  double xlo = xstart + (i-1)*dx;
  double xbar = x - xlo;
  double ylo = ystart + (j-1)*dy;
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

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dspline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = floor((x - xstart)/dx);
  if(i < 0){
    i = 0;
    x = xstart;
  }
  else if(i > coeff_size-2){
    i = coeff_size - 2;
    x = xstart + (coeff_size - 1)*dx;
  }
 
  double xlo = xstart + i*dx;
  double xbar = (x - xlo)/dx;

  return (coeff[i][1] + xbar*(2*coeff[i][2] + 3*xbar*coeff[i][3])) / dx;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dxspline(double x, double y, double *xstart, double ystart,
                double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  if(i > coeff_size-3){
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }

  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y

  double a0, a1, a2, a3;
  double p0, p1, p2, p3;

  p1 = dspline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = dspline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);

  a0 = p1;

  if(i == 0){
    p3 = dspline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = dspline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);

    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = dspline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = dspline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return a0 + ybar*(a1 + ybar*(a2 + a3*ybar));
}


/* ---------------------------------------------------------------------- */

double PairMesoCNT::dxspline(double x, double y, double xstart, double ystart,
		double dx, double dy, double ****coeff, int coeff_size)
{
  int j = ceil((y - ystart)/dy);
  if(j < 1){
    j = 1;
    y = ystart;
  }
  else if(j > coeff_size-1){
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  int i = ceil((x - xstart)/dx);
  if(i < 1){
    i = 1;
    x = xstart;
  }
  else if(i > coeff_size-1){
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }
  
  double xlo = xstart + (i-1)*dx;
  double xbar = x - xlo;
  double ylo = ystart + (j-1)*dy;
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

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dyspline(double x, double y, double *xstart, double ystart,
                double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  if(i > coeff_size-3){
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }

  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y

  double a0, a1, a2, a3;
  double p0, p1, p2, p3;

  p1 = spline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = spline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);

  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);

    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return (a1 + ybar*(2*a2 + 3*a3*ybar)) / dy;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dyspline(double x, double y, double xstart, double ystart,
		double dx, double dy, double ****coeff, int coeff_size)
{
  int j = ceil((y - ystart)/dy);
  if(j < 1){
    j = 1;
    y = ystart;
  }
  else if(j > coeff_size-1){
    j = coeff_size - 1;
    y = ystart + (coeff_size-1)*dy;
  }

  int i = ceil((x - xstart)/dx);
  if(i < 1){
    i = 1;
    x = xstart;
  }
  else if(i > coeff_size-1){
    i = coeff_size - 1;
    x = xstart + (coeff_size-1)*dx;
  }
  
  double xlo = xstart + (i-1)*dx;
  double xbar = x - xlo;
  double ylo = ystart + (j-1)*dy;
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


/* ---------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double *data, double **coeff, int data_size)
{
  for(int i = 0; i < data_size-1; i++){
    if(i == 0){
      coeff[i][0] = data[i];
      coeff[i][1] = 0.25*(-data[i+2] + 6*data[i+1] - 5*data[i]);
      coeff[i][2] = 0;
      coeff[i][3] = 0.25*(data[i+2] - 2*data[i+1] + data[i]);
    }
    else if(i == data_size-2){
      coeff[i][0] = data[i];
      coeff[i][1] = 0.5*(data[i+1] - data[i-1]);
      coeff[i][3] = 0.25*(-data[i+1] + 2*data[i] - data[i-1]);
      coeff[i][2] = -3*coeff[i][3];
    }
    else{
      coeff[i][0] = data[i];
      coeff[i][1] = 0.5*(data[i+1] - data[i-1]);
      coeff[i][2] = 0.5*(-data[i+2] + 4*data[i+1] - 5*data[i] + 2*data[i-1]);
      coeff[i][3] = 0.5*(data[i+2] - 3*data[i+1] + 3*data[i] - data[i-1]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ***coeff, int data_size)
{
  for(int i = 0; i < data_size; i++){
    for(int j = 0; j < data_size-1; j++){
      if(j == 0){
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = data[j+1][i] - data[j][i];
        coeff[i][j][3] = 0.5*(data[j+2][i] - 2*data[j+1][i] + data[j][i]);
        coeff[i][j][2] = -coeff[i][j][3];
      }
      else if(j == data_size-2){
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = 0.5*(data[j+1][i] - data[j-1][i]);
        coeff[i][j][3] = 0.5*(-data[j+1][i] + 2*data[j][i] - data[j-1][i]);
        coeff[i][j][2] = -2*coeff[i][j][3];
      }
      else{
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = 0.5*(data[j+1][i] - data[j-1][i]);
        coeff[i][j][2] = 0.5*(-data[j+2][i] + 4*data[j+1][i]
                        - 5*data[j][i] + 2*data[j-1][i]);
        coeff[i][j][3] = 0.5*(data[j+2][i] - 3*data[j+1][i]
                        + 3*data[j][i] - data[j-1][i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ****coeff, 
		double dx, double dy, int data_size)
{
  double **u = data;
  double ****gamma = coeff;
  int n = data_size;

  double **p, **q, **s, **b;
  memory->create(p,n,n,"pair:p");
  memory->create(q,n,n,"pair:q");
  memory->create(s,n,n,"pair:s");
  memory->create(b,n,n,"pair:b");

  double d;
  double *bprime, *dprime;
  memory->create(bprime,n,"pair:bprime");
  memory->create(dprime,n,"pair:dprime");

  double dxrec = 1.0 / dx;
  double dyrec = 1.0 / dy;
  double dxrecsq = dxrec * dxrec;
  double dyrecsq = dyrec * dyrec;
  double dxreccb = dxrecsq * dxrec;
  double dyreccb = dyrecsq * dyrec;

  double ax[4][4] =
  {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {-3*dxrecsq, -2*dxrec, 3*dxrecsq, -dxrec},
    {2*dxreccb, dxrecsq, -2*dxreccb, dxrecsq}
  };
  double ay[4][4] =
  {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {-3*dyrecsq, -2*dyrec, 3*dyrecsq, -dyrec},
    {2*dyreccb, dyrecsq, -2*dyreccb, dyrecsq}
  };

  double k_ay[4][4];
  
  //Compute finite difference derivatives at boundaries
  
  for(int j = 0; j < n; j++){
    p[0][j] = (u[1][j] - u[0][j]) * dxrec;
    p[n-1][j] = (u[n-1][j] - u[n-2][j]) * dxrec;
  }

  for(int i = 0; i < n; i++){
    q[i][0] = (u[i][1] - u[i][0]) * dyrec;
    q[i][n-1] = (u[i][n-1] - u[i][n-2]) * dyrec;
  }

  s[0][0] = (p[0][1] - p[0][0]) * dyrec;
  s[0][n-1] = (p[0][n-1] - p[0][n-2]) * dyrec;
  s[n-1][0] = (p[n-1][1] - p[n-1][0]) * dyrec;
  s[n-1][n-1] = (p[n-1][n-1] - p[n-1][n-2]) * dyrec;

  //Compute derivatives inside domain
 
  //Sweep in x

  for(int i = 1; i < n-1; i++){
    if(i > 1) b[i][i-1] = dx;
    b[i][i] = 4*dx;
    if(i < n-2) b[i][i+1] = dx;
  }
  bprime[1] = b[1][1];
  for(int i = 2; i < n-1; i++){
    bprime[i] = b[i][i] - b[i][i-1]*b[i-1][i]/bprime[i-1];
  }

  //Compute p

  for(int j = 0; j < n; j++){
    for(int i = 1; i < n-1; i++){
      d = 3 * (u[i+1][j] - u[i-1][j]);
      if(i == 1) d -= dx * p[i-1][j];
      if(i == n-2) d -= dx * p[i+1][j];
      dprime[i] = d;
      if(i != 1) dprime[i] -= b[i][i-1] * dprime[i-1] / bprime[i-1];
    }
    
    p[n-2][j] = dprime[n-2] / bprime[n-2];
    for(int i = n-3; i > 0; i--){
      p[i][j] = (dprime[i] - b[i][i+1]*p[i+1][j]) / bprime[i];
    }
  }

  //Compute s

  for(int j = 0; j < n; j += n-1){
    for(int i = 1; i < n-1; i++){
      d = 3 * (q[i+1][j] - q[i-1][j]);
      if(i == 1) d -= dx * s[i-1][j];
      if(i == n-2) d -= dx * s[i+1][j];
      dprime[i] = d;
      if(i != 1) dprime[i] -= b[i][i-1] * dprime[i-1] / bprime[i-1];
    }
    
    s[n-2][j] = dprime[n-2] / bprime[n-2];
    for(int i = n-3; i > 0; i--){
      s[i][j] = (dprime[i] - b[i][i+1]*s[i+1][j]) / bprime[i];
    }
  }
 
  //Sweep in y

  for(int i = 1; i < n-1; i++){
    if(i > 1) b[i][i-1] = dy;
    b[i][i] = 4*dy;
    if(i < n-2) b[i][i+1] = dy;
  }
  bprime[1] = b[1][1];
  for(int i = 2; i < n-1; i++){
    bprime[i] = b[i][i] - b[i][i-1]*b[i-1][i]/bprime[i-1];
  }

  //Compute q

  for(int i = 0; i < n; i++){
    for(int j = 1; j < n-1; j++){
      d = 3 * (u[i][j+1] - u[i][j-1]);
      if(j == 1) d -= dy * q[i][j-1];
      if(j == n-2) d -= dy * q[i][j+1];
      dprime[j] = d;
      if(j != 1) dprime[j] -= b[j][j-1] * dprime[j-1] / bprime[j-1];
    }
    
    q[i][n-2] = dprime[n-2] / bprime[n-2];
    for(int j = n-3; j > 0; j--){
      q[i][j] = (dprime[j] - b[j][j+1]*q[i][j+1]) / bprime[j];
    }
  }

  //Compute s

  for(int i = 0; i < n; i++){
    for(int j = 1; j < n-1; j++){
      d = 3 * (p[i][j+1] - p[i][j-1]);
      if(j == 1) d -= dy * s[i][j-1];
      if(j == n-2) d -= dy * s[i][j+1];
      dprime[j] = d;
      if(j != 1) dprime[j] -= b[j][j-1] * dprime[j-1] / bprime[j-1];
    }
    
    s[i][n-2] = dprime[n-2] / bprime[n-2];
    for(int j = n-3; j > 0; j--){
      s[i][j] = (dprime[j] - b[j][j+1]*s[i][j+1]) / bprime[j];
    }
  }

  //Compute spline coefficients

  for(int i = 1; i < n; i++){
    for(int j = 1; j < n; j++){
      
      for(int l = 0; l < 4; l++){
        for(int m = 0; m < 4; m++){
          for(int n = 0; n < 4; n++){
	    for(int o = 0; o < 4; o++){
              gamma[i][j][l][m] = 0;
	    }
	  }
	}
      }
	    
      double k[4][4] = 
      { 
        {u[i-1][j-1], q[i-1][j-1], u[i-1][j], q[i-1][j]},
        {p[i-1][j-1], s[i-1][j-1], p[i-1][j], s[i-1][j]},
	{u[i][j-1], q[i][j-1], u[i][j], q[i][j]},
        {p[i][j-1], s[i][j-1], p[i][j], s[i][j]}
      };

      for(int l = 0; l < 4; l++){
        for(int m = 0; m < 4; m++){
          for(int n = 0; n < 4; n++){
	    for(int o = 0; o < 4; o++){
              gamma[i][j][l][m] += ax[l][n] * k[n][o] * ay[m][o];
	    }
	  }
	}
      }
    }
  }

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(s);
  memory->destroy(b);

  memory->destroy(bprime);
  memory->destroy(dprime);
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double *data, 
		double *startx, double *dx, int ninput)
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
  double x,xtemp;

  for(int i = 0; i < ninput; i++){
    if(i > 0) xtemp = x;
    if(NULL == fgets(line,MAXLINE,fp)){
      std::string str("Premature end of file in pair table ");
      str += file;
      error->one(FLERR,str.c_str());
    }
    if(2 != sscanf(line,"%lg %lg",&x, &data[i])) ++cerror; 
    if(i == 0) *startx = x;
    if(i > 0){
      if(i == 1) *dx = x - xtemp;
      if((*dx -  x + xtemp) / *dx > SMALL) ++serror;
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror) {
    char str[128];
    sprintf(str,"%d of %d lines in table were incomplete\n"
            "  or could not be parsed completely",cerror,ninput);
    error->warning(FLERR,str);
  }

  // warn if spacing between data points is not constant
  
  if (serror) {
    char str[128];
    sprintf(str, "%d spacings were different\n"
	    "  from first entry",serror);
    error->warning(FLERR,str);
  }

  fclose(fp);
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double **data, 
		double *startx, double *starty, double *dx, double *dy,
		int ninput)
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
  double x,y,xtemp,ytemp;

  for(int i = 0; i < ninput; i++){ 
    if(i > 0) ytemp = y;
    for(int j = 0; j < ninput; j++){
      if(j > 0) xtemp = x;
      if(NULL == fgets(line,MAXLINE,fp)){
	std::string str("Premature end of file in pair table ");
	str += file;
        error->one(FLERR,str.c_str());
      }
      //if(3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[j][i])) ++cerror; 
      if(3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[j][i])) ++cerror; 
      if(j == 0) startx[i] = x;
      if(j > 0){
        if(j == 1) dx[i] = x - xtemp;
        if((dx[i] - x + xtemp)/dx[i] > SMALL) ++sxerror;
      }
    }
    if(i == 0) *starty = y;
    if(i > 0){
      if(i == 1) *dy = y - ytemp;
      if((*dy - y + ytemp)/ *dy > SMALL) ++syerror;
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror) {
    char str[128];
    sprintf(str,"%d of %d lines in table were incomplete\n"
            "  or could not be parsed completely",cerror,ninput);
    error->warning(FLERR,str);
  }
  
  // warn if spacing between data points is not constant
  
  if (sxerror) {
    char str[128];
    sprintf(str, "%d spacings in first column were different\n"
	    "  from first block entries",sxerror);
    error->warning(FLERR,str);
  }

  if (syerror) {
    char str[128];
    sprintf(str, "%d spacings in second column were different\n"
	    "  from first entry",syerror);
    error->warning(FLERR,str);
  }
  
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::uinf(double *param)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  if(sin_alphasq < SMALL){
    return (xi2 - xi1) * spline(h,start_uinf,del_uinf,uinf_coeff,uinf_points)
	    * qelectron;
  }
  else{
    double omega = 1.0 / (1.0 - comega*sin_alphasq);
    double a = omega * sin_alpha;
    double zeta1 = xi1 * a;
    double zeta2 = xi2 * a;

    double g = diameter_angstrom + DELTA2;
    double zetamin;
    if(h >= g) zetamin = 0;
    else zetamin = sqrt(g*g - h*h) 
      * s5((h - diameter_angstrom - DELTA1) * DELTAREC);
    double zetamax = sqrt(cutoffsq_angstrom - h*h);
    double diff_zetarec = 1.0 / (zetamax - zetamin);

    double psi1, psi2, phi1, phi2;
    if(zeta1 < 0){ 
      psi1 = -(zeta1 + zetamin) * diff_zetarec;
      phi1 = -spline(psi1,h,startzeta_phi[0],starth_phi,
        delzeta_phi[0],delh_phi,phi_coeff,phi_points);
    }
    else{
      psi1 = (zeta1 - zetamin) * diff_zetarec; 
      phi1 = spline(psi1,h,startzeta_phi[0],starth_phi,
        delzeta_phi[0],delh_phi,phi_coeff,phi_points);
    }
    if(zeta2 < 0){ 
      psi2 = -(zeta2 + zetamin) * diff_zetarec;
      phi2 = -spline(psi2,h,startzeta_phi[0],starth_phi,
        delzeta_phi[0],delh_phi,phi_coeff,phi_points);
    }
    else{
      psi2 = (zeta2 - zetamin) * diff_zetarec; 
      phi2 = spline(psi2,h,startzeta_phi[0],starth_phi,
        delzeta_phi[0],delh_phi,phi_coeff,phi_points);
    }
    
    double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
    double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;

    return gamma * (phi2 - phi1) * qelectron / a;
  }
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::usemi(double *param)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;
  double etae = param[6] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  double cos_alpha = cos(alpha);
  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double theta = 1.0 - ctheta*sin_alphasq; 

  double a1 = omega * sin_alpha;
  double a2 = theta * etae;

  int points = 100;
  double delxi = (xi2 - xi1) / (points - 1);
  
  double sum = 0;

  for(int i = 0; i < points; i++){
    double xibar = xi1 + i*delxi;
    double g = xibar * a1;
    double hbar = sqrt(h*h + g*g);
    double zetabar = xibar*cos_alpha - a2;
    
    double c = 1.0;
    if(i == 0 || i == points-1) c = 0.5;

    sum += c * spline(hbar,zetabar,starth_usemi,startxi_usemi,
		  delh_usemi,delxi_usemi,usemi_coeff,usemi_points);
  }

  double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;

  return delxi * gamma * sum * qelectron;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::finf(double *param, double *r1, double *r2, 
		double *p1, double *p2, double *m, double **f)
{
  using namespace MathExtra;
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;
  double eta1 = param[4];
  double eta2 = param[5];

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  
  if(sin_alphasq < HIGH){
    f[0][0] = 0.5 * (xi2 - xi1) * dspline(h,start_uinf,del_uinf,
        uinf_coeff,uinf_points) * forceunit;
    f[1][0] = f[0][0];
    f[0][1] = 0;
    f[1][1] = 0;
    f[0][2] = spline(h,start_uinf,del_uinf,uinf_coeff,uinf_points)
	    * forceunit;
    f[1][2] = -f[0][2];
  }

  if(sin_alphasq > LOW){ 
  double sin_alpharec = 1.0 / sin_alpha;
  double sin_alpharecsq = sin_alpharec * sin_alpharec;
  double cos_alpha = cos(alpha);
  double cot_alpha = cos_alpha * sin_alpharec;

  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double a1 = omega * sin_alpha;
  double a1rec = 1.0 / a1;
  double domega = 2 * comega * cos_alpha * a1 * omega;
    
  double gamma_orth = spline(h,start_gamma,del_gamma,
		  gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;
  double gammarec = 1.0 / gamma;
  double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
  double dh_gamma = dspline(h,start_gamma,del_gamma,
		  gamma_coeff,gamma_points) * sin_alphasq;

  double zeta1 = xi1 * a1;
  double zeta2 = xi2 * a1;
     
  double g = diameter_angstrom + DELTA2;
  double zetamin;
  double zetabarmin, s5arg, s5value;
  if(h >= g) zetamin = 0;
  else{ 
    zetabarmin = sqrt(g*g - h*h);
    s5arg = (h - diameter_angstrom - DELTA1) * DELTAREC;
    s5value = s5(s5arg);
    zetamin = zetabarmin * s5value;
  }
  double zetamax = sqrt(cutoffsq_angstrom - h*h);
  double dzetamin;
  if(h >= g) dzetamin = 0;
  else dzetamin = -h / zetabarmin * s5value 
    + zetabarmin * DELTAREC * ds5(s5arg);
  double dzetamax = -h / zetamax;
  double diff_zetarec = 1.0 / (zetamax - zetamin);
  double diff_zeta1 = fabs(zeta1) - zetamin;
  double diff_zeta2 = fabs(zeta2) - zetamin;

  double psi1 = diff_zeta1 * diff_zetarec;
  double psi2 = diff_zeta2 * diff_zetarec;
  
  double phi1 = spline(psi1,h,startzeta_phi[0],starth_phi,
   	delzeta_phi[0],delh_phi,phi_coeff,phi_points);
  double dpsi_phibar1 = dxspline(psi1,h,startzeta_phi[0],starth_phi,
	delzeta_phi[0],delh_phi,phi_coeff,phi_points);
  double dh_phibar1 = dyspline(psi1,h,startzeta_phi[0],starth_phi,
  	delzeta_phi[0],delh_phi,phi_coeff,phi_points);

  double phi2 = spline(psi2,h,startzeta_phi[0],starth_phi,
   	delzeta_phi[0],delh_phi,phi_coeff,phi_points);
  double dpsi_phibar2 = dxspline(psi2,h,startzeta_phi[0],starth_phi,
	delzeta_phi[0],delh_phi,phi_coeff,phi_points);
  double dh_phibar2 = dyspline(psi2,h,startzeta_phi[0],starth_phi,
  	delzeta_phi[0],delh_phi,phi_coeff,phi_points);

  double dzeta_phi1 = dpsi_phibar1 * diff_zetarec;
  double dzeta_phi2 = dpsi_phibar2 * diff_zetarec;

  double b1 = -diff_zetarec * dzetamin;
  double b2 = -diff_zetarec * (dzetamax - dzetamin);
  double dh_psi1 = b1 + b2*psi1;
  double dh_psi2 = b1 + b2*psi2;

  double dh_phi1 = dh_phibar1 + dpsi_phibar1*dh_psi1;
  double dh_phi2 = dh_phibar2 + dpsi_phibar2*dh_psi2;
  if(zeta1 < 0){ 
    phi1 *= -1;
    dh_phi1 *= -1;
  }
  if(zeta2 < 0){ 
    phi2 *= -1;
    dh_phi2 *= -1;
  }

  double diff_dzeta_phi = dzeta_phi2 - dzeta_phi1;
    
  double a2 = gamma * a1rec;
  double u = a2 * (phi2 - phi1);
  double a3 = u * gammarec;
  
  double dh_u = dh_gamma * a3 + a2 * (dh_phi2 - dh_phi1);
  double dalpha_u = dalpha_gamma * a3
	+ a1rec * (domega*sin_alpha + omega*cos_alpha)
	* (gamma*(xi2*dzeta_phi2 - xi1*dzeta_phi1) - u);

  double lrec = 1.0 / (xi2 - xi1);
  double cx = h * gamma * sin_alpharecsq * diff_dzeta_phi;
  double cy = gamma * cot_alpha * diff_dzeta_phi;

  if(sin_alphasq < HIGH && HIGH != LOW){
    double weight = s((sin_alphasq - LOW) / (HIGH - LOW));
    f[0][0] = weight * f[0][0] 
	    + (1-weight) * lrec * (xi2*dh_u - cx) * forceunit;
    f[1][0] = weight * f[1][0] 
	    + (1-weight) * lrec * (-xi1*dh_u + cx) * forceunit;
    f[0][1] = (1-weight) * lrec * (dalpha_u - xi2*cy) * forceunit;
    f[1][1] = (1-weight) * lrec * (-dalpha_u + xi1*cy) * forceunit;
    f[0][2] = weight * f[0][2] 
	    + (1-weight) * gamma * dzeta_phi1 * forceunit;
    f[1][2] = weight * f[1][2] 
	    - (1-weight) * gamma * dzeta_phi2 * forceunit;
  }
  else{
    f[0][0] = lrec * (xi2*dh_u - cx) * forceunit;
    f[1][0] = lrec * (-xi1*dh_u + cx) * forceunit;
    f[0][1] = lrec * (dalpha_u - xi2*cy) * forceunit;
    f[1][1] = lrec * (-dalpha_u + xi1*cy) * forceunit;
    f[0][2] = gamma * dzeta_phi1 * forceunit;
    f[1][2] = -gamma * dzeta_phi2 * forceunit;
  }
  }

  double p[3], delr1p[3], delr2p[3];
  double ftotal[3], fmoment[3], m1[3], m2[3], mtotal[3];
  add3(f[0],f[1],ftotal);
  scale3(-0.5,ftotal);

  add3(p1,p2,p);
  scale3(0.5,p);
  sub3(r1,p,delr1p);
  sub3(r2,p,delr2p);
  cross3(delr1p,f[0],m1);
  cross3(delr2p,f[1],m2);
  add3(m1,m2,mtotal);

  double letarec = 1.0 / (eta2 - eta1);
  cross3(mtotal,m,fmoment);
  scale3(letarec,fmoment);

  add3(ftotal,fmoment,f[2]);
  sub3(ftotal,fmoment,f[3]);
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::fsemi(double *param, double *r1, double *r2, 
		double *p1, double *p2, double *m, double **f)
{
  using namespace MathExtra;
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;
  double eta1 = param[4];
  double eta2 = param[5];
  double etae = param[6] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  double cos_alpha = cos(alpha);

  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double omegasq = omega * omega;
  double domega = 2 * comega * sin_alpha * cos_alpha * omegasq;

  double theta = 1.0 - ctheta*sin_alphasq; 
  double dtheta = -2 * ctheta * sin_alpha * cos_alpha;

  double a1 = omega * sin_alpha;
  double a1sq = a1 * a1;
  double a2 = theta * etae;

  double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;
  double gammarec = 1.0 / gamma;
  double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
  double dh_gamma = dspline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points) * sin_alphasq;
 
  int points = 100;
  double delxi = (xi2 - xi1) / (points - 1);
  double a3 = delxi * gamma;
  
  double jh = 0;
  double jh1 = 0;
  double jh2 = 0;
  double jxi = 0;
  double jxi1 = 0;
  double ubar = 0;

  for(int i = 0; i < points; i++){
    double xibar = xi1 + i*delxi;
    double g = xibar * a1;
    double hbar = sqrt(h*h + g*g);
    double zetabar = xibar*cos_alpha - a2;

    double c = 1.0;
    if(i == 0 || i == points-1) c = 0.5;

    double u = c * spline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,usemi_points);
    double uh;
    if(hbar < SMALL) uh = 0;
    else uh = c / hbar * dxspline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,usemi_points);
    double uxi = c * dyspline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,usemi_points);

    double uh1 = xibar * uh;
    jh += uh;
    jh1 += uh1;
    jh2 += xibar * uh1;
    jxi += uxi;
    jxi1 += xibar * uxi;
    ubar += u;
  }

  jh *= a3;
  jh1 *= a3;
  jh2 *= a3;
  jxi *= a3;
  jxi1 *= a3;
  ubar *= a3;

  double a4 = gammarec * ubar;
  double dh_ubar = dh_gamma*a4 + h*jh;
  double dalpha_ubar = dalpha_gamma*a4 
    + a1*(domega*sin_alpha + omega*cos_alpha)*jh2
    - sin_alpha * jxi1 - dtheta*etae*jxi;

  double cx = h * (omegasq*jh1 + cos_alpha*ctheta*jxi);
  double cy = sin_alpha * (cos_alpha*omegasq*jh1 + (ctheta-1)*jxi);
  double cz1 = a1sq*jh1 + cos_alpha*jxi;
  double cz2 = a1sq*jh2 + cos_alpha*jxi1;

  double lrec = 1.0 / (xi2 - xi1);
  f[0][0] = lrec * (xi2*dh_ubar - cx) * forceunit;
  f[1][0] = lrec * (cx - xi1*dh_ubar) * forceunit;
  f[0][1] = lrec * (dalpha_ubar - xi2*cy) * forceunit;
  f[1][1] = lrec * (xi1*cy - dalpha_ubar) * forceunit;
  f[0][2] = lrec * (cz2 + ubar - xi2*cz1) * forceunit;
  f[1][2] = lrec * (xi1*cz1 - cz2 - ubar) * forceunit;

  double p[3], delr1p[3], delr2p[3];
  double ftotal[3], fmoment[3], m1[3], m2[3], mtotal[3];
  add3(f[0],f[1],ftotal);
  scale3(-0.5,ftotal);

  add3(p1,p2,p);
  scale3(0.5,p);
  sub3(r1,p,delr1p);
  sub3(r2,p,delr2p);
  cross3(delr1p,f[0],m1);
  cross3(delr2p,f[1],m2);
  add3(m1,m2,mtotal);

  double letarec = 1.0 / (eta2 - eta1);
  cross3(mtotal,m,fmoment);
  scale3(letarec,fmoment);

  add3(ftotal,fmoment,f[2]);
  sub3(ftotal,fmoment,f[3]);

}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::geominf(const double *r1, const double *r2, 
		const double *p1, const double *p2, 
		double *param, double *m, double **basis)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3], l[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3];
  double *ex, *ey, *ez;
  double psi, denom, frac, taur, taup;
  double h, alpha, xi1, xi2, eta1, eta2;

  ex = basis[0];
  ey = basis[1];
  ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);
  
  sub3(p,r,delr);

  sub3(r2,r1,l);
  normalize3(l,l);
  sub3(p2,p1,m);
  normalize3(m,m);

  psi = dot3(l,m);
  if(psi > 1.0) psi = 1.0;
  else if(psi < -1.0) psi = -1.0;
  denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  if(denom < HIGH){
    taur = dot3(delr,l);
    taup = 0;
  }
  else{
    frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);
  
  h = len3(delrbar);
    
  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1/h,ex);
  cross3(ez,ex,ey);

  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  sub3(p1,pbar,delp1);
  sub3(p2,pbar,delp2);
  xi1 = dot3(delr1,l);
  xi2 = dot3(delr2,l);
  eta1 = dot3(delp1,m);
  eta2 = dot3(delp2,m);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[4] = eta1;
  param[5] = eta2;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::geomsemi(const double *r1, const double *r2, 
		const double *pin1, const double *pin2, const double *qe,
		double *param, double *m, double **basis)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3], l[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3], delpqe[3];
  const double *p1, *p2;
  double *ex, *ey, *ez;
  double dist1, dist2;
  double psi, denom, frac, taur, taup, rhoe;
  double h, alpha, xi1, xi2, eta1, eta2, etae;

  dist1 = distsq3(pin1,qe);
  dist2 = distsq3(pin2,qe);
  if(dist2 > dist1){
    p1 = pin1;
    p2 = pin2;
  }
  else{
    p1 = pin2;
    p2 = pin1;
  }
  
  ex = basis[0];
  ey = basis[1];
  ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);
  
  sub3(p,r,delr);

  sub3(r2,r1,l);
  normalize3(l,l);
  sub3(p2,p1,m);
  normalize3(m,m);

  psi = dot3(l,m);
  if(psi > 1.0) psi = 1.0;
  else if(psi < -1.0) psi = -1.0;
  denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  sub3(p,qe,delpqe);
  rhoe = dot3(delpqe,m); 

  if(denom < SMALL){
    taur = dot3(delr,l) - rhoe*psi;
    taup = -rhoe;
    etae = 0;
  }
  else{
    frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
    etae = -rhoe - taup;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);
  
  h = len3(delrbar);
    
  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1/h,ex);
  cross3(ez,ex,ey);

  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  sub3(p1,pbar,delp1);
  sub3(p2,pbar,delp2);
  xi1 = dot3(delr1,l);
  xi2 = dot3(delr2,l);
  eta1 = dot3(delp1,m);
  eta2 = dot3(delp2,m);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[4] = eta1;
  param[5] = eta2;
  param[6] = etae;
}


/* ---------------------------------------------------------------------- */

void PairMesoCNT::weight(const double *r1, const double *r2, 
		const double *p1, const double *p2, double *output)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3];
  double rho, rhoctemp, rhoc, rhomin;
  
  add3(r1,r2,r);
  add3(p1,p2,p);
  scale3(0.5,r);
  scale3(0.5,p);
  
  rhoctemp = sqrt(0.25*distsq3(r1,r2) + radiussq);
  rhoc = rhoctemp + sqrt(0.25*distsq3(p1,p2) + radiussq) + rc;
  rhomin = 2.0e-9;
  rho = sqrt(distsq3(r,p));
 
  if(rho > rhoc){ 
    output[0] = 0;
    output[1] = 0;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 0;
    output[6] = 0;
  }
  else if(rho < rhomin){
    output[0] = 1;
    output[1] = 0;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 0;
    output[6] = 0;
  }
  else{
    double denom = 1.0 / (rhoc - rhomin);
    double rhodiff = rho - rhomin;
    double rhobar = rhodiff * denom;
    double a1 = ds(rhobar) * denom;
    double a2 = 0.5 / rho;
    double a3 = a1 * a2;
    double a4 = 0.25 * rhobar / rhoctemp;
    double a5 = a1 * a4;

    sub3(r,p,delr);

    output[0] = s(rhobar);
    double temp1 = a3 * (r[0] - p[0]);
    double temp2 = a3 * (r[1] - p[1]);
    double temp3 = a3 * (r[2] - p[2]);
    double temp4 = a5 * (r1[0] - r2[0]);
    double temp5 = a5 * (r1[1] - r2[1]);
    double temp6 = a5 * (r1[2] - r2[2]);
    output[1] = temp1 - temp4;
    output[2] = temp2 - temp5;
    output[3] = temp3 - temp6;
    output[4] = temp1 + temp4;
    output[5] = temp2 + temp5;
    output[6] = temp3 + temp6;
  }
}

/* ---------------------------------------------------------------------- */

int PairMesoCNT::heaviside(double x)
{
  if(x < 0) return 0;
  else return 1;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::s(double x)
{
  return heaviside(-x) + heaviside(x)*heaviside(1-x)*(1 - x*x*(3 - 2*x));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::ds(double x)
{
  return 6 * heaviside(x) * heaviside(1-x) * x * (x-1);
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::s5(double x)
{
  double x2 = x*x;
  return heaviside(-x) 
	  + heaviside(x)*heaviside(1-x)*(1 - x2*x*(6*x2 - 15*x + 10));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::ds5(double x)
{
  double x2 = x*x;
  return -30 * heaviside(x) * heaviside(1-x) * x2 * (x2 - 2*x + 1);
}


/* ---------------------------------------------------------------------- */

void PairMesoCNT::sort(int *list, int size){
  int i,j,loc1,loc2;
  int *tag = atom->tag;
  for(i = 1; i < size; i++){
    j = i;
    loc1 = list[j-1];
    loc2 = list[j];
    while(j > 0 && tag[loc1] > tag[loc2]){
      list[j] = loc1;
      list[j-1] = loc2;
      j--;
      loc1 = list[j-1];
      loc2 = list[j];
    }
  } 
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::zero(double **a){
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      a[i][j] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::plus(const double * const *a, 
  const double * const *b, double **ans){
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      ans[i][j] = a[i][j] + b[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::scalar(double a, double **m){
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      m[i][j] *= a;
    }
  }
}


/* ---------------------------------------------------------------------- */

void PairMesoCNT::outer(const double *a, const double *b, double **c){
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      c[i][j] = a[i] * b[j];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::matrix_vector(
  const double * const *m, const double *v, double *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}


/* ---------------------------------------------------------------------- */

void PairMesoCNT::trans_matrix_vector(
  const double * const *m, const double *v, double *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}


