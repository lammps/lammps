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
#define SMALL 1.0e-6
#define SWITCH 1.0e-6
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
  
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  int ntypes = atom->ntypes;
  
  memory->create(cutsq,ntypes+1,ntypes+1,"pair:cutsq");
  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;

  memory->create(uinf_coeff,uinf_points,4,"pair:uinf_coeff");
  memory->create(gamma_coeff,gamma_points,4,"pair:gamma_coeff");
  memory->create(phi_coeff,phi_points,phi_points,4,4,"pair:phi_coeff");
  memory->create(usemi_coeff,usemi_points,usemi_points,4,4,"pair:usemi_coeff");

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
  if (narg != 10) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  // CNT constants
  n = force->inumeric(FLERR,arg[2]);
  sig = force->numeric(FLERR,arg[3]);
  eps = force->numeric(FLERR,arg[4]);
  nsig = force->numeric(FLERR,arg[5]);

  // file names
  uinf_file = arg[6];
  gamma_file = arg[7];
  phi_file = arg[8];
  usemi_file = arg[9];

  // units
  ang = force->angstrom;
  angrec = 1.0 / ang;
  e = force->qelectron;
  erec = 1.0 / e;
  funit = e * angrec;

  // potential variables
  r = 1.431 * 3 * n / MY_2PI * ang;
  rsq = r * r;
  d = 2 * r;
  d_ang = d * angrec;
  cutoff = rc + d;
  cutoffsq = cutoff * cutoff;
  cutoff_ang = cutoff * angrec;
  cutoffsq_ang = cutoff_ang * cutoff_ang;
  rc = 3.0 * sig;
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

  printf("Coefficients computed!\n");
  double *w,*dw,**w2,**dxw2,**dyw2;
  double **w_coeff,****w2_coeff;
  int points = 101;
  memory->create(w,points,"pair:w");
  memory->create(dw,points,"pair:dw");
  memory->create(w2,points,points,"pair:w2");
  memory->create(dxw2,points,points,"pair:dxw2");
  memory->create(dyw2,points,points,"pair:dyw2");
  memory->create(w_coeff,points,4,"pair:w_coeff");
  memory->create(w2_coeff,points,points,4,4,"pair:w2_coeff");
  double xstart = -1;
  double xend = 1;
  double ystart = -1;
  double yend = 1;
  double dx = (xend - xstart) / (points - 1);
  double dy = (yend - ystart) / (points - 1);

  std::ofstream wfile("w.dat");
  std::ofstream dwfile("dw.dat");
  std::ofstream w2file("w2.dat");
  std::ofstream dxw2file("dxw2.dat");
  std::ofstream dyw2file("dyw2.dat"); 

  for (int i = 0; i < points; i++) {
    double x = xstart + i*dx;
    w[i] = 1.0 / (1.0 + x*x);
    dw[i] = -2.0*x / pow(1.0 + x*x, 2);
    wfile << x << " " << w[i] << std::endl;
    dwfile << x << " " << dw[i] << std::endl;
    for (int j = 0; j < points; j++) {
      double y = ystart + j*dy;
      w2[i][j] = 1.0 / (1.0 + x*x + 2.0*y*y);
      dxw2[i][j] = -2.0*x / pow(1.0 + x*x + 2.0*y*y, 2);
      dyw2[i][j] = -4.0*y / pow(1.0 + x*x + 2.0*y*y, 2);
      w2file << x << " " << y << " " << w2[i][j] << std::endl;
      dxw2file << x << " " << y << " " << dxw2[i][j] << std::endl;
      dyw2file << x << " " << y << " " << dyw2[i][j] << std::endl;
    }
  }

  wfile.close();
  dwfile.close();
  w2file.close();
  dxw2file.close();
  dyw2file.close();

  spline_coeff(w,w_coeff,dx,points);
  spline_coeff(w2,w2_coeff,dx,dy,points);

  int points2 = 501;
  double xstart2 = -2;
  double xend2 = 2;
  double ystart2 = -2;
  double yend2 = 2;
  double dx2 = (xend2 - xstart2) / (points2 - 1);
  double dy2 = (yend2 - ystart2) / (points2 - 1);
  
  std::ofstream wifile("wi.dat");
  std::ofstream dwifile("dwi.dat");
  std::ofstream w2ifile("w2i.dat");
  std::ofstream dxw2ifile("dxw2i.dat");
  std::ofstream dyw2ifile("dyw2i.dat");

  for (int i = 0; i < points2; i++) {
    double x = xstart2 + i*dx2;
    wifile << x << " " << spline(x,xstart,dx,w_coeff,points) << std::endl;
    dwifile << x << " " << dspline(x,xstart,dx,w_coeff,points) << std::endl;
    for (int j = 0; j < points2; j++) {
      double y = ystart2 + j*dy2;
      w2ifile << x << " " << y << " " << spline(x,y,xstart,ystart,dx,dy,w2_coeff,points) << std::endl;
      dxw2ifile << x << " " << y << " " << dxspline(x,y,xstart,ystart,dx,dy,w2_coeff,points) << std::endl;
      dyw2ifile << x << " " << y << " " << dyspline(x,y,xstart,ystart,dx,dy,w2_coeff,points) << std::endl;
    }
  }

  wifile.close();
  dwifile.close();
  w2ifile.close();
  dxw2ifile.close();
  dyw2ifile.close();

  memory->destroy(w);
  memory->destroy(dw);
  memory->destroy(w2);
  memory->destroy(dxw2);
  memory->destroy(dyw2);
  memory->destroy(w_coeff);
  memory->destroy(w2_coeff);
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
		double &startx, double &dx, int ninput)
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

  for (int i = 0; i < ninput; i++){
    if (NULL == fgets(line,MAXLINE,fp)) {
      std::string str("Premature end of file in pair table ");
      str += file;
      error->one(FLERR,str.c_str());
    }
    if (i > 0) xtemp = x;
    if (2 != sscanf(line,"%lg %lg",&x,&data[i])) cerror++;
    if (i > 0) {
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
		double &startx, double &starty, 
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

  for (int i = 0; i < ninput; i++){
    if (i > 0) xtemp = x;
    for (int j = 0; j < ninput; j++){
      if (NULL == fgets(line,MAXLINE,fp)) {
        std::string str("Premature end of file in pair table ");
        str += file;
        error->one(FLERR,str.c_str());
      }
      if (j > 0) ytemp = y;
      if (3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[i][j])) cerror++;
      if (j > 0) {
	dytemp = y - ytemp;
    	if (j == 1) dy = dytemp;
	if (fabs(dytemp - dy)/dy > SMALL) syerror++;
      }
    }
    if (i == 0) startx = x;
    if (i > 0) {
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
  
}

/* ----------------------------------------------------------------------
   geometric parameters for semi-infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::geomsemi(const double *r1, const double *r2, 
		const double *p1, const double *p2, const double *qe,
		double *param, double **basis)
{
  
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
   potential energy for infinite CNT case
------------------------------------------------------------------------- */

double PairMesoCNT::uinf(const double *param)
{

}

/* ----------------------------------------------------------------------
   potential energy for semi-infinite CNT case
------------------------------------------------------------------------- */

double PairMesoCNT::usemi(const double *param)
{

}

/* ----------------------------------------------------------------------
   forces for infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::finf(const double *param, double **f)
{

}

/* ----------------------------------------------------------------------
   forces for semi-infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::fsemi(const double *param, double **f)
{

}
