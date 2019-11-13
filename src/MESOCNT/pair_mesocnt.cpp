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
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMesoCNT::~PairMesoCNT()
{
  if (allocated) {

  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{
  
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{

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

  return rc;
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

void PairMesoCNT::read_file(char *file, double *data, 
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
		    "  in pair table ",cerror,size);
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

void PairMesoCNT::read_file(char *file, double **data, 
		double &startx, double &starty, 
		double &dx, double &dy, int size)
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
      if (2 != sscanf(line,"%lg %lg %lg",&x,&y,&data[i][j])) cerror++;
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
		    "  in pair table ",cerror,size);
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
   compute cubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double *data, double **coeff,
		double dx, int size)
{

}

/* ----------------------------------------------------------------------
   compute bicubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ****coeff,
		double dx, double dy, int size)
{

}


/* ----------------------------------------------------------------------
   cubic spline evaluation
------------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  
}

/* ----------------------------------------------------------------------
   cubic spline derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dspline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  
}

/* ----------------------------------------------------------------------
   bicubic spline evaluation
------------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  
}

/* ----------------------------------------------------------------------
   bicubic spline partial x derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dxspline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  
}

/* ----------------------------------------------------------------------
   bicubic spline partial y derivative
------------------------------------------------------------------------- */

double PairMesoCNT::dyspline(double x, double y, 
		double xstart, double ystart, 
		double dx, double dy,
		double ****coeff, int coeff_size)
{
  
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
