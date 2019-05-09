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
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024

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
    memory->destroy(gamma_data);
    memory->destroy(u_inf_data);
    memory->destroy(delh_u_semi);
    memory->destroy(delzeta_phi);

    memory->destroy(u_semi_data);
    memory->destroy(phi_data);
    memory->destroy(gamma_coeff);
    memory->destroy(u_inf_coeff);
    
    memory->destroy(u_semi_coeff);
    memory->destroy(phi_coeff);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{

}

/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  
  memory->create(gamma_data,gamma_points,"pair:gamma_data");
  memory->create(u_inf_data,pot_points,"pair:u_inf_data");
  memory->create(delh_u_semi,pot_points,"pair:delh_u_semi");
  memory->create(delzeta_phi,pot_points,"pair:delzeta_phi");
  
  memory->create(u_semi_data,pot_points,pot_points,"pair:u_semi_data");
  memory->create(phi_data,pot_points,pot_points,"pair:phi_data");
  memory->create(gamma_coeff,gamma_points-1,4,"pair:gamma_coeff");
  memory->create(u_inf_coeff,pot_points-1,4,"pair:u_inf_coeff");
  
  memory->create(u_semi_coeff,pot_points,pot_points-1,4,"pair:u_semi_coeff");
  memory->create(phi_coeff,pot_points,pot_points-1,4,"pair:phi_coeff");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  gamma_points = force->inumeric(FLERR,arg[0]);
  pot_points = force->inumeric(FLERR,arg[1]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  sigma = force->numeric(FLERR,arg[2]);
  epsilon = force->numeric(FLERR,arg[3]);
  n_sigma = force->numeric(FLERR,arg[4]);

  gamma_file = arg[5];
  u_inf_file = arg[6];
  u_semi_file = arg[7];
  phi_file = arg[8];

  //Parse and bcast data
  int me;
  MPI_Comm_rank(world,&me);
  if (me == 0) { 
    read_file(gamma_file,gamma_data,&del_gamma,gamma_points);
    read_file(u_inf_file,u_inf_data,&del_u_inf,pot_points);
    read_file(u_semi_file,u_semi_data,delh_u_semi,&delxi_u_semi,pot_points);
    read_file(phi_file,phi_data,delzeta_phi,delh_phi,pot_points);
  }
  
  MPI_Bcast(gamma_data,gamma_points,MPI_DOUBLE,0,world);
  MPI_Bcast(u_inf_data,pot_points,MPI_DOUBLE,0,world); 
  MPI_Bcast(delh_u_semi,pot_points,MPI_DOUBLE,0,world);
  MPI_Bcast(delzeta_phi,pot_points,MPI_DOUBLE,0,world);
  for(int i = 0; i < pot_points; i++){
    MPI_Bcast(u_semi_data[i],pot_points,MPI_DOUBLE,0,world);
    MPI_Bcast(phi_data[i],pot_points,MPI_DOUBLE,0,world);
  }

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double *data, double *dx, int ninput)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    std::string str("Cannot open file ");
    str += file;
    error->one(FLERR,str.c_str());
  }

  // read values from file
  
  int cerror = 0;
  int serror = 0;
  double x,xtemp;

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for(int i = 0; i < ninput; i++){
    if(i > 0) xtemp = x;
    if(NULL == fgets(line,MAXLINE,fp))
      error->one(FLERR,"Premature end of file in pair table");
    if(2 != sscanf(line,"%lg %lg",&x, &data[i])) ++cerror; 
    if(i > 0){
      if(i == 1) *dx = x - xtemp;
      if(*dx != x - xtemp) ++serror;
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

}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double **data, 
		double *dx, double *dy, int ninput)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
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

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for(int i = 0; i < ninput; i++){
    if(i > 0) ytemp = y;
    for(int j = 0; j < ninput; j++){
      if(j > 0) xtemp = x;
      if(NULL == fgets(line,MAXLINE,fp))
        error->one(FLERR,"Premature end of file in pair table");
      if(3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[j][i])) ++cerror; 
      if(j > 0){
        if(j == 1) dx[i] = x - xtemp;
        if(dx[i] != x - xtemp) ++sxerror;
      }
    }
    if(i > 0){
      if(i == 1) *dy = y - ytemp;
      if(*dy != y - ytemp) ++syerror;
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

}
