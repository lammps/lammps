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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "fix_report_balance.h"
#include "timer.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "min.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "output.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixReportBalance::FixReportBalance(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
/*  if (narg < 5) error->all(FLERR,"Illegal fix adapt command");
  nevery = atoi(arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");

  // count # of adaptations

  nadapt = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 3;
    } else break;
  }

  if (nadapt == 0) error->all(FLERR,"Illegal fix adapt command");
  adapt = new ReportBalance[nadapt];

  // parse keywords

  nadapt = 0;
  diamflag = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      adapt[nadapt].pstyle = new char[n];
      strcpy(adapt[nadapt].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      adapt[nadapt].pparam = new char[n];
      strcpy(adapt[nadapt].pparam,arg[iarg+2]);
      force->bounds(arg[iarg+3],atom->ntypes,
		    adapt[nadapt].ilo,adapt[nadapt].ihi);
      force->bounds(arg[iarg+4],atom->ntypes,
		    adapt[nadapt].jlo,adapt[nadapt].jhi);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
	n = strlen(&arg[iarg+5][2]) + 1;
	adapt[nadapt].var = new char[n];
	strcpy(adapt[nadapt].var,&arg[iarg+5][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = KSPACE;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	int n = strlen(&arg[iarg+1][2]) + 1;
	adapt[nadapt].var = new char[n];
	strcpy(adapt[nadapt].var,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = ATOM;
      if (strcmp(arg[iarg+1],"diameter") == 0) {
	adapt[nadapt].aparam = DIAMETER;
	diamflag = 1;
      } else error->all(FLERR,"Illegal fix adapt command");
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
	int n = strlen(&arg[iarg+2][2]) + 1;
	adapt[nadapt].var = new char[n];
	strcpy(adapt[nadapt].var,&arg[iarg+2][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 3;
    } else break;
  }

  // optional keywords

  resetflag = 0;
  scaleflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"reset") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      if (strcmp(arg[iarg+1],"no") == 0) resetflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) resetflag = 1;
      else error->all(FLERR,"Illegal fix adapt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      if (strcmp(arg[iarg+1],"no") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix adapt command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adapt command");
  }

  // allocate pair style arrays

  int n = atom->ntypes;
  for (int m = 0; m < nadapt; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt:array_orig");
*/
}

/* ---------------------------------------------------------------------- */

FixReportBalance::~FixReportBalance()
{
}

/* ---------------------------------------------------------------------- */

int FixReportBalance::setmask()
{
  int mask = 0;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReportBalance::post_run()
{
  double time_other;
  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);


  time_other = timer->array[TIME_LOOP] -
    (timer->array[TIME_PAIR] + timer->array[TIME_BOND] + 
     timer->array[TIME_KSPACE] + timer->array[TIME_NEIGHBOR] +
     timer->array[TIME_COMM] + timer->array[TIME_OUTPUT]);

  if (me == 0 && screen) {
    fprintf(screen,"---------------------------------------------\n");
    fprintf(screen,"             Load Balance Report\n");
    fprintf(screen,"---------------------------------------------\n");
    fprintf(screen,"timer min avg max std\n");
    fprintf(screen,"---------------------------------------------\n");
  }  

  if (me == 0 && logfile) {
    fprintf(logfile,"---------------------------------------------\n");
    fprintf(logfile,"             Load Balance Report\n");
    fprintf(logfile,"---------------------------------------------\n");
    fprintf(logfile,"timer min avg max std\n");
    fprintf(logfile,"---------------------------------------------\n");
  }
  
  double **timer_mat = memory->create(timer_mat,6,nprocs,"report/balance");
  report_time(timer->array[TIME_PAIR],me,nprocs,"pair   ",timer_mat[0]);
  report_time(timer->array[TIME_NEIGHBOR],me,nprocs,"neigh  ",timer_mat[1]);
  report_time(timer->array[TIME_COMM],me,nprocs,"comm   ",timer_mat[2]);
  report_time(timer->array[TIME_BOND],me,nprocs,"bond   ",timer_mat[3]);
  report_time(timer->array[TIME_KSPACE],me,nprocs,"kspace ",timer_mat[4]);
  report_time(time_other,me,nprocs,"other  ",timer_mat[5]);
  memory->destroy(timer_mat);
  
  echo "plot 'mat.txt' using 1 title 'pair' with lines, 'mat.txt' using 2 title 'comm' with lines" | gnuplot -persist
  
  if (me == 0 && screen)
    fprintf(screen,"---------------------------------------------\n");
  if (me == 0 && logfile)
    fprintf(logfile,"---------------------------------------------\n");
}

/* ---------------------------------------------------------------------- */

void FixReportBalance::report_time(double ptime, int me, int nprocs,
                                   const char *name, double *ptimes)
{
  double min_time, max_time, avg_time, std_time;
  
  MPI_Gather(&ptime,1,MPI_DOUBLE,ptimes,1,MPI_DOUBLE,0,world);

  if (me == 0) {
    min_time = max_time = avg_time = ptimes[0];
    for (int i = 1; i < nprocs; i++) {
      min_time = MIN(min_time,ptimes[i]);
      max_time = MAX(max_time,ptimes[i]);
      avg_time += ptimes[i];
    }
    avg_time /= nprocs;
  
    std_time = 0.0;
    for (int i = 0; i < nprocs; i++) {
      double st = ptimes[i] - avg_time;
      std_time += st * st;
    }
    std_time = sqrt(std_time/nprocs);
    if (screen)
      fprintf(screen,"%s %.3f %.3f %.3f %.3f\n",name,min_time,avg_time,
                                                max_time,std_time);
    if (logfile)
      fprintf(logfile,"%s %.3f %.3f %.3f %.3f\n",name,min_time,avg_time,
                                                 max_time,std_time);
  }
}
