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
   Contributing author: Mike Brown (ORNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "fix_report_balance.h"
#include "timer.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "min.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "output.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReportBalance::FixReportBalance(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  writefile = 0;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (narg == 5) {
    if (strcmp(arg[3],"file") != 0) 
      error->all(FLERR,"Illegal fix report/balance command");
    writefile = 1;
    strncpy(outfile,arg[4],512);
  } else if (narg != 3)      
    error->all(FLERR,"Illegal fix report/balance command");
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
  report_time(timer->array[TIME_PAIR],"pair   ",timer_mat[0]);
  report_time(timer->array[TIME_NEIGHBOR],"neigh  ",timer_mat[1]);
  report_time(timer->array[TIME_COMM],"comm   ",timer_mat[2]);
  report_time(timer->array[TIME_BOND],"bond   ",timer_mat[3]);
  report_time(timer->array[TIME_KSPACE],"kspace ",timer_mat[4]);
  report_time(time_other,"other  ",timer_mat[5]);
  
  if (writefile) {
    if (me == 0) {
      FILE *fp = fopen(outfile,"w");
      if (fp == NULL) {
        char str[612];
        sprintf(str,"Cannot open fix report/balance file %s",outfile);
        error->one(FLERR,str);
      }
      output_per_process(fp,timer_mat);
      fclose(fp);
    } else
      output_per_process(NULL,timer_mat);
  }
      
  memory->destroy(timer_mat);
  
  if (me == 0 && screen)
    fprintf(screen,"---------------------------------------------\n");
  if (me == 0 && logfile)
    fprintf(logfile,"---------------------------------------------\n");
}

/* ---------------------------------------------------------------------- */

void FixReportBalance::report_time(double ptime, const char *name, 
                                   double *ptimes)
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

/* ---------------------------------------------------------------------- */

void FixReportBalance::output_per_process(FILE *fp, double **mat)
{
  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char node_names[MPI_MAX_PROCESSOR_NAME*nprocs];
  MPI_Get_processor_name(node_name,&name_length);
  MPI_Gather(&node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,&node_names,
             MPI_MAX_PROCESSOR_NAME,MPI_CHAR,0,world);

  if (fp) {
    fprintf(fp,"# plot '%s' using 1:3 title 'pair' with lines, ",outfile);
    fprintf(fp,"'%s' using 1:4 title 'neigh' with lines, ",outfile);
    fprintf(fp,"'%s' using 1:5 title 'comm' with lines, ",outfile);
    fprintf(fp,"'%s' using 1:6 title 'bond' with lines, ",outfile);
    fprintf(fp,"'%s' using 1:7 title 'kspace' with lines, ",outfile);
    fprintf(fp,"'%s' using 1:8 title 'other' with lines\n",outfile);
    fprintf(fp,"#\n#\n");
    fprintf(fp,"# rank node pair neigh comm bond kspace other\n");
    for (int i = 0; i < nprocs; i++)
      fprintf(fp,"%d %s %.3f %.3f %.3f %.3f %.3f %.3f\n",i,
              &node_names[i*MPI_MAX_PROCESSOR_NAME],mat[0][i],mat[1][i],
              mat[2][i],mat[3][i],mat[4][i],mat[5][i]);
  }
}

