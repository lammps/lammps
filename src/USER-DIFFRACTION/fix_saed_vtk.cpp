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
   Contributing author: Pieter in 't Veld (SNL)
   Incorporating SAED: Shawn Coleman (Arkansas)
------------------------------------------------------------------------- */

#include "fix_saed_vtk.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_saed.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE};
enum{ONE,RUNNING,WINDOW};
enum{FIRST,MULTI};

#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

FixSAEDVTK::FixSAEDVTK(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), ids(NULL), fp(NULL), vector(NULL),
  vector_total(NULL), vector_list(NULL), compute_saed(NULL), filename(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix saed/vtk command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;

  nvalues = 0;
  int iarg = 6;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0){
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues != 1) error->all(FLERR,"Illegal fix saed/vtk command");

  options(narg,arg);

  which = 0;
  ids = NULL;

  nvalues = 0;

  iarg = 6;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ) {

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) error->all(FLERR,"Illegal fix saed/vtk command");

      n = strlen(suffix) + 1;
      ids = new char[n];
      strcpy(ids,suffix);
      delete [] suffix;

      int icompute = modify->find_compute(ids);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix saed/vtk does not exist");

      Compute *compute = modify->compute[icompute];

      // Check that specified compute is for SAED
      compute_saed = (ComputeSAED*) modify->compute[icompute];
      if (strcmp(compute_saed->style,"saed") != 0)
        error->all(FLERR,"Fix saed/vtk has invalid compute assigned");

      // Gather variables from specified compute_saed
      double *saed_var = compute_saed->saed_var;
      lambda   = saed_var[0];
      Kmax     = saed_var[1];
      Zone[0]  = saed_var[2];
      Zone[1]  = saed_var[3];
      Zone[2]  = saed_var[4];
      c[0]     = saed_var[5];
      c[1]     = saed_var[6];
      c[2]     = saed_var[7];
      dR_Ewald = saed_var[8];
      double manual_double = saed_var[9];
      manual = false;
      if (manual_double == 1) manual = true;

      // Standard error check for fix/ave/time
      if (compute->vector_flag == 0)
        error->all(FLERR,"Fix saed/vtk compute does not calculate a vector");
      if (compute->extvector != 0)
        error->all(FLERR,"Illegal fix saed/vtk command");

      nrows = compute->size_vector;
      nvalues++;
      iarg++;
    } else break;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix saed/vtk command");
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Illegal fix saed/vtk command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix saed/vtk command");

  // allocate memory for averaging

  vector = vector_total = NULL;
  vector_list = NULL;

  if (ave == WINDOW)
    memory->create(vector_list,nwindow,nvalues,"saed/vtk:vector_list");

  memory->create(vector,nrows,"saed/vtk:vector");
  memory->create(vector_total,nrows,"saed/vtk:vector_total");

  extlist = NULL;

  vector_flag = 1;
  size_vector = nrows;

  if (nOutput == 0) {

    // SAED specific paramaters needed
    int *periodicity = domain->periodicity;
    // Zone flag to capture entire recrocal space volume
    if (  (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    } else {
        R_Ewald = (1 / lambda);
        double Rnorm = R_Ewald/ sqrt(Zone[0] * Zone[0] +
                       Zone[1] * Zone[1] +  Zone[2]* Zone[2]);
        Zone[0] = Zone[0] * Rnorm;
        Zone[1] = Zone[1] * Rnorm;
        Zone[2] = Zone[2] * Rnorm;
    }

    double *prd;
    double ave_inv = 0.0;
    prd = domain->prd;

    if (periodicity[0]){
      prd_inv[0] = 1 / prd[0];
      ave_inv += prd_inv[0];
    }
    if (periodicity[1]){
      prd_inv[1] = 1 / prd[1];
      ave_inv += prd_inv[1];
    }
    if (periodicity[2]){
      prd_inv[2] = 1 / prd[2];
      ave_inv += prd_inv[2];
    }

    // Using the average inverse dimensions for non-periodic direction
    ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
    if (!periodicity[0]){
      prd_inv[0] = ave_inv;
    }
    if (!periodicity[1]){
      prd_inv[1] = ave_inv;
    }
    if (!periodicity[2]){
      prd_inv[2] = ave_inv;
    }

    // Use manual mapping of reciprocal lattice
    if (manual) {
      for (int i=0; i<3; i++) {
        prd_inv[i] = 1.0;
      }
    }

    // Find integer dimensions of the reciprocal lattice box bounds
    if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
      for (int i=0; i<3; i++) {
        dK[i] = prd_inv[i]*c[i];
        Knmax[i] = ceil(Kmax / dK[i]);
        Knmin[i] = -Knmax[i];
      }
    } else {

      for (int i=0; i<3; i++) {
        Knmax[i] = -10000;
        Knmin[i] =  10000;
      }
      double dinv2 = 0.0;
      double r = 0.0;
      double K[3];
      int Ksearch[3];
      for (int i=0; i<3; i++) {
        dK[i] = prd_inv[i]*c[i];
        Ksearch[i] = ceil(Kmax / dK[i]);
      }

      for (int k = -Ksearch[2]; k <= Ksearch[2]; k++) {
        for (int j = -Ksearch[1]; j <= Ksearch[1]; j++) {
          for (int i = -Ksearch[0]; i <= Ksearch[0]; i++) {
            K[0] = i * dK[0];
            K[1] = j * dK[1];
            K[2] = k * dK[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < Kmax * Kmax) {

              r=0.0;
              for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
              r = sqrt(r);
              if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){

                if ( i < Knmin[0] ) Knmin[0] = i;
                if ( j < Knmin[1] ) Knmin[1] = j;
                if ( k < Knmin[2] ) Knmin[2] = k;

                if ( i > Knmax[0] ) Knmax[0] = i;
                if ( j > Knmax[1] ) Knmax[1] = j;
                if ( k > Knmax[2] ) Knmax[2] = k;
              }
            }
          }
        }
      }
    }

   // Finding dimensions for vtk files
    for (int i=0; i<3; i++) {
      if ( ( (Knmin[i] > 0) && (Knmax[i] > 0) ) || ( (Knmin[i] < 0) && (Knmax[i] < 0) ) ){
        Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] );
      } else Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] ) + 1;
    }
  }

  // initialization

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  for (int i = 0; i < nrows; i++)
     vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

}

/* ---------------------------------------------------------------------- */

FixSAEDVTK::~FixSAEDVTK()
{
  delete [] extlist;
  delete [] filename;
  delete [] ids;
  memory->destroy(vector);
  memory->destroy(vector_total);
  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixSAEDVTK::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSAEDVTK::init()
{
  // set current indices for all computes,fixes,variables


  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix saed/vtk does not exist");

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixSAEDVTK::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixSAEDVTK::end_of_step()
{
  // skip if not step which requires doing something
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixSAEDVTK::invoke_vector(bigint ntimestep)
{
  // zero if first step
  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix saed/vtk does not exist");

  if (irepeat == 0)
    for (int i = 0; i < nrows; i++)
       vector[i] = 0.0;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // invoke compute if not previously invoked

  Compute *compute = modify->compute[icompute];

  if (!(compute->invoked_flag & INVOKED_VECTOR)) {
    compute->compute_vector();
    compute->invoked_flag |= INVOKED_VECTOR;
  }

  double *vector = compute->vector;

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for ( int i = 0; i < nrows; i++)
    vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (int i = 0; i < nrows; i++) vector_total[i] = vector[i];
    norm = 1;

  } else if (ave == RUNNING) {
    for (int i = 0; i < nrows; i++) vector_total[i] += vector[i];
    norm++;

  } else if (ave == WINDOW) {
    for (int i = 0; i < nrows; i++) {
      vector_total[i] += vector[i];
      if (window_limit) vector_total[i] -= vector_list[iwindow][i];
        vector_list[iwindow][i] = vector[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // output result to file

  if (fp && me == 0) {
    if (nOutput > 0) {
      fclose(fp);

      char nName [128];
      snprintf(nName,128,"%s.%d.vtk",filename,nOutput);
      fp = fopen(nName,"w");

      if (fp == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix saed/vtk file %s",nName);
        error->one(FLERR,str);
      }
    }

    fprintf(fp,"# vtk DataFile Version 3.0 c_%s\n",ids);
    fprintf(fp,"Image data set\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS %d %d %d\n", Dim[0],  Dim[1], Dim[2]);
    fprintf(fp,"ASPECT_RATIO %g %g %g\n", dK[0], dK[1], dK[2]);
    fprintf(fp,"ORIGIN %g %g %g\n", Knmin[0] * dK[0],  Knmin[1] * dK[1], Knmin[2] * dK[2]);
    fprintf(fp,"POINT_DATA %d\n",  Dim[0] *  Dim[1] * Dim[2] );
    fprintf(fp,"SCALARS intensity float\n");
    fprintf(fp,"LOOKUP_TABLE default\n");

    filepos = ftell(fp);

    if (overwrite) fseek(fp,filepos,SEEK_SET);

     // Finding the intersection of the reciprical space and Ewald sphere
      int NROW1 = 0;
      int NROW2 = 0;
      double dinv2 = 0.0;
      double r = 0.0;
      double K[3];

      // Zone flag to capture entire recrocal space volume
      if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
        for (int k = Knmin[2]; k <= Knmax[2]; k++) {
          for (int j = Knmin[1]; j <= Knmax[1]; j++) {
            for (int i = Knmin[0]; i <= Knmax[0]; i++) {
              K[0] = i * dK[0];
              K[1] = j * dK[1];
              K[2] = k * dK[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax * Kmax) {
                 fprintf(fp,"%g\n",vector_total[NROW1]/norm);
                 fflush(fp);
                 NROW1++;
                 NROW2++;
              } else {
              fprintf(fp,"%d\n",-1);
              fflush(fp);
              NROW2++;
              }
            }
          }
        }
      } else {
        for (int k = Knmin[2]; k <= Knmax[2]; k++) {
          for (int j = Knmin[1]; j <= Knmax[1]; j++) {
            for (int i = Knmin[0]; i <= Knmax[0]; i++) {
              K[0] = i * dK[0];
              K[1] = j * dK[1];
              K[2] = k * dK[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax * Kmax) {
                r=0.0;
                for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
                r = sqrt(r);
                if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                 fprintf(fp,"%g\n",vector_total[NROW1]/norm);
                 fflush(fp);
                 NROW2++;
                 NROW1++;
                } else {
                  fprintf(fp,"%d\n",-1);
                  fflush(fp);
                  NROW2++;
                }
              } else {
              fprintf(fp,"%d\n",-1);
              fflush(fp);
              NROW2++;
             }
            }
          }
        }
      }
    }
  nOutput++;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixSAEDVTK::compute_vector(int i)
{
  if (norm) {
    return vector_total[i]/norm;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixSAEDVTK::options(int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  overwrite = 0;

  // optional args
  int iarg = 6 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix saed/vtk command");
      if (me == 0) {

         nOutput = 0;
         int n = strlen(arg[iarg+1]) + 1;
         filename = new char[n];
         strcpy(filename,arg[iarg+1]);

         char nName [128];
         snprintf(nName,128,"%s.%d.vtk",filename,nOutput);
         fp = fopen(nName,"w");

        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix saed/vtk file %s",nName);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix saed/vtk command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix saed/vtk command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix saed/vtk command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix saed/vtk command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix saed/vtk command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix saed/vtk command");
  }
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixSAEDVTK::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixSAEDVTK::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix saed/vtk missed timestep");
}
