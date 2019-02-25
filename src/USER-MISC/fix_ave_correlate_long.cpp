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
     Jorge Ramirez (jorge.ramirez@upm.es, Universidad Politecnica de Madrid),
     Alexei Likhtman (University of Reading)
   Structure and syntax of fix inspired by fix_ave_correlate
   Scalar Correlator f(tau)=<A(t)A(t+tau)> and
     Cross-correlator f(tau)=<A(t)B(t+tau)>
   see J. Chem. Phys. 133, 154103 (2010)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "fix_ave_correlate_long.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE,FIX,VARIABLE};
enum{AUTO,UPPER,LOWER,AUTOUPPER,AUTOLOWER,FULL};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

static const char cite_fix_ave_correlate_long[] =
"fix ave/correlate/long command:\n\n"
"@Article{Ramirez10,\n"
" author = {Jorge Rami{\'}rez and Sathish K. Sukumaran and Bart Vorselaars and Alexei E. Likhtman},\n"
" title =   {Efficient on the fly calculation of time correlation functions in computer simulations},"
" journal = {J.~Chem.~Phys.},\n"
" year =    2010,\n"
" volume =  133,\n"
" pages =   {154103}\n"
"}\n\n";

/* ---------------------------------------------------------------------- */

FixAveCorrelateLong::FixAveCorrelateLong(LAMMPS * lmp, int narg, char **arg):
  Fix (lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ave_correlate_long);

  // At least nevery nfrez and one value are needed
  if (narg < 6) error->all(FLERR,"Illegal fix ave/correlate/long command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nfreq = force->inumeric(FLERR,arg[4]);

  restart_global = 1;
  global_freq = nfreq;

  // parse values until one isn't recognized

  which = new int[narg-5];
  argindex = new int[narg-5];
  ids = new char*[narg-5];
  value2index = new int[narg-5];
  nvalues = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0 ||
        strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/correlate/long command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      delete [] suffix;

      nvalues++;
      iarg++;
    } else break;
  }

  // optional args

  type = AUTO;
  startstep = 0;
  fp = NULL;
  overwrite = 0;
  numcorrelators=20;
  p = 16;
  m = 2;
  char *title1 = NULL;
  char *title2 = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      if (strcmp(arg[iarg+1],"auto") == 0) type = AUTO;
      else if (strcmp(arg[iarg+1],"upper") == 0) type = UPPER;
      else if (strcmp(arg[iarg+1],"lower") == 0) type = LOWER;
      else if (strcmp(arg[iarg+1],"auto/upper") == 0) type = AUTOUPPER;
      else if (strcmp(arg[iarg+1],"auto/lower") == 0) type = AUTOLOWER;
      else if (strcmp(arg[iarg+1],"full") == 0) type = FULL;
      else error->all(FLERR,"Illegal fix ave/correlate/long command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ncorr") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      numcorrelators = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nlen") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      p = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ncount") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      m = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix ave/correlate/long file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/correlate/long command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/correlate/long command");
  }

  if (p % m != 0) error->all(FLERR,"fix_correlator: p mod m must be 0");
  dmin = p/m;
  length = numcorrelators*p;
  npcorr = 0;
  kmax = 0;

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/correlate/long command");
  if (nfreq % nevery)
    error->all(FLERR,"Illegal fix ave/correlate/long command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/correlate/long does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/correlate/long compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/correlate/long compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
        error->all(FLERR,"Fix ave/correlate/long compute vector "
                   "is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/correlate/long does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,"Fix ave/correlate/long fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix ave/correlate/long fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,
                   "Fix ave/correlate/long fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,"Fix for fix ave/correlate/long "
                   "not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/correlate/long does not exist");
      if (input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR,
                   "Fix ave/correlate/long variable is not equal-style variable");
    }
  }

  // npair = # of correlation pairs to calculate
  if (type == AUTO) npair = nvalues;
  if (type == UPPER || type == LOWER) npair = nvalues*(nvalues-1)/2;
  if (type == AUTOUPPER || type == AUTOLOWER) npair = nvalues*(nvalues+1)/2;
  if (type == FULL) npair = nvalues*nvalues;

  // print file comment lines
  if (fp && me == 0) {
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-correlated data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else {
      fprintf(fp,"# Time");
      if (type == AUTO)
        for (int i = 0; i < nvalues; i++)
          fprintf(fp," %s*%s",arg[5+i],arg[5+i]);
      else if (type == UPPER)
        for (int i = 0; i < nvalues; i++)
          for (int j = i+1; j < nvalues; j++)
            fprintf(fp," %s*%s",arg[5+i],arg[5+j]);
      else if (type == LOWER)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < i-1; j++)
            fprintf(fp," %s*%s",arg[5+i],arg[5+j]);
      else if (type == AUTOUPPER)
        for (int i = 0; i < nvalues; i++)
          for (int j = i; j < nvalues; j++)
            fprintf(fp," %s*%s",arg[5+i],arg[5+j]);
      else if (type == AUTOLOWER)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < i; j++)
            fprintf(fp," %s*%s",arg[5+i],arg[5+j]);
      else if (type == FULL)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < nvalues; j++)
            fprintf(fp," %s*%s",arg[5+i],arg[5+j]);
      fprintf(fp,"\n");
    }
    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;

  // allocate and initialize memory for calculated values and correlators

  memory->create(values,nvalues,"correlator:values");
  memory->create(shift,npair,numcorrelators,p,"correlator:shift");
  memory->create(shift2,npair,numcorrelators,p,"correlator:shift2"); //NOT OPTMAL
  memory->create(correlation,npair,numcorrelators,p,"correlator:correlation");
  memory->create(accumulator,npair,numcorrelators,"correlator:accumulator");
  memory->create(accumulator2,npair,numcorrelators,"correlator:accumulator2"); // NOT OPTIMAL

  memory->create(ncorrelation,numcorrelators,p,"correlator:ncorrelation");
  memory->create(naccumulator,numcorrelators,"correlator:naccumulator");
  memory->create(insertindex,numcorrelators,"correlator:insertindex");
  memory->create(t,length,"correlator:t");
  memory->create(f,npair,length,"correlator:f");

  for (int i=0;i<npair;i++)
    for (int j=0;j<numcorrelators;j++) {
      for (unsigned int k=0;k<p;k++) {
        shift[i][j][k]=-2E10;
        shift2[i][j][k]=0.0;
        correlation[i][j][k]=0.0;
      }
      accumulator[i][j]=0.0;
      accumulator2[i][j]=0.0;
    }

  for (int i=0;i<numcorrelators;i++) {
    for (unsigned int j=0;j<p;j++) ncorrelation[i][j]=0;
    naccumulator[i]=0;
    insertindex[i]=0;
  }

  for (int i=0;i<length;i++) t[i]=0.0;
  for (int i=0;i<npair;i++)
    for (int j=0;j<length;j++) f[i][j]=0.0;


  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
  last_accumulated_step = -1;

}

/* ---------------------------------------------------------------------- */

FixAveCorrelateLong::~FixAveCorrelateLong()
{
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  memory->destroy(values);
  memory->destroy(shift);
  memory->destroy(shift2);
  memory->destroy(correlation);
  memory->destroy(accumulator);
  memory->destroy(accumulator2);
  memory->destroy(ncorrelation);
  memory->destroy(naccumulator);
  memory->destroy(insertindex);
  memory->destroy(t);
  memory->destroy(f);

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixAveCorrelateLong::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveCorrelateLong::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/correlate/long does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/correlate/long does not exist");
      value2index[i] = ifix;

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/correlate/long does not exist");
      value2index[i] = ivariable;
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveCorrelateLong::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveCorrelateLong::end_of_step()
{
  int i,m;
  double scalar;

  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/correlate/long");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // accumulate results of computes,fixes,variables to origin
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];
    scalar = 0.0;

    // invoke compute if not previously invoked

    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        scalar = compute->scalar;
      } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        scalar = compute->vector[argindex[i]-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {
      if (argindex[i] == 0)
        scalar = modify->fix[m]->compute_scalar();
      else
        scalar = modify->fix[m]->compute_vector(argindex[i]-1);

    // evaluate equal-style variable

    } else if (which[i] == VARIABLE)
      scalar = input->variable->compute_equal(m);

    values[i] = scalar;
  }

  // fistindex = index in values ring of earliest time sample

  nvalid += nevery;
  modify->addstep_compute(nvalid);

  // calculate all Cij() enabled by latest values

  accumulate();
  if (ntimestep % nfreq) return;

  // output result to file
  evaluate();

  if (fp && me == 0) {
    if(overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,"# Timestep: " BIGINT_FORMAT "\n", ntimestep);
    for (unsigned int i=0;i<npcorr;++i) {
      fprintf(fp, "%lg ", t[i]*update->dt);
      for (int j=0;j<npair;++j) {
        fprintf(fp, "%lg ", f[j][i]);
      }
    fprintf(fp, "\n");
    }
    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) ftruncate(fileno(fp),fileend);
    }
  }

  return;

}

void FixAveCorrelateLong::evaluate() {
  unsigned int jm=0;

  // First correlator
  for (unsigned int j=0;j<p;++j) {
    if (ncorrelation[0][j] > 0) {
      t[jm] = j;
      for (int i=0;i<npair;++i)
        f[i][jm] = correlation[i][0][j]/ncorrelation[0][j];
      ++jm;
    }
  }

  // Subsequent correlators
  for (int k=1;k<kmax;++k) {
    for (unsigned int j=dmin;j<p;++j) {
      if (ncorrelation[k][j]>0) {
        t[jm] = j * pow((double)m, k);
        for (int i=0;i<npair;++i)
          f[i][jm] = correlation[i][k][j] / ncorrelation[k][j];
        ++jm;
      }
    }
  }

  npcorr = jm;
}


/* ----------------------------------------------------------------------
   accumulate correlation data using more recently added values
------------------------------------------------------------------------- */

void FixAveCorrelateLong::accumulate()
{
  int i,j,ipair;

  if (update->ntimestep <= last_accumulated_step) return;

  if (type == AUTO) {
    for (i=0; i<nvalues;i++) add(i,values[i]);
  } else if (type == UPPER) {
    ipair = 0;
    for (i=0;i<nvalues;i++)
      for (j=i+1;j<nvalues;j++) add(ipair++,values[i],values[j]);
  } else if (type == LOWER) {
    ipair = 0;
    for (i=0;i<nvalues;i++)
      for (j=0;j<i;j++) add(ipair++,values[i],values[j]);
  } else if (type == AUTOUPPER) {
    ipair = 0;
    for (i=0;i<nvalues;i++)
      for (j=i;j<nvalues;j++) {
        if (i==j) add(ipair++,values[i]);
        else add(ipair++,values[i],values[j]);
      }
  } else if (type == AUTOLOWER) {
    ipair = 0;
    for (i=0;i<nvalues;i++)
      for (j=0;j<=i;j++) {
        if (i==j) add(ipair++,values[i]);
        else add(ipair++,values[i],values[j]);
      }
  } else if (type == FULL) {
    ipair = 0;
    for (i=0;i<nvalues;i++)
      for (j=0;j<nvalues;j++) {
        if (i==j) add(ipair++,values[i]);
        else add(ipair++,values[i],values[j]);
      }
  }
  last_accumulated_step = update->ntimestep;
}


/* ----------------------------------------------------------------------
   Add a scalar value to the autocorrelator k of pair i
------------------------------------------------------------------------- */
void FixAveCorrelateLong::add(const int i, const double w, const int k){
  // If we exceed the correlator side, the value is discarded
  if (k == numcorrelators) return;
  if (k > kmax) kmax=k;

  // Insert new value in shift array
  shift[i][k][insertindex[k]] = w;

  // Add to accumulator and, if needed, add to next correlator
  accumulator[i][k] += w;
  if (i==0) ++naccumulator[k];
  if (naccumulator[k]==m) {
    add(i,accumulator[i][k]/m, k+1);
    accumulator[i][k]=0;
    if (i==npair-1) naccumulator[k]=0;
  }

  // Calculate correlation function
  unsigned int ind1=insertindex[k];
  if (k==0) { // First correlator is different
    int ind2=ind1;
    for (unsigned int j=0;j<p;++j) {
      if (shift[i][k][ind2] > -1e10) {
        correlation[i][k][j]+= shift[i][k][ind1]*shift[i][k][ind2];
        if (i==0) ++ncorrelation[k][j];
      }
      --ind2;
      if (ind2<0) ind2+=p;
    }
  } else {
    int ind2=ind1-dmin;
    for (unsigned int j=dmin;j<p;++j) {
      if (ind2<0) ind2+=p;
      if (shift[i][k][ind2] > -1e10) {
        correlation[i][k][j]+= shift[i][k][ind1]*shift[i][k][ind2];
        if (i==0) ++ncorrelation[k][j];
      }
      --ind2;
    }
  }

  if (i==npair-1) {
    ++insertindex[k];
    if (insertindex[k]==p) insertindex[k]=0;
  }
}


/* ----------------------------------------------------------------------
   Add 2 scalar values to the cross-correlator k of pair i
------------------------------------------------------------------------- */
void FixAveCorrelateLong::add(const int i, const double wA, const double wB,
                        const int k) {
  if (k == numcorrelators) return;
  if (k > kmax) kmax=k;

  shift[i][k][insertindex[k]] = wA;
  shift2[i][k][insertindex[k]] = wB;

  accumulator[i][k] += wA;
  accumulator2[i][k] += wB;
  if (i==0) ++naccumulator[k];
  if (naccumulator[k]==m) {
    add(i,accumulator[i][k]/m, accumulator2[i][k]/m,k+1);
    accumulator[i][k]=0;
    accumulator2[i][k]=0;
    if (i==npair-1) naccumulator[k]=0;
  }

  unsigned int ind1=insertindex[k];
  if (k==0) {
    int ind2=ind1;
    for (unsigned int j=0;j<p;++j) {
      if (shift[i][k][ind2] > -1e10) {
        correlation[i][k][j]+= shift[i][k][ind1]*shift2[i][k][ind2];
        if (i==0) ++ncorrelation[k][j];
      }
      --ind2;
      if (ind2<0) ind2+=p;
    }
  }
  else {
    int ind2=ind1-dmin;
    for (unsigned int j=dmin;j<p;++j) {
      if (ind2<0) ind2+=p;
      if (shift[i][k][ind2] > -1e10) {
        correlation[i][k][j]+= shift[i][k][ind1]*shift2[i][k][ind2];
        if (i==0) ++ncorrelation[k][j];
      }
      --ind2;
    }
  }

  if (i==npair-1) {
    ++insertindex[k];
    if (insertindex[k]==p) insertindex[k]=0;
  }
}


/* ----------------------------------------------------------------------
   nvalid = next step on which end_of_step does something
   this step if multiple of nevery, else next multiple
   startstep is lower bound
------------------------------------------------------------------------- */

bigint FixAveCorrelateLong::nextvalid()
{
  bigint nvalid = update->ntimestep;
  if (startstep > nvalid) nvalid = startstep;
  if (nvalid % nevery) nvalid = (nvalid/nevery)*nevery + nevery;
  return nvalid;
}


/* ----------------------------------------------------------------------
   memory_usage
------------------------------------------------------------------------- */
double FixAveCorrelateLong::memory_usage() {
  //    shift:            npair x numcorrelators x p
  //    shift2:           npair x numcorrelators x p
  //    correlation:      npair x numcorrelators x p
  //    accumulator:      npair x numcorrelators
  //    accumulator2:     npair x numcorrelators
  //    ncorrelation:     numcorrelators x p
  //    naccumulator:     numcorrelators
  //    insertindex:      numcorrelators
  //    t:              numcorrelators x p
  //    f:              npair x numcorrelators x p
  double bytes = (4*npair*numcorrelators*p + 2*npair*numcorrelators
                  + numcorrelators*p)*sizeof(double)
    + numcorrelators*p*sizeof(unsigned long int)
    + 2*numcorrelators*sizeof(unsigned int);
  return bytes;
}

/* ----------------------------------------------------------------------
   Write Restart data to restart file
------------------------------------------------------------------------- */
// Save everything except t and f
void FixAveCorrelateLong::write_restart(FILE *fp) {
  if (me == 0) {
    int nsize = 3*npair*numcorrelators*p + 2*npair*numcorrelators
                + numcorrelators*p + 2*numcorrelators + 6;
    int n=0;
    double *list;
    memory->create(list,nsize,"correlator:list");
    list[n++]=npair;
    list[n++]=numcorrelators;
    list[n++]=p;
    list[n++]=m;
    list[n++] = last_accumulated_step;
    for (int i=0;i<npair;i++)
      for (int j=0;j<numcorrelators;j++) {
        for (unsigned int k=0;k<p;k++) {
          list[n++]=shift[i][j][k];
          list[n++]=shift2[i][j][k];
          list[n++]=correlation[i][j][k];
        }
        list[n++]=accumulator[i][j];
        list[n++]=accumulator2[i][j];
      }
    for (int i=0;i<numcorrelators;i++) {
      for (unsigned int j=0;j<p;j++) list[n++]=ncorrelation[i][j];
      list[n++]=naccumulator[i];
      list[n++]=insertindex[i];
    }

    int size = n*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
    memory->destroy(list);
  }
}


/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixAveCorrelateLong::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int npairin = static_cast<int> (list[n++]);
  int numcorrelatorsin = static_cast<int> (list[n++]);
  int pin = static_cast<int> (list[n++]);
  int min = static_cast<int> (list[n++]);
  last_accumulated_step = static_cast<int> (list[n++]);

  if ((npairin!=npair) || (numcorrelatorsin!=numcorrelators)
      || (pin!=(int)p) || (min!=(int)m))
    error->all(FLERR,"Fix ave/correlate/long: restart and input data are different");

  for (int i=0;i<npair;i++)
    for (int j=0;j<numcorrelators;j++) {
      for (unsigned int k=0;k<p;k++) {
        shift[i][j][k] = list[n++];
        shift2[i][j][k] = list[n++];
        correlation[i][j][k] = list[n++];
      }
      accumulator[i][j] = list[n++];
      accumulator2[i][j] = list[n++];
    }
  for (int i=0;i<numcorrelators;i++) {
    for (unsigned int j=0;j<p;j++)
      ncorrelation[i][j] = static_cast<unsigned long int>(list[n++]);
    naccumulator[i] = static_cast<unsigned int> (list[n++]);
    insertindex[i] = static_cast<unsigned int> (list[n++]);
  }
}
