// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
     Benoit Leblanc, Dave Rigby, Paul Saxe (Materials Design)
     Reese Jones (Sandia)
------------------------------------------------------------------------- */

#include "fix_ave_correlate.h"

#include "arg_info.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ONE,RUNNING};
enum{AUTO,UPPER,LOWER,AUTOUPPER,AUTOLOWER,FULL};

/* ---------------------------------------------------------------------- */

FixAveCorrelate::FixAveCorrelate(LAMMPS * lmp, int narg, char **arg):
  Fix (lmp, narg, arg),
  nvalues(0), which(nullptr), argindex(nullptr), value2index(nullptr), ids(nullptr), fp(nullptr),
  count(nullptr), values(nullptr), corr(nullptr), save_count(nullptr), save_corr(nullptr)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/correlate command");

  MPI_Comm_rank(world,&me);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  nfreq = utils::inumeric(FLERR,arg[5],false,lmp);

  global_freq = nfreq;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-6,&arg[6],0,earg,lmp);

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  nvalues = 0;

  int iarg = 0;
  while (iarg < nargnew) {
    ArgInfo argi(arg[iarg]);

    if (argi.get_type() == ArgInfo::NONE) break;
    if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR,"Invalid fix ave/correlate command");

    which[nvalues] = argi.get_type();
    argindex[nvalues] = argi.get_index1();
    ids[nvalues] = argi.copy_name();

    nvalues++;
    iarg++;
  }

  // optional args

  type = AUTO;
  ave = ONE;
  startstep = 0;
  prefactor = 1.0;
  fp = nullptr;
  overwrite = 0;
  char *title1 = nullptr;
  char *title2 = nullptr;
  char *title3 = nullptr;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      if (strcmp(arg[iarg+1],"auto") == 0) type = AUTO;
      else if (strcmp(arg[iarg+1],"upper") == 0) type = UPPER;
      else if (strcmp(arg[iarg+1],"lower") == 0) type = LOWER;
      else if (strcmp(arg[iarg+1],"auto/upper") == 0) type = AUTOUPPER;
      else if (strcmp(arg[iarg+1],"auto/lower") == 0) type = AUTOLOWER;
      else if (strcmp(arg[iarg+1],"full") == 0) type = FULL;
      else error->all(FLERR,"Illegal fix ave/correlate command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else error->all(FLERR,"Illegal fix ave/correlate command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"prefactor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      prefactor = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open fix ave/correlate file {}:"" {}",
                                       arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      delete [] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      delete [] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/correlate command");
      delete [] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/correlate command");
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/correlate command");
  if (nfreq % nevery)
    error->all(FLERR,"Illegal fix ave/correlate command");
  if (ave == ONE && nfreq < (nrepeat-1)*nevery)
    error->all(FLERR,"Illegal fix ave/correlate command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/correlate command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/correlate does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/correlate compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/correlate compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
        error->all(FLERR,"Fix ave/correlate compute vector "
                   "is accessed out-of-range");

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/correlate does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,"Fix ave/correlate fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix ave/correlate fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,
                   "Fix ave/correlate fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,"Fix for fix ave/correlate "
                   "not computed at compatible time");

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/correlate does not exist");
      if (argindex[i] == 0 && input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR,
                   "Fix ave/correlate variable is not equal-style variable");
      if (argindex[i] && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,
                   "Fix ave/correlate variable is not vector-style variable");
    }
  }

  // npair = # of correlation pairs to calculate

  if (type == AUTO) npair = nvalues;
  if (type == UPPER || type == LOWER) npair = nvalues*(nvalues-1)/2;
  if (type == AUTOUPPER || type == AUTOLOWER) npair = nvalues*(nvalues+1)/2;
  if (type == FULL) npair = nvalues*nvalues;

  // print file comment lines

  if (fp && me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-correlated data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# Timestep Number-of-time-windows\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else {
      fprintf(fp,"# Index TimeDelta Ncount");
      if (type == AUTO)
        for (int i = 0; i < nvalues; i++)
          fprintf(fp," %s*%s",earg[i],earg[i]);
      else if (type == UPPER)
        for (int i = 0; i < nvalues; i++)
          for (int j = i+1; j < nvalues; j++)
            fprintf(fp," %s*%s",earg[i],earg[j]);
      else if (type == LOWER)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < i-1; j++)
            fprintf(fp," %s*%s",earg[i],earg[j]);
      else if (type == AUTOUPPER)
        for (int i = 0; i < nvalues; i++)
          for (int j = i; j < nvalues; j++)
            fprintf(fp," %s*%s",earg[i],earg[j]);
      else if (type == AUTOLOWER)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < i; j++)
            fprintf(fp," %s*%s",earg[i],earg[j]);
      else if (type == FULL)
        for (int i = 0; i < nvalues; i++)
          for (int j = 0; j < nvalues; j++)
            fprintf(fp," %s*%s",earg[i],earg[j]);
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing file header");

    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // allocate and initialize memory for averaging
  // set count and corr to zero since they accumulate
  // also set save versions to zero in case accessed via compute_array()

  memory->create(values,nrepeat,nvalues,"ave/correlate:values");
  memory->create(count,nrepeat,"ave/correlate:count");
  memory->create(save_count,nrepeat,"ave/correlate:save_count");
  memory->create(corr,nrepeat,npair,"ave/correlate:corr");
  memory->create(save_corr,nrepeat,npair,"ave/correlate:save_corr");

  int i,j;
  for (i = 0; i < nrepeat; i++) {
    save_count[i] = count[i] = 0;
    for (j = 0; j < npair; j++)
      save_corr[i][j] = corr[i][j] = 0.0;
  }

  // this fix produces a global array

  array_flag = 1;
  size_array_rows = nrepeat;
  size_array_cols = npair+2;
  extarray = 0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  lastindex = -1;
  firstindex = 0;
  nsample = 0;
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveCorrelate::~FixAveCorrelate()
{
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  memory->destroy(values);
  memory->destroy(count);
  memory->destroy(save_count);
  memory->destroy(corr);
  memory->destroy(save_corr);

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixAveCorrelate::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveCorrelate::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/correlate does not exist");
      value2index[i] = icompute;

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/correlate does not exist");
      value2index[i] = ifix;

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/correlate does not exist");
      value2index[i] = ivariable;
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    lastindex = -1;
    firstindex = 0;
    nsample = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveCorrelate::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveCorrelate::end_of_step()
{
  int i,j,m;
  double scalar;

  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/correlate");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // accumulate results of computes,fixes,variables to origin
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // lastindex = index in values ring of latest time sample

  lastindex++;
  if (lastindex == nrepeat) lastindex = 0;

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];

    // invoke compute if not previously invoked

    if (which[i] == ArgInfo::COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
        if (!(compute->invoked_flag & Compute::INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        scalar = compute->scalar;
      } else {
        if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        scalar = compute->vector[argindex[i]-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == ArgInfo::FIX) {
      if (argindex[i] == 0)
        scalar = modify->fix[m]->compute_scalar();
      else
        scalar = modify->fix[m]->compute_vector(argindex[i]-1);

    // evaluate equal-style or vector-style variable

    } else if (which[i] == ArgInfo::VARIABLE) {
      if (argindex[i] == 0)
        scalar = input->variable->compute_equal(m);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(m,&varvec);
        int index = argindex[i];
        if (nvec < index) scalar = 0.0;
        else scalar = varvec[index-1];
      }
    }

    values[lastindex][i] = scalar;
  }

  // fistindex = index in values ring of earliest time sample
  // nsample = number of time samples in values ring

  if (nsample < nrepeat) nsample++;
  else {
    firstindex++;
    if (firstindex == nrepeat) firstindex = 0;
  }

  nvalid += nevery;
  modify->addstep_compute(nvalid);

  // calculate all Cij() enabled by latest values

  accumulate();
  if (ntimestep % nfreq) return;

  // save results in save_count and save_corr

  for (i = 0; i < nrepeat; i++) {
    save_count[i] = count[i];
    if (count[i])
      for (j = 0; j < npair; j++)
        save_corr[i][j] = prefactor*corr[i][j]/count[i];
    else
      for (j = 0; j < npair; j++)
        save_corr[i][j] = 0.0;
  }

  // output result to file

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d\n",ntimestep,nrepeat);
    for (i = 0; i < nrepeat; i++) {
      fprintf(fp,"%d %d %d",i+1,i*nevery,count[i]);
      if (count[i])
        for (j = 0; j < npair; j++)
          fprintf(fp," %g",prefactor*corr[i][j]/count[i]);
      else
        for (j = 0; j < npair; j++)
          fprintf(fp," 0.0");
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing out correlation data");

    fflush(fp);

    if (overwrite) {
      long fileend = ftell(fp);
      if ((fileend > 0) && (ftruncate(fileno(fp),fileend)))
        perror("Error while tuncating output");
    }
  }

  // zero accumulation if requested
  // recalculate Cij(0)

  if (ave == ONE) {
    for (i = 0; i < nrepeat; i++) {
      count[i] = 0;
      for (j = 0; j < npair; j++)
        corr[i][j] = 0.0;
    }
    nsample = 1;
    accumulate();
  }
}

/* ----------------------------------------------------------------------
   accumulate correlation data using more recently added values
------------------------------------------------------------------------- */

void FixAveCorrelate::accumulate()
{
  int i,j,k,m,n,ipair;

  for (k = 0; k < nsample; k++) count[k]++;

  if (type == AUTO) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++) {
        corr[k][ipair++] += values[m][i]*values[n][i];
      }
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == UPPER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = i+1; j < nvalues; j++)
          corr[k][ipair++] += values[m][i]*values[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == LOWER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j < i; j++)
          corr[k][ipair++] += values[m][i]*values[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == AUTOUPPER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = i; j < nvalues; j++)
          corr[k][ipair++] += values[m][i]*values[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == AUTOLOWER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j <= i; j++)
          corr[k][ipair++] += values[m][i]*values[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == FULL) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j < nvalues; j++)
          corr[k][ipair++] += values[m][i]*values[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  }
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveCorrelate::compute_array(int i, int j)
{
  if (j == 0) return 1.0*i*nevery;
  else if (j == 1) return 1.0*save_count[i];
  else if (save_count[i]) return save_corr[i][j-2];
  return 0.0;
}

/* ----------------------------------------------------------------------
   nvalid = next step on which end_of_step does something
   this step if multiple of nevery, else next multiple
   startstep is lower bound
------------------------------------------------------------------------- */

bigint FixAveCorrelate::nextvalid()
{
  bigint nvalid = update->ntimestep;
  if (startstep > nvalid) nvalid = startstep;
  if (nvalid % nevery) nvalid = (nvalid/nevery)*nevery + nevery;
  return nvalid;
}
