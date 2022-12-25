// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { ONE, RUNNING };
enum { AUTO, UPPER, LOWER, AUTOUPPER, AUTOLOWER, FULL };

/* ---------------------------------------------------------------------- */

FixAveCorrelate::FixAveCorrelate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), fp(nullptr), count(nullptr), cvalues(nullptr), corr(nullptr),
    save_count(nullptr), save_corr(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ave/correlate", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);

  time_depend = 1;
  global_freq = nfreq;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR, narg - 6, &arg[6], 0, earg, lmp);

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values

  int iarg = 0;
  while (iarg < nargnew) {
    ArgInfo argi(arg[iarg]);
    value_t val;

    if (argi.get_type() == ArgInfo::NONE) break;
    if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR, "Unknown fix ave/correlate data type: {}", arg[iarg]);

    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    values.push_back(val);
    iarg++;
  }
  nvalues = values.size();

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
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate type", error);
      if (strcmp(arg[iarg+1],"auto") == 0) type = AUTO;
      else if (strcmp(arg[iarg+1],"upper") == 0) type = UPPER;
      else if (strcmp(arg[iarg+1],"lower") == 0) type = LOWER;
      else if (strcmp(arg[iarg+1],"auto/upper") == 0) type = AUTOUPPER;
      else if (strcmp(arg[iarg+1],"auto/lower") == 0) type = AUTOLOWER;
      else if (strcmp(arg[iarg+1],"full") == 0) type = FULL;
      else error->all(FLERR,"Unknown fix ave/correlate type: {}");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate ave", error);
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else error->all(FLERR,"Unknown fix ave/correlate ave mode: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate start", error);
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"prefactor") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate prefactor", error);
      prefactor = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate file", error);
      if (comm->me == 0) {
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
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate title1", error);
      delete[] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate title2", error);
      delete[] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/correlate title3", error);
      delete[] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unkown fix ave/correlate keyword: {}", arg[iarg]);
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0) error->all(FLERR,"Illegal fix ave/correlate nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR,"Illegal fix ave/correlate nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR,"Illegal fix ave/correlate nfreq value: {}", nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Inconsistent fix ave/correlate nevery/nrepeat/nfreq values");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Fix ave/correlate overwrite keyword requires ave running setting");

  for (auto &val : values) {

    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, "Compute ID {} for fix ave/correlate does not exist", val.id);
      if (val.argindex == 0 && val.val.c->scalar_flag == 0)
        error->all(FLERR, "Fix ave/correlate compute {} does not calculate a scalar", val.id);
      if (val.argindex && val.val.c->vector_flag == 0)
        error->all(FLERR, "Fix ave/correlate compute {} does not calculate a vector", val.id);
      if (val.argindex && val.argindex > val.val.c->size_vector)
        error->all(FLERR, "Fix ave/correlate compute {} vector is accessed out-of-range", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for fix ave/correlate does not exist", val.id);
      if (val.argindex == 0 && val.val.f->scalar_flag == 0)
        error->all(FLERR, "Fix ave/correlate fix {} does not calculate a scalar", val.id);
      if (val.argindex && val.val.f->vector_flag == 0)
        error->all(FLERR, "Fix ave/correlate fix {} does not calculate a vector", val.id);
      if (val.argindex && val.argindex > val.val.f->size_vector)
        error->all(FLERR, "Fix ave/correlate fix {} vector is accessed out-of-range", val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for fix ave/correlate not computed at compatible time", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for fix ave/correlate does not exist", val.id);
      if (val.argindex == 0 && input->variable->equalstyle(val.val.v) == 0)
        error->all(FLERR, "Fix ave/correlate variable {} is not equal-style variable", val.id);
      if (val.argindex && input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR, "Fix ave/correlate variable {} is not vector-style variable", val.id);
    }
  }

  // npair = # of correlation pairs to calculate

  if (type == AUTO) npair = nvalues;
  if (type == UPPER || type == LOWER) npair = nvalues*(nvalues-1)/2;
  if (type == AUTOUPPER || type == AUTOLOWER) npair = nvalues*(nvalues+1)/2;
  if (type == FULL) npair = nvalues*nvalues;

  // print file comment lines

  if (fp && comm->me == 0) {
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
      error->one(FLERR,"Error writing ave/correlate header: {}", utils::getsyserror());

    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // allocate and initialize memory for averaging
  // set count and corr to zero since they accumulate
  // also set save versions to zero in case accessed via compute_array()

  memory->create(cvalues,nrepeat,nvalues,"ave/correlate:cvalues");
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
  memory->destroy(cvalues);
  memory->destroy(count);
  memory->destroy(save_count);
  memory->destroy(corr);
  memory->destroy(save_corr);

  if (fp && comm->me == 0) fclose(fp);
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

  for (auto &val : values) {

    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, "Compute ID {} for fix ave/correlate does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for fix ave/correlate does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/correlate does not exist", val.id);
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
  int i,j;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // accumulate results of computes,fixes,variables to origin
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // lastindex = index in values ring of latest time sample

  lastindex++;
  if (lastindex == nrepeat) lastindex = 0;

  i = 0;
  for (auto &val : values) {
    double scalar = 0.0;

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {

      if (val.argindex == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
          val.val.c->compute_scalar();
          val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        scalar = val.val.c->scalar;
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        scalar = val.val.c->vector[val.argindex-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0)
        scalar = val.val.f->compute_scalar();
      else
        scalar = val.val.f->compute_vector(val.argindex-1);

    // evaluate equal-style or vector-style variable

    } else if (val.which == ArgInfo::VARIABLE) {
      if (val.argindex == 0)
        scalar = input->variable->compute_equal(val.val.v);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(val.val.v,&varvec);
        int index = val.argindex;
        if (nvec < index) scalar = 0.0;
        else scalar = varvec[index-1];
      }
    }

    cvalues[lastindex][i++] = scalar;
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

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    fmt::print(fp,"{} {}\n",ntimestep,nrepeat);
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
      error->one(FLERR,"Error writing out fix ave/correlate data: {}", utils::getsyserror());

    fflush(fp);

    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}", utils::getsyserror());
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
        corr[k][ipair++] += cvalues[m][i]*cvalues[n][i];
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
          corr[k][ipair++] += cvalues[m][i]*cvalues[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == LOWER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j < i; j++)
          corr[k][ipair++] += cvalues[m][i]*cvalues[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == AUTOUPPER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = i; j < nvalues; j++)
          corr[k][ipair++] += cvalues[m][i]*cvalues[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == AUTOLOWER) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j <= i; j++)
          corr[k][ipair++] += cvalues[m][i]*cvalues[n][j];
      m--;
      if (m < 0) m = nrepeat-1;
    }
  } else if (type == FULL) {
    m = n = lastindex;
    for (k = 0; k < nsample; k++) {
      ipair = 0;
      for (i = 0; i < nvalues; i++)
        for (j = 0; j < nvalues; j++)
          corr[k][ipair++] += cvalues[m][i]*cvalues[n][j];
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
