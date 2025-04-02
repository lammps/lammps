/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_ave_gauss.h"

#include "arg_info.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAveGauss::FixAveGauss(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), window_list(nullptr), result_list(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal fix ave/gauss command");

  // this fix's data is always accessible (but might be meaningless)

  nevery = 1;
  global_freq = 1;

  dynamic_group_allow = 1;

  // fixed arguments

  nfreq = utils::inumeric(FLERR,arg[3],false,lmp);
  nwindow = utils::inumeric(FLERR,arg[4],false,lmp);

  if (nfreq <= 0) error->all(FLERR, 3, "Illegal fix ave/gauss nfreq value: {}", nfreq);
  if (nwindow  <= 0) error->all(FLERR, 4, "Illegal fix ave/gauss nwindow value: {}", nwindow);

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  // first 5 args are fixed
  const int ioffset = 5;
  int iarg = ioffset;
  while (iarg < narg) {
    if (utils::strmatch(arg[iarg],"^v_")) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/gauss command");

  values.clear();
  delays.clear();

  options(iarg,narg,arg);

  // parse values

  for (int i = 0; i < nvalues; i++) {
    ArgInfo argi(arg[i + ioffset]);

    value_t val;
    val.iarg = i + ioffset;
    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();

    if ((val.which != ArgInfo::VARIABLE) || (argi.get_dim() > 1))
      error->all(FLERR, val.iarg, "Invalid fix ave/gauss argument: {}", arg[i]);

    values.push_back(val);
  }
  if (nvalues != (int)values.size())
    error->all(FLERR, Error::NOPOINTER,
               "Could not parse value data consistently for fix ave/gauss");

  // setup and error check
  for (auto &val : values) {
    int ivariable = input->variable->find(val.id.c_str());
    if (ivariable < 0)
      error->all(FLERR, val.iarg, "Variable name {} for fix ave/gauss does not exist", val.id);
    if ((val.argindex == 0) && (input->variable->equalstyle(ivariable) == 0))
      error->all(FLERR, val.iarg, "Fix ave/gauss variable {} is not equal-style variable", val.id);
    if ((val.argindex) && (input->variable->vectorstyle(ivariable) == 0))
      error->all(FLERR, val.iarg, "Fix ave/gauss variable {} is not vector-style variable", val.id);
  }

  // allocate memory for averaging

  window_list = nullptr;
  result_list = nullptr;

  // one window of length nwindow
  memory->create(window_list, nwindow, nvalues, "ave/gauss:window_list");
  for (int i = 0; i < nwindow; i++)
    for (int j = 0; j < nvalues; j++)
      window_list[i][j] = 0.0;


  // for lookback results, keep the longest delay rows, produce 2 outputs per value
  int delaymax = *std::max_element(delays.begin(), delays.end());
  nresult = delaymax + 1;
  memory->create(result_list, nresult, nvalues*2, "ave/gauss:result_list");
  for (int i = 0; i < nresult; i++)
    for (int j = 0; j < nvalues*2; j++)
      result_list[i][j] = 0.0;

  // this fix produces a global vector and array
  vector_flag = 1;
  size_vector = nvalues*2;
  extvector = 0;
  array_flag = 1;
  size_array_rows = nvalues*2;
  size_array_cols = delays.size();
  extarray = 0;

  // initializations

  iwindow = window_filled = 0;
  iresult = 0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = update->ntimestep;
  nvalid_last = -1;
  nfull_next = update->ntimestep;
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveGauss::~FixAveGauss()
{
  values.clear();
  delays.clear();

  memory->destroy(window_list);
  memory->destroy(result_list);
}

/* ---------------------------------------------------------------------- */

int FixAveGauss::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::init()
{
  // update indices/pointers for all computes,fixes,variables

  for (auto &val : values) {
    val.v = input->variable->find(val.id.c_str());
    if (val.v < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix ave/gauss does not exist", val.id);
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    nvalid = update->ntimestep;
    nfull_next = update->ntimestep;
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveGauss::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  append_values(ntimestep);

  if (ntimestep == nfull_next) {
    update_results(ntimestep);
    nfull_next = ntimestep + nfreq;
    // skip ahead until we need to fill the window again
    nvalid = std::max(ntimestep, nfull_next - 1 - nwindow);
  }

  nvalid += 1;
  modify->addstep_compute(nvalid);
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::append_values(bigint ntimestep)
{
  // accumulate results of variable to local copy

  modify->clearstep_compute();

  int i;
  double scalar;
  double *row = window_list[iwindow];

  for (i = 0; i < values.size(); i++) {
    auto &val = values[i];

    if (val.argindex == 0)
      scalar = input->variable->compute_equal(val.v);
    else {
      double *varvec;
      int nvec = input->variable->compute_vector(val.v,&varvec);
      if (val.argindex > nvec) scalar = 0.0;
      else scalar = varvec[val.argindex-1];
    }

    row[i] = scalar;
  }

  iwindow += 1;
  if (iwindow >= nwindow) {
    window_filled = 1;
    iwindow = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::update_results(bigint ntimestep)
{
  int count = window_filled ? nwindow : iwindow;
  if (count<1) return;

  double invcount = 1.0 / count;
  double *result = result_list[iresult];
  for (int i = 0; i < nvalues; i++)
    result[i*2] = result[i*2+1] = 0.0;

  // first pass: mean
  for (int j = 0; j < count; j++) {
    double *row = window_list[j];
    for (int i = 0; i < nvalues; i++)
      result[i*2] += row[i] * invcount;
  }

  // second pass: de-biased variance
  for (int j = 0; j < count; j++) {
    double *row = window_list[j];
    for (int i = 0; i < nvalues; i++) {
      double x = row[i] - result[i*2];
      result[i*2+1] += x * x * invcount;
    }
  }

  // return as stddev
  for (int i = 0; i < nvalues; i++)
    result[i*2+1] = sqrt(result[i*2+1]);

  iresult += 1;
  if (iresult >= nresult)
    iresult = 0;
}


/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveGauss::compute_scalar()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveGauss::compute_vector(int i)
{
  return compute_array(i, 0);
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveGauss::compute_array(int i, int j)
{
  if (i >= nvalues*2) return 0.0;
  if (j >= delays.size()) return 0.0;
  int row = (iresult - 1 - delays[j] + nresult) % nresult;
  if (row >= nresult) return 0.0;
  return result_list[row][i];
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveGauss::options(int iarg, int narg, char **arg)
{
  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/gauss delay", error);
      int delay  = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (delay < 0)
        error->all(FLERR, iarg+1, "Illegal fix ave/gauss delay argument {}; must be > 0", delay);
      delays.push_back(delay);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix ave/gauss keyword {}", arg[iarg]);
  }

  if (delays.empty()) {
    delays.push_back(0);
  }
}
