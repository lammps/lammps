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
   Contributing author: Sebastian Huetter (OvGU)
------------------------------------------------------------------------- */

#include "fix_ave_moments.h"

#include "arg_info.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathSpecial::square;
using MathSpecial::cube;

enum { MEAN, STDDEV, VARIANCE, SKEW, KURTOSIS };

/* ---------------------------------------------------------------------- */

FixAveMoments::FixAveMoments(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), result_list(nullptr), window_list(nullptr)
{
  // this fix's data is always accessible (but might be meaningless)
  global_freq = 1;
  dynamic_group_allow = 1;
  time_depend = 1;

  // EXAMPLE:
  //    fix ID group-ID ave/moments Nevery Nrepeat Nfreq value1 ... valueN moment1 ... momentM keyword value ...

  // the first six arguments are fixed & need at least one input and moment
  const int nfixedargs = 6;
  if (narg < nfixedargs + 2) utils::missing_cmd_args(FLERR, "fix ave/moments", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  nfreq = utils::inumeric(FLERR,arg[5],false,lmp);

  // scan values to count them

  nvalues = 0;
  // first input name is position after the fixed args
  int iarg = nfixedargs;
  while (iarg < narg) {
    if (utils::strmatch(arg[iarg],"^[cfv]_")) {
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues == 0)
    error->all(FLERR, nfixedargs,
               "No values from computes, fixes, or variables used in fix ave/moments command");

  // next, the moments
  iarg = consume_moments(iarg, narg, arg);
  if (moments.empty())
    error->all(FLERR, nfixedargs,
               "No values from computes, fixes, or variables used in fix ave/moments command");

  // parse optional keywords which must follow the data

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  int *amap = nullptr;
  nvalues = utils::expand_args(FLERR, nvalues, &arg[nfixedargs], /* mode=scalar */ 0, earg, lmp, &amap);

  if (earg != &arg[nfixedargs]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  for (int i = 0; i < nvalues; i++) {
    ArgInfo argi(arg[i]);

    value_t val;
    val.keyword = arg[i];
    val.which = argi.get_type();

    val.argindex = argi.get_index1();
    val.iarg = (expand ? amap[i] : i) + nfixedargs;
    val.varlen = 0;
    val.id = argi.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::NONE) || (val.which == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR, val.iarg, "Invalid fix ave/moments argument: {}", arg[i]);

    values.push_back(std::move(val));
  }
  if (nvalues != (int)values.size())
    error->all(FLERR, Error::NOPOINTER,
               "Could not parse value data consistently for fix ave/moments");

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix ave/moments nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR, 4, "Illegal fix ave/moments nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR, 5, "Illegal fix ave/moments nfreq value: {}", nfreq);

  for (auto &val : values) {
    switch (val.which) {
      case ArgInfo::COMPUTE:
        val.val.c = modify->get_compute_by_id(val.id);
        if (!val.val.c)
          error->all(FLERR, val.iarg, "Compute ID {} for fix ave/moments does not exist", val.id);
        if (val.argindex == 0 && (val.val.c->scalar_flag == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments compute {} does not calculate a scalar", val.id);
        if (val.argindex && (val.val.c->vector_flag == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments compute {} does not calculate a vector", val.id);
        if (val.argindex && (val.argindex > val.val.c->size_vector) &&
            (val.val.c->size_vector_variable == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments compute {} vector is accessed out-of-range{}",
                     val.id, utils::errorurl(20));
        if (val.argindex && val.val.c->size_vector_variable) val.varlen = 1;
        break;

      case ArgInfo::FIX:
        val.val.f = modify->get_fix_by_id(val.id);
        if (!val.val.f) error->all(FLERR,"Fix ID {} for fix ave/moments does not exist", val.id);
        if ((val.argindex == 0) && (val.val.f->scalar_flag == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments fix {} does not calculate a scalar", val.id);
        if (val.argindex && (val.val.f->vector_flag == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments fix {} does not calculate a vector", val.id);
        if (val.argindex && (val.val.f->size_vector_variable))
          error->all(FLERR, val.iarg, "Fix ave/moments fix {} vector cannot be variable length", val.id);
        if (val.argindex && (val.argindex > val.val.f->size_vector))
          error->all(FLERR, val.iarg, "Fix ave/moments fix {} vector is accessed out-of-range{}",
                     val.id, utils::errorurl(20));
        if (nevery % val.val.f->global_freq)
          error->all(FLERR, val.iarg, "Fix {} for fix ave/moments not computed at compatible time{}",
                     val.id, utils::errorurl(7));
        break;

      case ArgInfo::VARIABLE:
        int ivariable = input->variable->find(val.id.c_str());
        if (ivariable < 0)
          error->all(FLERR, val.iarg, "Variable name {} for fix ave/moments does not exist", val.id);
        if ((val.argindex == 0) && (input->variable->equalstyle(ivariable) == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments variable {} is not equal-style variable", val.id);
        if ((val.argindex) && (input->variable->vectorstyle(ivariable) == 0))
          error->all(FLERR, val.iarg, "Fix ave/moments variable {} is not vector-style variable",
                     val.id);
        break;
    }
  }

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
    memory->sfree(amap);
  }

  // allocate memory for averaging

  window_list = nullptr;
  result_list = nullptr;

  // one window of nvalues columns and nrepeat rows (=all scalars of one value are consecutive)
  memory->create(window_list, nvalues, nrepeat, "ave/moments:window_list");
  for (int i = 0; i < nvalues; i++)
    for (int j = 0; j < nrepeat; j++)
      window_list[i][j] = 0.0;

  // this fix produces a global vector and array

  vector_flag = 1;
  size_vector = nvalues * moments.size();
  array_flag = 1;
  size_array_rows = size_vector;
  size_array_cols = nhistory;

  // produce nmoments outputs per value with nhistory depth
  memory->create(result_list, nhistory, size_vector, "ave/moments:result_list");
  for (int i = 0; i < nhistory; i++)
    for (int j = 0; j < size_vector; j++)
      result_list[i][j] = 0.0;

  // intensive/extensive flags set by compute,fix,variable that produces value

  extvector = -1;
  extarray = -2;
  extlist = new int[size_vector];
  int extvalue = 0;
  int i = 0;
  for (auto &val : values) {
    switch (val.which) {
      case ArgInfo::COMPUTE:
        if (val.argindex == 0) extvalue = val.val.c->extscalar;
        else if (val.val.f->extvector >= 0) extvalue = val.val.c->extvector;
        else extvalue = val.val.c->extlist[val.argindex-1];
        break;

      case ArgInfo::FIX:
        if (val.argindex == 0) extvalue = val.val.f->extscalar;
        else if (val.val.f->extvector >= 0) extvalue = val.val.f->extvector;
        else extvalue = val.val.f->extlist[val.argindex-1];
        break;

      case ArgInfo::VARIABLE:
        extvalue = 0;
        break;
    }
    if (extvalue == -1)
      error->all(FLERR, Error::NOLASTLINE, "Fix ave/moments cannot set output array "
                 "intensive/extensive from these inputs");
    if (extarray < -1) extarray = extvalue;
    else if (extvalue != extarray)
      error->all(FLERR, Error::NOLASTLINE, "Fix ave/moments cannot set output array "
                 "intensive/extensive from these inputs");
    for (int j=0; j < (int)moments.size(); j++)
      extlist[i + j] = extvalue;
    i += moments.size();
  }

  // initializations

  iwindow = window_filled = 0;
  iresult = 0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_comp_next = -1;
  nvalid = -1;
  setnextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveMoments::~FixAveMoments()
{
  values.clear();
  moments.clear();
  delete[] extlist;

  memory->destroy(window_list);
  memory->destroy(result_list);
}

/* ---------------------------------------------------------------------- */

int FixAveMoments::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveMoments::init()
{
  // update indices/pointers for all computes,fixes,variables

  for (auto &val : values) {
    switch (val.which) {
      case ArgInfo::COMPUTE:
        val.val.c = modify->get_compute_by_id(val.id);
        if (!val.val.c)
          error->all(FLERR, Error::NOLASTLINE, "Compute ID {} for fix ave/moments does not exist",
                     val.id);
        break;

      case ArgInfo::FIX:
        val.val.f = modify->get_fix_by_id(val.id);
        if (!val.val.f)
          error->all(FLERR, Error::NOLASTLINE, "Fix ID {} for fix ave/moments does not exist", val.id);
        break;

      case ArgInfo::VARIABLE:
        val.val.v = input->variable->find(val.id.c_str());
        if (val.val.v < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix ave/moments does not exist",
                     val.id);
        break;
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    setnextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveMoments::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveMoments::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // always take new values
  append_values();

  // if window boundary reached, do a compute, otherwise just schedule next take
  if (ntimestep == nvalid_comp_next) {
    update_results();
    setnextvalid();
  } else {
    nvalid += nevery;
  }

  modify->addstep_compute(nvalid);
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveMoments::compute_scalar()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveMoments::compute_vector(int i)
{
  return compute_array(i, 0);
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveMoments::compute_array(int i, int j)
{
  if (i >= size_vector) return 0.0;
  if (j >= nhistory) return 0.0;
  // locate the j'th previous result in the ring buffer, relative to the
  // row before iresult (the current insert cursor)
  int row = (iresult - 1 - j + nhistory) % nhistory;
  return result_list[row][i];
}

/* ----------------------------------------------------------------------
   parse moment names
------------------------------------------------------------------------- */

int FixAveMoments::consume_moments(int iarg, int narg, char **arg)
{
  moments.clear();

  while (iarg < narg) {
    if (strcmp(arg[iarg],"mean") == 0)
      moments.push_back(MEAN);
    else if (strcmp(arg[iarg],"stddev") == 0)
      moments.push_back(STDDEV);
    else if (strcmp(arg[iarg],"variance") == 0)
      moments.push_back(VARIANCE);
    else if (strcmp(arg[iarg],"skew") == 0)
      moments.push_back(SKEW);
    else if (strcmp(arg[iarg],"kurtosis") == 0)
      moments.push_back(KURTOSIS);
    else
      break;
    iarg++;
  }
  return iarg;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveMoments::options(int iarg, int narg, char **arg)
{
  // option defaults

  nhistory = 1;
  startstep = 0;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"history") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/moments history", error);
      nhistory = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nhistory <= 0)
        error->all(FLERR, iarg+2, "Illegal ave/moments history argument {}; must be > 0",
                   nhistory);
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/moments start", error);
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix ave/moments keyword {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   return next timestep no earlier than `after`, rounded to next
   multiple of freq
------------------------------------------------------------------------- */

bigint next_after(const bigint ts, const bigint after, const int freq)
{
  if (ts >= after) return ts;
  return ts + ((after - ts) / freq + 1) * freq;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   this is either a step to take data
    or a step to take and compute the values (nfreq multiple)
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

void FixAveMoments::setnextvalid()
{
  bigint ntimestep = update->ntimestep;

  if (nvalid_comp_next > ntimestep) {
    // next window end boundary is still in the future, just increment
    nvalid = ntimestep + nevery;
    return;
  }

  // get next window end first
  bigint next_comp = (ntimestep/nfreq)*nfreq + nfreq;
  nvalid_comp_next = next_after(next_comp, startstep, nfreq);

  // from there, calculate the first time we have to take a value
  bigint ntake = nvalid_comp_next - static_cast<bigint>(nrepeat-1)*nevery;
  nvalid = next_after(ntake, ntimestep, nevery);
}

/* ---------------------------------------------------------------------- */

void FixAveMoments::get_values(std::vector<double>& scalars)
{
  // accumulate results of computes,fixes,variables to local copy
  int i = 0;
  double scalar = 0.0;
  for (auto &val : values) {
    switch (val.which) {
      case ArgInfo::COMPUTE:
        // invoke compute if not previously invoked
        // ensure no out-of-range access to variable-length compute vector
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
          if (val.varlen && (val.val.c->size_vector < val.argindex)) scalar = 0.0;
          else scalar = val.val.c->vector[val.argindex-1];
        }
        break;

      case ArgInfo::FIX:
        // access fix fields, guaranteed to be ready
        if (val.argindex == 0)
          scalar = val.val.f->compute_scalar();
        else
          scalar = val.val.f->compute_vector(val.argindex-1);
        break;

      case ArgInfo::VARIABLE:
        // evaluate equal-style or vector-style variable
        // if index exceeds vector length, use a zero value
        //   this can be useful if vector length is not known a priori
        if (val.argindex == 0)
          scalar = input->variable->compute_equal(val.val.v);
        else {
          double *varvec;
          int nvec = input->variable->compute_vector(val.val.v,&varvec);
          if (val.argindex > nvec) scalar = 0.0;
          else scalar = varvec[val.argindex-1];
        }
        break;
    }

    scalars[i] = scalar;
    ++i;
  }
}

/* ---------------------------------------------------------------------- */

void FixAveMoments::append_values()
{
  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  std::vector<double> scalars(nvalues);
  get_values(scalars);

  // transpose for faster access later
  for (int i=0; i<nvalues; i++) {
    window_list[i][iwindow] = scalars[i];
  }

  if (++iwindow >= nrepeat) {
    window_filled = 1;
    iwindow = 0;
  }
}

void FixAveMoments::update_results()
{
  const int count = window_filled ? nrepeat : iwindow;
  // Delay until we can safely do all moments. Avoids branching in the hot loop.
  if (count<3) return;

  double *result = result_list[iresult];

  // zero out previous values

  for (int i = 0; i < size_vector; i++)
    result[i] = 0.0;

  const double inv_n = 1.0 / count;
  const double fk2 = (double)count / (count - 1);
  const double fk3 = square((double)count) / ((count - 1) * (count - 2));
  const double np1_nm3 = (count+1.0)/(count-3.0);
  const double _3_nm1_nm3 = 3.0 * (count-1.0)/(count-3.0);

  // Each value is a series that can be processed individually
  for (int i = 0; i < nvalues; i++) {
    const double* series = window_list[i];

    // first pass: mean
    double mean = 0.0;
    for (int j = 0; j<count; j++)
      mean += series[j] * inv_n;

    // second pass: calculate biased sample moments
    double m2 = 0.0;
    double m3 = 0.0;
    double m4 = 0.0;
    for (int j = 0; j<count; j++) {
      const double dx = series[j] - mean;
      double y = square(dx) * inv_n;
      m2 += y;
      y *= dx;
      m3 += y;
      y *= dx;
      m4 += y;
    }
    // obtain unbiased cumulants as defined by Cramér and
    // reported i.e. in https://doi.org/10.1111/1467-9884.00122
    const double k2 = fk2 * m2;
    const double k3 = fk3 * m3;
    const double k4 = fk3 * (np1_nm3 * m4 - _3_nm1_nm3 * square(m2));
    // corrected sample standard deviation
    const double stddev = sqrt(k2);
    // adjusted Fisher-Pearson standardized moment coefficient G1
    const double G1 = k3 / cube(stddev);
    // adjusted Fisher-Pearson standardized moment coefficient G2 (from unbiased cumulant)
    const double G2 = k4 / square(k2);

    // map to result array, starting at value interleave offset
    double* rfirst = &result[i * moments.size()];
    for (int j = 0; j < (int)moments.size(); j++) {
      switch(moments[j]) {
        case MEAN:
          rfirst[j] = mean;
          break;
        case STDDEV:
          rfirst[j] = stddev;
          break;
        case VARIANCE:
          rfirst[j] = k2;
          break;
        case SKEW:
          rfirst[j] = G1;
          break;
        case KURTOSIS:
          rfirst[j] = G2;
          break;
      }
    }
  }

  if (++iresult >= nhistory)
    iresult = 0;
}
