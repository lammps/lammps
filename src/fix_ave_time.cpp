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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "fix_ave_time.h"

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
#include <stdexcept>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ ONE, RUNNING, WINDOW };
enum{ SCALAR, VECTOR };

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), fp(nullptr), offlist(nullptr), format(nullptr), format_user(nullptr),
  vector(nullptr), vector_total(nullptr), vector_list(nullptr),
  column(nullptr), array(nullptr), array_total(nullptr), array_list(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ave/time", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  nfreq = utils::inumeric(FLERR,arg[5],false,lmp);

  global_freq = nfreq;

  dynamic_group_allow = 1;
  time_depend = 1;

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;
  int iarg = 6;
  while (iarg < narg) {
    if (utils::strmatch(arg[iarg],"^[cfv]_")) {
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues == 0)
    error->all(FLERR,"No values from computes, fixes, or variables used in fix ave/time command");

  // parse optional keywords

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = utils::expand_args(FLERR,nvalues,&arg[6],mode,earg,lmp);
  key2col.clear();

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  for (int i = 0; i < nvalues; i++) {
    ArgInfo argi(arg[i]);

    value_t val;
    val.keyword = arg[i];
    val.which = argi.get_type();
    key2col[arg[i]] = i;

    if ((val.which == ArgInfo::NONE) || (val.which == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR,"Invalid fix ave/time argument: {}", arg[i]);

    val.argindex = argi.get_index1();
    val.varlen = 0;
    val.offcol = 0;
    val.id = argi.get_name();
    val.val.c = nullptr;

    values.push_back(val);
  }
  if (nvalues != (int)values.size())
    error->all(FLERR, "Could not parse value data consistently for fix ave/time");

  // set off columns now that nvalues is finalized

  for (int i = 0; i < noff; i++) {
    if (offlist[i] < 1 || offlist[i] > nvalues)
      error->all(FLERR,"Invalid fix ave/time off column: {}", offlist[i]);
    values[offlist[i]-1].offcol = 1;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable
  // set variable_length if any compute is variable length

  if (nevery <= 0) error->all(FLERR,"Illegal fix ave/time nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR,"Illegal fix ave/time nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR,"Illegal fix ave/time nfreq value: {}", nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Inconsistent fix ave/time nevery/nrepeat/nfreq values");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Fix ave/time overwrite keyword requires ave running setting");

  for (auto &val : values) {

    if ((val.which == ArgInfo::COMPUTE) && (mode == SCALAR)) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for fix ave/time does not exist", val.id);
      if (val.argindex == 0 && (val.val.c->scalar_flag == 0))
        error->all(FLERR,"Fix ave/time compute {} does not calculate a scalar", val.id);
      if (val.argindex && (val.val.c->vector_flag == 0))
        error->all(FLERR,"Fix ave/time compute {} does not calculate a vector", val.id);
      if (val.argindex && (val.argindex > val.val.c->size_vector) &&
          (val.val.c->size_vector_variable == 0))
        error->all(FLERR, "Fix ave/time compute {} vector is accessed out-of-range", val.id);
      if (val.argindex && val.val.c->size_vector_variable) val.varlen = 1;

    } else if ((val.which == ArgInfo::COMPUTE) && (mode == VECTOR)) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for fix ave/time does not exist", val.id);
      if ((val.argindex == 0) && (val.val.c->vector_flag == 0))
        error->all(FLERR,"Fix ave/time compute {} does not calculate a vector", val.id);
      if (val.argindex && (val.val.c->array_flag == 0))
        error->all(FLERR,"Fix ave/time compute {} does not calculate an array", val.id);
      if (val.argindex && (val.argindex > val.val.c->size_array_cols))
        error->all(FLERR,"Fix ave/time compute {} array is accessed out-of-range", val.id);
      if ((val.argindex == 0) && (val.val.c->size_vector_variable)) val.varlen = 1;
      if (val.argindex && (val.val.c->size_array_rows_variable)) val.varlen = 1;

    } else if ((val.which == ArgInfo::FIX) && (mode == SCALAR)) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for fix ave/time does not exist", val.id);
      if ((val.argindex == 0) && (val.val.f->scalar_flag == 0))
        error->all(FLERR,"Fix ave/time fix {} does not calculate a scalar", val.id);
      if (val.argindex && (val.val.f->vector_flag == 0))
        error->all(FLERR,"Fix ave/time fix {} does not calculate a vector", val.id);
      if (val.argindex && (val.val.f->size_vector_variable))
        error->all(FLERR,"Fix ave/time fix {} vector cannot be variable length", val.id);
      if (val.argindex && (val.argindex > val.val.f->size_vector))
        error->all(FLERR,"Fix ave/time fix {} vector is accessed out-of-range", val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for fix ave/time not computed at compatible time", val.id);

    } else if ((val.which == ArgInfo::FIX) && (mode == VECTOR)) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for fix ave/time does not exist", val.id);
      if ((val.argindex == 0) && (val.val.f->vector_flag == 0))
        error->all(FLERR,"Fix ave/time fix {} does not calculate a vector", val.id);
      if (val.argindex && (val.val.f->array_flag == 0))
        error->all(FLERR,"Fix ave/time fix {} does not calculate an array", val.id);
      if (val.argindex && (val.val.f->size_array_rows_variable))
        error->all(FLERR,"Fix ave/time fix {} array cannot be variable length", val.id);
      if (val.argindex && (val.argindex > val.val.f->size_array_cols))
        error->all(FLERR,"Fix ave/time fix {} array is accessed out-of-range", val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for fix ave/time not computed at compatible time", val.id);

    } else if ((val.which == ArgInfo::VARIABLE) && (mode == SCALAR)) {
      int ivariable = input->variable->find(val.id.c_str());
      if (ivariable < 0)
        error->all(FLERR,"Variable name {} for fix ave/time does not exist", val.id);
      if ((val.argindex == 0) && (input->variable->equalstyle(ivariable) == 0))
        error->all(FLERR,"Fix ave/time variable {} is not equal-style variable", val.id);
      if ((val.argindex) && (input->variable->vectorstyle(ivariable) == 0))
        error->all(FLERR,"Fix ave/time variable {} is not vector-style variable", val.id);

    } else if ((val.which == ArgInfo::VARIABLE) && (mode == VECTOR)) {
      int ivariable = input->variable->find(val.id.c_str());
      if (ivariable < 0)
        error->all(FLERR,"Variable name {} for fix ave/time does not exist", val.id);
      if ((val.argindex == 0) && (input->variable->vectorstyle(ivariable) == 0))
        error->all(FLERR,"Fix ave/time variable {} is not vector-style variable", val.id);
      if (val.argindex)
        error->all(FLERR,"Fix ave/time mode vector variable {} cannot be indexed", val.id);
      val.varlen = 1;
    }
  }

  // all_variable_length = 1 if all values are variable length
  // any_variable_length = 1 if any values are variable length

  all_variable_length = 1;
  any_variable_length = 0;
  for (auto &val : values) {
    if (val.varlen == 0) all_variable_length = 0;
    if (val.varlen) any_variable_length = 1;
  }

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array
  // if all columns are variable length, just set nrows = 1 for now

  column = nullptr;
  if (mode == VECTOR) {
    if (all_variable_length == 0) nrows = column_length(0);
    else nrows = 1;
    memory->create(column,nrows,"ave/time:column");
  }

  // enable locking of row count by this fix for computes of variable length
  // only if nrepeat > 1 or ave = RUNNING/WINDOW,
  //   so that locking spans multiple timesteps

  if (any_variable_length && ((nrepeat > 1) || (ave == RUNNING) || (ave == WINDOW))) {
    for (auto &val : values) {
      if (val.varlen && val.which == ArgInfo::COMPUTE) val.val.c->lock_enable();
      lockforever = 0;
    }
  }

  // print file comment lines
  // for mode = VECTOR, cannot use arg to print
  // since array args may have been expanded to multiple vectors

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-averaged data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else if (mode == SCALAR) {
      fprintf(fp,"# TimeStep");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    } else fprintf(fp,"# TimeStep Number-of-rows\n");
    if (title3 && mode == VECTOR) fprintf(fp,"%s\n",title3);
    else if (mode == VECTOR) {
      fprintf(fp,"# Row");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    }
    if (yaml_flag) fputs("---\n",fp);
    if (ferror(fp)) error->one(FLERR,"Error writing file header: {}", utils::getsyserror());
    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // allocate memory for averaging

  vector = vector_total = nullptr;
  vector_list = nullptr;
  array = array_total = nullptr;
  array_list = nullptr;

  if (mode == SCALAR) {
    vector = new double[nvalues];
    vector_total = new double[nvalues];
    if (ave == WINDOW) memory->create(vector_list,nwindow,nvalues,"ave/time:vector_list");
  } else allocate_arrays();

  // this fix produces either a global scalar or vector or array
  // SCALAR mode produces either a scalar or vector
  // VECTOR mode produces either a vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value

  extlist = nullptr;

  if (mode == SCALAR) {
    if (nvalues == 1) {
      scalar_flag = 1;
      auto &val = values[0];
      if (val.which == ArgInfo::COMPUTE) {
        if (val.argindex == 0) extscalar = val.val.c->extscalar;
        else if (val.val.c->extvector >= 0) extscalar = val.val.c->extvector;
        else extscalar = val.val.c->extlist[val.argindex-1];
      } else if (val.which == ArgInfo::FIX) {
        if (val.argindex == 0) extscalar = val.val.f->extscalar;
        else if (val.val.f->extvector >= 0) extscalar = val.val.f->extvector;
        else extscalar = val.val.f->extlist[val.argindex-1];
      } else if (val.which == ArgInfo::VARIABLE) {
        extscalar = 0;
      }

    } else {
      vector_flag = 1;
      size_vector = nrows = nvalues;
      extvector = -1;
      extlist = new int[nvalues];
      int i = 0;
      for (auto &val : values) {
        if (val.which == ArgInfo::COMPUTE) {
          if (val.argindex == 0) extlist[i] = val.val.c->extscalar;
          else if (val.val.c->extvector >= 0) extlist[i] = val.val.c->extvector;
          else extlist[i] = val.val.c->extlist[val.argindex-1];
        } else if (val.which == ArgInfo::FIX) {
          if (val.argindex == 0) extlist[i] = val.val.f->extscalar;
          else if (val.val.f->extvector >= 0) extlist[i] = val.val.f->extvector;
          else extlist[i] = val.val.f->extlist[val.argindex-1];
        } else if (val.which == ArgInfo::VARIABLE) {
          extlist[i] = 0;
        }
        ++i;
      }
    }

  } else {
    if (nvalues == 1) {
      auto &val = values[0];
      vector_flag = 1;
      size_vector = nrows;
      if (all_variable_length) size_vector_variable = 1;
      if (val.which == ArgInfo::COMPUTE) {
        if (val.argindex == 0) {
          extvector = val.val.c->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = val.val.c->extlist[i];
          }
        } else extvector = val.val.c->extarray;
      } else if (val.which == ArgInfo::FIX) {
        if (val.argindex == 0) {
          extvector = val.val.f->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = val.val.f->extlist[i];
          }
        } else extvector = val.val.f->extarray;
      } else if (val.which == ArgInfo::VARIABLE) {
        extlist = new int[nrows];
        for (int i = 0; i < nrows; i++) extlist[i] = 0;
      }

    } else {
      array_flag = 1;
      size_array_rows = nrows;
      size_array_cols = nvalues;
      if (all_variable_length) size_array_rows_variable = 1;
      int extvalue = 0;
      extarray = -2;
      for (auto &val : values) {
        if (val.which == ArgInfo::COMPUTE) {
          if (val.argindex == 0) extvalue = val.val.c->extvector;
          else extvalue = val.val.c->extarray;
        } else if (val.which == ArgInfo::FIX) {
          if (val.argindex == 0) extvalue = val.val.f->extvector;
          else extvalue = val.val.f->extarray;
        } else if (val.which == ArgInfo::VARIABLE) {
          extvalue = 0;
        }
        if (extvalue == -1)
          error->all(FLERR,"Fix ave/time cannot set output array intensive/extensive "
                     "from these inputs");
        if (extarray < -1) extarray = extvalue;
        else if (extvalue != extarray)
          error->all(FLERR,"Fix ave/time cannot set output array intensive/extensive "
                     "from these inputs");
      }
    }
  }

  // initializations
  // set vector_total to zero since it accumulates
  // array_total already zeroed in allocate_arrays

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  if (mode == SCALAR)
    for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  // decrement lock counter in compute chunk/atom, it if still exists

  if (any_variable_length && ((nrepeat > 1) || (ave == RUNNING) || (ave == WINDOW))) {
    for (auto &val : values) {
      if (val.varlen) {
        auto icompute = modify->get_compute_by_id(val.id);
        if (icompute) {
          if ((ave == RUNNING) || (ave == WINDOW))
            icompute->unlock(this);
          icompute->lock_disable();
        }
      }
    }
  }

  delete[] format_user;
  delete[] extlist;

  if (fp && comm->me == 0) {
    if (yaml_flag) fputs("...\n", fp);
    fclose(fp);
  }
  memory->destroy(column);

  delete[] vector;
  delete[] vector_total;
  memory->destroy(array);
  memory->destroy(array_total);
  memory->destroy(array_list);
}

/* ---------------------------------------------------------------------- */

int FixAveTime::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveTime::init()
{
  // update indices/pointers for all computes,fixes,variables

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for fix ave/time does not exist", val.id);
    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for fix ave/time does not exist", val.id);
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/time does not exist", val.id);
    }
  }

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

void FixAveTime::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  if (mode == SCALAR) invoke_scalar(ntimestep);
  else invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_scalar(bigint ntimestep)
{
  // zero if first sample within single Nfreq epoch
  // if any input is variable length, initialize current length
  // check for exceeding length is done below

  if (irepeat == 0) {
    if (any_variable_length) {
      modify->clearstep_compute();
      column_length(1);
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);
    }
    for (int i = 0; i < nvalues; i++) vector[i] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  int i = 0;
  double scalar = 0.0;
  for (auto &val : values) {

    // invoke compute if not previously invoked
    // ensure no out-of-range access to variable-length compute vector

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
        if (val.varlen && (val.val.c->size_vector < val.argindex)) scalar = 0.0;
        else scalar = val.val.c->vector[val.argindex-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0)
        scalar = val.val.f->compute_scalar();
      else
        scalar = val.val.f->compute_vector(val.argindex-1);

    // evaluate equal-style or vector-style variable
    // ensure no out-of-range access to vector-style variable

    } else if (val.which == ArgInfo::VARIABLE) {
      if (val.argindex == 0)
        scalar = input->variable->compute_equal(val.val.v);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(val.val.v,&varvec);
        if (nvec < val.argindex) scalar = 0.0;
        else scalar = varvec[val.argindex-1];
      }
    }

    // add value to vector or just set directly if offcol is set

    if (val.offcol) vector[i] = scalar;
    else vector[i] += scalar;
    ++i;
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep + nfreq - static_cast<bigint>(nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nvalues; i++)
    if (values[i].offcol == 0) vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (i = 0; i < nvalues; i++) vector_total[i] = vector[i];
    norm = 1;

  } else if (ave == RUNNING) {
    for (i = 0; i < nvalues; i++) vector_total[i] += vector[i];
    norm++;

  } else if (ave == WINDOW) {
    for (i = 0; i < nvalues; i++) {
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

  // ensure any columns with offcol set are effectively set to last value

  for (i = 0; i < nvalues; i++)
    if (values[i].offcol) vector_total[i] = norm*vector[i];

  // output result to file

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    if (yaml_flag) {
      if (!yaml_header || overwrite) {
        yaml_header = true;
        fputs("keywords: ['Step', ", fp);
        for (const auto &val : values) fmt::print(fp, "'{}', ", val.keyword);
        fputs("]\ndata:\n", fp);
      }
      fmt::print(fp, "  - [{}, ", ntimestep);
      for (i = 0; i < nvalues; i++) fmt::print(fp,"{}, ",vector_total[i]/norm);
      fputs("]\n", fp);
    } else {
      fmt::print(fp,"{}",ntimestep);
      for (i = 0; i < nvalues; i++) fprintf(fp,format,vector_total[i]/norm);
      fprintf(fp,"\n");
      if (ferror(fp)) error->one(FLERR,"Error writing out time averaged data");
    }
    fflush(fp);

    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}", utils::getsyserror());
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_vector(bigint ntimestep)
{
  // first sample within single Nfreq epoch
  // zero out arrays that accumulate over many samples, but not across epochs
  // invoke setup_chunks() to determine current nchunk
  //   re-allocate per-chunk arrays if needed
  // invoke lock() in two cases:
  //   if nrepeat > 1: so nchunk cannot change until Nfreq epoch is over,
  //     will be unlocked on last repeat of this Nfreq
  //   if ave = RUNNING/WINDOW and not yet locked:
  //     set forever, will be unlocked in fix destructor
  // wrap setup_chunks in clearstep/addstep b/c it may invoke computes
  //   both nevery and nfreq are future steps,
  //   since call below to cchunk->ichunk()
  //     does not re-invoke internal cchunk compute on this same step

  if (irepeat == 0) {
    if (any_variable_length) {
      modify->clearstep_compute();
      int nrows_new = column_length(1);
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);

      if (all_variable_length && nrows_new != nrows) {
        nrows = nrows_new;
        memory->destroy(column);
        memory->create(column,nrows,"ave/time:column");
        allocate_arrays();
      }

      int lockforever_flag = 0;
      for (auto &val : values) {
        if (!val.varlen || (val.which != ArgInfo::COMPUTE)) continue;
        if ((nrepeat > 1) && (ave == ONE)) {
          val.val.c->lock(this,ntimestep,ntimestep+static_cast<bigint>(nrepeat-1)*nevery);
        } else if (((ave == RUNNING) || (ave == WINDOW)) && !lockforever) {
          val.val.c->lock(this,update->ntimestep,-1);
          lockforever_flag = 1;
        }
      }
      if (lockforever_flag) lockforever = 1;
    }

    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++) array[i][j] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  int j = 0;
  for (auto &val : values) {

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {
      if (val.argindex == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        double *cvector = val.val.c->vector;
        for (int i = 0; i < nrows; i++)
          column[i] = cvector[i];

      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
          val.val.c->compute_array();
          val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
        }
        double **carray = val.val.c->array;
        int icol = val.argindex-1;
        for (int i = 0; i < nrows; i++)
          column[i] = carray[i][icol];
      }

    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0)
        for (int i = 0; i < nrows; i++)
          column[i] = val.val.f->compute_vector(i);
      else {
        int icol = val.argindex-1;
        for (int i = 0; i < nrows; i++)
          column[i] = val.val.f->compute_array(i,icol);
      }

    // evaluate vector-style variable
    // ensure nvec = nrows, else error
    // could be different on this timestep than when column_length(1) set nrows

    } else if (val.which == ArgInfo::VARIABLE) {
      double *varvec;
      int nvec = input->variable->compute_vector(val.val.v,&varvec);
      if (nvec != nrows)
        error->all(FLERR,"Fix ave/time vector-style variable {} changed length", val.id);
      for (int i = 0; i < nrows; i++)
        column[i] = varvec[i];
    }

    // add columns of values to array or just set directly if offcol is set

    if (val.offcol) {
      for (int i = 0; i < nrows; i++)
        array[i][j] = column[i];
    } else {
      for (int i = 0; i < nrows; i++)
        array[i][j] += column[i];
    }
    ++j;
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - static_cast<bigint>(nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // unlock any variable length computes at end of Nfreq epoch
  // do not unlock if ave = RUNNING or WINDOW

  if (any_variable_length && (nrepeat > 1) && (ave == ONE)) {
    for (auto &val : values) {
      if (!val.varlen) continue;
      if ((val.which == ArgInfo::COMPUTE) && val.val.c) val.val.c->unlock(this);
    }
  }

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++)
      if (values[j].offcol == 0) array[i][j] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++) array_total[i][j] = array[i][j];
    norm = 1;

  } else if (ave == RUNNING) {
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++) array_total[i][j] += array[i][j];
    norm++;

  } else if (ave == WINDOW) {
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++) {
        array_total[i][j] += array[i][j];
        if (window_limit) array_total[i][j] -= array_list[iwindow][i][j];
        array_list[iwindow][i][j] = array[i][j];
      }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // ensure any columns with offcol set are effectively set to last value

  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++)
      if (values[j].offcol) array_total[i][j] = norm*array[i][j];

  // output result to file

  if (fp && comm->me == 0) {
    if (overwrite) platform::fseek(fp,filepos);
    if (yaml_flag) {
      if (!yaml_header || overwrite) {
        yaml_header = true;
        fputs("keywords: [", fp);
        for (const auto &val : values) fmt::print(fp, "'{}', ", val.keyword);
        fputs("]\ndata:\n", fp);
      }
      fmt::print(fp, "  {}:\n", ntimestep);
      for (int i = 0; i < nrows; i++) {
        fputs("  - [", fp);
        for (int j = 0; j < nvalues; j++) fmt::print(fp,"{}, ",array_total[i][j]/norm);
        fputs("]\n", fp);
      }
    } else {
      fmt::print(fp,"{} {}\n",ntimestep,nrows);
      for (int i = 0; i < nrows; i++) {
        fprintf(fp,"%d",i+1);
        for (int j = 0; j < nvalues; j++) fprintf(fp,format,array_total[i][j]/norm);
        fprintf(fp,"\n");
      }
    }
    fflush(fp);
    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}", utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

int FixAveTime::column_length(int dynamic)
{
  int length,lengthone;

  // determine nrows for static values

  if (!dynamic) {
    length = 0;
    for (auto &val : values) {
      if (val.varlen) continue;
      if (val.which == ArgInfo::COMPUTE) {
        if (val.argindex == 0)
          lengthone = val.val.c->size_vector;
        else lengthone = val.val.c->size_array_rows;
      } else if (val.which == ArgInfo::FIX) {
        if (val.argindex == 0) lengthone = val.val.f->size_vector;
        else lengthone = val.val.f->size_array_rows;
      } else if (val.which == ArgInfo::VARIABLE) {
        // variables are always varlen = 1, so dynamic
      }
      if (length == 0) length = lengthone;
      else if (lengthone != length)
        error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
    }
  }

  // determine new nrows for dynamic values
  // either all must be the same
  // or must match other static values
  // don't need to check if not MODE = VECTOR, just invoke lock_length()

  if (dynamic) {
    length = 0;
    for (auto &val : values) {
      if (val.varlen == 0) continue;
      if (val.which == ArgInfo::COMPUTE) {
        lengthone = val.val.c->lock_length();
      } else if (val.which == ArgInfo::VARIABLE) {
        double *varvec;
        lengthone = input->variable->compute_vector(val.val.v,&varvec);
      }
      if (mode == SCALAR) continue;
      if (all_variable_length) {
        if (length == 0) length = lengthone;
        else if (lengthone != length)
          error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
      } else {
        if (lengthone != nrows)
          error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
      }
    }
  }

  return length;
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveTime::compute_scalar()
{
  if (norm) return vector_total[0]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveTime::compute_vector(int i)
{
  if (i >= nrows) return 0.0;
  if (norm) {
    if (mode == SCALAR) return vector_total[i]/norm;
    if (mode == VECTOR) return array_total[i][0]/norm;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveTime::compute_array(int i, int j)
{
  if (i >= nrows) return 0.0;
  if (norm) return array_total[i][j]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   modify settings
------------------------------------------------------------------------- */

int FixAveTime::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "colname") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "fix_modify colname", error);
    int icol = -1;
    if (utils::is_integer(arg[1])) {
      icol = utils::inumeric(FLERR, arg[1], false, lmp);
      if (icol < 0) icol = values.size() + icol + 1;
      icol--;
    } else {
      try {
        icol = key2col.at(arg[1]);
      } catch (std::out_of_range &) {
        icol = -1;
      }
    }
    if ((icol < 0) || (icol >= (int) values.size()))
      error->all(FLERR, "Thermo_modify colname column {} invalid", arg[1]);
    values[icol].keyword = arg[2];
    return 3;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveTime::options(int iarg, int narg, char **arg)
{
  // option defaults

  fp = nullptr;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  noff = 0;
  offlist = nullptr;
  overwrite = 0;
  yaml_flag = yaml_header = false;
  format_user = nullptr;
  format = (char *) " %g";
  title1 = nullptr;
  title2 = nullptr;
  title3 = nullptr;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      yaml_flag = utils::strmatch(arg[iarg+1],"\\.[yY][aA]?[mM][lL]$");
      if (comm->me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open fix ave/time file {}: {}",
                     arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/time command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/time command");
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/time command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"off") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      memory->grow(offlist,noff+1,"ave/time:offlist");
      offlist[noff++] = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      delete[] format_user;
      format_user = utils::strdup(arg[iarg+1]);
      format = format_user;
      iarg += 2;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      delete[] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      delete[] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      delete[] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix ave/time command option {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays for mode = VECTOR of size Nrows x Nvalues
------------------------------------------------------------------------- */

void FixAveTime::allocate_arrays()
{
  memory->destroy(array);
  memory->destroy(array_total);
  memory->create(array,nrows,nvalues,"ave/time:array");
  memory->create(array_total,nrows,nvalues,"ave/time:array_total");
  if (ave == WINDOW) {
    memory->destroy(array_list);
    memory->create(array_list,nwindow,nrows,nvalues,"ave/time:array_list");
  }

  // reinitialize regrown array_total since it accumulates

  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++) array_total[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveTime::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= static_cast<bigint>(nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}
