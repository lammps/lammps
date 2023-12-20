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
   Contributing author: Shawn Coleman (ARL)
------------------------------------------------------------------------- */
#include "fix_ave_histo_weight.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { ONE, RUNNING };
enum { SCALAR, VECTOR, WINDOW };
enum { DEFAULT, GLOBAL, PERATOM, LOCAL };
enum { IGNORE, END, EXTRA };
enum { SINGLE, VALUE };

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixAveHistoWeight::FixAveHistoWeight(LAMMPS *lmp, int narg, char **arg) :
    FixAveHisto(lmp, narg, arg)
{
  // nvalues = 2 required for histo/weight

  if (nvalues != 2)
    error->all(FLERR, "Illegal fix ave/histo/weight command: must have two data arguments");

  // check that length of 2 values is the same

  int size[2] = {0, 0};

  for (int i = 0; i < nvalues; i++) {
    auto &val = values[i];
    if (val.which == ArgInfo::X || val.which == ArgInfo::V || val.which == ArgInfo::F) {
      size[i] = atom->nlocal;
    } else if (val.which == ArgInfo::COMPUTE && kind == GLOBAL && mode == SCALAR) {
      size[i] = val.val.c->size_vector;
    } else if (val.which == ArgInfo::COMPUTE && kind == GLOBAL && mode == VECTOR) {
      size[i] = val.val.c->size_array_rows;
    } else if (val.which == ArgInfo::COMPUTE && kind == PERATOM) {
      size[i] = atom->nlocal;
    } else if (val.which == ArgInfo::COMPUTE && kind == LOCAL) {
      size[i] = val.val.c->size_local_rows;
    } else if (val.which == ArgInfo::FIX && kind == GLOBAL && mode == SCALAR) {
      size[i] = val.val.f->size_vector;
    } else if (val.which == ArgInfo::FIX && kind == GLOBAL && mode == VECTOR) {
      size[i]= val.val.f->size_array_rows;
    } else if (val.which == ArgInfo::FIX && kind == PERATOM) {
      size[i] = atom->nlocal;
    } else if (val.which == ArgInfo::FIX && kind == LOCAL) {
      size[i] = val.val.f->size_local_rows;
    } else if (val.which == ArgInfo::VARIABLE && kind == PERATOM) {
      size[i] = atom->nlocal;
    }
  }

  if (size[0] != size[1])
    error->all(FLERR,"Fix ave/histo/weight value and weight vector lengths do not match");
}

/* ---------------------------------------------------------------------- */

void FixAveHistoWeight::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0) {
    stats[0] = stats[1] = 0.0;
    stats[2] = BIG;
    stats[3] = -BIG;
    for (int i = 0; i < nbins; i++) bin[i] = 0.0;
  }

  // first calculate weight factors, then bin single value
  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // calculate weight factors which are 2nd value (i = 1)

  double weight = 0.0;
  double *weights = nullptr;
  int stride = 0;
  auto &val = values[1];
  int j = val.argindex;

  // atom attributes

  if (val.which == ArgInfo::X) {
    weights = &atom->x[0][j];
    stride = 3;
  } else if (val.which == ArgInfo::V) {
    weights = &atom->v[0][j];
    stride = 3;
    bin_atoms(&atom->v[0][j],3);
  } else if (val.which == ArgInfo::F) {
    weights = &atom->f[0][j];
    stride = 3;
  }

  // invoke compute if not previously invoked

  if (val.which == ArgInfo::COMPUTE) {

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
          val.val.c->compute_scalar();
          val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        weight = val.val.c->scalar;
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        weight = val.val.c->vector[j-1];
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        weights = val.val.c->vector;
        stride = 1;
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
          val.val.c->compute_array();
          val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
        }
        if (val.val.c->array) weights = &val.val.c->array[0][j-1];
        stride = val.val.c->size_array_cols;
      }
    } else if (kind == PERATOM) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
        val.val.c->compute_peratom();
        val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
      }
      if (j == 0) {
        weights = val.val.c->vector_atom;
        stride = 1;
      } else if (val.val.c->array_atom) {
        weights = &val.val.c->array_atom[0][j-1];
        stride = val.val.c->size_peratom_cols;
      }
    } else if (kind == LOCAL) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_LOCAL)) {
        val.val.c->compute_local();
        val.val.c->invoked_flag |= Compute::INVOKED_LOCAL;
      }
      if (j == 0) {
        weights = val.val.c->vector_local;
        stride = 1;
      } else if (val.val.c->array_local) {
        weights = &val.val.c->array_local[0][j-1];
        stride = val.val.c->size_local_cols;
      }
    }

  // access fix fields, guaranteed to be ready

  } else if (val.which == ArgInfo::FIX) {

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) weight = val.val.f->compute_scalar();
      else weight = val.val.f->compute_vector(j-1);
    } else if (kind == GLOBAL && mode == VECTOR) {
      error->all(FLERR,"Fix ave/histo/weight option not yet supported");
      // NOTE: need to allocate local storage
      if (j == 0) {
        int n = val.val.f->size_vector;
        for (int i = 0; i < n; i++) weights[n] = val.val.f->compute_vector(i);
      } else {
        int n = val.val.f->size_vector;
        for (int i = 0; i < n; i++) weights[n] = val.val.f->compute_array(i,j-1);
      }
    } else if (kind == PERATOM) {
      if (j == 0) {
        weights = val.val.f->vector_atom;
        stride = 1;
      } else if (val.val.f->array_atom) {
        weights = &val.val.f->array_atom[0][j-1];
        stride = val.val.f->size_peratom_cols;
      }
    } else if (kind == LOCAL) {
      if (j == 0) {
        weights = val.val.f->vector_local;
        stride = 1;
      } else if (val.val.f->array_local) {
        weights = &val.val.f->array_local[0][j-1];
        stride = val.val.f->size_local_cols;
      }
    }

  // evaluate equal-style variable

  } else if (val.which == ArgInfo::VARIABLE && kind == GLOBAL) {
    weight = input->variable->compute_equal(val.val.v);

  } else if (val.which == ArgInfo::VARIABLE && kind == PERATOM) {
    if (atom->nmax > maxatom) {
      memory->destroy(vector);
      maxatom = atom->nmax;
      memory->create(vector,maxatom,"ave/histo/weight:vector");
    }
    input->variable->compute_atom(val.val.v,igroup,vector,1,0);
    weights = vector;
    stride = 1;
  }

  // bin values using weights, values are 1st value (i = 0)

  val = values[0];
  j = val.argindex;

  // atom attributes

  if (val.which == ArgInfo::X && weights != nullptr)
    bin_atoms_weights(&atom->x[0][j],3,weights,stride);
  else if (val.which == ArgInfo::V && weights != nullptr)
    bin_atoms_weights(&atom->v[0][j],3,weights,stride);
  else if (val.which == ArgInfo::F && weights != nullptr)
    bin_atoms_weights(&atom->f[0][j],3,weights,stride);

  // invoke compute if not previously invoked

  if (val.which == ArgInfo::COMPUTE) {

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
          val.val.c->compute_scalar();
          val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        bin_one_weights(val.val.c->scalar,weight);
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        bin_one_weights(val.val.c->vector[j-1],weight);
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        bin_vector_weights(val.val.c->size_vector,val.val.c->vector,1,
                           weights,stride);
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
          val.val.c->compute_array();
          val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
        }
        if (val.val.c->array)
          bin_vector_weights(val.val.c->size_array_rows,&val.val.c->array[0][j-1],
                             val.val.c->size_array_cols,weights,stride);
      }

    } else if (kind == PERATOM) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
        val.val.c->compute_peratom();
        val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
      }
      if (j == 0)
        bin_atoms_weights(val.val.c->vector_atom,1,weights, stride);
      else if (val.val.c->array_atom)
        bin_atoms_weights(&val.val.c->array_atom[0][j-1],
                          val.val.c->size_peratom_cols,weights,stride);

    } else if (kind == LOCAL) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_LOCAL)) {
        val.val.c->compute_local();
        val.val.c->invoked_flag |= Compute::INVOKED_LOCAL;
      }
      if (j == 0)
        bin_vector_weights(val.val.c->size_local_rows,
                           val.val.c->vector_local,1,weights,stride);
      else if (val.val.c->array_local)
        bin_vector_weights(val.val.c->size_local_rows,
                           &val.val.c->array_local[0][j-1],
                           val.val.c->size_local_cols,weights,stride);
    }

    // access fix fields, guaranteed to be ready

  } else if (val.which == ArgInfo::FIX) {

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) bin_one_weights(val.val.f->compute_scalar(),weight);
      else bin_one_weights(val.val.f->compute_vector(j-1),weight);

    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        int n = val.val.f->size_vector;
        for (int i = 0; i < n; i++)
          bin_one_weights(val.val.f->compute_vector(i),weights[i*stride]);
      } else {
        int n = val.val.f->size_vector;
        for (int i = 0; i < n; i++)
          bin_one_weights(val.val.f->compute_array(i,j-1),weights[i*stride]);
      }

    } else if (kind == PERATOM) {
      if (j == 0)
        bin_atoms_weights(val.val.f->vector_atom,1,weights,stride);
      else if (val.val.f->array_atom)
        bin_atoms_weights(&val.val.f->array_atom[0][j-1],val.val.f->size_peratom_cols,
                          weights,stride);


    } else if (kind == LOCAL) {
      if (j == 0) bin_vector_weights(val.val.f->size_local_rows,val.val.f->vector_local,1,
                                     weights,stride);
      else if (val.val.f->array_local)
        bin_vector_weights(val.val.f->size_local_rows,&val.val.f->array_local[0][j-1],
                           val.val.f->size_local_cols,weights,stride);
    }

    // evaluate equal-style variable

  } else if (val.which == ArgInfo::VARIABLE && kind == GLOBAL) {
    bin_one_weights(input->variable->compute_equal(val.val.v),weight);

  } else if (val.which == ArgInfo::VARIABLE && kind == PERATOM) {
    if (atom->nmax > maxatom) {
      memory->destroy(vector);
      maxatom = atom->nmax;
      memory->create(vector,maxatom,"ave/histo/weight:vector");
    }
    input->variable->compute_atom(val.val.v,igroup,vector,1,0);
    bin_atoms_weights(vector,1,weights,stride);
  }

  // code beyond this point is identical to FixAveHisto

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

  // merge histogram stats across procs if necessary

  if (kind == PERATOM || kind == LOCAL) {
    MPI_Allreduce(stats,stats_all,2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&stats[2],&stats_all[2],1,MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(&stats[3],&stats_all[3],1,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(bin,bin_all,nbins,MPI_DOUBLE,MPI_SUM,world);

    stats[0] = stats_all[0];
    stats[1] = stats_all[1];
    stats[2] = stats_all[2];
    stats[3] = stats_all[3];
    for (int i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (int i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (int i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
    stats_total[0] += stats[0];
    if (window_limit) stats_total[0] -= stats_list[iwindow][0];
    stats_list[iwindow][0] = stats[0];
    stats_total[1] += stats[1];
    if (window_limit) stats_total[1] -= stats_list[iwindow][1];
    stats_list[iwindow][1] = stats[1];

    int m;
    if (window_limit) m = nwindow;
    else m = iwindow+1;

    stats_list[iwindow][2] = stats[2];
    stats_total[2] = stats_list[0][2];
    for (int i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (int i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (int i = 0; i < nbins; i++) {
      bin_total[i] += bin[i];
      if (window_limit) bin_total[i] -= bin_list[iwindow][i];
      bin_list[iwindow][i] = bin[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
  }

  // output result to file

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    fmt::print(fp,"{} {} {} {} {} {}\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0)
      for (int i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
    else
      for (int i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);

    if (ferror(fp))
      error->one(FLERR,"Error writing out histogram data");

    fflush(fp);
    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}", utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   bin a single value with weight)
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_one_weights(double value, double weight)
{
  stats[2] = MIN(stats[2],value);
  stats[3] = MAX(stats[3],value);

  if (value < lo) {
    if (beyond == IGNORE) {
      stats[1] += weight;
      return;
    } else bin[0] += weight;
  } else if (value > hi) {
    if (beyond == IGNORE) {
      stats[1] += weight;
      return;
    } else bin[nbins-1] += weight;
  } else {
    int ibin = static_cast<int> ((value-lo)*bininv);
    ibin = MIN(ibin,nbins-1);
    if (beyond == EXTRA) ibin++;
    bin[ibin] += weight;
  }

  stats[0] += weight;
}

/* ----------------------------------------------------------------------
   bin a vector of values with weights
   values and weights each have a stride
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_vector_weights(int n, double *values,
                                           int stride, double *weights,
                                           int stridewt)
{
  int m = 0;
  int m2 = 0;
  for (int i = 0; i < n; i++) {
    bin_one_weights(values[m],weights[m2]);
    m += stride;
    m2 += stridewt;
  }
}

/* ----------------------------------------------------------------------
   bin a per-atom vector of values with weights
   values and weights each have a stride
   only bin if atom is in group
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_atoms_weights(double *values, int stride,
                                          double *weights, int stridewt)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  int m2 = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) bin_one_weights(values[m],weights[m2]);
    m += stride;
    m2 += stridewt;
  }
}
