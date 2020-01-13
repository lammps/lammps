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
   Contributing author: Shawn Coleman (ARL)
------------------------------------------------------------------------- */
#include "fix_ave_histo_weight.h"
#include <mpi.h>
#include <unistd.h>
#include "fix.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};
enum{SCALAR,VECTOR,WINDOW};
enum{DEFAULT,GLOBAL,PERATOM,LOCAL};
enum{IGNORE,END,EXTRA};
enum{SINGLE,VALUE};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixAveHistoWeight::FixAveHistoWeight(LAMMPS *lmp, int narg, char **arg) :
  FixAveHisto(lmp, narg, arg)
{
  // nvalues = 2 required for histo/weight

  if (nvalues != 2) error->all(FLERR,"Illegal fix ave/histo/weight command");

  // check that length of 2 values is the same

  int size[2];

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == X || which[i] == V || which[i] == F) {
      size[i] = atom->nlocal;
    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      size[i] = modify->compute[icompute]->size_vector;
    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      size[i] = modify->compute[icompute]->size_array_rows;
    } else if (which[i] == COMPUTE && kind == PERATOM) {
      size[i] = atom->nlocal;
    } else if (which[i] == COMPUTE && kind == LOCAL) {
      int icompute = modify->find_compute(ids[i]);
      size[i] = modify->compute[icompute]->size_local_rows;
    } else if (which[i] == FIX && kind == GLOBAL && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      size[i] = modify->fix[ifix]->size_vector;
    } else if (which[i] == FIX && kind == GLOBAL && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      size[i]= modify->fix[ifix]->size_array_rows;
    } else if (which[i] == FIX && kind == PERATOM) {
      size[i] = atom->nlocal;
    } else if (which[i] == FIX && kind == LOCAL) {
      int ifix = modify->find_fix(ids[i]);
      size[i] = modify->fix[ifix]->size_local_rows;
    } else if (which[i] == VARIABLE && kind == PERATOM) {
      size[i] = atom->nlocal;
    }
  }

  if (size[0] != size[1])
    error->all(FLERR,"Fix ave/histo/weight value and weight vector "
               "lengths do not match");
}

/* ---------------------------------------------------------------------- */

void FixAveHistoWeight::end_of_step()
{
  int i,j,m;

  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/histo");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0) {
    stats[0] = stats[1] = 0.0;
    stats[2] = BIG;
    stats[3] = -BIG;
    for (i = 0; i < nbins; i++) bin[i] = 0.0;
  }

  // first calculate weight factors, then bin single value
  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // calculate weight factors which are 2nd value (i = 1)

  double weight = 0.0;
  double *weights = NULL;
  int stride = 0;
  i = 1;

  m = value2index[i];
  j = argindex[i];

  // atom attributes

  if (which[i] == X) {
    weights = &atom->x[0][j];
    stride = 3;
  } else if (which[i] == V){
    weights = &atom->v[0][j];
    stride = 3;
    bin_atoms(&atom->v[0][j],3);
  } else if (which[i] == F) {
    weights = &atom->f[0][j];
    stride = 3;
  }

  // invoke compute if not previously invoked

  if (which[i] == COMPUTE) {

    Compute *compute = modify->compute[m];

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        weight = compute->scalar;
      } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        weight = compute->vector[j-1];
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        weights = compute->vector;
        stride = 1;
      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        if (compute->array) weights = &compute->array[0][j-1];
        stride = compute->size_array_cols;
      }
    } else if (kind == PERATOM) {
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }
      if (j == 0) {
        weights = compute->vector_atom;
        stride = 1;
      } else if (compute->array_atom) {
        weights = &compute->array_atom[0][j-1];
        stride = compute->size_peratom_cols;
      }
    } else if (kind == LOCAL) {
      if (!(compute->invoked_flag & INVOKED_LOCAL)) {
        compute->compute_local();
        compute->invoked_flag |= INVOKED_LOCAL;
      }
      if (j == 0) {
        weights = compute->vector_local;
        stride = 1;
      } else if (compute->array_local) {
        weights = &compute->array_local[0][j-1];
        stride = compute->size_local_cols;
      }
    }

  // access fix fields, guaranteed to be ready

  } else if (which[i] == FIX) {

    Fix *fix = modify->fix[m];

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) weight = fix->compute_scalar();
      else weight = fix->compute_vector(j-1);
    } else if (kind == GLOBAL && mode == VECTOR) {
      error->all(FLERR,"Fix ave/histo/weight option not yet supported");
      // NOTE: need to allocate local storage
      if (j == 0) {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_vector(i);
      } else {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_array(i,j-1);
      }
    } else if (kind == PERATOM) {
      if (j == 0) {
        weights = fix->vector_atom;
        stride = 1;
      } else if (fix->array_atom) {
        weights = fix->array_atom[j-1];
        stride = fix->size_peratom_cols;
      }
    } else if (kind == LOCAL) {
      if (j == 0) {
        weights = fix->vector_local;
        stride = 1;
      } else if (fix->array_local) {
        weights = &fix->array_local[0][j-1];
        stride = fix->size_local_cols;
      }
    }

  // evaluate equal-style variable

  } else if (which[i] == VARIABLE && kind == GLOBAL) {
    weight = input->variable->compute_equal(m);

  } else if (which[i] == VARIABLE && kind == PERATOM) {
    if (atom->nmax > maxatom) {
      memory->destroy(vector);
      maxatom = atom->nmax;
      memory->create(vector,maxatom,"ave/histo/weight:vector");
    }
    input->variable->compute_atom(m,igroup,vector,1,0);
    weights = vector;
    stride = 1;
  }

  // bin values using weights, values are 1st value (i = 0)

  i = 0;
  m = value2index[i];
  j = argindex[i];

  // atom attributes

  if (which[i] == X && weights != NULL)
    bin_atoms_weights(&atom->x[0][j],3,weights,stride);
  else if (which[i] == V && weights != NULL)
    bin_atoms_weights(&atom->v[0][j],3,weights,stride);
  else if (which[i] == F && weights != NULL)
    bin_atoms_weights(&atom->f[0][j],3,weights,stride);

  // invoke compute if not previously invoked

  if (which[i] == COMPUTE) {
    Compute *compute = modify->compute[m];
    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        bin_one_weights(compute->scalar,weight);
      } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        bin_one_weights(compute->vector[j-1],weight);
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        bin_vector_weights(compute->size_vector,compute->vector,1,
                           weights,stride);
      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        if (compute->array)
          bin_vector_weights(compute->size_array_rows,&compute->array[0][j-1],
                             compute->size_array_cols,weights,stride);
      }

    } else if (kind == PERATOM) {
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }
      if (j == 0)
        bin_atoms_weights(compute->vector_atom,1,weights, stride);
      else if (compute->array_atom)
        bin_atoms_weights(&compute->array_atom[0][j-1],
                          compute->size_peratom_cols,weights,stride);

    } else if (kind == LOCAL) {
      if (!(compute->invoked_flag & INVOKED_LOCAL)) {
        compute->compute_local();
        compute->invoked_flag |= INVOKED_LOCAL;
      }
      if (j == 0)
        bin_vector_weights(compute->size_local_rows,
                           compute->vector_local,1,weights,stride);
      else if (compute->array_local)
        bin_vector_weights(compute->size_local_rows,
                           &compute->array_local[0][j-1],
                           compute->size_local_cols,weights,stride);
    }

    // access fix fields, guaranteed to be ready

  } else if (which[i] == FIX) {

    Fix *fix = modify->fix[m];

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) bin_one_weights(fix->compute_scalar(),weight);
      else bin_one_weights(fix->compute_vector(j-1),weight);

    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        int n = fix->size_vector;
        for (i = 0; i < n; i++)
          bin_one_weights(fix->compute_vector(i),weights[i*stride]);
      } else {
        int n = fix->size_vector;
        for (i = 0; i < n; i++)
          bin_one_weights(fix->compute_array(i,j-1),weights[i*stride]);
      }

    } else if (kind == PERATOM) {
      if (j == 0)
        bin_atoms_weights(fix->vector_atom,1,weights,stride);
      else if (fix->array_atom)
        bin_atoms_weights(fix->array_atom[j-1],fix->size_peratom_cols,
                          weights,stride);


    } else if (kind == LOCAL) {
      if (j == 0) bin_vector_weights(fix->size_local_rows,fix->vector_local,1,
                                     weights,stride);
      else if (fix->array_local)
        bin_vector_weights(fix->size_local_rows,&fix->array_local[0][j-1],
                           fix->size_local_cols,weights,stride);
    }

    // evaluate equal-style variable

  } else if (which[i] == VARIABLE && kind == GLOBAL) {
    bin_one_weights(input->variable->compute_equal(m),weight);

  } else if (which[i] == VARIABLE && kind == PERATOM) {
    if (atom->nmax > maxatom) {
      memory->destroy(vector);
      maxatom = atom->nmax;
      memory->create(vector,maxatom,"ave/histo/weight:vector");
    }
    input->variable->compute_atom(m,igroup,vector,1,0);
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
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
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
    for (i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
    stats_total[0] += stats[0];
    if (window_limit) stats_total[0] -= stats_list[iwindow][0];
    stats_list[iwindow][0] = stats[0];
    stats_total[1] += stats[1];
    if (window_limit) stats_total[1] -= stats_list[iwindow][1];
    stats_list[iwindow][1] = stats[1];

    if (window_limit) m = nwindow;
    else m = iwindow+1;

    stats_list[iwindow][2] = stats[2];
    stats_total[2] = stats_list[0][2];
    for (i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (i = 0; i < nbins; i++) {
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

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d %g %g %g %g\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0)
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
    else
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);

    if (ferror(fp))
      error->one(FLERR,"Error writing out histogram data");

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) ftruncate(fileno(fp),fileend);
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
