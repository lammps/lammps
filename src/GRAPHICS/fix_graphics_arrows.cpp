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

#include "fix_graphics_arrows.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "group.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, DIPOLE, FORCE, VELOCITY, VARIABLE, CHUNK };

/* ---------------------------------------------------------------------- */

FixGraphicsArrows::FixGraphicsArrows(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), vec{0.0, 0.0, 0.0}, xstr(nullptr), ystr(nullptr), zstr(nullptr),
    cchunk(nullptr), cpos(nullptr), cvec(nullptr), id_chunk(nullptr), id_pos(nullptr),
    id_vec(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix graphics/arrows", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/arrows nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;
  scalar_flag = 1;
  extscalar = 0;

  mode = NONE;
  numobjs = 0;
  varflag = 0;
  scale = 1.0;
  autoscale = false;
  autovalue = 1.0;
  xvar = yvar = zvar = -1;

  int iarg = 7;
  if (strcmp(arg[4], "dipole") == 0) {
    mode = DIPOLE;
  } else if (strcmp(arg[4], "force") == 0) {
    mode = FORCE;
  } else if (strcmp(arg[4], "velocity") == 0) {
    mode = VELOCITY;
  } else if (strcmp(arg[4], "variable") == 0) {
    mode = VARIABLE;
    if (narg < 9) utils::missing_cmd_args(FLERR, "fix graphics/arrows variable", error);
    iarg = 9;
    if (strstr(arg[5], "v_") == arg[5]) {
      varflag = 1;
      xstr = utils::strdup(arg[5] + 2);
    } else {
      vec[0] = utils::numeric(FLERR, arg[5], false, lmp);
    }
    if (strstr(arg[6], "v_") == arg[6]) {
      varflag = 1;
      ystr = utils::strdup(arg[6] + 2);
    } else {
      vec[1] = utils::numeric(FLERR, arg[6], false, lmp);
    }
    if (strstr(arg[7], "v_") == arg[7]) {
      varflag = 1;
      zstr = utils::strdup(arg[7] + 2);
    } else {
      vec[2] = utils::numeric(FLERR, arg[7], false, lmp);
    }
    radius = utils::numeric(FLERR, arg[8], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 6, "Arrow radius must be > 0");
  } else if (strcmp(arg[4], "chunk") == 0) {
    mode = CHUNK;
    if (narg < 10) utils::missing_cmd_args(FLERR, "fix graphics/arrows chunk", error);
    iarg = 10;

    id_chunk = utils::strdup(arg[5]);
    cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
    if (!cchunk)
      error->all(FLERR, 5, "Chunk/atom compute {} does not exist or is incorrect style for fix {}",
                 id_chunk, style);

    id_pos = utils::strdup(arg[6]);
    cpos = modify->get_compute_by_id(id_pos);
    if (!cpos)
      error->all(FLERR, 6, "Per-chunk compute {} does not exist for fix graphics/arrows", id_pos);

    id_vec = utils::strdup(arg[7]);
    cvec = modify->get_compute_by_id(id_vec);
    if (!cvec)
      error->all(FLERR, 7, "Per-chunk compute {} does not exist for fix graphics/arrows", id_vec);

    scale = utils::numeric(FLERR, arg[8], false, lmp);
    radius = utils::numeric(FLERR, arg[9], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 9, "Arrow radius must be > 0");
  } else {
    error->all(FLERR, 4, "Unknown fix graphics/arrows keyword: {}", arg[4]);
  }

  // we have the same arguments for these modes
  if ((mode == DIPOLE) || (mode == FORCE) || (mode == VELOCITY)) {
    scale = utils::numeric(FLERR, arg[5], false, lmp);
    radius = utils::numeric(FLERR, arg[6], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 6, "Arrow radius must be > 0");
  }

  // parse optional keywords

  while (iarg < narg) {
    if (strcmp(arg[iarg], "autoscale") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/arrows autoscale", error);
      autoscale = true;
      autovalue = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (autovalue <= 0.0)
        error->all(FLERR, iarg + 1, "Fix graphics/arrow autoscale value must be > 0");
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/arrows keyword: {}", arg[iarg]);
    }
  }

  // checks

  if ((mode == DIPOLE) && !atom->mu_flag)
    error->all(FLERR, 4, "Fix graphics/arrows dipole mode requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

FixGraphicsArrows::~FixGraphicsArrows()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] id_chunk;
  delete[] id_pos;
  delete[] id_vec;
}

/* ---------------------------------------------------------------------- */

int FixGraphicsArrows::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsArrows::init()
{
  if (xstr) {
    int ivar = input->variable->find(xstr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows x value does not exist", xstr);
    if ((input->variable->atomstyle(ivar) == 0) && (input->variable->equalstyle(ivar) == 0))
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom- or equal-style variable", xstr);
    xvar = ivar;
  }
  if (ystr) {
    int ivar = input->variable->find(ystr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows y value does not exist", ystr);
    if ((input->variable->atomstyle(ivar) == 0) && (input->variable->equalstyle(ivar) == 0))
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom- or equal-style variable", ystr);
    yvar = ivar;
  }
  if (zstr) {
    int ivar = input->variable->find(zstr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows z value does not exist", zstr);
    if ((input->variable->atomstyle(ivar) == 0) && (input->variable->equalstyle(ivar) == 0))
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom- or equal-style variable", zstr);
    zvar = ivar;
  }

  // refresh and check per-chunk computes
  if (mode == CHUNK) {
    cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
    if (!cchunk)
      error->all(FLERR, Error::NOLASTLINE,
                 "Chunk/atom compute {} does not exist or is incorrect style for fix {}", id_chunk,
                 style);
    cpos = modify->get_compute_by_id(id_pos);
    if (!cpos)
      error->all(FLERR, Error::NOLASTLINE,
                 "Per-chunk compute {} does not exist for fix graphics/arrows", id_pos);
    if (!cpos->array_flag || (cpos->size_array_cols < 3))
      error->all(FLERR, Error::NOLASTLINE,
                 "Per-chunk compute {} is not compatible with fix graphics/arrows", id_pos);
    cvec = modify->get_compute_by_id(id_vec);
    if (!cvec)
      error->all(FLERR, Error::NOLASTLINE,
                 "Per-chunk compute {} does not exist for fix graphics/arrows", id_vec);
    if (!cvec->array_flag || (cvec->size_array_cols < 3))
      error->all(FLERR, Error::NOLASTLINE,
                 "Per-chunk compute {} is not compatible with fix graphics/arrows", id_vec);
  }

  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsArrows::end_of_step()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  // evaluate variable or per-chunk computes if necessary, wrap with clear/add

  if (varflag || (mode == CHUNK)) modify->clearstep_compute();

  const auto *const *const x = atom->x;
  const auto *const *const f = atom->f;
  const auto *const *const v = atom->v;
  const auto *const *const mu = atom->mu;
  const auto *const mask = atom->mask;
  const auto *const type = atom->type;
  const auto nlocal = atom->nlocal;

  // count (local) number of arrows

  if (mode != CHUNK) {
    int n = 0;
    for (int i = 0; i < nlocal; ++i)
      if (mask[i] & groupbit) ++n;

    numobjs = n;
    memory->create(imgobjs, numobjs, "fix_graphics_arrows:imgobjs");
    memory->create(imgparms, numobjs, 10, "fix_graphics_arrows:imgparms");
    double *xdata = nullptr;
    double *ydata = nullptr;
    double *zdata = nullptr;

    // update variable data, if needed
    if (mode == VARIABLE) {
      if (xstr) {
        if (input->variable->atomstyle(xvar)) {
          memory->create(xdata, nlocal, "fix_graphics_arrows:xdata");
          input->variable->compute_atom(xvar, igroup, xdata, 1, 0);
        } else {
          vec[0] = input->variable->compute_equal(xvar);
        }
      }
      if (ystr) {
        if (input->variable->atomstyle(yvar)) {
          memory->create(ydata, nlocal, "fix_graphics_arrows:ydata");
          input->variable->compute_atom(yvar, igroup, ydata, 1, 0);
        } else {
          vec[1] = input->variable->compute_equal(yvar);
        }
      }
      if (zstr) {
        if (input->variable->atomstyle(zvar)) {
          memory->create(zdata, nlocal, "fix_graphics_arrows:zdata");
          input->variable->compute_atom(zvar, igroup, zdata, 1, 0);
        } else {
          vec[2] = input->variable->compute_equal(zvar);
        }
      }
    }

    // determine automatic scale value for enabled automatic scaling of vectors, if needed

    auto num = group->count(igroup);
    if (autoscale && (num > 0)) {
      double mysum = 0.0;
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          if (mode == DIPOLE) {
            mysum += MathExtra::len3(mu[i]);
          } else if (mode == FORCE) {
            mysum += MathExtra::len3(f[i]);
          } else if (mode == VELOCITY) {
            mysum += MathExtra::len3(v[i]);
          } else if (mode == VARIABLE) {
            double myvec[3];
            if (xdata)
              myvec[0] = xdata[i];
            else
              myvec[0] = vec[0];
            if (ydata)
              myvec[1] = ydata[i];
            else
              myvec[1] = vec[1];
            if (zdata)
              myvec[2] = zdata[i];
            else
              myvec[2] = vec[2];
            mysum += MathExtra::len3(myvec);
          }
        }
      }

      // compute average across all processors and invert (num was checked for != 0 above)
      MPI_Allreduce(&mysum, &scale, 1, MPI_DOUBLE, MPI_SUM, world);
      if (scale != 0.0)
        scale = autovalue / (scale / static_cast<double>(num));
      else
        scale = 1.0;
    }

    n = 0;
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        imgobjs[n] = Graphics::ARROW;
        imgparms[n][0] = type[i];
        imgparms[n][1] = x[i][0];
        imgparms[n][2] = x[i][1];
        imgparms[n][3] = x[i][2];
        if (mode == DIPOLE) {
          imgparms[n][4] = mu[i][0];
          imgparms[n][5] = mu[i][1];
          imgparms[n][6] = mu[i][2];
        } else if (mode == FORCE) {
          imgparms[n][4] = f[i][0];
          imgparms[n][5] = f[i][1];
          imgparms[n][6] = f[i][2];
        } else if (mode == VELOCITY) {
          imgparms[n][4] = v[i][0];
          imgparms[n][5] = v[i][1];
          imgparms[n][6] = v[i][2];
        } else if (mode == VARIABLE) {
          if (xdata)
            imgparms[n][4] = xdata[i];
          else
            imgparms[n][4] = vec[0];
          if (ydata)
            imgparms[n][5] = ydata[i];
          else
            imgparms[n][5] = vec[1];
          if (zdata)
            imgparms[n][6] = zdata[i];
          else
            imgparms[n][6] = vec[2];
        }
        imgparms[n][7] = scale;
        imgparms[n][8] = 2.0 * radius;
        imgparms[n][9] = 0.2;
        ++n;
      }
    }
    memory->destroy(xdata);
    memory->destroy(ydata);
    memory->destroy(zdata);

  } else {    // mode == CHUNK
    numobjs = cchunk->setup_chunks();
    memory->create(imgobjs, numobjs, "fix_graphics_arrows:imgobjs");
    memory->create(imgparms, numobjs, 10, "fix_graphics_arrows:imgparms");

    // invoke per-chunk computes
    if (!(cpos->invoked_flag & Compute::INVOKED_ARRAY)) {
      cpos->compute_array();
      cpos->invoked_flag |= Compute::INVOKED_ARRAY;
    }
    if (!(cvec->invoked_flag & Compute::INVOKED_ARRAY)) {
      cvec->compute_array();
      cvec->invoked_flag |= Compute::INVOKED_ARRAY;
    }

    const double *const *const pdata = cpos->array;
    const double *const *const vdata = cvec->array;

    // determine automatic scale value for enabled automatic scaling of vectors, if needed

    double mysum = 0.0;
    if (autoscale) {
      for (int n = 0; n < numobjs; ++n) mysum += MathExtra::len3(vdata[n]);

      if (mysum != 0.0)
        scale = autovalue / (mysum / static_cast<double>(numobjs));
      else
        scale = 1.0;
    }
    for (int n = 0; n < numobjs; ++n) {
      imgobjs[n] = Graphics::ARROW;
      imgparms[n][0] = 1;
      imgparms[n][1] = pdata[n][0];
      imgparms[n][2] = pdata[n][1];
      imgparms[n][3] = pdata[n][2];
      imgparms[n][4] = vdata[n][0];
      imgparms[n][5] = vdata[n][1];
      imgparms[n][6] = vdata[n][2];
      imgparms[n][7] = scale;
      imgparms[n][8] = 2.0 * radius;
      imgparms[n][9] = 0.2;
    }
  }

  if (varflag || (mode == CHUNK))
    modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   current scale value determined for autoscalar or set statically
------------------------------------------------------------------------- */

double FixGraphicsArrows::compute_scalar()
{
  return scale;
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsArrows::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
