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

#include "fix_ave_grid.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid2d.h"
#include "grid3d.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ALL,SAMPLE,NONORM};
enum{ONE,RUNNING,WINDOW};

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.5;

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  which(nullptr), argindex(nullptr), ids(nullptr),
  value2index(nullptr), value2grid(nullptr), value2data(nullptr),
  grid2d(nullptr), grid3d(nullptr),
  grid_buf1(nullptr), grid_buf2(nullptr),
  vec2d(nullptr), vec3d(nullptr), array2d(nullptr), array3d(nullptr),
  count2d(nullptr), count3d(nullptr),
  vec2d_sample(nullptr), vec3d_sample(nullptr), array2d_sample(nullptr), array3d_sample(nullptr),
  count2d_sample(nullptr), count3d_sample(nullptr),
  vec2d_epoch(nullptr), vec3d_epoch(nullptr), array2d_epoch(nullptr), array3d_epoch(nullptr),
  count2d_epoch(nullptr), count3d_epoch(nullptr),
  vec2d_running(nullptr), vec3d_running(nullptr), array2d_running(nullptr), array3d_running(nullptr),
  count2d_running(nullptr), count3d_running(nullptr)
{
  if (narg < 10) error->all(FLERR,"Illegal fix ave/grid command");

  pergrid_flag = 1;
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  pergrid_freq = utils::inumeric(FLERR,arg[5],false,lmp);
  time_depend = 1;

  // NOTE: allow Dxyz as well at some point

  nxgrid = utils::inumeric(FLERR,arg[6],false,lmp);
  nygrid = utils::inumeric(FLERR,arg[7],false,lmp);
  nzgrid = utils::inumeric(FLERR,arg[8],false,lmp);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-9,&arg[9],1,earg,lmp);

  if (earg != &arg[9]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  value2grid = new int[nargnew];
  value2data = new int[nargnew];

  modeatom = modegrid = 0;
  nvalues = 0;

  int iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = nullptr;

    if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues] = 1;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues] = 2;
      modeatom = 1;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues] = 1;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues] = 2;
      modeatom = 1;

    } else if (strcmp(arg[iarg],"density/number") == 0) {
      which[nvalues] = ArgInfo::DENSITY_NUMBER;
      argindex[nvalues] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"density/mass") == 0) {
      which[nvalues] = ArgInfo::DENSITY_MASS;
      argindex[nvalues] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      which[nvalues] = ArgInfo::MASS;
      argindex[nvalues] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      which[nvalues] = ArgInfo::TEMPERATURE;
      argindex[nvalues] = 0;
      modeatom = 1;

    } else {
      ArgInfo argi(arg[iarg]);

      if (argi.get_type() == ArgInfo::NONE) break;
      if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR,"Invalid fix ave/grid command");

      which[nvalues] = argi.get_type();
      argindex[nvalues] = argi.get_index1();
      ids[nvalues] = argi.copy_name();

      if (strchr(ids[nvalues],':')) modegrid = 1;
      else modeatom = 1;

      if (modegrid && which[nvalues] == ArgInfo::VARIABLE)
        error->all(FLERR,"Fix ave/grid cannot use variable for grid info");
    }

    nvalues++;
    iarg++;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/grid command");

  if (modeatom && modegrid)
    error->all(FLERR,"Fix ave/grid cannot operate on per-atom and "
               "per-grid values");

  // optional args

  normflag = ALL;
  aveflag = ONE;
  nwindow = 0;
  biasflag = 0;
  id_bias = nullptr;
  adof = domain->dimension;
  cdof = 0.0;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"all") == 0) normflag = ALL;
      else if (strcmp(arg[iarg+1],"sample") == 0) normflag = SAMPLE;
      else if (strcmp(arg[iarg+1],"none") == 0) normflag = NONORM;
      else error->all(FLERR,"Illegal fix ave/grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"one") == 0) aveflag = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) aveflag = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) aveflag = WINDOW;
      else error->all(FLERR,"Illegal fix ave/grid command");
      if (aveflag == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/grid command");
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/grid command");
        iarg++;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/grid command");
      biasflag = 1;
      id_bias = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"adof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/grid command");
      adof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cdof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ave/grid command");
      cdof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else error->all(FLERR,"Illegal fix ave/grid command");
  }

  // if wildcard expansion occurred, free earg memory from exapnd_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  dimension = domain->dimension;

  if (nevery <= 0 || nrepeat <= 0 || pergrid_freq <= 0)
    error->all(FLERR,"Illegal fix ave/grid command");
  if (pergrid_freq % nevery || nrepeat*nevery > pergrid_freq)
    error->all(FLERR,"Illegal fix ave/grid command");

  if (nxgrid < 1 || nygrid < 1 || nzgrid < 1)
    error->all(FLERR,"Invalid fix ave/grid grid size");
  if (dimension == 2 && nzgrid != 1)
    error->all(FLERR,"Fix ave/grid grid Nz must be 1 for 2d simulation");

  if (biasflag) {
    tbias = modify->get_compute_by_id(id_bias);
    if (!tbias)
      error->all(FLERR,"Could not find compute ID for temperature bias");
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute does not calculate a velocity bias");
  }

  // error checks for ATOM mode

  if (modeatom) {
    for (int i = 0; i < nvalues; i++) {
      if (which[i] == ArgInfo::COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        if (icompute < 0)
          error->all(FLERR,"Compute ID for fix ave/grid does not exist");
        if (modify->compute[icompute]->peratom_flag == 0)
          error->all(FLERR,
                     "Fix ave/atom compute does not calculate per-atom values");
        if (argindex[i] == 0 &&
            modify->compute[icompute]->size_peratom_cols != 0)
          error->all(FLERR,"Fix ave/atom compute does not "
                     "calculate a per-atom vector");
        if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
          error->all(FLERR,"Fix ave/atom compute does not "
                     "calculate a per-atom array");
        if (argindex[i] &&
            argindex[i] > modify->compute[icompute]->size_peratom_cols)
          error->all(FLERR,"Fix ave/atom compute array is accessed out-of-range");

      } else if (which[i] == ArgInfo::FIX) {
        int ifix = modify->find_fix(ids[i]);
        if (ifix < 0)
          error->all(FLERR,"Fix ID for fix ave/atom does not exist");
        if (modify->fix[ifix]->peratom_flag == 0)
          error->all(FLERR,"Fix ave/atom fix does not calculate per-atom values");
        if (argindex[i] == 0 && modify->fix[ifix]->size_peratom_cols != 0)
          error->all(FLERR,
                     "Fix ave/atom fix does not calculate a per-atom vector");
        if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
          error->all(FLERR,
                     "Fix ave/atom fix does not calculate a per-atom array");
        if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
          error->all(FLERR,"Fix ave/atom fix array is accessed out-of-range");
        if (nevery % modify->fix[ifix]->peratom_freq)
          error->all(FLERR,
                     "Fix for fix ave/atom not computed at compatible time");

      } else if (which[i] == ArgInfo::VARIABLE) {
        int ivariable = input->variable->find(ids[i]);
        if (ivariable < 0)
          error->all(FLERR,"Variable name for fix ave/atom does not exist");
        if (input->variable->atomstyle(ivariable) == 0)
          error->all(FLERR,"Fix ave/atom variable is not atom-style variable");
      }
    }
  }

  // setup and error checks for GRID mode

  if (modegrid) {
    for (int i = 0; i < nvalues; i++) {
      if (which[i] == ArgInfo::COMPUTE) {

        auto words = utils::parse_gridid(FLERR,ids[i],error);
        const auto &idcompute = words[0];
        const auto &gname = words[1];
        const auto &dname = words[2];

        delete[] ids[i];
        ids[i] = utils::strdup(idcompute);

        auto icompute = modify->get_compute_by_id(idcompute);
        if (!icompute) error->all(FLERR,"Could not find fix ave/grid compute ID: {}",idcompute);
        if (icompute->pergrid_flag == 0)
          error->all(FLERR,"Fix ave/grid compute {} does not compute per-grid info",idcompute);

        int dim;
        int igrid = icompute->get_grid_by_name(gname,dim);
        if (igrid < 0)
          error->all(FLERR,"Fix ave/grid compute {} does not recognize grid name {}",
                     idcompute,gname);

        int ncol;
        int idata = icompute->get_griddata_by_name(igrid,dname,ncol);
        if (idata < 0)
          error->all(FLERR,"Fix ave/grid compute {} does not recognize data name {}",
                     idcompute,dname);

        if (argindex[i] == 0 && ncol)
          error->all(FLERR,"Fix ave/grid compute {} data {} is not per-grid vector",idcompute,dname);
        if (argindex[i] && ncol == 0)
          error->all(FLERR,"Fix ave/grid compute {} data {} is not per-grid array",idcompute,dname);
        if (argindex[i] && argindex[i] > ncol)
          error->all(FLERR,"Fix ave/grid compute {} array {} is accessed out-of-range",
                     idcompute,dname);

        value2grid[i] = igrid;
        value2data[i] = idata;

      } else if (which[i] == ArgInfo::FIX) {

        auto words = utils::parse_gridid(FLERR,ids[i],error);
        const auto &idfix = words[0];
        const auto &gname = words[1];
        const auto &dname = words[2];

        delete[] ids[i];
        ids[i] = utils::strdup(idfix);

        Fix *ifix = modify->get_fix_by_id(idfix);
        if (!ifix) error->all(FLERR,"Could not find fix ave/grid fix ID: {}",idfix);
        if (ifix->pergrid_flag == 0)
          error->all(FLERR,"Fix ave/grid fix {} does not compute per-grid info",idfix);
        if (nevery % ifix->pergrid_freq)
          error->all(FLERR, "Fix ID {} for fix grid/atom not computed at compatible time",idfix);

        int dim;
        int igrid = ifix->get_grid_by_name(gname,dim);
        if (igrid < 0)
          error->all(FLERR,"Fix ave/grid fix {} does not recognize grid name {}",idfix,gname);

        int ncol;
        int idata = ifix->get_griddata_by_name(igrid,dname,ncol);
        if (idata < 0)
          error->all(FLERR,"Fix ave/grid fix {} does not recognize data name {}",idfix,dname);

        if (argindex[i] == 0 && ncol)
          error->all(FLERR, "Fix ave/grid fix {} data {} is not per-grid vector",idfix,dname);
        if (argindex[i] && ncol == 0)
          error->all(FLERR,"Fix ave/grid fix {} data {} is not per-grid array",idfix,dname);
        if (argindex[i] && argindex[i] > ncol)
          error->all(FLERR,"Fix ave/grid fix {} array {} is accessed out-of-range",idfix,dname);

        value2grid[i] = igrid;
        value2data[i] = idata;
      }
    }
  }

  // instantiate the Grid class and allocate per-grid memory

  double maxdist;

  if (modeatom) {
    maxdist = 0.5 * neighbor->skin;
    shift = OFFSET + SHIFT;
  } else if (modegrid) {
    maxdist = 0.0;
    shift = 0.0;
  }

  if (dimension == 2) {
    grid2d = new Grid2d(lmp, world, nxgrid, nygrid, maxdist, 0, SHIFT,
                        nxlo_in, nxhi_in, nylo_in, nyhi_in,
                        nxlo_out, nxhi_out, nylo_out, nyhi_out);

    // ngrid_buf12 converted to nvalues + count

    grid2d->setup(ngrid_buf1, ngrid_buf2);
    ngrid_buf1 *= nvalues + 1;
    ngrid_buf2 *= nvalues + 1;

    memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1);

  } else {
    grid3d = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid, maxdist, 0, SHIFT,
                        nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                        nxlo_out, nxhi_out, nylo_out, nyhi_out,
                        nzlo_out, nzhi_out);

    // ngrid_buf12 converted to nvalues + count

    grid3d->setup(ngrid_buf1, ngrid_buf2);
    ngrid_buf1 *= nvalues + 1;
    ngrid_buf2 *= nvalues + 1;

    memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
      (nzhi_out - nzlo_out + 1);
  }

  // allocate all per-grid data

  allocate_grid();

  // zero the output grid and counts since dump may access it on timestep 0
  // zero running grid if used by ave

  running_count = 0;
  window_count = 0;
  window_oldest = -1;
  window_newest = 0;

  zero_grid(count2d,vec2d,array2d,count3d,vec3d,array3d);
  if (aveflag == RUNNING || aveflag == WINDOW)
    zero_grid(count2d_running,vec2d_running,array2d_running,
              count3d_running,vec3d_running,array3d_running);

  // bin indices and skip flags for ATOM mode
  // vresult for per-atom variable evaluation

  maxatom = 0;
  bin = nullptr;
  skip = nullptr;

  maxvar = 0;
  vresult = nullptr;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  irepeat = 0;
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveGrid::~FixAveGrid()
{
  delete[] which;
  delete[] argindex;
  for (int m = 0; m < nvalues; m++) delete[] ids[m];
  delete[] ids;
  delete[] value2index;
  delete[] value2grid;
  delete[] value2data;

  delete grid2d;
  delete grid3d;

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  // deallocate all per-grid data

  deallocate_grid();

  if (modeatom) {
    memory->destroy(bin);
    memory->destroy(skip);
    memory->destroy(vresult);
  }
}

/* ---------------------------------------------------------------------- */

int FixAveGrid::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::init()
{
  if (biasflag) {
    tbias = modify->get_compute_by_id(id_bias);
    if (!tbias)
      error->all(FLERR,"Could not find compute ID for temperature bias");
  }

  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/grid does not exist");
      value2index[m] = icompute;

    } else if (which[m] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/grid does not exist");
      value2index[m] = ifix;

    } else if (which[m] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/grid does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }

  // check that grid sizes for all fields match grid size for this fix

  if (modegrid) {
    Compute *compute;
    Fix *fix;
    Grid2d *grid2d;
    Grid3d *grid3d;

    int nxtmp,nytmp,nztmp;

    for (int m = 0; m < nvalues; m++) {
      if (dimension == 2) {
        if (which[m] == ArgInfo::COMPUTE) {
          compute = modify->compute[value2index[m]];
          grid2d = (Grid2d *) compute->get_grid_by_index(value2grid[m]);
        } else {
          fix = modify->fix[value2index[m]];
          grid2d = (Grid2d *) fix->get_grid_by_index(value2grid[m]);
        }
        grid2d->get_size(nxtmp,nytmp);
        if (nxtmp != nxgrid || nytmp != nygrid)
          error->all(FLERR,"Fix ave/grid value grid sizes do not match");

      } else {
        if (which[m] == ArgInfo::COMPUTE) {
          compute = modify->compute[value2index[m]];
          grid3d = (Grid3d *) compute->get_grid_by_index(value2grid[m]);
        } else {
          fix = modify->fix[value2index[m]];
          grid3d = (Grid3d *) fix->get_grid_by_index(value2grid[m]);
        }
        grid3d->get_size(nxtmp,nytmp,nztmp);
        if (nxtmp != nxgrid || nytmp != nygrid || nztmp != nzgrid)
          error->all(FLERR,"Fix ave/grid value grid sizes do not match");
      }
    }
  }

  // set triclinic flag

  triclinic = domain->triclinic;

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

void FixAveGrid::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::end_of_step()
{
  int m,ix,iy,iz;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // if first sample in epoch, zero owned and ghost grid points
  
  if (irepeat == 0) zero_grid(count2d_sample,vec2d_sample,array2d_sample,
                              count3d_sample,vec3d_sample,array3d_sample);
  if (irepeat == 0 && normflag == SAMPLE)
    zero_grid(count2d_epoch,vec2d_epoch,array2d_epoch,
              count3d_epoch,vec3d_epoch,array3d_epoch);

  // accumulate per-grid values for one sample for either ATOM or GRID mode
  // per-atom compute/fix/variable may invoke computes so wrap with clear/add

  if (modeatom) {
    modify->clearstep_compute();
    atom2grid();
  } else {
    grid2grid();
  }

  // return if irepeat < nrepeat and norm = SAMPLE or NONORM

  irepeat++;
  if (irepeat < nrepeat && (normflag == ALL || normflag == NONORM)) {
    nvalid += nevery;
    if (modeatom) modify->addstep_compute(nvalid);
    return;
  }

  // for ATOM mode, perform ghost to owned grid communication
  // done every sample for norm = SAMPLE, else once per Nfreq epoch
  // nvalues+1 includes atom count

  if (modeatom) {
    if (dimension == 2)
      grid2d->reverse_comm(Grid2d::FIX,this,nvalues+1,sizeof(double),0,
                           grid_buf1,grid_buf2,MPI_DOUBLE);
    else
      grid3d->reverse_comm(Grid3d::FIX,this,nvalues+1,sizeof(double),0,
                           grid_buf1,grid_buf2,MPI_DOUBLE);
  }

  // for norm = SAMPLE, sum per-sample values and count to per-epoch
  // for ATOM mode, also need to normalize for single sample
  // then return if irepeat < nrepeat

  if (normflag == SAMPLE) {
    if (modeatom) normalize_atom(1);
    sum_sample_to_epoch();

    if (irepeat < nrepeat && (normflag == ALL || normflag == NONORM)) {
      nvalid += nevery;
      if (modeatom) modify->addstep_compute(nvalid);
      return;
    }
  }

  // reset irepeat and nvalid

  irepeat = 0;
  nvalid = ntimestep+pergrid_freq - ((bigint) nrepeat-1)*nevery;
  if (modeatom) modify->addstep_compute(nvalid);

  // just return if this proc owns no grid points

  if (ngridout == 0) return;

  // average the final results for the entire epoch
  // for norm = ALL, normalize values_sample by Nrepeat and counts
  // for norm = SAMPLE, normalize values_sample_sum by Nrepeat

  if (modeatom) {
    if (normflag == ALL || normflag == NONORM) normalize_atom(nrepeat);
    else if (normflag == SAMPLE) {
      normalize_grid(nrepeat,vec2d_epoch,array2d_epoch,vec3d_epoch,array3d_epoch);
      copy_epoch_to_sample();
    }
  }

  if (modegrid) 
    normalize_grid(nrepeat,vec2d_sample,array2d_sample,
                   vec3d_sample,array3d_sample);

  // create Nfreq output

  if (aveflag == ONE) {
    copy_sample_to_output();

  } else if (aveflag == RUNNING) {
    running_count++;
    sum_sample_to_running();
    copy_running_to_output();
    normalize_grid(running_count,vec2d_running,array2d_running,
                   vec3d_running,array3d_running);
      
  } else if (aveflag == WINDOW) {
    sum_sample_to_running();

    if (window_oldest >= 0) {
      subtract_window_from_running();
      window_count--;
    }
    window_oldest++;
    if (window_oldest == nwindow) window_oldest = 0;

    copy_sample_to_window(window_newest);
    window_count++;
    window_newest++;
    if (window_newest == nwindow) window_newest = 0;

    copy_running_to_output();
    normalize_grid(window_count,vec2d_running,array2d_running,
                   vec3d_running,array3d_running);
  }
}

/* ----------------------------------------------------------------------
   sum per-atom contributions to owned+ghost grid cells
   sets one of vec2d,array2d,vec3d,array3d sample
   also set count2d or count3d for atom count per bin
------------------------------------------------------------------------- */

void FixAveGrid::atom2grid()
{
  int i,j,m,n,ix,iy,iz;

  // bin[i][dim] = indices of bin each atom is in
  // not set if group mask does not match
  // also count atoms contributing to each bin
  // check if any atom is out of bounds for my local grid

  double *boxlo,*prd;
  int *periodicity = domain->periodicity;

  if (triclinic) {
    boxlo = domain->boxlo_lamda;
    prd = domain->prd_lamda;
  } else {
    boxlo = domain->boxlo;
    prd = domain->prd;
  }

  double dxinv = nxgrid/prd[0];
  double dyinv = nygrid/prd[1];
  double dzinv = nzgrid/prd[2];

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxatom) {
    memory->destroy(bin);
    memory->destroy(skip);
    maxatom = atom->nmax;
    memory->create(bin,maxatom,dimension,"ave/grid:bin");
    memory->create(skip,maxatom,"ave/grid:skip");
  }

  if (triclinic) domain->x2lamda(nlocal);

  int flag = 0;

  if (dimension == 2) {
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) {
        skip[i] = 1;
        continue;
      }

      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;

      if (ix < nxlo_out || ix > nxhi_out) {
        if (periodicity[0]) flag = 1;
        else skip[i] = 1;
        continue;
      }
      if (iy < nylo_out || iy > nyhi_out) {
        if (periodicity[1]) flag = 1;
        else skip[i] = 1;
        continue;
      }

      skip[i] = 0;
      count2d_sample[iy][ix] += 1.0;
      bin[i][0] = iy;
      bin[i][1] = ix;
    }

  } else {
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) {
        skip[i] = 1;
        continue;
      }

      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;

      if (ix < nxlo_out || ix > nxhi_out) {
        if (periodicity[0]) flag = 1;
        else skip[i] = 1;
        continue;
      }
      if (iy < nylo_out || iy > nyhi_out) {
        if (periodicity[1]) flag = 1;
        else skip[i] = 1;
        continue;
      }
      if (iz < nzlo_out || iz > nyhi_out) {
        if (periodicity[2]) flag = 1;
        else skip[i] = 1;
        continue;
      }

      skip[i] = 0;
      count3d_sample[iz][iy][ix] += 1.0;
      bin[i][0] = iz;
      bin[i][1] = iy;
      bin[i][2] = ix;
    }
  }

  if (flag) error->one(FLERR,"Out of range fix ave/grid atoms");

  if (triclinic) domain->lamda2x(nlocal);

  // loop over user-specified values

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    // V,F adds velocity,force to value

    if (which[m] == ArgInfo::V || which[m] == ArgInfo::F) {

      double **attribute;
      if (which[m] == ArgInfo::V) attribute = atom->v;
      else if (which[m] == ArgInfo::F) attribute = atom->f;

      if (dimension == 2) {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              vec2d_sample[bin[i][0]][bin[i][1]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              array2d_sample[bin[i][0]][bin[i][1]][m] += attribute[i][j];
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              vec3d_sample[bin[i][0]][bin[i][1]][bin[i][2]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              array3d_sample[bin[i][0]][bin[i][1]][bin[i][2]][m] += attribute[i][j];
          }
      }

    // DENSITY_NUMBER adds 1 to value
    // DENSITY_MASS or MASS adds mass to value

    } else if ((which[m] == ArgInfo::DENSITY_NUMBER) ||
               (which[m] == ArgInfo::DENSITY_MASS) ||
               (which[m] == ArgInfo::MASS)) {

      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;
      double one;

      if (dimension == 2) {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec2d_sample[bin[i][0]][bin[i][1]] += one;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array2d_sample[bin[i][0]][bin[i][1]][m] += one;
            }
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec3d_sample[bin[i][0]][bin[i][1]][bin[i][2]] += one;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array3d_sample[bin[i][0]][bin[i][1]][bin[i][2]][m] += one;
            }
          }
      }

    // TEMPERATURE adds KE to values
    // subtract and restore velocity bias if requested

    } else if (which[m] == ArgInfo::TEMPERATURE) {

      if (biasflag) {
        if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
        tbias->remove_bias_all();
      }

      double **v = atom->v;
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;
      double vsq,one;

      if (dimension == 2) {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec2d_sample[bin[i][0]][bin[i][1]] += one*vsq;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array2d_sample[bin[i][0]][bin[i][1]][m] += one*vsq;;
            }
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec3d_sample[bin[i][0]][bin[i][1]][bin[i][2]] += one*vsq;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array3d_sample[bin[i][0]][bin[i][1]][bin[i][2]][m] += one*vsq;
            }
          }
      }

      if (biasflag) tbias->restore_bias_all();

    // per-atom compute or fix or variable
    // invoke compute if not previously invoked
    // evaluate atom-style variable

    } else if (which[m] == ArgInfo::COMPUTE || which[m] == ArgInfo::FIX ||
               which[m] == ArgInfo::VARIABLE) {
      double *ovector,**oarray;

      if (which[m] == ArgInfo::COMPUTE) {
        Compute *compute = modify->compute[n];
        if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
          compute->compute_peratom();
          compute->invoked_flag |= Compute::INVOKED_PERATOM;
        }
        if (j == 0) ovector = compute->vector_atom;
        else oarray = compute->array_atom;

      } else if (which[m] == ArgInfo::FIX) {
        Fix *fix = modify->fix[n];
        if (j == 0) ovector = fix->vector_atom;
        else oarray = fix->array_atom;
      } else if (which[m] == ArgInfo::VARIABLE) {
        if (nlocal > maxvar) {
          memory->destroy(vresult);
          maxvar = atom->nmax;
          memory->create(vresult,maxvar,"ave/grid:vresult");
        }
        input->variable->compute_atom(n,igroup,vresult,1,0);
        ovector = vresult;
      }

      if (dimension == 2) {
        if (nvalues == 1) {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec2d_sample[bin[i][0]][bin[i][1]] += ovector[i];
            }
          } else {
            int jm1 = j = 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec2d_sample[bin[i][0]][bin[i][1]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array2d_sample[bin[i][0]][bin[i][1]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array2d_sample[bin[i][0]][bin[i][1]][m] += oarray[i][jm1];
            }
          }
        }

      } else {
        if (nvalues == 1) {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec3d_sample[bin[i][0]][bin[i][1]][bin[i][2]] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec3d_sample[bin[i][0]][bin[i][1]][bin[i][2]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array3d_sample[bin[i][0]][bin[i][1]][bin[i][2]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array3d_sample[bin[i][0]][bin[i][1]][bin[i][2]][m] += oarray[i][jm1];
            }
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   copy per-grid values from other computes/fixes to owned grid cells
   sets one of vec2d,array2d,vec3d,array3d
------------------------------------------------------------------------- */

void FixAveGrid::grid2grid()
{
  int j,m,n,ix,iy,iz;

  // loop over user-specified values

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];
    int idata = value2data[m];

    Compute *compute;
    Fix *fix;

    if (which[m] == ArgInfo::COMPUTE) {
      compute = modify->compute[n];
      if (!(compute->invoked_flag & Compute::INVOKED_PERGRID)) {
        compute->compute_pergrid();
        compute->invoked_flag |= Compute::INVOKED_PERGRID;
      }
    } else if (which[m] == ArgInfo::FIX) fix = modify->fix[n];

    if (dimension == 2) {
      double **ovec2d,***oarray2d;
      if (which[m] == ArgInfo::COMPUTE) {
        if (j == 0)
          ovec2d = (double **) compute->get_griddata_by_index(idata);
        else
          oarray2d = (double ***) compute->get_griddata_by_index(idata);
      } else {
        if (j == 0)
          ovec2d = (double **) fix->get_griddata_by_index(idata);
        else
          oarray2d = (double ***) fix->get_griddata_by_index(idata);
      }

      if (nvalues == 1) {
        if (j == 0) {
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              vec2d_sample[iy][ix] += ovec2d[iy][ix];
        } else {
          int jm1 = j - 1;
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              vec2d_sample[iy][ix] += oarray2d[iy][ix][jm1];
        }
      } else {
        if (j == 0) {
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              array2d_sample[iy][ix][m] += ovec2d[iy][ix];
        } else {
          int jm1 = j - 1;
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              array2d_sample[iy][ix][m] += oarray2d[iy][ix][jm1];
        }
      }

    } else {
      double ***ovec3d,****oarray3d;
      if (which[m] == ArgInfo::COMPUTE) {
        if (j == 0)
          ovec3d = (double ***) compute->get_griddata_by_index(idata);
        else
          oarray3d = (double ****) compute->get_griddata_by_index(idata);
      } else {
        if (j == 0) {
          ovec3d = (double ***) fix->get_griddata_by_index(idata);
        } else
          oarray3d = (double ****) fix->get_griddata_by_index(idata);
      }

      if (nvalues == 1) {
        if (j == 0) {
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                vec3d_sample[iz][iy][ix] += ovec3d[iz][iy][ix];
        } else {
          int jm1 = j - 1;
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                vec3d_sample[iz][iy][ix] += oarray3d[iz][iy][ix][jm1];
        }
      } else {
        if (j == 0) {
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                array3d_sample[iz][iy][ix][m] += ovec3d[iz][iy][ix];
        } else {
          int jm1 = j - 1;
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                array3d_sample[iz][iy][ix][m] += oarray3d[iz][iy][ix][jm1];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all grids which will be used
   if ATOM mode, also allocate per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::allocate_grid()
{
  if (dimension == 2) {
    
    if (nvalues == 1) {
      memory->create2d_offset(vec2d, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:vec2d");
      memory->create2d_offset(vec2d_sample, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:vec2d_sample");
      if (normflag == SAMPLE)
        memory->create2d_offset(vec2d_epoch, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:vec2d_epoch");
      if (aveflag == RUNNING) 
        memory->create2d_offset(vec2d_running, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:vec2d_running");
    } else {
      memory->create3d_offset_last(array2d, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                   nvalues, "ave/grid:array3d");
      memory->create3d_offset_last(array2d_sample, nylo_out, nyhi_out, nxlo_out, nxhi_out, 
                                   nvalues, "ave/grid:array3d_sample");
      if (normflag == SAMPLE)
        memory->create3d_offset_last(array2d_epoch, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                     nvalues, "ave/grid:array3d_epoch");
      if (aveflag == RUNNING)
        memory->create3d_offset_last(array2d_running, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                     nvalues, "ave/grid:array3d_running");
    }

    if (modeatom) {
      memory->create2d_offset(count2d, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:count2d");
      memory->create2d_offset(count2d_sample, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:count2d_sample");
      if (normflag == SAMPLE)
        memory->create2d_offset(count2d_epoch, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:count2d_epoch");
      if (normflag == RUNNING) {
        memory->create2d_offset(count2d_running, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:count2d_running");
      }
    }

  } else if (dimension == 3) {

    if (nvalues == 1) {
      memory->create3d_offset(vec3d, nzlo_out, nzhi_out, nylo_out,
                              nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:vec3d");
      memory->create3d_offset(vec3d_sample, nzlo_out, nzhi_out, nylo_out, 
                              nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:vec3d_sample");
      if (normflag == SAMPLE)
        memory->create3d_offset(vec3d_epoch, nzlo_out, nzhi_out, nylo_out, 
                                nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:vec3d_epoch");
      if (aveflag == RUNNING) 
        memory->create3d_offset(vec3d_running, nzlo_out, nzhi_out, nylo_out, 
                                nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:vec3d_running");
    } else {
      memory->create4d_offset_last(array3d, nzlo_out, nzhi_out, nylo_out,
                                   nyhi_out, nxlo_out, nxhi_out, nvalues,
                                   "ave/grid:array3d");
      memory->create4d_offset_last(array3d_sample, nzlo_out, nzhi_out, nylo_out, 
                                   nyhi_out, nxlo_out, nxhi_out, 
                                   nvalues, "ave/grid:array3d_sample");
      if (normflag == SAMPLE)
        memory->create4d_offset_last(array3d_epoch, nzlo_out, nzhi_out, nylo_out, 
                                     nyhi_out, nxlo_out, nxhi_out,
                                     nvalues, "ave/grid:array3d_epoch");
      if (aveflag == RUNNING)
        memory->create4d_offset_last(array3d_running, nzlo_out, nzhi_out, nylo_out, 
                                     nyhi_out, nxlo_out, nxhi_out,
                                     nvalues, "ave/grid:array3d_running");
    }
    
    if (modeatom) {
      memory->create3d_offset(count3d, nzlo_out, nzhi_out, nylo_out,
                              nyhi_out, nxlo_out, nxhi_out, "ave/grid:count3d");
      memory->create3d_offset(count3d_sample, nzlo_out, nzhi_out, nylo_out, 
                              nyhi_out, nxlo_out, nxhi_out,
                              "ave/grid:count3d_sample");
      if (normflag == SAMPLE)
        memory->create3d_offset(count3d_epoch, nzlo_out, nzhi_out, nylo_out, 
                                nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:count3d_epoch");
      if (normflag == RUNNING)
        memory->create3d_offset(count3d_running, nzlo_out, nzhi_out, nylo_out, 
                                nyhi_out, nxlo_out, nxhi_out,
                                "ave/grid:count3d_running");
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate all grids and counts
------------------------------------------------------------------------- */

void FixAveGrid::deallocate_grid()
{
  memory->destroy2d_offset(vec2d,nylo_out,nxlo_out);
  memory->destroy3d_offset_last(array2d,nylo_out,nxlo_out);
  memory->destroy2d_offset(count2d,nylo_out,nxlo_out);

  memory->destroy2d_offset(vec2d_sample,nylo_out,nxlo_out);
  memory->destroy3d_offset_last(array2d_sample,nylo_out,nxlo_out);
  memory->destroy2d_offset(count2d_sample,nylo_out,nxlo_out);

  memory->destroy2d_offset(vec2d_epoch,nylo_out,nxlo_out);
  memory->destroy3d_offset_last(array2d_epoch,nylo_out,nxlo_out);
  memory->destroy2d_offset(count2d_epoch,nylo_out,nxlo_out);

  memory->destroy2d_offset(vec2d_running,nylo_out,nxlo_out);
  memory->destroy3d_offset_last(array2d_running,nylo_out,nxlo_out);
  memory->destroy2d_offset(count2d_running,nylo_out,nxlo_out);

  memory->destroy3d_offset(vec3d,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(count3d,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(vec3d_sample,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d_sample,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(count3d_sample,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(vec3d_epoch,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d_epoch,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(count3d_epoch,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(vec3d_running,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d_running,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(count3d_running,nzlo_out,nylo_out,nxlo_out);
}

/* ----------------------------------------------------------------------
   zero values for a grid incluing ghost cells
   if ATOM mode, also zero per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::zero_grid(double **c2d, double **v2d, double ***a2d, 
                           double ***c3d, double ***v3d, double ****a3d)
{
  if (!ngridout) return;

  if (dimension == 2) {
    if (modeatom) 
      memset(&c2d[nylo_out][nxlo_out],0,ngridout*sizeof(double));
    if (nvalues == 1)
      memset(&v2d[nylo_out][nxlo_out],0,ngridout*sizeof(double));
    else
      memset(&a2d[nylo_out][nxlo_out][0],0,ngridout*nvalues*sizeof(double));
  } else {
    if (modeatom)
      memset(&c3d[nzlo_out][nylo_out][nxlo_out],0,ngridout*sizeof(double));
    if (nvalues == 1)
      memset(&v3d[nzlo_out][nylo_out][nxlo_out],0,ngridout*sizeof(double));
    else
      memset(&a3d[nzlo_out][nylo_out][nxlo_out][0],0,
             ngridout*nvalues*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   sum sample grid values to Nfreq epoch grid, just for owned grid cells
   if ATOM mode, also sum per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::sum_sample_to_epoch()
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d_epoch[iy][ix] += vec2d_sample[iy][ix];
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            array2d_epoch[iy][ix][m] += array2d_sample[iy][ix][m];
    }
    if (modeatom) 
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count2d_epoch[iy][ix] += count2d_sample[iy][ix];

  } else {
    if (nvalues == 1) {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d_epoch[iz][iy][ix] += vec3d_sample[iz][iy][ix];
    } else {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              array3d_epoch[iz][iy][ix][m] += array3d_sample[iz][iy][ix][m];
    }
    if (modeatom) 
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            count3d_epoch[iz][iy][ix] += count3d_sample[iz][iy][ix];
  }
}

/* ----------------------------------------------------------------------
   copy Nfreq epoch grid values to sample grid, just for owned grid cells
   if ATOM mode, also copy per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::copy_epoch_to_sample()
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d_sample[iy][ix] = vec2d_epoch[iy][ix];
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            array2d_sample[iy][ix][m] = array2d_epoch[iy][ix][m];
    }
    if (modeatom) 
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count2d_sample[iy][ix] = count2d_epoch[iy][ix];

  } else {
    if (nvalues == 1) {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d_sample[iz][iy][ix] = vec3d_epoch[iz][iy][ix];
    } else {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              array3d_sample[iz][iy][ix][m] = array3d_epoch[iz][iy][ix][m];
    }
    if (modeatom) 
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            count3d_sample[iz][iy][ix] = count3d_epoch[iz][iy][ix];
  }
}

/* ----------------------------------------------------------------------
   copy sample grid values to output grid, just for owned grid cells
   if ATOM mode, also copy per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::copy_sample_to_output()
{
  int ix,iy,iz,m;

  if (!ngridout) return;
        
  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = vec2d_sample[iy][ix];
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            array2d[iy][ix][m] = array2d_sample[iy][ix][m];
    }
    if (modeatom) 
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count2d[iy][ix] = count2d_sample[iy][ix];

  } else {
    if (nvalues == 1) {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = vec3d_sample[iz][iy][ix];
    } else {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              array3d[iz][iy][ix][m] = array3d_sample[iz][iy][ix][m];
    }
    if (modeatom) 
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            count3d[iz][iy][ix] = count3d_sample[iz][iy][ix];
  }
}

/* ----------------------------------------------------------------------
   copy running grid values to output grid, just for owned grid cells
   if ATOM mode, also copy per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::copy_running_to_output()
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = vec2d_running[iy][ix];
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            array2d[iy][ix][m] = array2d_running[iy][ix][m];
    }
    if (modeatom) 
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count2d[iy][ix] = count2d_running[iy][ix];

  } else {
    if (nvalues == 1) {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = vec3d_running[iz][iy][ix];
    } else {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              array3d[iz][iy][ix][m] = array3d_running[iz][iy][ix][m];
    }
    if (modeatom) 
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            count3d[iz][iy][ix] = count3d_running[iz][iy][ix];
  }
}

/* ----------------------------------------------------------------------
   sum sample grid values to running grid, just for owned grid cells
   if ATOM mode, also copy per-cell counts
------------------------------------------------------------------------- */

void FixAveGrid::sum_sample_to_running()
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d_running[iy][ix] += vec2d_sample[iy][ix];
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            array2d_running[iy][ix][m] += array2d_sample[iy][ix][m];
    }
    if (modeatom) 
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count2d_running[iy][ix] += count2d_sample[iy][ix];

  } else {
    if (nvalues == 1) {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d_running[iz][iy][ix] += vec3d_sample[iz][iy][ix];
    } else {
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              array3d_running[iz][iy][ix][m] += array3d_sample[iz][iy][ix][m];
    }
    if (modeatom) 
      for (iz = nylo_in; iz <= nyhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            count3d_running[iz][iy][ix] += count3d_sample[iz][iy][ix];
  }
}

/* ----------------------------------------------------------------------
   normalize grid values and per-cell counts for ATOM mode for owned cells
   numsamples = normalization factor
   value_sample & count are either for one sample or summed over epoch
   result = sample_value / count
     exception is DENSITY_NUMBER:
       result = value / (current binvol * Nrepeat)
     exception is DENSITY_MASS:
       result = (value * mv2d) / (current binvol * Nrepeat)
     exception is TEMPERATURE:
       result = (value * mvv2e) / (Nrepeat*cdof + adof*count) * boltz)
------------------------------------------------------------------------- */

void FixAveGrid::normalize_atom(int numsamples)
{
  int ix,iy,iz,m;
  double count,norm;

  double mvv2e = force->mvv2e;
  double mv2d = force->mv2d;
  double boltz = force->boltz;

  double *prd = domain->prd;
  double dx = prd[0]/nxgrid;
  double dy = prd[1]/nygrid;
  double dz = prd[2]/nzgrid;

  double binvol;
  if (dimension == 2) binvol = dx*dy;
  else binvol = dx*dy*dz;

  double repeat = numsamples;
  double invrepeat = 1.0/repeat;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          count = count2d_sample[iy][ix];
          if (count) {
            if (which[0] == ArgInfo::DENSITY_NUMBER)
              norm = 1.0 / (binvol * repeat);
            else if (which[0] == ArgInfo::DENSITY_MASS)
              norm = mv2d / (binvol * nrepeat);
            else if (which[0] == ArgInfo::TEMPERATURE)
              norm = mvv2e /((repeat*cdof + adof*count) * boltz);
            else if (normflag == NONORM)
              norm = 1.0/invrepeat;
            else
              norm = 1.0/count;
            vec2d_sample[iy][ix] *= norm;
            count2d_sample[iy][ix] *= invrepeat;
          }
        }
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          count = count2d_sample[iy][ix];
          if (count) {
            for (m = 0; m < nvalues; m++) {
              if (which[m] == ArgInfo::DENSITY_NUMBER)
                norm = 1.0 / (binvol * repeat);
              else if (which[m] == ArgInfo::DENSITY_MASS)
                norm = mv2d / (binvol * nrepeat);
              else if (which[m] == ArgInfo::TEMPERATURE)
                norm = mvv2e /((repeat*cdof + adof*count) * boltz);
              else if (normflag == NONORM)
                norm = 1.0/invrepeat;
              else
                norm = 1.0/count;
              array2d_sample[iy][ix][m] *= norm;
            }
            count2d_sample[iy][ix] *= invrepeat;
          }
        }
    }
  } else {
    if (nvalues == 1) {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++) {
            count = count3d_sample[iz][iy][ix];
            if (count) {
              if (which[0] == ArgInfo::DENSITY_NUMBER)
                norm = 1.0 / (binvol * repeat);
              else if (which[0] == ArgInfo::DENSITY_MASS)
                norm = mv2d / (binvol * nrepeat);
              else if (which[0] == ArgInfo::TEMPERATURE)
                norm = mvv2e /((repeat*cdof + adof*count) * boltz);
              else if (normflag == NONORM)
                norm = 1.0/invrepeat;
              else
                norm = 1.0/count;
              vec3d_sample[iz][iy][ix] *= norm;
              count3d_sample[iz][iy][ix] *= invrepeat;
            }
          }
    } else {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++) {
            count = count3d_sample[iz][iy][ix];
            if (count) {
              for (m = 0; m < nvalues; m++) {
                if (which[m] == ArgInfo::DENSITY_NUMBER)
                  norm = 1.0 / (binvol * repeat);
                else if (which[m] == ArgInfo::DENSITY_MASS)
                  norm = mv2d / (binvol * nrepeat);
                else if (which[m] == ArgInfo::TEMPERATURE)
                  norm = mvv2e /((repeat*cdof + adof*count) * boltz);
                else if (normflag == NONORM)
                  norm = 1.0/invrepeat;
                else
                  norm = 1.0/count;
                array3d_sample[iz][iy][ix][m] *= norm;
              }
              count3d_sample[iz][iy][ix] *= invrepeat;
            }
          }
    }
  }
}

/* ----------------------------------------------------------------------
   normalize grid values by Nrepeat
   used for GRID mode
   used for ATOM MODE with norm = SAMPLE
------------------------------------------------------------------------- */

void FixAveGrid::normalize_grid(int numsamples, double **v2d, double ***a2d, 
                                double ***v3d, double ****a3d)
{
  int ix,iy,iz,m;

  double invrepeat = 1.0/numsamples;
  
  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          v2d[iy][ix] *= invrepeat;
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m <= nvalues; m++)
            a2d[iy][ix][m] *= invrepeat;
    }
  } else {
    if (nvalues == 1) {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            v3d[iz][iy][ix] *= invrepeat;
    } else {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m <= nvalues; m++)
              a3d[iz][iy][ix][m] *= invrepeat;
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
   nvalues per grid point + count
------------------------------------------------------------------------- */

void FixAveGrid::pack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  int i,j,m;

  auto buf = (double *) vbuf;
  double *count,*data,*values;
  m = 0;

  if (dimension == 2) {
    count = &count2d_sample[nylo_out][nxlo_out];
    if (nvalues == 1) data = &vec2d_sample[nylo_out][nxlo_out];
    else data = &array2d_sample[nylo_out][nxlo_out][0];
  } else if (dimension == 3) {
    count = &count3d_sample[nzlo_out][nylo_out][nxlo_out];
    if (nvalues == 1) data = &vec3d_sample[nzlo_out][nylo_out][nxlo_out];
    else data = &array3d_sample[nzlo_out][nylo_out][nxlo_out][0];
  }

  if (nvalues == 1) {
    for (i = 0; i < nlist; i++) {
      buf[m++] = count[list[i]];
      buf[m++] = data[list[i]];
    }
  } else {
    for (i = 0; i < nlist; i++) {
      buf[m++] = count[list[i]];
      values = &data[nvalues*list[i]];
      for (j = 0; j < nvalues; j++)
        buf[m++] = values[j];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
   nvalues per grid point + count
------------------------------------------------------------------------- */

void FixAveGrid::unpack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  int i,j,m;

  auto buf = (double *) vbuf;
  double *count,*data,*values;
  m = 0;

  if (dimension == 2) {
    count = &count2d_sample[nylo_out][nxlo_out];
    if (nvalues == 1) data = &vec2d_sample[nylo_out][nxlo_out];
    else data = &array2d_sample[nylo_out][nxlo_out][0];
  } else if (dimension == 3) {
    count = &count3d_sample[nzlo_out][nylo_out][nxlo_out];
    if (nvalues == 1) data = &vec3d_sample[nzlo_out][nylo_out][nxlo_out];
    else data = &array3d_sample[nzlo_out][nylo_out][nxlo_out][0];
  }

  if (nvalues == 1) {
    for (i = 0; i < nlist; i++) {
      count[list[i]] += buf[m++];
      data[list[i]] += buf[m++];
    }
  } else {
    for (i = 0; i < nlist; i++) {
      count[list[i]] = buf[m++];
      values = &data[nvalues*list[i]];
      for (j = 0; j < nvalues; j++)
        values[j] = buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   called by load balancer when proc subdomains are adjusted
   not supported for now, b/c requires per-grid values to persist, i.e. a remap()
------------------------------------------------------------------------- */

void FixAveGrid::reset_grid()
{
  error->all(FLERR,"Fix ave/grid does not support load balancing (yet)");
}


/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
   return -1 if grid name not found
------------------------------------------------------------------------- */

int FixAveGrid::get_grid_by_name(const std::string &name, int &dim)
{
  if (name == "grid") {
    dim = dimension;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to Grid data struct for grid with index
   this class can store M named grids, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixAveGrid::get_grid_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2) return grid2d;
    else return grid3d;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return index of data associated with name in grid with index igrid
   this class can store M named grids, indexed 0 to M-1
   each grid can store G named data sets, indexed 0 to G-1
     a data set name can be associated with multiple grids
   set ncol for data set, 0 = vector, 1-N for array with N columns
     vector = single value per grid pt, array = N values per grid pt
   return -1 if data name not found
------------------------------------------------------------------------- */

int FixAveGrid::get_griddata_by_name(int igrid, const std::string &name, int &ncol)
{
  if ((igrid == 0) && (name == "data")) {
    if (nvalues == 1) ncol = 0;
    else ncol = nvalues;
    return 0;
  }

  // count is only produced for ATOM mode

  if (modeatom && (igrid == 0) && (name == "count")) {
    ncol = 0;
    return 1;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to multidim data array associated with index
   this class can store G named data sets, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixAveGrid::get_griddata_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2) {
      if (nvalues == 1) return vec2d;
      else return array2d;
    } else {
      if (nvalues == 1) return vec3d;
      else return array3d;
    }
  }
  if (index == 1) {
    if (dimension == 2) return count2d;
    else return count3d;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = (double) ngridout * nvalues * sizeof(double); // vec/array 2d/3d
  if (modeatom) bytes += (double) ngridout * sizeof(double);   // count2d/3d
  return bytes;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveGrid::nextvalid()
{
  bigint nvalid = (update->ntimestep/pergrid_freq)*pergrid_freq + pergrid_freq;
  if (nvalid-pergrid_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= ((bigint)nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += pergrid_freq;
  return nvalid;
}
