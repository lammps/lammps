// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

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
enum{DISCARD,KEEP};

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly

static constexpr int OFFSET = 16384;

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_bias(nullptr), which(nullptr), argindex(nullptr), ids(nullptr),
  value2index(nullptr), value2grid(nullptr), value2data(nullptr), grid2d(nullptr), grid3d(nullptr),
  grid_buf1(nullptr), grid_buf2(nullptr), grid_output(nullptr), grid_sample(nullptr),
  grid_nfreq(nullptr), grid_running(nullptr), grid_window(nullptr), grid2d_previous(nullptr),
  grid3d_previous(nullptr), grid_sample_previous(nullptr), grid_nfreq_previous(nullptr),
  grid_running_previous(nullptr), grid_window_previous(nullptr), bin(nullptr), skip(nullptr),
  vresult(nullptr)
{
  if (narg < 10) utils::missing_cmd_args(FLERR,"fix ave/grid", error);

  pergrid_flag = 1;
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  pergrid_freq = utils::inumeric(FLERR,arg[5],false,lmp);
  time_depend = 1;

  if (nevery <= 0 || nrepeat <= 0 || pergrid_freq <= 0)
    error->all(FLERR,"Illegal fix ave/grid command");
  if (pergrid_freq % nevery || nrepeat*nevery > pergrid_freq)
    error->all(FLERR,"Illegal fix ave/grid command");

  // NOTE: allow Dxyz as well at some point ?

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

  // parse values until one isn't recognized (optional keyword)

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

      // if arg is not a per-atom or per-grid value
      // then it's an optional arg after the values

      ArgInfo argi(arg[iarg]);
      if (argi.get_type() == ArgInfo::NONE || argi.get_type() == ArgInfo::UNKNOWN) break;
      if (argi.get_dim() > 1) error->all(FLERR,"Invalid fix ave/grid command");

      // atom value has no colon

      if (!strchr(arg[iarg],':')) {
        modeatom = 1;
        ids[nvalues] = argi.copy_name();
        which[nvalues] = argi.get_type();
        argindex[nvalues] = argi.get_index1();

      // per-grid value has colons

      } else {
        modegrid = 1;

        int igrid,idata,index;
        int iflag =
          utils::check_grid_reference((char *) "Fix ave/grid",
                                      arg[iarg],nevery,ids[nvalues],igrid,idata,index,lmp);
        if (iflag < 0) error->all(FLERR,"Invalid grid reference in fix ave/grid command");

        which[nvalues] = iflag;
        value2grid[nvalues] = igrid;
        value2data[nvalues] = idata;
        argindex[nvalues] = index;
      }
    }

    nvalues++;
    iarg++;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/grid command");

  if (modeatom && modegrid)
    error->all(FLERR,"Fix ave/grid cannot operate on per-atom and "
               "per-grid values");

  // optional args

  discardflag = DISCARD;
  normflag = ALL;
  aveflag = ONE;
  nwindow = 0;
  biasflag = 0;
  id_bias = nullptr;
  adof = domain->dimension;
  cdof = 0.0;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"discard") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) discardflag = DISCARD;
      else if (strcmp(arg[iarg+1],"no") == 0) discardflag = KEEP;
      else error->all(FLERR,"Illegal fix ave/grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"norm") == 0) {
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

  // more error checks
  // for fix inputs, check that fix frequency is acceptable

  dimension = domain->dimension;

  if ((nxgrid < 1) || (nygrid < 1) || (nzgrid < 1))
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
          error->all(FLERR, "Fix ave/atom compute does not calculate per-atom values");
        if (argindex[i] == 0 &&
            modify->compute[icompute]->size_peratom_cols != 0)
          error->all(FLERR,"Fix ave/atom compute does not calculate a per-atom vector");
        if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
          error->all(FLERR,"Fix ave/atom compute does not calculate a per-atom array");
        if (argindex[i] && argindex[i] > modify->compute[icompute]->size_peratom_cols)
          error->all(FLERR,"Fix ave/atom compute array is accessed out-of-range");

      } else if (which[i] == ArgInfo::FIX) {
        int ifix = modify->find_fix(ids[i]);
        if (ifix < 0)
          error->all(FLERR,"Fix ID for fix ave/atom does not exist");
        if (modify->fix[ifix]->peratom_flag == 0)
          error->all(FLERR,"Fix ave/atom fix does not calculate per-atom values");
        if (argindex[i] == 0 && modify->fix[ifix]->size_peratom_cols != 0)
          error->all(FLERR, "Fix ave/atom fix does not calculate a per-atom vector");
        if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
          error->all(FLERR, "Fix ave/atom fix does not calculate a per-atom array");
        if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
          error->all(FLERR,"Fix ave/atom fix array is accessed out-of-range");
        if (nevery % modify->fix[ifix]->peratom_freq)
          error->all(FLERR, "Fix for fix ave/atom not computed at compatible time");

      } else if (which[i] == ArgInfo::VARIABLE) {
        int ivariable = input->variable->find(ids[i]);
        if (ivariable < 0)
          error->all(FLERR,"Variable name for fix ave/atom does not exist");
        if (input->variable->atomstyle(ivariable) == 0)
          error->all(FLERR,"Fix ave/atom variable is not atom-style variable");
      }
    }
  }

  // instantiate Grid class and buffers
  // allocate/zero per-grid data

  allocate_grid();

  grid_sample = allocate_one_grid();
  grid_nfreq = allocate_one_grid();
  if (aveflag == RUNNING || aveflag == WINDOW) grid_running = allocate_one_grid();
  if (aveflag == WINDOW) {
    grid_window = new GridData*[nwindow];
    for (int i = 0; i < nwindow; i++)
      grid_window[i] = allocate_one_grid();
  }

  // output may occur via dump on timestep 0

  grid_output = new GridData();
  output_grid(grid_nfreq);

  // initialize running and window values

  running_count = 0;
  window_count = 0;
  window_oldest = -1;
  window_newest = 0;

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

  // deallocate Grid class and buffers

  if (dimension == 2) delete grid2d;
  else delete grid3d;

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  // deallocate per-grid data

  deallocate_one_grid(grid_sample,nxlo_out,nylo_out,nzlo_out);
  deallocate_one_grid(grid_nfreq,nxlo_out,nylo_out,nzlo_out);
  if (aveflag == RUNNING || aveflag == WINDOW)
    deallocate_one_grid(grid_running,nxlo_out,nylo_out,nzlo_out);
  if (aveflag == WINDOW) {
    for (int i = 0; i < nwindow; i++)
      deallocate_one_grid(grid_window[i],nxlo_out,nylo_out,nzlo_out);
    delete [] grid_window;
  }

  delete grid_output;

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
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // if first sample in nfreq, zero owned and ghost grid points

  if (irepeat == 0) {
    zero_grid(grid_sample);
    if (normflag == SAMPLE) zero_grid(grid_nfreq);
  }

  // accumulate per-grid values for one sample for either ATOM or GRID mode
  // per-atom compute/fix/variable may invoke computes so wrap with clear/add

  if (modeatom) {
    modify->clearstep_compute();
    atom2grid();
  } else {
    grid2grid();
  }

  // return if irepeat < nrepeat
  // unless ATOM mode and norm = SAMPLE

  irepeat++;
  if (irepeat < nrepeat && (modegrid || normflag != SAMPLE)) {
    nvalid += nevery;
    if (modeatom) modify->addstep_compute(nvalid);
    return;
  }

  // for ATOM mode, perform ghost to owned grid communication
  //   done once per Nfreq for norm = ONE
  //   done every sample for norm = SAMPLE
  // nvalues+1 includes atom count

  if (modeatom) {
    if (dimension == 2)
      grid2d->reverse_comm(Grid2d::FIX,this,0,nvalues+1,sizeof(double),
                           grid_buf1,grid_buf2,MPI_DOUBLE);
    else
      grid3d->reverse_comm(Grid3d::FIX,this,0,nvalues+1,sizeof(double),
                           grid_buf1,grid_buf2,MPI_DOUBLE);
  }

  // if ATOM mode and norm = SAMPLE:
  //   normalize the single sample
  //   sum sample grid to Nfreq grid
  //   return if irepeat < nrepeat

  if (modeatom && normflag == SAMPLE) {
    normalize_atom(1,grid_sample);
    add_grid(grid_sample,grid_nfreq);
    zero_grid(grid_sample);

    if (irepeat < nrepeat) {
      nvalid += nevery;
      modify->addstep_compute(nvalid);
      return;
    }
  }

  // this is an Nfreq timestep
  // reset irepeat and nvalid

  irepeat = 0;
  nvalid = ntimestep+pergrid_freq - ((bigint) nrepeat-1)*nevery;
  if (modeatom) modify->addstep_compute(nvalid);

  // just return if this proc owns no grid points

  if (ngridout == 0) return;

  // average final results over the entire Nfreq of samples
  // store final result in Nfreq grid
  // for ATOM mode:
  //   for norm = ALL, normalize sample grid by counts over all samples
  //   for norm = SAMPLE, normalize Nfreq grid by Nrepeat
  //   for norm = NONORM, normalize sample grid by Nrepeat, not by counts
  //              this check is made inside normalize_atom()
  // for GRID mode:
  //   normalize sample grid by Nrepeat

  if (modeatom) {
    if (normflag == ALL) {
      normalize_atom(nrepeat,grid_sample);
      normalize_count(nrepeat,grid_sample);
      copy_grid(grid_sample,grid_nfreq);

    } else if (normflag == SAMPLE) {
      normalize_grid(nrepeat,grid_nfreq);
      normalize_count(nrepeat,grid_nfreq);

    } else if (normflag == NONORM) {
      normalize_atom(nrepeat,grid_sample);
      normalize_count(nrepeat,grid_sample);
      copy_grid(grid_sample,grid_nfreq);
    }
  }

  if (modegrid) {
    normalize_grid(nrepeat,grid_sample);
    copy_grid(grid_sample,grid_nfreq);
  }

  // create Nfreq output
  // for aveflag == RUNNING
  //   running_count = # of Nfreq entries in grid_running
  // for aveflag == WINDOW
  //   window_count = # of Nfreq entries in grid_window

  if (aveflag == ONE) {
    output_grid(grid_nfreq);

  } else if (aveflag == RUNNING) {
    running_count++;
    add_grid(grid_nfreq,grid_running);
    copy_grid(grid_running,grid_nfreq);
    normalize_grid(running_count,grid_nfreq);
    if (modeatom) normalize_count(running_count,grid_nfreq);
    output_grid(grid_nfreq);

  } else if (aveflag == WINDOW) {

    // update grid_running to be sum over grid_window entries
    //   add grid_nfreq to grid_running
    //   if window is full, subtract oldest window entry from grid_running
    // copy grid_nfreq into window

    add_grid(grid_nfreq,grid_running);
    if (window_count == nwindow)
      subtract_grid(grid_window[window_oldest],grid_running);
    copy_grid(grid_nfreq,grid_window[window_newest]);

    // update status of window
    // window_count = # of entries in window
    // window_oldest = index of oldest entry in grid_window
    // window_newest = index in window where next grid_nfreq will be copied

    if (window_count < nwindow) window_count++;
    if (window_count == nwindow) window_oldest++;
    if (window_oldest == nwindow) window_oldest = 0;
    window_newest++;
    if (window_newest == nwindow) window_newest = 0;

    // copy grid running to grid_nfreq and perform normalization

    copy_grid(grid_running,grid_nfreq);
    normalize_grid(window_count,grid_nfreq);
    output_grid(grid_nfreq);
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

  double **count2d = grid_sample->count2d;
  double **vec2d = grid_sample->vec2d;
  double ***array2d = grid_sample->array2d;
  double ***count3d = grid_sample->count3d;
  double ***vec3d = grid_sample->vec3d;
  double ****array3d = grid_sample->array3d;

  // bin[i][dim] = indices of bin each atom is in
  // skip atom if group mask does not match
  // check if any atom is out of bounds for my local grid
  // for nonperiodic dim, remap atom to first/last bin if out of bounds
  // count atoms contributing to each bin

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

      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;

      if (ix < nxlo_out || ix > nxhi_out) {
        if (periodicity[0]) flag = 1;
        else if (discardflag == KEEP) {
          if (ix < nxlo_out && nxlo_out == 0) ix = 0;
          else if (ix > nxhi_out && nxhi_out == nxgrid-1) ix = nxgrid-1;
          else flag = 1;
        } else skip[i] = 1;
        continue;
      }
      if (iy < nylo_out || iy > nyhi_out) {
        if (periodicity[1]) flag = 1;
        else if (discardflag == KEEP) {
          if (iy < nylo_out && nylo_out == 0) iy = 0;
          else if (iy > nyhi_out && nyhi_out == nygrid-1) iy = nygrid-1;
          else flag = 1;
        } else skip[i] = 1;
        continue;
      }

      skip[i] = 0;
      count2d[iy][ix] += 1.0;
      bin[i][0] = iy;
      bin[i][1] = ix;
    }

  } else {
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) {
        skip[i] = 1;
        continue;
      }

      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + OFFSET) - OFFSET;

      if (ix < nxlo_out || ix > nxhi_out) {
        if (periodicity[0]) flag = 1;
        else if (discardflag == KEEP) {
          if (ix < nxlo_out && nxlo_out == 0) ix = 0;
          else if (ix > nxhi_out && nxhi_out == nxgrid-1) ix = nxgrid-1;
          else flag = 1;
        } else skip[i] = 1;
      }

      if (iy < nylo_out || iy > nyhi_out) {
        if (periodicity[1]) flag = 1;
        else if (discardflag == KEEP) {
          if (iy < nylo_out && nylo_out == 0) iy = 0;
          else if (iy > nyhi_out && nyhi_out == nygrid-1) iy = nygrid-1;
          else flag = 1;
        } else skip[i] = 1;
        continue;
      }

      if (iz < nzlo_out || iz > nzhi_out) {
        if (periodicity[2]) flag = 1;
        else if (discardflag == KEEP) {
          if (iz < nzlo_out && nzlo_out == 0) iz = 0;
          else if (iz > nzhi_out && nzhi_out == nzgrid-1) iz = nzgrid-1;
          else flag = 1;
        } else skip[i] = 1;
        continue;
      }

      skip[i] = 0;
      count3d[iz][iy][ix] += 1.0;
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
              vec2d[bin[i][0]][bin[i][1]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              array2d[bin[i][0]][bin[i][1]][m] += attribute[i][j];
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i])
              array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += attribute[i][j];
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
              vec2d[bin[i][0]][bin[i][1]] += one;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array2d[bin[i][0]][bin[i][1]][m] += one;
            }
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += one;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              if (which[m] == ArgInfo::DENSITY_NUMBER) one = 1.0;
              else if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += one;
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
              vec2d[bin[i][0]][bin[i][1]] += one*vsq;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array2d[bin[i][0]][bin[i][1]][m] += one*vsq;;
            }
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += one*vsq;
            }
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (!skip[i]) {
              vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
              if (rmass) one = rmass[i];
              else one = mass[type[i]];
              array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += one*vsq;
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
                vec2d[bin[i][0]][bin[i][1]] += ovector[i];
            }
          } else {
            int jm1 = j = 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec2d[bin[i][0]][bin[i][1]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array2d[bin[i][0]][bin[i][1]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array2d[bin[i][0]][bin[i][1]][m] += oarray[i][jm1];
            }
          }
        }

      } else {
        if (nvalues == 1) {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (!skip[i])
                array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += oarray[i][jm1];
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

  double **vec2d = grid_sample->vec2d;
  double ***array2d = grid_sample->array2d;
  double ***vec3d = grid_sample->vec3d;
  double ****array3d = grid_sample->array3d;

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
              vec2d[iy][ix] += ovec2d[iy][ix];
        } else {
          int jm1 = j - 1;
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              vec2d[iy][ix] += oarray2d[iy][ix][jm1];
        }
      } else {
        if (j == 0) {
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              array2d[iy][ix][m] += ovec2d[iy][ix];
        } else {
          int jm1 = j - 1;
          for (iy = nylo_in; iy <= nyhi_in; iy++)
            for (ix = nxlo_in; ix <= nxhi_in; ix++)
              array2d[iy][ix][m] += oarray2d[iy][ix][jm1];
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
                vec3d[iz][iy][ix] += ovec3d[iz][iy][ix];
        } else {
          int jm1 = j - 1;
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                vec3d[iz][iy][ix] += oarray3d[iz][iy][ix][jm1];
        }
      } else {
        if (j == 0) {
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                array3d[iz][iy][ix][m] += ovec3d[iz][iy][ix];
        } else {
          int jm1 = j - 1;
          for (iz = nzlo_in; iz <= nzhi_in; iz++)
            for (iy = nylo_in; iy <= nyhi_in; iy++)
              for (ix = nxlo_in; ix <= nxhi_in; ix++)
                array3d[iz][iy][ix][m] += oarray3d[iz][iy][ix][jm1];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   normalize grid values for ATOM mode for owned cells
   grid values are summed over numsamples (can be 1 or Nrepeat)
   result = sample_value / count
     exception is DENSITY_NUMBER:
       result = value / (current binvol * Nrepeat)
     exception is DENSITY_MASS:
       result = (value * mv2d) / (current binvol * Nrepeat)
     exception is TEMPERATURE:
       result = (value * mvv2e) / (Nrepeat*cdof + adof*count) * boltz)
   exception normalization is same for norm = ALL, SAMPLE, NONORM
     so NONORM if test is after exception if tests
------------------------------------------------------------------------- */

void FixAveGrid::normalize_atom(int numsamples, GridData *grid)
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

  double repeat = numsamples;
  double invrepeat = 1.0/repeat;

  double binvol;
  if (dimension == 2) binvol = dx*dy;
  else binvol = dx*dy*dz;

  double density_number_norm = 1.0 / (binvol * repeat);
  double density_mass_norm = mv2d / (binvol * repeat);

  if (dimension == 2) {
    double **count2d = grid->count2d;

    if (nvalues == 1) {
      double **vec2d = grid->vec2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          count = count2d[iy][ix];
          if (count) {
            if (which[0] == ArgInfo::DENSITY_NUMBER)
              norm = density_number_norm;
            else if (which[0] == ArgInfo::DENSITY_MASS)
              norm = density_mass_norm;
            else if (which[0] == ArgInfo::TEMPERATURE)
              norm = mvv2e / ((repeat*cdof + adof*count) * boltz);
            else if (normflag == NONORM)
              norm = invrepeat;
            else
              norm = 1.0/count;
            vec2d[iy][ix] *= norm;
          }
        }

    } else {
      double ***array2d = grid->array2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          count = count2d[iy][ix];
          if (count) {
            for (m = 0; m < nvalues; m++) {
              if (which[m] == ArgInfo::DENSITY_NUMBER)
                norm = density_number_norm;
              else if (which[m] == ArgInfo::DENSITY_MASS)
                norm = density_mass_norm;
              else if (which[m] == ArgInfo::TEMPERATURE)
                norm = mvv2e / ((repeat*cdof + adof*count) * boltz);
              else if (normflag == NONORM)
                norm = invrepeat;
              else
                norm = 1.0/count;
              array2d[iy][ix][m] *= norm;
            }
          }
        }
    }

  } else if (dimension == 3) {
    double ***count3d = grid->count3d;

    if (nvalues == 1) {
      double ***vec3d = grid->vec3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++) {
            count = count3d[iz][iy][ix];
            if (count) {
              if (which[0] == ArgInfo::DENSITY_NUMBER)
                norm = density_number_norm;
              else if (which[0] == ArgInfo::DENSITY_MASS)
                norm = density_mass_norm;
              else if (which[0] == ArgInfo::TEMPERATURE)
                norm = mvv2e / ((repeat*cdof + adof*count) * boltz);
              else if (normflag == NONORM)
                norm = invrepeat;
              else
                norm = 1.0/count;
              vec3d[iz][iy][ix] *= norm;
            }
          }

    } else {
      double ****array3d = grid->array3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++) {
            count = count3d[iz][iy][ix];
            if (count) {
              for (m = 0; m < nvalues; m++) {
                if (which[m] == ArgInfo::DENSITY_NUMBER)
                  norm = density_number_norm;
                else if (which[m] == ArgInfo::DENSITY_MASS)
                norm = density_mass_norm;
                else if (which[m] == ArgInfo::TEMPERATURE)
                  norm = mvv2e / ((repeat*cdof + adof*count) * boltz);
                else if (normflag == NONORM)
                  norm = invrepeat;
                else
                  norm = 1.0/count;
                array3d[iz][iy][ix][m] *= norm;
              }
            }
          }
    }
  }
}

/* ----------------------------------------------------------------------
   normalize grid values by numsamples
   used for ATOM MODE when norm = SAMPLE
   used for GRID mode
------------------------------------------------------------------------- */

void FixAveGrid::normalize_grid(int numsamples, GridData *grid)
{
  int ix,iy,iz,m;

  double invrepeat = 1.0/numsamples;

  if (dimension == 2) {
    if (nvalues == 1) {
      double **vec2d = grid->vec2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] *= invrepeat;
    } else {
      double ***array2d = grid->array2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m < nvalues; m++)
            array2d[iy][ix][m] *= invrepeat;
    }

  } else if (dimension == 3) {
    if (nvalues == 1) {
      double ***vec3d = grid->vec3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] *= invrepeat;
    } else {
      double ****array3d = grid->array3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m < nvalues; m++)
              array3d[iz][iy][ix][m] *= invrepeat;
    }
  }
}

/* ----------------------------------------------------------------------
   normalize grid counts by numsamples
   used for ATOM mode
------------------------------------------------------------------------- */

void FixAveGrid::normalize_count(int numsamples, GridData *grid)
{
  int ix,iy,iz;

  double invrepeat = 1.0/numsamples;

  if (dimension == 2) {
    double **count2d = grid->count2d;
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        count2d[iy][ix] *= invrepeat;

  } else if (dimension == 3) {
    double ***count3d = grid->count3d;
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          count3d[iz][iy][ix] *= invrepeat;
  }
}

/* ----------------------------------------------------------------------
   allocate instance of Grid2d or Grid3d
------------------------------------------------------------------------- */

void FixAveGrid::allocate_grid()
{
  if (modeatom) maxdist = 0.5 * neighbor->skin;
  else if (modegrid) maxdist = 0.0;

  if (dimension == 2) {
    grid2d = new Grid2d(lmp, world, nxgrid, nygrid);
    grid2d->set_distance(maxdist);
    grid2d->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in,
                       nxlo_out, nxhi_out, nylo_out, nyhi_out);

    // ngrid_buf12 converted to nvalues + count

    grid2d->setup_comm(ngrid_buf1, ngrid_buf2);
    ngrid_buf1 *= nvalues + 1;
    ngrid_buf2 *= nvalues + 1;

    grid_buf1 = grid_buf2 = nullptr;
    if (ngrid_buf1) memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    if (ngrid_buf2) memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1);

  } else {
    grid3d = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
    grid3d->set_distance(maxdist);
    grid3d->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                       nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

    // ngrid_buf12 converted to nvalues + count

    grid3d->setup_comm(ngrid_buf1, ngrid_buf2);
    ngrid_buf1 *= nvalues + 1;
    ngrid_buf2 *= nvalues + 1;

    grid_buf1 = grid_buf2 = nullptr;
    if (ngrid_buf1) memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    if (ngrid_buf2) memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
      (nzhi_out - nzlo_out + 1);
  }
}

/* ----------------------------------------------------------------------
   allocate a data grid and zero its values
   if ATOM mode, also allocate per-grid count
------------------------------------------------------------------------- */

FixAveGrid::GridData *FixAveGrid::allocate_one_grid()
{
  GridData *grid = new GridData();

  grid->vec2d = nullptr;
  grid->array2d = nullptr;
  grid->count2d = nullptr;
  grid->vec3d = nullptr;
  grid->array3d = nullptr;
  grid->count3d = nullptr;

  if (dimension == 2) {
    if (nvalues == 1)
      memory->create2d_offset(grid->vec2d, nylo_out, nyhi_out,
                              nxlo_out, nxhi_out, "ave/grid:vec2d");
    else
      memory->create3d_offset_last(grid->array2d, nylo_out, nyhi_out,
                                   nxlo_out, nxhi_out,
                                   nvalues, "ave/grid:array3d");
    if (modeatom)
      memory->create2d_offset(grid->count2d, nylo_out, nyhi_out,
                              nxlo_out, nxhi_out, "ave/grid:count2d");

  } else if (dimension == 3) {
    if (nvalues == 1)
      memory->create3d_offset(grid->vec3d, nzlo_out, nzhi_out, nylo_out,
                              nyhi_out, nxlo_out, nxhi_out, "ave/grid:vec3d");
    else
      memory->create4d_offset_last(grid->array3d, nzlo_out, nzhi_out, nylo_out,
                                   nyhi_out, nxlo_out, nxhi_out, nvalues,
                                   "ave/grid:array3d");
    if (modeatom)
      memory->create3d_offset(grid->count3d, nzlo_out, nzhi_out, nylo_out,
                              nyhi_out, nxlo_out, nxhi_out, "ave/grid:count3d");
  }

  zero_grid(grid);

  return grid;
}


/* ----------------------------------------------------------------------
   create clone of a data grid
   allocate a new grid and copy only the pointers from the source grid
   used by reset_grid() to keep old data values until remap is complete
------------------------------------------------------------------------- */

FixAveGrid::GridData *FixAveGrid::clone_one_grid(GridData *src)
{
  GridData *grid = new GridData();

  grid->vec2d = src->vec2d;
  grid->array2d = src->array2d;
  grid->count2d = src->count2d;
  grid->vec3d = src->vec3d;
  grid->array3d = src->array3d;
  grid->count3d = src->count3d;

  return grid;
}

/* ----------------------------------------------------------------------
   deallocate a data grid and all its memory
   if ATOM mode, also deallocate per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::deallocate_one_grid(GridData *grid,
                                     int xoffset, int yoffset, int zoffset)
{
  if (dimension == 2) {
    if (nvalues == 1)
      memory->destroy2d_offset(grid->vec2d,yoffset,xoffset);
    else
      memory->destroy3d_offset_last(grid->array2d,yoffset,xoffset);
    if (modeatom)
      memory->destroy2d_offset(grid->count2d,yoffset,xoffset);

  } else if (dimension == 3) {
    if (nvalues == 1)
      memory->destroy3d_offset(grid->vec3d,zoffset,yoffset,xoffset);
    else
      memory->destroy4d_offset_last(grid->array3d,zoffset,yoffset,xoffset);
    if (modeatom)
      memory->destroy3d_offset(grid->count3d,zoffset,yoffset,xoffset);
  }

  delete grid;
}

/* ----------------------------------------------------------------------
   size of a data grid
   if ATOM mode, also include per-grid count
------------------------------------------------------------------------- */

double FixAveGrid::size_grid(GridData * /*grid*/)
{
  int nper = nvalues;
  if (modeatom) nper++;

  double bytes;
  if (dimension == 2)
    bytes = sizeof(double) * nper * (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1);
  else
    bytes = sizeof(double) * nper * (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
      (nzhi_out - nzlo_out + 1);

  return bytes;
}

/* ----------------------------------------------------------------------
   zero values for a data grid including ghost cells
   for ATOM mode, also zero per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::zero_grid(GridData *grid)
{
  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1)
      memset(&grid->vec2d[nylo_out][nxlo_out], 0, sizeof(double) * ngridout);
    else
      memset(&grid->array2d[nylo_out][nxlo_out][0], 0, sizeof(double) * ngridout * nvalues);
    if (modeatom) memset(&grid->count2d[nylo_out][nxlo_out], 0, sizeof(double) * ngridout);

  } else {
    if (nvalues == 1)
      memset(&grid->vec3d[nzlo_out][nylo_out][nxlo_out], 0, sizeof(double) * ngridout);
    else
      memset(&grid->array3d[nzlo_out][nylo_out][nxlo_out][0], 0, sizeof(double)* ngridout * nvalues);
    if (modeatom)
      memset(&grid->count3d[nzlo_out][nylo_out][nxlo_out], 0, sizeof(double) * ngridout);
  }
}

/* ----------------------------------------------------------------------
   copy src grid values to result grid, just for owned grid cells
   for ATOM mode, also copy per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::copy_grid(GridData *src, GridData *result)
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      double **vsrc = src->vec2d;
      double **vresult = result->vec2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vresult[iy][ix] = vsrc[iy][ix];
    } else {
      double ***asrc = src->array2d;
      double ***aresult = result->array2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m < nvalues; m++)
            aresult[iy][ix][m] = asrc[iy][ix][m];
    }
    if (modeatom) {
      double **csrc = src->count2d;
      double **cresult = result->count2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          cresult[iy][ix] = csrc[iy][ix];
    }

  } else if (dimension == 3) {
    if (nvalues == 1) {
      double ***vsrc = src->vec3d;
      double ***vresult = result->vec3d;
        for (iz = nzlo_in; iz <= nzhi_in; iz++)
          for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vresult[iz][iy][ix] = vsrc[iz][iy][ix];
    } else {
      double ****asrc = src->array3d;
      double ****aresult = result->array3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m < nvalues; m++)
              aresult[iz][iy][ix][m] = asrc[iz][iy][ix][m];
    }
    if (modeatom) {
      double ***csrc = src->count3d;
      double ***cresult = result->count3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            cresult[iz][iy][ix] = csrc[iz][iy][ix];
    }
  }
}

/* ----------------------------------------------------------------------
   add values in src grid to result grid, just for owned grid cells
   for ATOM mode, also sum per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::add_grid(GridData *src, GridData *result)
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      double **vsrc = src->vec2d;
      double **vresult = result->vec2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vresult[iy][ix] += vsrc[iy][ix];
    } else {
      double ***asrc = src->array2d;
      double ***aresult = result->array2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m < nvalues; m++)
            aresult[iy][ix][m] += asrc[iy][ix][m];
    }
    if (modeatom) {
      double **csrc = src->count2d;
      double **cresult = result->count2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          cresult[iy][ix] += csrc[iy][ix];
    }

  } else if (dimension == 3) {
    if (nvalues == 1) {
      double ***vsrc = src->vec3d;
      double ***vresult = result->vec3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vresult[iz][iy][ix] += vsrc[iz][iy][ix];
    } else {
      double ****asrc = src->array3d;
      double ****aresult = result->array3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m < nvalues; m++)
              aresult[iz][iy][ix][m] += asrc[iz][iy][ix][m];
    }
    if (modeatom) {
      double ***csrc = src->count3d;
      double ***cresult = result->count3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            cresult[iz][iy][ix] += csrc[iz][iy][ix];
    }
  }
}

/* ----------------------------------------------------------------------
   subtact values in src grid from result grid, just for owned grid cells
   for ATOM mode, also sum per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::subtract_grid(GridData *src, GridData *result)
{
  int ix,iy,iz,m;

  if (!ngridout) return;

  if (dimension == 2) {
    if (nvalues == 1) {
      double **vsrc = src->vec2d;
      double **vresult = result->vec2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          vresult[iy][ix] -= vsrc[iy][ix];
    } else {
      double ***asrc = src->array2d;
      double ***aresult = result->array2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          for (m = 0; m < nvalues; m++)
            aresult[iy][ix][m] -= asrc[iy][ix][m];
    }
    if (modeatom) {
      double **csrc = src->count2d;
      double **cresult = result->count2d;
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          cresult[iy][ix] -= csrc[iy][ix];
    }

  } else if (dimension == 3) {
    if (nvalues == 1) {
      double ***vsrc = src->vec3d;
      double ***vresult = result->vec3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            vresult[iz][iy][ix] -= vsrc[iz][iy][ix];
    } else {
      double ****asrc = src->array3d;
      double ****aresult = result->array3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            for (m = 0; m < nvalues; m++)
              aresult[iz][iy][ix][m] -= asrc[iz][iy][ix][m];
    }
    if (modeatom) {
      double ***csrc = src->count3d;
      double ***cresult = result->count3d;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            cresult[iz][iy][ix] -= csrc[iz][iy][ix];
    }
  }
}

/* ----------------------------------------------------------------------
   set output grid pointers to src grid data
   for ATOM mode, also set pointers to per-grid count
------------------------------------------------------------------------- */

void FixAveGrid::output_grid(GridData *src)
{
  if (dimension == 2) {
    if (nvalues == 1)
      grid_output->vec2d = src->vec2d;
    else
      grid_output->array2d = src->array2d;
    if (modeatom)
      grid_output->count2d = src->count2d;

  } else if (dimension == 3) {
    if (nvalues == 1)
      grid_output->vec3d = src->vec3d;
    else
      grid_output->array3d = src->array3d;
    if (modeatom)
      grid_output->count3d = src->count3d;
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
   nvalues per grid point + count-
   only invoked for ATOM mode
------------------------------------------------------------------------ */

void FixAveGrid::pack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  int i,j,m;

  auto buf = (double *) vbuf;
  double *count,*data,*values;
  m = 0;

  if (dimension == 2) {
    count = &grid_sample->count2d[nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid_sample->vec2d[nylo_out][nxlo_out];
    else data = &grid_sample->array2d[nylo_out][nxlo_out][0];
  } else {
    count = &grid_sample->count3d[nzlo_out][nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid_sample->vec3d[nzlo_out][nylo_out][nxlo_out];
    else data = &grid_sample->array3d[nzlo_out][nylo_out][nxlo_out][0];
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
   only invoked for ATOM mode
------------------------------------------------------------------------- */

void FixAveGrid::unpack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  int i,j,m;

  auto buf = (double *) vbuf;
  double *count,*data,*values;

  if (dimension == 2) {
    count = &grid_sample->count2d[nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid_sample->vec2d[nylo_out][nxlo_out];
    else data = &grid_sample->array2d[nylo_out][nxlo_out][0];
  } else {
    count = &grid_sample->count3d[nzlo_out][nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid_sample->vec3d[nzlo_out][nylo_out][nxlo_out];
    else data = &grid_sample->array3d[nzlo_out][nylo_out][nxlo_out][0];
  }

  m = 0;
  if (nvalues == 1) {
    for (i = 0; i < nlist; i++) {
      count[list[i]] += buf[m++];
      data[list[i]] += buf[m++];
    }
  } else {
    for (i = 0; i < nlist; i++) {
      count[list[i]] += buf[m++];
      values = &data[nvalues*list[i]];
      for (j = 0; j < nvalues; j++)
        values[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack old grid values to buf to send to another proc
   invoked for both GRID and ATOM mode
------------------------------------------------------------------------- */

void FixAveGrid::pack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;

  int running_flag = 0;
  if (aveflag == RUNNING || aveflag == WINDOW) running_flag = 1;
  int window_flag = 0;
  if (aveflag == WINDOW) window_flag = 1;

  int m = 0;
  for (int i = 0; i < nlist; i++) {
    m += pack_one_grid(grid_sample_previous,list[i],&buf[m]);
    m += pack_one_grid(grid_nfreq_previous,list[i],&buf[m]);
    if (running_flag) m += pack_one_grid(grid_running_previous,list[i],&buf[m]);
    if (window_flag)
      for (int iwindow = 0; iwindow < nwindow; iwindow++)
        m += pack_one_grid(grid_window_previous[iwindow],list[i],&buf[m]);
  }
}

/* ----------------------------------------------------------------------
   unpack received owned values from buf into new grid
   invoked for both GRID and ATOM mode
------------------------------------------------------------------------- */

void FixAveGrid::unpack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
   auto buf = (double *) vbuf;

  int running_flag = 0;
  if (aveflag == RUNNING || aveflag == WINDOW) running_flag = 1;
  int window_flag = 0;
  if (aveflag == WINDOW) window_flag = 1;

  int m = 0;
  for (int i = 0; i < nlist; i++) {
    m += unpack_one_grid(&buf[m],grid_sample,list[i]);
    m += unpack_one_grid(&buf[m],grid_nfreq,list[i]);
    if (running_flag) m += unpack_one_grid(&buf[m],grid_running,list[i]);
    if (window_flag)
      for (int iwindow = 0; iwindow < nwindow; iwindow++)
        m += unpack_one_grid(&buf[m],grid_window[iwindow],list[i]);
  }
}

/* ----------------------------------------------------------------------
   pack values for a single grid cell to buf
   return number of values packed
------------------------------------------------------------------------- */

int FixAveGrid::pack_one_grid(GridData *grid, int index, double *buf)
{
  double *count,*data,*values;

  if (dimension == 2) {
    count = &grid->count2d[nylo_out_previous][nxlo_out_previous];
    if (nvalues == 1) data = &grid->vec2d[nylo_out_previous][nxlo_out_previous];
    else data = &grid->array2d[nylo_out_previous][nxlo_out_previous][0];
  } else {
    count = &grid->count3d[nzlo_out_previous][nylo_out_previous][nxlo_out_previous];
    if (nvalues == 1) data = &grid->vec3d[nzlo_out_previous][nylo_out_previous][nxlo_out_previous];
    else data = &grid->array3d[nzlo_out_previous][nylo_out_previous][nxlo_out_previous][0];
  }

  int m = 0;
  if (modeatom) buf[m++] = count[index];
  if (nvalues == 1) buf[m++] = data[index];
  else {
    values = &data[nvalues*index];
    for (int j = 0; j < nvalues; j++)
      buf[m++] = values[j];
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for a single grid cell from buf
   return number of values unpacked
------------------------------------------------------------------------- */

int FixAveGrid::unpack_one_grid(double *buf, GridData *grid, int index)
{
  double *count,*data,*values;

  if (dimension == 2) {
    count = &grid->count2d[nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid->vec2d[nylo_out][nxlo_out];
    else data = &grid->array2d[nylo_out][nxlo_out][0];
  } else {
    count = &grid->count3d[nzlo_out][nylo_out][nxlo_out];
    if (nvalues == 1) data = &grid->vec3d[nzlo_out][nylo_out][nxlo_out];
    else data = &grid->array3d[nzlo_out][nylo_out][nxlo_out][0];
  }

  int m = 0;
  if (modeatom) count[index] = buf[m++];
  if (nvalues == 1) data[index] = buf[m++];
  else {
    values = &data[nvalues*index];
    for (int j = 0; j < nvalues; j++)
      values[j] = buf[m++];
  }

  return m;
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   called by load balancer when proc subdomains are adjusted
   persist per-grid data by performing a grid remap
------------------------------------------------------------------------- */

void FixAveGrid::reset_grid()
{
  // check if new grid partitioning is different on any proc
  // if not, just return

  if (dimension == 2) {
    int tmp[8];
    Grid2d *gridnew = new Grid2d(lmp, world, nxgrid, nygrid);
    gridnew->set_distance(maxdist);
    gridnew->setup_grid(tmp[0], tmp[1], tmp[2], tmp[3],
                        tmp[4], tmp[5], tmp[6], tmp[7]);

    if (grid2d->identical(gridnew)) {
      delete gridnew;
      return;
    } else delete gridnew;

  } else {

    int tmp[12];
    Grid3d *gridnew = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
    gridnew->set_distance(maxdist);
    gridnew->setup_grid(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5],
                        tmp[6], tmp[7], tmp[8], tmp[9], tmp[10], tmp[11]);

    if (grid3d->identical(gridnew)) {
      delete gridnew;
      return;
    } else delete gridnew;
  }

  // for ATOM mode, perform ghost to owned grid comm for grid_sample
  // necessary b/c remap will only communicate owned grid cell data
  //   so can't lose ghost data not yet summed to owned cells
  // nvalues+1 includes atom count

  if (modeatom) {
    if (dimension == 2)
      grid2d->reverse_comm(Grid2d::FIX,this,0,nvalues+1,sizeof(double),
                           grid_buf1,grid_buf2,MPI_DOUBLE);
    else
      grid3d->reverse_comm(Grid3d::FIX,this,0,nvalues+1,sizeof(double),
                           grid_buf1,grid_buf2,MPI_DOUBLE);
  }

  // deallocate local comm buffers b/c new ones will be allocated

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  // make copy of ptrs to grid data which needs to persist

  if (dimension == 2) grid2d_previous = grid2d;
  else grid3d_previous = grid3d;

  nxlo_out_previous = nxlo_out;
  nylo_out_previous = nylo_out;
  nzlo_out_previous = nzlo_out;

  grid_sample_previous = clone_one_grid(grid_sample);
  grid_nfreq_previous = clone_one_grid(grid_nfreq);
  if (aveflag == RUNNING || aveflag == WINDOW)
    grid_running_previous = clone_one_grid(grid_running);
  if (aveflag == WINDOW) {
    grid_window_previous = new GridData*[nwindow];
    for (int i = 0; i < nwindow; i++)
      grid_window_previous[i] = clone_one_grid(grid_window[i]);
  }

  delete grid_sample;
  delete grid_nfreq;
  if (aveflag == RUNNING || aveflag == WINDOW) delete grid_running;
  if (aveflag == WINDOW) {
    for (int i = 0; i < nwindow; i++)
      delete grid_window[i];
    delete [] grid_window;
  }

  // allocate grid instance and grid data for new decomposition

  allocate_grid();

  grid_sample = allocate_one_grid();
  grid_nfreq = allocate_one_grid();
  if (aveflag == RUNNING || aveflag == WINDOW) grid_running = allocate_one_grid();
  if (aveflag == WINDOW) {
    grid_window = new GridData*[nwindow];
    for (int i = 0; i < nwindow; i++)
      grid_window[i] = allocate_one_grid();
  }

  // perform remap from previous decomp to new decomp
  // nper = # of remapped values per grid cell
  //   depends on atom vs grid mode, and running/window options

  int nremap_buf1,nremap_buf2;
  if (dimension == 2)
    grid2d->setup_remap(grid2d_previous,nremap_buf1,nremap_buf2);
  else
    grid3d->setup_remap(grid3d_previous,nremap_buf1,nremap_buf2);

  int n = 2;                                          // grid_sample & grid_nfreq
  if (aveflag == RUNNING || aveflag == WINDOW) n++;   // grid_running
  if (aveflag == WINDOW) n += nwindow;                // grid_window
  int nper = n*nvalues;
  if (modeatom) nper += n;

  double *remap_buf1 = nullptr;
  double *remap_buf2 = nullptr;
  if (nremap_buf1) memory->create(remap_buf1, nper*nremap_buf1, "ave/grid:remap_buf1");
  if (nremap_buf2) memory->create(remap_buf2, nper*nremap_buf2, "ave/grid:remap_buf2");

  if (dimension == 2)
    grid2d->remap(Grid2d::FIX,this,0,nper,sizeof(double),remap_buf1,remap_buf2,MPI_DOUBLE);
  else
    grid3d->remap(Grid3d::FIX,this,0,nper,sizeof(double),remap_buf1,remap_buf2,MPI_DOUBLE);

  memory->destroy(remap_buf1);
  memory->destroy(remap_buf2);

  // delete grid instance and grid data for previous decomposition

  if (dimension == 2) delete grid2d_previous;
  else delete grid3d_previous;

  deallocate_one_grid(grid_sample_previous,nxlo_out_previous,nylo_out_previous,nzlo_out_previous);
  deallocate_one_grid(grid_nfreq_previous,nxlo_out_previous,nylo_out_previous,nzlo_out_previous);
  if (aveflag == RUNNING || aveflag == WINDOW)
    deallocate_one_grid(grid_running_previous,nxlo_out_previous,nylo_out_previous,nzlo_out_previous);
  if (aveflag == WINDOW) {
    for (int i = 0; i < nwindow; i++)
      deallocate_one_grid(grid_window_previous[i],nxlo_out_previous,nylo_out_previous,nzlo_out_previous);
    delete [] grid_window_previous;
  }

  // set output data in case load balance fix comes after fix ave/grid

  output_grid(grid_nfreq);
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
      if (nvalues == 1) return grid_output->vec2d;
      else return grid_output->array2d;
    } else {
      if (nvalues == 1) return grid_output->vec3d;
      else return grid_output->array3d;
    }
  }
  if (index == 1) {
    if (dimension == 2) return grid_output->count2d;
    else return grid_output->count3d;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of local per-grid data
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = 0.0;

  bytes += size_grid(grid_sample);
  bytes += size_grid(grid_nfreq);
  if (aveflag == RUNNING || aveflag == WINDOW)
    bytes += size_grid(grid_running);
  if (aveflag == WINDOW)
    bytes += nwindow * size_grid(grid_window[0]);

  if (modeatom) {
    bytes += sizeof(int) * maxatom * dimension;   // bin array
    bytes += sizeof(int) * maxatom;               // skip vector
    bytes += sizeof(double) * maxvar;             // vresult for per-atom variable
  }

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
