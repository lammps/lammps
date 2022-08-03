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
#include "grid2d.h"
#include "grid3d.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SAMPLE,ALL};
enum{NOSCALE,ATOM};
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
  vec2d(nullptr), array2d(nullptr), vec3d(nullptr), array3d(nullptr),
  count2d(nullptr), count3d(nullptr)
{
  if (narg < 10) error->all(FLERR,"Illegal fix ave/grid command");

  pergrid_flag = 1;
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  pergrid_freq = utils::inumeric(FLERR,arg[5],false,lmp);
  time_depend = 1;

  // NOTE: allow Dxyz as well

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
      argindex[nvalues++] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues++] = 1;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues++] = 2;
      modeatom = 1;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 1;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 2;
      modeatom = 1;

    } else if (strcmp(arg[iarg],"density/number") == 0) {
      which[nvalues] = ArgInfo::DENSITY_NUMBER;
      argindex[nvalues++] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"density/mass") == 0) {
      which[nvalues] = ArgInfo::DENSITY_MASS;
      argindex[nvalues++] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      which[nvalues] = ArgInfo::MASS;
      argindex[nvalues++] = 0;
      modeatom = 1;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      which[nvalues] = ArgInfo::TEMPERATURE;
      argindex[nvalues++] = 0;
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

      nvalues++;

    // unrecognized arg (option)

    } else break;

    iarg++;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/grid command");

  if (modeatom && modegrid) 
    error->all(FLERR,"Fix ave/grid cannot operate on per-atom and "
               "per-grid values");

  // optional args

  normflag = ALL;
  scaleflag = ATOM;
  ave = ONE;
  nwindow = 0;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"all") == 0) {
        normflag = ALL;
        scaleflag = ATOM;
      } else if (strcmp(arg[iarg+1],"sample") == 0) {
        normflag = SAMPLE;
        scaleflag = ATOM;
      } else if (strcmp(arg[iarg+1],"none") == 0) {
        normflag = SAMPLE;
        scaleflag = NOSCALE;
      } else error->all(FLERR,"Illegal fix ave/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/grid command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/grid command");
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/grid command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;

    } else error->all(FLERR,"Illegal fix ave/grid command");
  }

  // if wildcard expansion occurred, free earg memory from exapnd_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
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

        char *idcompute,*gname,*dname;
        utils::grid_parse(FLERR,ids[i],idcompute,gname,dname,error);

        Compute *icompute = modify->get_compute_by_id(idcompute);
        if (!icompute) 
          error->all(FLERR,"Could not find fix ave/grid compute ID: {}",
                     idcompute);
        if (icompute->pergrid_flag == 0)
          error->all(FLERR,
                     "Fix ave/grid compute {} does not compute per-grid info",
                     idcompute);

        int dim;
        int igrid = icompute->get_grid_by_name(gname,dim);
        if (igrid < 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} does not recognize grid name {}",
                     idcompute,gname);
      
        int ncol;
        int idata = icompute->get_griddata_by_name(igrid,dname,ncol);
        if (idata < 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} does not recognize data name {}",
                     idcompute,dname);

        if (argindex[i] == 0 && ncol)
          error->all(FLERR,
                     "Fix ave/grid compute {} data {} is not per-grid vector",
                     idcompute,dname);
        if (argindex[i] && ncol == 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} data {} is not per-grid array",
                     idcompute,dname);
        if (argindex[i] && argindex[i] > ncol)
          error->all(FLERR,
                     "Fix ave/grid compute {} array {} is accessed out-of-range",
                     idcompute,dname);
      
	value2grid[iarg] = igrid;
	value2data[iarg] = idata;

        delete [] idcompute;
        delete [] gname;
        delete [] dname;

      } else if (which[i] == ArgInfo::FIX) {

        char *idfix,*gname,*dname;
        utils::grid_parse(FLERR,ids[i],idfix,gname,dname,error);

        Fix *ifix = modify->get_fix_by_id(idfix);
        if (!ifix) error->all(FLERR,"Could not find fix ave/grid fix ID: {}",
                              idfix);
        if (ifix->pergrid_flag == 0)
          error->all(FLERR,"Fix ave/grid fix {} does not compute per-grid info",
                     idfix);
        if (nevery % ifix->pergrid_freq)
          error->all(FLERR,
                     "Fix for fix grid/atom not computed at compatible time");

        int dim;
        int igrid = ifix->get_grid_by_name(gname,dim);
        if (igrid < 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} does not recognize grid name {}",
                     idfix,gname);
      
        int ncol;
        int idata = ifix->get_griddata_by_name(igrid,dname,ncol);
        if (idata < 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} does not recognize data name {}",
                     idfix,dname);

        if (argindex[i] == 0 && ncol)
          error->all(FLERR,
                     "Fix ave/grid compute {} data {} is not per-grid vector",
                     idfix,dname);
        if (argindex[i] && ncol == 0) 
          error->all(FLERR,
                     "Fix ave/grid compute {} data {} is not per-grid array",
                     idfix,dname);
        if (argindex[i] && argindex[i] > ncol)
          error->all(FLERR,
                     "Fix ave/grid compute {} array {} is accessed out-of-range",
                     idfix,dname);
      
	value2grid[iarg] = igrid;
	value2data[iarg] = idata;

        delete [] idfix;
        delete [] gname;
        delete [] dname;
      }
    }
  }

  // instantiate the Grid class and allocate per-grid memory
  // NOTE: need to extend ghost grid for ATOM mode ?

  if (modeatom) shift = OFFSET + SHIFT;
  else shift = 0.0;

  if (dimension == 2) {
    if (modeatom) 
    grid2d = new Grid2d(lmp, world, nxgrid, nygrid, 0, 0.0, shift,
                        nxlo_in, nxhi_in, nylo_in, nyhi_in,
                        nxlo_out, nxhi_out, nylo_out, nyhi_out);

    grid2d->setup(ngrid_buf1, ngrid_buf2);
    memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1);

    if (nvalues == 1)
      memory->create2d_offset(vec2d, nylo_out, nyhi_out, nxlo_out, nxhi_out, 
                              "ave/grid:vec2d");
    else
      memory->create3d_offset_last(array2d, nylo_out, nyhi_out, nxlo_out,
                                   nxhi_out, nvalues, "ave/grid:array2d");

    if (modeatom) 
      memory->create2d_offset(count2d, nylo_out, nyhi_out, nxlo_out, nxhi_out, 
                              "ave/grid:count2d");
    
  } else {
    grid3d = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid, 0, 0.0, shift,
                        nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in, 
                        nxlo_out, nxhi_out, nylo_out, nyhi_out, 
                        nzlo_out, nzhi_out);

    grid3d->setup(ngrid_buf1, ngrid_buf2);
    memory->create(grid_buf1, ngrid_buf1, "ave/grid:grid_buf1");
    memory->create(grid_buf2, ngrid_buf2, "ave/grid:grid_buf2");

    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) * 
      (nzhi_out - nzlo_out + 1);

    if (nvalues == 1)
      memory->create3d_offset(vec3d, nzlo_out, nzhi_out, nylo_out, 
                              nyhi_out, nxlo_out, nxhi_out, 
                              "ave/grid:vec3d");
    else
      memory->create4d_offset_last(array3d, nzlo_out, nzhi_out, nylo_out, 
                                   nyhi_out, nxlo_out, nxhi_out, nvalues, 
                                   "ave/grid:array3d");

    if (modeatom)
      memory->create3d_offset(count3d, nzlo_out, nzhi_out, nylo_out, 
                              nyhi_out, nxlo_out, nxhi_out, 
                              "ave/grid:vec3d");
  }

  // zero the grid since dump may access it on timestep 0

  zero_grid();

  // bin indices for ATOM mode
  // vresult for per-atom variable evaluation

  maxatom = 0;
  bin = nullptr;

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
  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;
  delete [] value2grid;
  delete [] value2data;

  delete grid2d;
  delete grid3d;
  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  memory->destroy2d_offset(vec2d,nylo_out,nxlo_out);
  memory->destroy2d_offset(array2d,nylo_out,nxlo_out);
  memory->destroy2d_offset(count2d,nylo_out,nxlo_out);
  memory->destroy3d_offset(vec3d,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(count3d,nzlo_out,nylo_out,nxlo_out);

  memory->destroy(bin);
  memory->destroy(vresult);
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
  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/atom does not exist");
      value2index[m] = icompute;

    } else if (which[m] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/atom does not exist");
      value2index[m] = ifix;

    } else if (which[m] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/atom does not exist");
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

  // zero owned and ghost grid points if first step
  // zero atom count per bin for ATOM mode

  if (irepeat == 0) zero_grid();

  if (modeatom) {
    if (dimension == 2) {
      if (ngridout) memcpy(&count2d[0][0],0,ngridout*sizeof(int));
    } else {
      if (ngridout) memcpy(&count3d[0][0][0],0,ngridout*sizeof(int));
    }
  }

  // ATOM mode
  // accumulate per-atom attributes,computes,fixes,variables to local grid
  
  // set local per-grid values for either ATOM or GRID mode
  // per-atom compute/fix/variable may invoke computes so wrap with clear/add

  if (modeatom) {
    modify->clearstep_compute();
    atom2grid();
  } else {
    grid2grid();
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    if (modeatom) modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+peratom_freq - ((bigint) nrepeat-1)*nevery;
  if (modeatom) modify->addstep_compute(nvalid);

  // ghost to owned grid communication for atom mode
  // NOTE: still need to implement pack/unpack methods

  if (modeatom) {
    if (dimension == 2) 
      grid2d->reverse_comm(Grid2d::FIX,this,1,sizeof(double),0,
                           grid_buf1,grid_buf2,MPI_DOUBLE);
    else 
      grid3d->reverse_comm(Grid3d::FIX,this,1,sizeof(double),0,
                           grid_buf1,grid_buf2,MPI_DOUBLE);
  }

  // just return if this proc owns no grid points
  
  if (ngridout == 0) return;
  
  // average the final result for the Nfreq timestep
  // just loop over owned grid points
  
  double invrepeat = 1.0/nrepeat;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
	for (ix = nxlo_in; ix <= nxhi_in; ix++)
	  vec2d[iy][ix] *= invrepeat;
    } else {
      for (iy = nylo_in; iy <= nyhi_in; iy++)
	for (ix = nxlo_in; ix <= nxhi_in; ix++)
	  for (m = 0; m <= nvalues; m++)
	    array2d[iy][ix][m] *= invrepeat;
    }
  } else {
    if (nvalues == 1) {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
	for (iy = nylo_in; iy <= nyhi_in; iy++)
	  for (ix = nxlo_in; ix <= nxhi_in; ix++)
	    vec3d[iz][iy][ix] *= invrepeat;
    } else {
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
	for (iy = nylo_in; iy <= nyhi_in; iy++)
	  for (ix = nxlo_in; ix <= nxhi_in; ix++)
	    for (m = 0; m <= nvalues; m++)
	    array3d[iz][iy][ix][m] *= invrepeat;
    }
  }
}

/* ----------------------------------------------------------------------
   sum per-atom contributions to owned+ghost grid cells
   sets one of vec2d,array2d,vec3d,array3d
   also set count2d or count3d for atom count per bin
------------------------------------------------------------------------- */

void FixAveGrid::atom2grid()
{
  int i,j,k,m,n,ix,iy,iz;

  // bin[i][dim] = indices of bin each atom is in
  // not set if group mask does not match
  // also count atoms contributing to each bin

  // NOTE: error check if any atom out of grid bounds?
  
  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;
  
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxatom) {
    memory->destroy(bin);
    maxatom = atom->nmax;
    memory->create(bin,maxatom,dimension,"ave/grid:bin");
  }

  if (dimension == 2) {
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      count2d[iy][ix]++;
      bin[i][0] = iy;
      bin[i][1] = ix;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;
      count3d[iz][iy][ix]++;
      bin[i][0] = iz;
      bin[i][1] = iy;
      bin[i][2] = ix;
    }
  }

  // loop over user-specified values

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];
      
    if (which[m] == ArgInfo::X || which[m] == ArgInfo::V || 
        which[m] == ArgInfo::F) {
      double **attribute;
      if (which[m] == ArgInfo::X) attribute = atom->x;
      else if (which[m] == ArgInfo::V) attribute = atom->v;
      else if (which[m] == ArgInfo::F) attribute = atom->f;
      
      if (dimension == 2) {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit)
              vec2d[bin[i][0]][bin[i][1]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit)
              array2d[bin[i][0]][bin[i][1]][m] += attribute[i][j];
          }
      } else {
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit)
              vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += attribute[i][j];
          }
        } else
          for (i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit)
              array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += attribute[i][j];
          }
      }
	
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
              if (mask[i] & groupbit)
                vec2d[bin[i][0]][bin[i][1]] += ovector[i];
            }
          } else {
            int jm1 = j = 1;
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                vec2d[bin[i][0]][bin[i][1]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                array2d[bin[i][0]][bin[i][1]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                array2d[bin[i][0]][bin[i][1]][m] += oarray[i][jm1];
            }
          }
        }

      } else {
        if (nvalues == 1) {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += ovector[i];
            }
          } else { 
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                vec3d[bin[i][0]][bin[i][1]][bin[i][2]] += oarray[i][jm1];
            }
          }
        } else {
          if (j == 0) {
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
                array3d[bin[i][0]][bin[i][1]][bin[i][2]][m] += ovector[i];
            }
          } else {
            int jm1 = j - 1;
            for (i = 0; i < nlocal; i++) {
              if (mask[i] & groupbit)
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
        if (j == 0) 
          ovec3d = (double ***) fix->get_griddata_by_index(idata);
        else 
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
   zero grid values incluing ghost cells
------------------------------------------------------------------------- */

void FixAveGrid::zero_grid()
{
  if (dimension == 2) {
    if (nvalues == 1) {
      if (ngridout) memcpy(&vec2d[0][0],0,ngridout*sizeof(double));
    } else {
      if (ngridout) memcpy(&array2d[0][0][0],0,ngridout*nvalues*sizeof(double));
    }
  } else {
    if (nvalues == 1) {
      if (ngridout) memcpy(&vec3d[0][0][0],0,ngridout*sizeof(double));
    } else {
      if (ngridout) memcpy(&array3d[0][0][0][0],0,ngridout*nvalues*sizeof(double));
    }
  }
}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
   return -1 if grid name not found
------------------------------------------------------------------------- */

int FixAveGrid::get_grid_by_name(char *name, int &dim)
{
  if (strcmp(name,"grid") == 0) {
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

int FixAveGrid::get_griddata_by_name(int igrid, char *name, int &ncol)
{
  if (igrid == 0 && strcmp(name,"data") == 0) {
    if (nvalues == 1) ncol = 0;
    else ncol = nvalues;
    return 0;
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
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   NOTE: add more memory tallying
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = (double) ngridout * nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveGrid::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= ((bigint)nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}
