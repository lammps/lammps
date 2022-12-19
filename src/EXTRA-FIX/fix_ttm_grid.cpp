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
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm_grid.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "grid3d.h"
#include "memory.h"
#include "neighbor.h"
#include "random_mars.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int MAXLINE = 256;
static constexpr int CHUNK = 1024;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly

static constexpr int OFFSET = 16384;

/* ---------------------------------------------------------------------- */

FixTTMGrid::FixTTMGrid(LAMMPS *lmp, int narg, char **arg) :
  FixTTM(lmp, narg, arg)
{
  pergrid_flag = 1;
  pergrid_freq = 1;
  restart_file = 1;

  if (outfile) error->all(FLERR,"Fix ttm/grid does not support outfile option - "
                          "use dump grid command or restart files instead");

  skin_original = neighbor->skin;
}

/* ---------------------------------------------------------------------- */

FixTTMGrid::~FixTTMGrid()
{
  FixTTMGrid::deallocate_grid();
  deallocate_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::post_constructor()
{
  // allocate global grid on each proc
  // needs to be done in post_contructor() beccause is virtual method

  allocate_grid();

  // initialize electron temperatures on grid

  int ix,iy,iz;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix <= nxhi_out; ix++)
        T_electron[iz][iy][ix] = tinit;

  // zero net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  // set initial electron temperatures from user input file
  // communicate new T_electron values to ghost grid points

  if (infile) {
    read_electron_temperatures(infile);
    grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::init()
{
  FixTTM::init();

  if (neighbor->skin > skin_original)
    error->all(FLERR,"Cannot extend neighbor skin after fix ttm/grid defined");
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::post_force(int /*vflag*/)
{
  int ix,iy,iz;
  double gamma1,gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;

  // apply damping and thermostat to all atoms in fix group

  int flag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + OFFSET) - OFFSET;

      // flag if ix,iy,iz is not within my ghost cell range

      if (ix < nxlo_out || ix > nxhi_out ||
          iy < nylo_out || iy > nyhi_out ||
          iz < nzlo_out || iz > nzhi_out) {
        flag = 1;
        continue;
      }

      if (T_electron[iz][iy][ix] < 0)
        error->one(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[iz][iy][ix]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }

  if (flag) error->one(FLERR,"Out of range fix ttm/grid atoms");
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::end_of_step()
{
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;
  double volgrid = 1.0 / (dxinv*dyinv*dzinv);

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + OFFSET) - OFFSET;

      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  grid->reverse_comm(Grid3d::FIX,this,0,1,sizeof(double),
                     grid_buf1,grid_buf2,MPI_DOUBLE);

  // clang-format off

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  double stability_criterion = 1.0 -
    2.0*inner_dt/(electronic_specific_heat*electronic_density) *
    (electronic_thermal_conductivity*(dxinv*dxinv + dyinv*dyinv + dzinv*dzinv));

  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) /
      (electronic_thermal_conductivity*(dxinv*dxinv + dyinv*dyinv + dzinv*dzinv));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/grid");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    memcpy(&T_electron_old[nzlo_out][nylo_out][nxlo_out],
           &T_electron[nzlo_out][nylo_out][nxlo_out],ngridout*sizeof(double));

    // compute new electron T profile

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          T_electron[iz][iy][ix] =
            T_electron_old[iz][iy][ix] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *

             ((T_electron_old[iz][iy][ix-1] + T_electron_old[iz][iy][ix+1] -
               2.0*T_electron_old[iz][iy][ix])*dxinv*dxinv +
              (T_electron_old[iz][iy-1][ix] + T_electron_old[iz][iy+1][ix] -
               2.0*T_electron_old[iz][iy][ix])*dyinv*dyinv +
              (T_electron_old[iz-1][iy][ix] + T_electron_old[iz+1][iy][ix] -
               2.0*T_electron_old[iz][iy][ix])*dzinv*dzinv) -

             net_energy_transfer[iz][iy][ix]/volgrid);
        }

    // communicate new T_electron values to ghost grid points

    grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ----------------------------------------------------------------------
   read electron temperatures on grid from a user-specified file
------------------------------------------------------------------------- */

void FixTTMGrid::read_electron_temperatures(const std::string &filename)
{
  memory->create3d_offset(T_electron_read, nzlo_in, nzhi_in, nylo_in, nyhi_in, nxlo_in, nxhi_in,
                          "ttm/grid:T_electron_read");
  memset(&T_electron_read[nzlo_in][nylo_in][nxlo_in], 0, ngridown * sizeof(int));

  // proc 0 opens file

  FILE *fp = nullptr;
  if (comm->me == 0) {
    fp = utils::open_potential(filename, lmp, nullptr);
    if (!fp) error->one(FLERR, "Cannot open grid file: {}: {}", filename, utils::getsyserror());
  }

  // read the file
  // Grid3d::read_file() calls back to read_grid_lines() with chunks of lines

  grid->read_file(Grid3d::FIX,this,fp,CHUNK,MAXLINE);

  // close file

  if (comm->me == 0) fclose(fp);

  // check completeness of input data

  int flag = 0;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++)
    for (int iy = nylo_in; iy <= nyhi_in; iy++)
      for (int ix = nxlo_in; ix <= nxhi_in; ix++)
        if (T_electron_read[iz][iy][ix] == 0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall) error->all(FLERR, "Fix ttm/grid infile did not set all temperatures");

  memory->destroy3d_offset(T_electron_read, nzlo_in, nylo_in, nxlo_in);
}

/* ----------------------------------------------------------------------
   process a chunk of lines in buffer
   each proc stores values for grid points it owns
   called back to from Grid3d::read_file()
------------------------------------------------------------------------- */

int FixTTMGrid::unpack_read_grid(int /*nlines*/, char *buffer)
{
  // loop over chunk of lines of grid point values
  // skip comment lines
  // tokenize the line into ix,iy,iz grid index plus temperature value
  // if I own grid point, store the value

  int nread = 0;

  for (const auto &line : utils::split_lines(buffer)) {
    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() == 0) {
        ;    // ignore comment only or blank lines
      } else if (values.count() == 4) {
        ++nread;

        int ix = values.next_int() - 1;
        int iy = values.next_int() - 1;
        int iz = values.next_int() - 1;

        if (ix < 0 || ix >= nxgrid || iy < 0 || iy >= nygrid || iz < 0 || iz >= nzgrid)
          throw TokenizerException("Fix ttm/grid invalid grid index in input", "");

        if (ix >= nxlo_in && ix <= nxhi_in && iy >= nylo_in && iy <= nyhi_in && iz >= nzlo_in &&
            iz <= nzhi_in) {
          T_electron[iz][iy][ix] = values.next_double();
          T_electron_read[iz][iy][ix] = 1;
        }
      } else {
        throw TokenizerException("Incorrect format in fix ttm electron grid file", "");
      }
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }
  }

  return nread;
}

/* ----------------------------------------------------------------------
   pack state of Fix into one write, but not per-grid values
------------------------------------------------------------------------- */

void FixTTMGrid::write_restart(FILE *fp)
{
  double rlist[4];

  rlist[0] = nxgrid;
  rlist[1] = nygrid;
  rlist[2] = nzgrid;
  rlist[3] = seed;

  if (comm->me == 0) {
    int size = 4 * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),4,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTMGrid::restart(char *buf)
{
  auto rlist = (double *) buf;

  // check that restart grid size is same as current grid size

  int nxgrid_old = static_cast<int> (rlist[0]);
  int nygrid_old = static_cast<int> (rlist[1]);
  int nzgrid_old = static_cast<int> (rlist[2]);

  if (nxgrid_old != nxgrid || nygrid_old != nygrid || nzgrid_old != nzgrid)
    error->all(FLERR,"Must restart fix ttm with same grid size");

  // change RN seed from initial seed, to avoid same Langevin factors
  // just increment by 1, since for RanMars that is a new RN stream

  seed = static_cast<int> (rlist[3]) + 1;
  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   write electron temperatures on grid to file
   identical format to infile option, so info can be read in when restarting
   each proc contributes info for its portion of grid
------------------------------------------------------------------------- */

void FixTTMGrid::write_restart_file(const char *file)
{
  // proc 0 opens file and writes header

  if (comm->me == 0) {
    auto outfile = std::string(file) + ".ttm";
    fpout = fopen(outfile.c_str(),"w");
    if (fpout == nullptr)
      error->one(FLERR,"Cannot open fix ttm/grid restart file {}: {}",outfile,utils::getsyserror());

    fmt::print(fpout,"# DATE: {} UNITS: {} COMMENT: "
               "Electron temperature on {}x{}x{} grid at step {} - "
               "created by fix {}\n",
               utils::current_date(),update->unit_style,
               nxgrid,nygrid,nzgrid,update->ntimestep,style);
  }

  // write file
  // Grid3d::write_file() calls back to pack_write_file() and unpack_write_file()

  grid->write_file(Grid3d::FIX,this,0,1,sizeof(double), MPI_DOUBLE);

  // close file

  if (comm->me == 0) fclose(fpout);
}

/* ----------------------------------------------------------------------
   pack values from local grid into buf
------------------------------------------------------------------------- */

void FixTTMGrid::pack_write_grid(int /*which*/, void *vbuf)
{
  int ix, iy, iz;

  auto buf = (double *) vbuf;

  int m = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf[m++] = T_electron[iz][iy][ix];
}

/* ----------------------------------------------------------------------
   unpack values from buf and write them to restart file
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_write_grid(int /*which*/, void *vbuf, int *bounds)
{
  int ix, iy, iz;

  int xlo = bounds[0];
  int xhi = bounds[1];
  int ylo = bounds[2];
  int yhi = bounds[3];
  int zlo = bounds[4];
  int zhi = bounds[5];

  auto buf = (double *) vbuf;
  double value;

  int m = 0;
  for (iz = zlo; iz <= zhi; iz++)
    for (iy = ylo; iy <= yhi; iy++)
      for (ix = xlo; ix <= xhi; ix++) {
        value = buf[m++];
        fprintf(fpout, "%d %d %d %20.16g\n", ix+1, iy+1, iz+1, value);
      }
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   called by load balancer when proc subdomains are adjusted
------------------------------------------------------------------------- */

void FixTTMGrid::reset_grid()
{
  // check if new grid partitioning is different on any proc
  // if not, just return

  int tmp[12];
  double maxdist = 0.5 * neighbor->skin;
  Grid3d *gridnew = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
  gridnew->set_distance(maxdist);
  gridnew->set_stencil_grid(1,1);
  gridnew->setup_grid(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],
                      tmp[6],tmp[7],tmp[8],tmp[9],tmp[10],tmp[11]);

  if (grid->identical(gridnew)) {
    delete gridnew;
    return;
  } else delete gridnew;

  // delete grid data which doesn't need to persist from previous to new decomp

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);
  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);

  // make copy of ptrs to grid data which does need to persist

  grid_previous = grid;
  T_electron_previous = T_electron;
  nxlo_out_previous = nxlo_out;
  nylo_out_previous = nylo_out;
  nzlo_out_previous = nzlo_out;

  // allocate new per-grid data for new decomposition

  allocate_grid();

  // perform remap from previous decomp to new decomp

  int nremap_buf1,nremap_buf2;
  grid->setup_remap(grid_previous,nremap_buf1,nremap_buf2);

  double *remap_buf1,*remap_buf2;
  memory->create(remap_buf1, nremap_buf1, "ttm/grid:remap_buf1");
  memory->create(remap_buf2, nremap_buf2, "ttm/grid:remap_buf2");

  grid->remap(Grid3d::FIX,this,0,1,sizeof(double),remap_buf1,remap_buf2,MPI_DOUBLE);

  memory->destroy(remap_buf1);
  memory->destroy(remap_buf2);

  // delete grid data and grid for previous decomposition

  memory->destroy3d_offset(T_electron_previous,
                           nzlo_out_previous, nylo_out_previous,
                           nxlo_out_previous);
  delete grid_previous;

  // communicate temperatures to ghost cells on new grid

  grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                     grid_buf1,grid_buf2,MPI_DOUBLE);

  // zero new net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_forward_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_forward_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] += buf[i];
}

/* ----------------------------------------------------------------------
   pack old grid values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src =
    &T_electron_previous[nzlo_out_previous][nylo_out_previous][nxlo_out_previous];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::allocate_grid()
{
  double maxdist = 0.5 * neighbor->skin;

  grid = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
  grid->set_distance(maxdist);
  grid->set_stencil_grid(1,1);
  grid->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                   nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  ngridown = (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) *
    (nzhi_in - nzlo_in + 1);
  ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
    (nzhi_out - nzlo_out + 1);

  // setup grid communication and allocate grid data structs

  grid->setup_comm(ngrid_buf1, ngrid_buf2);

  memory->create(grid_buf1, ngrid_buf1, "ttm/grid:grid_buf1");
  memory->create(grid_buf2, ngrid_buf2, "ttm/grid:grid_buf2");

  memory->create3d_offset(T_electron_old, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/grid:T_electron_old");
  memory->create3d_offset(T_electron, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "ttm/grid:T_electron");
  memory->create3d_offset(net_energy_transfer, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/grid:net_energy_transfer");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::deallocate_grid()
{
  delete grid;
  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(T_electron, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);
}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
   return -1 if grid name not found
------------------------------------------------------------------------- */

int FixTTMGrid::get_grid_by_name(const std::string &name, int &dim)
{
  if (name == "grid") {
    dim = 3;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to Grid data struct for grid with index
   this class can store M named grids, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixTTMGrid::get_grid_by_index(int index)
{
  if (index == 0) return grid;

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

int FixTTMGrid::get_griddata_by_name(int igrid, const std::string &name, int &ncol)
{
  if ((igrid == 0) && (name == "data")) {
    ncol = 0;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to multidim data array associated with index
   this class can store G named data sets, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixTTMGrid::get_griddata_by_index(int index)
{
  if (index == 0) return T_electron;

  return nullptr;
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem
   or the net_energy transfer between the subsystems
------------------------------------------------------------------------- */

double FixTTMGrid::compute_vector(int n)
{
  int ix, iy, iz;

  if (outflag == 0) {
    double dx = domain->xprd / nxgrid;
    double dy = domain->yprd / nygrid;
    double dz = domain->zprd / nzgrid;
    double volgrid = dx * dy * dz;

    double e_energy_me = 0.0;
    double transfer_energy_me = 0.0;

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          e_energy_me +=
              T_electron[iz][iy][ix] * electronic_specific_heat * electronic_density * volgrid;
          transfer_energy_me += net_energy_transfer[iz][iy][ix] * update->dt;
        }

    MPI_Allreduce(&e_energy_me, &e_energy, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&transfer_energy_me, &transfer_energy, 1, MPI_DOUBLE, MPI_SUM, world);
    outflag = 1;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grid
------------------------------------------------------------------------- */

double FixTTMGrid::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) 3 * atom->nmax * sizeof(double);
  bytes += (double) 3 * ngridout * sizeof(double);
  return bytes;
}
