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
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt
// use SHIFT = 0.0 for now since it allows fix ave/chunk
//   to spatially average consistent with the TTM grid

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.5;

/* ---------------------------------------------------------------------- */

FixTTMGrid::FixTTMGrid(LAMMPS *lmp, int narg, char **arg) :
  FixTTM(lmp, narg, arg)
{
  pergrid_flag = 1;
  pergrid_freq = 1;

  if (outfile) error->all(FLERR,"Fix ttm/grid does not support outfile option - "
                          "use dump grid command instead");

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
    grid->forward_comm(Grid3d::FIX,this,1,sizeof(double),0,
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::init()
{
  FixTTM::init();

  if (neighbor->skin > skin_original)
    error->all(FLERR,"Cannot extend neighbor skin after fix ttm/griddefined");
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
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;

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
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;

      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  grid->reverse_comm(Grid3d::FIX,this,1,sizeof(double),0,
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

    grid->forward_comm(Grid3d::FIX,this,1,sizeof(double),0,
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }

  // clang-format on

  // output of grid temperatures to file

  if (outfile && (update->ntimestep % outevery == 0))
    write_electron_temperatures(fmt::format("{}.{}", outfile, update->ntimestep));
}

/* ----------------------------------------------------------------------
   read electron temperatures on grid from a user-specified file
   proc 0 reads one chunk at a time, broadcasts to other procs
   each proc stores values for grid points it owns
------------------------------------------------------------------------- */

void FixTTMGrid::read_electron_temperatures(const std::string &filename)
{
  int ***T_initial_set;
  memory->create3d_offset(T_initial_set, nzlo_in, nzhi_in, nylo_in, nyhi_in, nxlo_in, nxhi_in,
                          "ttm/grid:T_initial_set");
  memset(&T_initial_set[nzlo_in][nylo_in][nxlo_in], 0, ngridown * sizeof(int));

  // proc 0 opens file

  FILE *fp = nullptr;
  if (comm->me == 0) {
    fp = utils::open_potential(filename, lmp, nullptr);
    if (!fp) error->one(FLERR, "Cannot open grid file: {}: {}", filename, utils::getsyserror());
  }

  // read electron temperature values from file, one chunk at a time

  auto buffer = new char[CHUNK * MAXLINE];
  bigint ntotal = (bigint) nxgrid * nygrid * nzgrid;
  bigint nread = 0;

  while (nread < ntotal) {
    int nchunk = MIN(ntotal - nread, CHUNK);
    int eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, comm->me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");

    // loop over lines of grid point values
    // tokenize the line into ix,iy,iz grid index plus temperature value
    // if I own grid point, store the value

    for (const auto &line : utils::split_lines(buffer)) {
      try {
        ValueTokenizer values(utils::trim_comment(line));
        if (values.count() == 0) {
          ;    // ignore comment only lines
        } else if (values.count() == 4) {
          ++nread;

          int ix = values.next_int();
          int iy = values.next_int();
          int iz = values.next_int();

          if (ix < 0 || ix >= nxgrid || iy < 0 || iy >= nygrid || iz < 0 || iz >= nzgrid)
            throw TokenizerException("Fix ttm/grid invalid grid index in input", "");

          if (ix >= nxlo_in && ix <= nxhi_in && iy >= nylo_in && iy <= nyhi_in && iz >= nzlo_in &&
              iz <= nzhi_in) {
            T_electron[iz][iy][ix] = values.next_double();
            T_initial_set[iz][iy][ix] = 1;
          }
        } else {
          throw TokenizerException("Incorrect format in fix ttm electron grid file", "");
        }
      } catch (std::exception &e) {
        error->one(FLERR, e.what());
      }
    }
  }

  // close file

  if (comm->me == 0) fclose(fp);

  // clean up

  delete[] buffer;

  // check completeness of input data

  int flag = 0;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++)
    for (int iy = nylo_in; iy <= nyhi_in; iy++)
      for (int ix = nxlo_in; ix <= nxhi_in; ix++)
        if (T_initial_set[iz][iy][ix] == 0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall) error->all(FLERR, "Fix ttm/grid infile did not set all temperatures");

  memory->destroy3d_offset(T_initial_set, nzlo_in, nylo_in, nxlo_in);
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   called by load balancer when proc subdomains are adjusted
------------------------------------------------------------------------- */

void FixTTMGrid::reset_grid()
{
  // delete grid data which doesn't need to persist from previous to new decomp

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);
  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);

  // make copy of ptrs to grid data which does need to persist

  grid_previous = grid;
  T_electron_previous = T_electron;

  // allocate new per-grid data for new decomposition

  allocate_grid();

  // perform remap from previous decomp to new decomp

  int nremap_buf1,nremap_buf2;
  grid->remap_setup(grid_previous,nremap_buf1,nremap_buf2);

  double *remap_buf1,*remap_buf2;
  memory->create(remap_buf1, nremap_buf1, "ttm/grid:remap_buf1");
  memory->create(remap_buf2, nremap_buf2, "ttm/grid:remap_buf2");

  grid->remap(Grid3d::FIX,this,1,sizeof(double),remap_buf1,remap_buf2,MPI_DOUBLE);

  memory->destroy(remap_buf1);
  memory->destroy(remap_buf2);

  // delete grid data and grid for previous decomposition
  // NOTE: need to set offsets

  int nxlo_out_prev,nylo_out_prev,nzlo_out_prev;
  memory->destroy3d_offset(T_electron_previous,
                           nzlo_out_prev, nylo_out_prev, nxlo_out_prev);
  delete grid_previous;
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_forward_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_forward_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] += buf[i];
}

/* ----------------------------------------------------------------------
   pack old grid  values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_remap_grid(void *vbuf, int nlist, int *list)
{
  auto buf = (double *) vbuf;
  double *src = &T_electron_previous[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_remap_grid(void *vbuf, int nlist, int *list)
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

  grid = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid, maxdist, 1, SHIFT,
                    nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                    nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  ngridown = (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) *
    (nzhi_in - nzlo_in + 1);
  ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
    (nzhi_out - nzlo_out + 1);

  // setup grid communication and allocate grid data structs

  grid->setup(ngrid_buf1, ngrid_buf2);

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
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTMGrid::write_restart(FILE *fp)
{
  error->all(FLERR,"Fix ttm/grid command does not yet support "
             "writing a distributed grid to a restart file");
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTMGrid::restart(char *buf)
{
  error->all(FLERR,"Fix ttm/grid command does not yet support "
             "reading a distributed grid from a restart file");
}

/* ----------------------------------------------------------------------
   pack values from local grid into buf
   used by which = 0 and 1
   NOTE: remove this function when ready to release
------------------------------------------------------------------------- */

void FixTTMGrid::pack_gather_grid(int /*which*/, void *vbuf)
{
  int ix, iy, iz;

  auto buf = (double *) vbuf;

  int m = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) buf[m++] = T_electron[iz][iy][ix];
}

/* ----------------------------------------------------------------------
   which = 0: unpack values from buf into global gbuf based on their indices
   which = 1: print values from buf to FPout file
   NOTE: remove this function when ready to release
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_gather_grid(int which, void *vbuf, void *vgbuf, int xlo, int xhi, int ylo,
                                    int yhi, int zlo, int zhi)
{
  int ix, iy, iz;

  auto buf = (double *) vbuf;
  auto gbuf = (double *) vgbuf;

  if (which == 0) {
    int iglobal;
    int ilocal = 0;

    for (iz = zlo; iz <= zhi; iz++)
      for (iy = ylo; iy <= yhi; iy++)
        for (ix = xlo; ix <= xhi; ix++) {
          iglobal = nygrid * nxgrid * iz + nxgrid * iy + ix;
          gbuf[iglobal] = buf[ilocal++];
        }

  } else if (which == 1) {
    int ilocal = 0;
    double value;

    for (iz = zlo; iz <= zhi; iz++)
      for (iy = ylo; iy <= yhi; iy++)
        for (ix = xlo; ix <= xhi; ix++) {
          value = buf[ilocal++];
          fprintf(FPout, "%d %d %d %20.16g\n", ix, iy, iz, value);
        }
  }
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
