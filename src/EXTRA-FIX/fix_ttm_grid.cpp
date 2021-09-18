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
#include "gridcomm.h"
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
static constexpr int OFFSET = 16384;

/* ---------------------------------------------------------------------- */

FixTTMGrid::FixTTMGrid(LAMMPS *lmp, int narg, char **arg) :
  FixTTM(lmp, narg, arg)
{
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
    gc->forward_comm(GridComm::FIX,this,1,sizeof(double),0,gc_buf1,gc_buf2,MPI_DOUBLE);
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

  gc->reverse_comm(GridComm::FIX,this,1,sizeof(double),0,
                   gc_buf1,gc_buf2,MPI_DOUBLE);

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

    gc->forward_comm(GridComm::FIX,this,1,sizeof(double),0,gc_buf1,gc_buf2,MPI_DOUBLE);
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
  memset(&T_initial_set[nzlo_in][nylo_in][nxlo_in], 0, ngridmine * sizeof(int));

  // proc 0 opens file

  FILE *fp = nullptr;
  if (comm->me == 0) {
    fp = utils::open_potential(filename, lmp, nullptr);
    if (!fp) error->one(FLERR, "Cannot open grid file: {}: {}", filename, utils::getsyserror());
  }

  // read electron temperature values from file, one chunk at a time

  char *buffer = new char[CHUNK * MAXLINE];
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
          ; // ignore comment only lines
        } else if (values.count() == 4) {
          ++nread;

          int ix = values.next_int();
          int iy = values.next_int();
          int iz = values.next_int();

          if (ix < 0 || ix >= nxgrid || iy < 0 || iy >= nygrid || iz < 0 || iz >= nzgrid)
            throw parser_error("Fix ttm/grid invalid grid index in input");

          if (ix >= nxlo_in && ix <= nxhi_in && iy >= nylo_in && iy <= nyhi_in
              && iz >= nzlo_in && iz <= nzhi_in) {
            T_electron[iz][iy][ix] = values.next_double();
            T_initial_set[iz][iy][ix] = 1;
          }
        } else {
          throw parser_error("Incorrect format in fix ttm electron grid file");
        }
      } catch (std::exception &e) {
        error->one(FLERR,e.what());
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
   write out current electron temperatures to user-specified file
   only written by proc 0
------------------------------------------------------------------------- */

void FixTTMGrid::write_electron_temperatures(const std::string &filename)
{
  if (comm->me == 0) {
    FPout = fopen(filename.c_str(), "w");
    if (!FPout) error->one(FLERR, "Fix ttm/grid could not open output file");

    fmt::print(FPout,"# DATE: {} UNITS: {} COMMENT: Electron temperature "
               "{}x{}x{} grid at step {}. Created by fix {}\n", utils::current_date(),
               update->unit_style, nxgrid, nygrid, nzgrid, update->ntimestep, style);
  }

  gc->gather(GridComm::FIX, this, 1, sizeof(double), 1, nullptr, MPI_DOUBLE);

  if (comm->me == 0) fclose(FPout);
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_forward_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_forward_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_reverse_grid(int /*flag*/, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] += buf[i];
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::allocate_grid()
{
  // partition global grid across procs
  // n xyz lo/hi in = lower/upper bounds of global grid this proc owns
  // indices range from 0 to N-1 inclusive in each dim

  comm->partition_grid(nxgrid, nygrid, nzgrid, 0.0, nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in,
                       nzhi_in);

  // nlo,nhi = min/max index of global grid pt my owned atoms can be mapped to
  // finite difference stencil requires extra grid pt around my owned grid pts
  // max of these 2 quantities is the ghost cells needed in each dim
  // nlo_out,nhi_out = nlo_in,nhi_in + ghost cells

  double *boxlo = domain->boxlo;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  double dxinv = nxgrid / domain->xprd;
  double dyinv = nxgrid / domain->yprd;
  double dzinv = nxgrid / domain->zprd;

  int nlo, nhi;
  double cuthalf = 0.5 * neighbor->skin;

  nlo = static_cast<int>((sublo[0] - cuthalf - boxlo[0]) * dxinv + shift) - OFFSET;
  nhi = static_cast<int>((subhi[0] + cuthalf - boxlo[0]) * dxinv + shift) - OFFSET;
  nxlo_out = MIN(nlo, nxlo_in - 1);
  nxhi_out = MAX(nhi, nxhi_in + 1);

  nlo = static_cast<int>((sublo[1] - cuthalf - boxlo[1]) * dyinv + shift) - OFFSET;
  nhi = static_cast<int>((subhi[1] + cuthalf - boxlo[1]) * dyinv + shift) - OFFSET;
  nylo_out = MIN(nlo, nylo_in - 1);
  nyhi_out = MAX(nhi, nyhi_in + 1);

  nlo = static_cast<int>((sublo[2] - cuthalf - boxlo[2]) * dzinv + shift) - OFFSET;
  nhi = static_cast<int>((subhi[2] + cuthalf - boxlo[2]) * dzinv + shift) - OFFSET;
  nzlo_out = MIN(nlo, nzlo_in - 1);
  nzhi_out = MAX(nhi, nzhi_in + 1);

  bigint totalmine;
  totalmine =
      (bigint) (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) * (nzhi_out - nzlo_out + 1);
  if (totalmine > MAXSMALLINT) error->one(FLERR, "Too many owned+ghost grid points in fix ttm");
  ngridout = totalmine;

  totalmine = (bigint) (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) * (nzhi_in - nzlo_in + 1);
  ngridmine = totalmine;

  gc = new GridComm(lmp, world, nxgrid, nygrid, nzgrid, nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in,
                    nzhi_in, nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  gc->setup(ngc_buf1, ngc_buf2);

  memory->create(gc_buf1, ngc_buf1, "ttm/grid:gc_buf1");
  memory->create(gc_buf2, ngc_buf2, "ttm/grid:gc_buf2");

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
  delete gc;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);

  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(T_electron, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTMGrid::write_restart(FILE *fp)
{
  int rsize = nxgrid * nygrid * nzgrid + 4;
  double *rlist;
  memory->create(rlist, rsize, "ttm/grid:rlist");

  int n = 0;
  rlist[n++] = nxgrid;
  rlist[n++] = nygrid;
  rlist[n++] = nzgrid;
  rlist[n++] = seed;

  // gather rest of rlist on proc 0 as global grid values

  gc->gather(GridComm::FIX, this, 1, sizeof(double), 0, &rlist[4], MPI_DOUBLE);

  if (comm->me == 0) {
    int size = rsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(rlist, sizeof(double), rsize, fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTMGrid::restart(char *buf)
{
  int ix, iy, iz;

  int n = 0;
  double *rlist = (double *) buf;

  // check that restart grid size is same as current grid size

  int nxgrid_old = static_cast<int>(rlist[n++]);
  int nygrid_old = static_cast<int>(rlist[n++]);
  int nzgrid_old = static_cast<int>(rlist[n++]);

  if (nxgrid_old != nxgrid || nygrid_old != nygrid || nzgrid_old != nzgrid)
    error->all(FLERR, "Must restart fix ttm/grid with same grid size");

  // change RN seed from initial seed, to avoid same Langevin factors
  // just increment by 1, since for RanMars that is a new RN stream

  seed = static_cast<int>(rlist[n++]) + 1;
  delete random;
  random = new RanMars(lmp, seed + comm->me);

  // extract this proc's local grid values from global grid in rlist

  int iglobal;

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        iglobal = nygrid * nxgrid * iz + nxgrid * iy + ix;
        T_electron[iz][iy][ix] = rlist[n + iglobal];
      }

  // communicate new T_electron values to ghost grid points

  gc->forward_comm(GridComm::FIX, this, 1, sizeof(double), 0, gc_buf1, gc_buf2, MPI_DOUBLE);
}

/* ----------------------------------------------------------------------
   pack values from local grid into buf
   used by which = 0 and 1
------------------------------------------------------------------------- */

void FixTTMGrid::pack_gather_grid(int /*which*/, void *vbuf)
{
  int ix, iy, iz;

  double *buf = (double *) vbuf;

  int m = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) buf[m++] = T_electron[iz][iy][ix];
}

/* ----------------------------------------------------------------------
   which = 0: unpack values from buf into global gbuf based on their indices
   which = 1: print values from buf to FPout file
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_gather_grid(int which, void *vbuf, void *vgbuf, int xlo, int xhi, int ylo,
                                    int yhi, int zlo, int zhi)
{
  int ix, iy, iz;

  double *buf = (double *) vbuf;
  double *gbuf = (double *) vgbuf;

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
