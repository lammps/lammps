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
   Contributing authors: Original fix ttm
                         Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)

                         ttm/thermal
                         Bradly Baer (Vanderbilt University)
                         D. Greg Walker (Vanderbilt University)

------------------------------------------------------------------------- */

// clang-format on

#include "fix_ttm_thermal.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "potential_file_reader.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt
// use SHIFT = 0.0 for now since it allows fix ave/chunk
//   to spatially average consistent with the TTM grid

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.0;

/* ---------------------------------------------------------------------- */

FixTTMThermal::FixTTMThermal(LAMMPS *lmp, int narg, char **arg) :
    FixTTM(lmp, narg, arg), gamma_p_grid(nullptr), inductive_response_grid(nullptr),
    c_e_grid(nullptr), k_e_grid(nullptr)

{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix ttm/thermal", error);

  seed = utils::inumeric(FLERR, arg[3], false, lmp);
  e_property_file = arg[4];
  nxgrid = utils::inumeric(FLERR, arg[5], false, lmp);
  nygrid = utils::inumeric(FLERR, arg[6], false, lmp);
  nzgrid = utils::inumeric(FLERR, arg[7], false, lmp);

  inductive_power = 0.0;
  tinit = 0.0;
  gamma_p = 0.0;    // to avoid uninitialized data access when calling FixTTM::init()

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ttm/thermal set", error);
      tinit = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (tinit <= 0.0)
        error->all(FLERR, iarg + 1, "Fix ttm/thermal initial temperature must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "source") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ttm/thermal source", error);
      inductive_power = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "infile") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ttm/thermal infile", error);
      infile = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "outfile") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix ttm/thermal outfile", error);
      outevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      outfile = arg[iarg + 2];
      iarg += 3;
    } else
      error->all(FLERR, iarg, "Unknown fix ttm/thermal keyword {}", arg[iarg]);
  }

  // error check

  if (seed <= 0) error->all(FLERR, 3, "Invalid random number seed in fix ttm/thermal command");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR, Error::NOPOINTER, "Fix ttm grid sizes must be > 0");

  // grid OFFSET to perform
  // SHIFT to map atom to nearest or lower-left grid point

  shift = OFFSET + SHIFT;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp, seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes + 1];
  gfactor2 = new double[atom->ntypes + 1];

  // check for allowed maximum number of total grid points

  bigint totalgrid = (bigint) nxgrid * nygrid * nzgrid;
  if (totalgrid > MAXSMALLINT)
    error->all(FLERR, Error::NOPOINTER, "Too many grid points in fix ttm/thermal");
  ngridtotal = totalgrid;

  // allocate per-atom flangevin and zero it

  FixTTM::grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }

  // set 2 callbacks

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);
}

/* ---------------------------------------------------------------------- */

FixTTMThermal::~FixTTMThermal()
{
  FixTTMThermal::deallocate_grid();
}

/* ---------------------------------------------------------------------- */
static inline double safe_effective_kappa(double a, double b)
{
  if (a == 0 || b == 0) return 0;
  return 2.0 * a * b / (a + b);
}
/* ---------------------------------------------------------------------- */

void FixTTMThermal::post_constructor()
{
  FixTTM::post_constructor();

  // set electron grid properties from file
  read_electron_properties(e_property_file);
}

/* ---------------------------------------------------------------------- */

void FixTTMThermal::init()
{
  FixTTM::init();

  // to allow this, would have to reset grid bounds dynamically
  // for RCB balancing would have to reassign grid pts to procs
  //   and create a new GridComm, and pass old GC data to new GC

  if (domain->box_change)
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot use fix ttm/thermal with changing box shape, size, or sub-domains");
}

// clang-format off

/* ---------------------------------------------------------------------- */

void FixTTMThermal::post_force(int /*vflag*/)
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

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;

      if (T_electron[iz][iy][ix] < 0)
        error->one(FLERR, Error::NOLASTLINE, "Electronic temperature dropped below zero");

      for (int itype = 1; itype <= atom->ntypes; itype++) {
        gfactor1[itype] = - gamma_p_grid[iz][iy][ix] / force->ftm2v;
        gfactor2[itype] = sqrt(24.0*force->boltz*gamma_p_grid[iz][iy][ix]/update->dt/force->mvv2e)
          / force->ftm2v;
      }

      double tsqrt = sqrt(T_electron[iz][iy][ix]);
      gamma1 = gfactor1[type[i]];
      gamma2 = gfactor2[type[i]] * tsqrt;
      if (T_electron[iz][iy][ix] > 0) {
        flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
        flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
        flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

        f[i][0] += flangevin[i][0];
        f[i][1] += flangevin[i][1];
        f[i][2] += flangevin[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMThermal::end_of_step()
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

  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        net_energy_transfer[iz][iy][ix] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;

      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  outflag = 0;
  MPI_Allreduce(&net_energy_transfer[0][0][0],&net_energy_transfer_all[0][0][0],
                ngridtotal,MPI_DOUBLE,MPI_SUM,world);

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve
  // This could be moved out of the loop with an appropriate trigger
  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double voxel_coeff =(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz);

  std::vector<double> grid_fourier(nzgrid * nygrid * nxgrid);
  int index = 0;  // Location unimportant, only max value
    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++)
        grid_fourier[index++] = 2.0/c_e_grid[iz][iy][ix]*(k_e_grid[iz][iy][ix]*voxel_coeff);

  double fourier_max = *std::max_element(grid_fourier.begin(), grid_fourier.end());

  double stability_criterion = 1.0 - fourier_max*inner_dt;

  if (stability_criterion < 0.0) {
    inner_dt = 1/fourier_max;
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++)
          T_electron_old[iz][iy][ix] = T_electron[iz][iy][ix];

    // compute new electron T profile

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          int xright = ix + 1;
          int yright = iy + 1;
          int zright = iz + 1;
          if (xright == nxgrid) xright = 0;
          if (yright == nygrid) yright = 0;
          if (zright == nzgrid) zright = 0;
          int xleft = ix - 1;
          int yleft = iy - 1;
          int zleft = iz - 1;
          if (xleft == -1) xleft = nxgrid - 1;
          if (yleft == -1) yleft = nygrid - 1;
          if (zleft == -1) zleft = nzgrid - 1;

          // Initialize flags for vacuum
          int left = 1;
          int right =1;
          int in = 1;
          int out = 1;
          int up = 1;
          int down = 1;

          // Set flags to 0 if vacuum
          if (T_electron[iz][iy][xleft] == 0) left = 0;
          if (T_electron[iz][iy][xright] == 0) right = 0;
          if (T_electron[iz][yright][ix] == 0) in = 0;
          if (T_electron[iz][yleft][ix] == 0) out = 0;
          if (T_electron[zright][iy][ix] == 0) up = 0;
          if (T_electron[zleft][iy][ix] == 0) down = 0;

          if (T_electron[iz][iy][ix] > 0) {
            T_electron[iz][iy][ix] =
              T_electron_old[iz][iy][ix] +
              inner_dt/c_e_grid[iz][iy][ix] *
              ((safe_effective_kappa(k_e_grid[iz][iy][xleft],k_e_grid[iz][iy][ix]))*
               (T_electron_old[iz][iy][xleft]-T_electron_old[iz][iy][ix])/dx/dx*left +

               (safe_effective_kappa(k_e_grid[iz][iy][xright],k_e_grid[iz][iy][ix]))*
               (T_electron_old[iz][iy][xright]-T_electron_old[iz][iy][ix])/dx/dx*right +

               (safe_effective_kappa(k_e_grid[iz][yleft][ix],k_e_grid[iz][iy][ix]))*
               (T_electron_old[iz][yleft][ix]-T_electron_old[iz][iy][ix])/dy/dy*out +

               (safe_effective_kappa(k_e_grid[iz][yright][ix],k_e_grid[iz][iy][ix]))*
               (T_electron_old[iz][yright][ix]-T_electron_old[iz][iy][ix])/dy/dy*in +

               (safe_effective_kappa(k_e_grid[zleft][iy][ix],k_e_grid[iz][iy][ix]))*
               (T_electron_old[zleft][iy][ix]-T_electron_old[iz][iy][ix])/dz/dz*down +

               (safe_effective_kappa(k_e_grid[zright][iy][ix],k_e_grid[iz][iy][ix]))*
               (T_electron_old[zright][iy][ix]-T_electron_old[iz][iy][ix])/dz/dz*up

               -(net_energy_transfer_all[iz][iy][ix])/(del_vol)
               +(inductive_power*inductive_response_grid[iz][iy][ix]));
          }
        }
  }

  // output of grid electron temperatures to file
  if (!outfile.empty() && (update->ntimestep % outevery == 0))
    write_electron_temperatures(fmt::format("{}.{}",outfile,update->ntimestep));
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only read by proc 0, grid values are Bcast to other procs
------------------------------------------------------------------------- */

void FixTTMThermal::read_electron_properties(const std::string &filename)
{
  if (comm->me == 0) {

    int ***prop_initial_set;
    memory->create(prop_initial_set,nzgrid,nygrid,nxgrid,"ttm:prop_initial_set");
    memset(&prop_initial_set[0][0][0],0,ngridtotal*sizeof(int));

    // read electron properties from file
    bigint nread = 0;

    try {
      PotentialFileReader reader(lmp, filename, "electron property grid");

      while (nread < ngridtotal) {
        // reader will skip over comment-only lines
        auto values = reader.next_values(7);
        ++nread;

        int ix = values.next_int() - 1;
        int iy = values.next_int() - 1;
        int iz = values.next_int() - 1;
        double c_e_tmp  = values.next_double();
        double k_e_tmp  = values.next_double();
        double gamma_p_tmp  = values.next_double();
        double ind_tmp  = values.next_double();

        // check correctness of input data

        if ((ix < 0) || (ix >= nxgrid) || (iy < 0) || (iy >= nygrid) || (iz < 0) || (iz >= nzgrid))
          throw TokenizerException("Fix ttm/thermal invalid grid index in grid file","");

        if (c_e_tmp < 0.0)
          throw TokenizerException("Fix ttm/thermal electron specific heat must be > 0.0","");

        if (k_e_tmp < 0.0)
          throw TokenizerException("Fix ttm/thermal electron conductivity must be > 0.0","");

        if (gamma_p_tmp < 0.0)
          throw TokenizerException("Fix ttm/thermal electron coupling must be > 0.0","");

        if (ind_tmp < 0.0)
          throw TokenizerException("Fix ttm/thermal electron inductive response must be >= 0.0","");

        c_e_grid[iz][iy][ix] = c_e_tmp;
        k_e_grid[iz][iy][ix] = k_e_tmp;
        gamma_p_grid[iz][iy][ix] = gamma_p_tmp;
        inductive_response_grid[iz][iy][ix] = ind_tmp;
        prop_initial_set[iz][iy][ix] = 1;
      }
    } catch (std::exception &e) {
      error->one(FLERR, Error::NOLASTLINE, e.what());
    }

    // check completeness of input data

    for (int iz = 0; iz < nzgrid; iz++)
      for (int iy = 0; iy < nygrid; iy++)
        for (int ix = 0; ix < nxgrid; ix++)
          if (prop_initial_set[iz][iy][ix] == 0)
            error->all(FLERR, Error::NOLASTLINE,
                       "Fix ttm/thermal infile did not set all properties");

    memory->destroy(prop_initial_set);
  }
  MPI_Bcast(&c_e_grid[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
  MPI_Bcast(&k_e_grid[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma_p_grid[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
  MPI_Bcast(&inductive_response_grid[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTMThermal::compute_vector(int n)
{
  if (outflag == 0) {
    e_energy = 0.0;
    transfer_energy = 0.0;

    int ix,iy,iz;

    double dx = domain->xprd/nxgrid;
    double dy = domain->yprd/nygrid;
    double dz = domain->zprd/nzgrid;
    double del_vol = dx*dy*dz;

    for (int iz = 0; iz < nzgrid; iz++) {
      for (int iy = 0; iy < nygrid; iy++) {
        for (int ix = 0; ix < nxgrid; ix++) {
          e_energy += T_electron[iz][iy][ix]*c_e_grid[iz][iy][ix]*del_vol;
          transfer_energy += net_energy_transfer_all[iz][iy][ix]*update->dt;
        }
      }
    }
    outflag = 1;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grids
------------------------------------------------------------------------- */

double FixTTMThermal::memory_usage()
{
  double bytes = FixTTM::memory_usage();
  bytes += (double) 4*ngridtotal * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMThermal::allocate_grid()
{
  FixTTM::allocate_grid();
  memory->create(c_e_grid,nzgrid,nygrid,nxgrid,"ttm:c_e_grid");
  memory->create(k_e_grid,nzgrid,nygrid,nxgrid,"ttm:k_e_grid");
  memory->create(gamma_p_grid,nzgrid,nygrid,nxgrid,"ttm:gamma_p_grid");
  memory->create(inductive_response_grid,nzgrid,nygrid,nxgrid,"ttm:inductive_response_grid");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMThermal::deallocate_grid()
{
  memory->destroy(c_e_grid);
  memory->destroy(k_e_grid);
  memory->destroy(gamma_p_grid);
  memory->destroy(inductive_response_grid);
}
