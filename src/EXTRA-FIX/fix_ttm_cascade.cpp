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
   Contributing authors: Original fix ttm
                         Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)

                         ttm/cascade
                         Andrés Rojano (Lancaster University)
                         Samuel T. Murphy (Lancaster University)

------------------------------------------------------------------------- */

#include "fix_ttm_cascade.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "grid3d.h"
#include "memory.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "random_mars.h"
#include "safe_pointers.h"
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

FixTTMCascade::FixTTMCascade(LAMMPS *lmp, int narg, char **arg)
    : FixTTMGrid(lmp, 13, arg) // 13 is to pass non-keyword arguments to FixTTM
{
  cutoff_active = false;
  offset_active = false;
  cetable_active = false;
  ketable_active = false;
  tinit = 0.0;

  // parsing all keywords, standard (set and infile) and custom (cutoff, offset, cetab, and ketab)
  int iarg = 13;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      tinit = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (tinit <= 0.0)
        error->all(FLERR, "Fix ttm/cascade initial temperature must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "infile") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      infile = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      cutoff_active = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "offset") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      offset_active = true;
      time_offset = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cetab") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      cetable_active = true;
      tableinterpreader(arg[iarg + 1], "ce");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ketab") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm/cascade command");
      ketable_active = true;
      tableinterpreader(arg[iarg + 1], "ke");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix ttm/cascade command");
    }
  }

  // error check

  if (seed <= 0)
    error->all(FLERR, "Invalid random number seed in fix ttm/cascade command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR, "Fix ttm/cascade electronic_specific_heat must be > 0.0");

  if (cetable_active) {
    for (int i = 0; i < ce_values.size(); i++) {
      if (ce_values[i] <= 0.0)
        error->all(FLERR, "Fix ttm/cascade all electronic_specific_heat entries "
                          "in the file must be > 0.0");
    }
  }

  if (ketable_active) {
    for (int i = 0; i < ke_values.size(); i++) {
      if (ke_values[i] <= 0.0)
        error->all(FLERR, "Fix ttm/cascade all electronic_thermal_conductivity "
                          "entries in the file must be > 0.0");
    }
  }

  if (offset_active) {
    if (time_offset < 0)
      error->all(FLERR, "Fix ttm/cascade time_offset must be >= 0");
  }

  if (offset_active && cutoff_active)
    error->all(
        FLERR,
        "Fix ttm/cascade cannot have cutoff and offset active at the same time");

  if (electronic_density <= 0.0)
    error->all(FLERR, "Fix ttm/cascade electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,
               "Fix ttm/cascade electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0)
    error->all(FLERR, "Fix ttm/cascade gamma_p must be > 0.0");
  if (gamma_s < 0.0)
    error->all(FLERR, "Fix ttm/cascade gamma_s must be >= 0.0");
  if (v_0 < 0.0)
    error->all(FLERR, "Fix ttm/cascade v_0 must be >= 0.0");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR, "Fix ttm/cascade grid sizes must be > 0");
}

/* ---------------------------------------------------------------------- */

FixTTMCascade::~FixTTMCascade()
{
  //FixTTMGrid::~FixTTMGrid(); // Deals with this
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::post_force(int /*vflag*/)
{
  int ix,iy,iz;
  double gamma1,gamma2,gamma_cutoff,gamma_offset;

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

  update->update_time();
  double current_time = update->atime; // look for the current time
  gamma_offset = (offset_active && (current_time < time_offset)) ? 0 : 1;

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

      gamma_cutoff = 1.0;

      // evaluates whether the cutoff approach applies in the computation
      if (vsq > v_0_sq) {
        if (cutoff_active) {
          gamma1 *=
              (gamma_s / gamma_p); // avoids gamma_p for fast-moving atoms
          gamma_cutoff = 0;
        } else {
          gamma1 *= ((gamma_p * gamma_offset) + gamma_s) /
                    gamma_p; // standard ttm approach
        }
      } else {
        gamma1 *= gamma_offset; // decouples gamma_p when time_offset applies
      }

      gamma2 = gfactor2[type[i]] * tsqrt * gamma_cutoff * gamma_offset; // modifies gamma2 to disable stochastic forces if needed

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }

  if (flag) error->one(FLERR,"Out of range fix ttm/cascade atoms");
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::end_of_step(){
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

  double variable_electronic_specific_heat;
  double variable_electronic_thermal_conductivity;
  double el_th_diffusivity;
  double el_th_diffusivity_max = 0.0;
  double el_th_diffusivity_global_max = 0.0;

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

  if (cetable_active || ketable_active) {

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          variable_electronic_specific_heat = cetable_active ? linearinterpolation(T_electron[iz][iy][ix], "ce") : electronic_specific_heat;
          variable_electronic_thermal_conductivity = ketable_active ? linearinterpolation(T_electron[iz][iy][ix], "ke") : electronic_thermal_conductivity;
          el_th_diffusivity = variable_electronic_thermal_conductivity/(variable_electronic_specific_heat*electronic_density);
          if (el_th_diffusivity > el_th_diffusivity_max) el_th_diffusivity_max = el_th_diffusivity;
        }
    MPI_Allreduce(&el_th_diffusivity_max, &el_th_diffusivity_global_max, 1, MPI_DOUBLE, MPI_MAX, world);
  }


  else {
    el_th_diffusivity_global_max = electronic_thermal_conductivity/(electronic_specific_heat*electronic_density);
  }

  double stability_criterion = 1.0 -
    2.0*inner_dt* el_th_diffusivity_global_max * (dxinv*dxinv + dyinv*dyinv + dzinv*dzinv);

  if (stability_criterion < 0.0) {
    inner_dt = 0.5/ el_th_diffusivity_global_max / (dxinv*dxinv + dyinv*dyinv + dzinv*dzinv);
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/cascade");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    memcpy(&T_electron_old[nzlo_out][nylo_out][nxlo_out],
           &T_electron[nzlo_out][nylo_out][nxlo_out],ngridout*sizeof(double));

  // store thermal conductivity in a grid for rapid access

  if(ketable_active){

    memset(&thermal_conductivity_grid[nzlo_out][nylo_out][nxlo_out],0, ngridout*sizeof(double));

    for (iz = nzlo_out; iz <= nzhi_out; iz++)
      for (iy = nylo_out; iy <= nyhi_out; iy++)
        for (ix = nxlo_out; ix <= nxhi_out; ix++){
          thermal_conductivity_grid[iz][iy][ix] = linearinterpolation(T_electron_old[iz][iy][ix], "ke");
        }

  }


    // compute new electron T profile

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          variable_electronic_specific_heat = cetable_active ? linearinterpolation(T_electron_old[iz][iy][ix], "ce") : electronic_specific_heat;

          T_electron[iz][iy][ix] =
            T_electron_old[iz][iy][ix] +
            inner_dt/(variable_electronic_specific_heat*electronic_density) *
            (heat_flux_gradient(ix, iy, iz, dxinv, dyinv, dzinv) -
             net_energy_transfer[iz][iy][ix]/volgrid);
        }

    // communicate new T_electron values to ghost grid points

    grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMCascade::allocate_grid()
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

  memory->create(grid_buf1, ngrid_buf1, "ttm/cascade:grid_buf1");
  memory->create(grid_buf2, ngrid_buf2, "ttm/cascade:grid_buf2");

  memory->create3d_offset(T_electron_old, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:T_electron_old");
  memory->create3d_offset(T_electron, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "ttm/cascade:T_electron");
  memory->create3d_offset(net_energy_transfer, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:net_energy_transfer");
  memory->create3d_offset(thermal_conductivity_grid, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:thermal_conductivity_grid");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMCascade::deallocate_grid()
{
  delete grid;
  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(T_electron, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(thermal_conductivity_grid, nzlo_out, nylo_out, nxlo_out);
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem
   or the net_energy transfer between the subsystems
------------------------------------------------------------------------- */

double FixTTMCascade::compute_vector(int n)
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
              integrated_ce(T_electron[iz][iy][ix]) * electronic_density * volgrid;
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
   read the specific heat and thermal conductivity table files
------------------------------------------------------------------------- */

void FixTTMCascade::tableinterpreader(const std::string &filename,
                                   const std::string &keyword) {

  std::vector<double> &temp_vals =
      (keyword == "ce") ? temp_ce_values : temp_ke_values;
  std::vector<double> &dtemp_vals =
      (keyword == "ce") ? dtemp_ce_values : dtemp_ke_values;
  std::vector<double> &y_vals = (keyword == "ce") ? ce_values : ke_values;
  std::vector<double> &dy_vals = (keyword == "ce") ? dce_values : dke_values;

  int table_count = 0;

  if (comm->me == 0) {
    std::string table_label = (keyword == "ce") ? "specific heat table"
                                                : "thermal conductivity table";
    PotentialFileReader reader(lmp, filename, table_label);
    while (char *line = reader.next_line()) {
      double temp_value, y_value;
      if (sscanf(line, "%lg %lg", &temp_value, &y_value) == 2) {
        temp_vals.push_back(temp_value);
        y_vals.push_back(y_value);
      }
    }

    int nsize_table = static_cast<int>(temp_vals.size());
    dtemp_vals.resize(nsize_table);
    dy_vals.resize(nsize_table);
    if (keyword == "ce") {
      ce_integral_values.resize(nsize_table);
      ce_integral_values[0] = y_vals[0] * temp_vals[0];
    }
    for (int i = 0; i < nsize_table - 1; i++) {
      dtemp_vals[i] = temp_vals[i + 1] - temp_vals[i];
      dy_vals[i] = y_vals[i + 1] - y_vals[i];
      if (keyword == "ce") {
        ce_integral_values[i + 1] = ce_integral_values[i] + 0.5 * (y_vals[i + 1] + y_vals[i]) * (temp_vals[i + 1] - temp_vals[i]);
      }
      if (temp_vals[i] >= temp_vals[i + 1]) table_count += 1;
    }
  }


  int nsize_table = static_cast<int>(temp_vals.size());
  MPI_Bcast(&nsize_table, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    temp_vals.resize(nsize_table);
    y_vals.resize(nsize_table);
    dtemp_vals.resize(nsize_table);
    dy_vals.resize(nsize_table);
    if (keyword == "ce")
      ce_integral_values.resize(nsize_table);
  }

  MPI_Bcast(temp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(y_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dtemp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dy_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(&table_count, 1, MPI_INT, 0, world);
  if (keyword == "ce") MPI_Bcast(ce_integral_values.data(), nsize_table, MPI_DOUBLE, 0, world);

  // error check

  if (table_count) error->all(FLERR, "There are {} non-increasing temperature intervals in the {} table", table_count, keyword);

  if (temp_vals.size() < 2) error->all(FLERR, "The {} table has less than 2 rows, interpolation is invalid", keyword);


  }

/* ----------------------------------------------------------------------
   interpolates the specific heat and thermal conductivity values
------------------------------------------------------------------------- */

double FixTTMCascade::linearinterpolation(double temp, const std::string &keyword) {

  const std::vector<double> &temp_vals =
      (keyword == "ce") ? temp_ce_values : temp_ke_values;
  const std::vector<double> &dtemp_vals =
      (keyword == "ce") ? dtemp_ce_values : dtemp_ke_values;
  const std::vector<double> &y_vals = (keyword == "ce") ? ce_values : ke_values;
  const std::vector<double> &dy_vals =
      (keyword == "ce") ? dce_values : dke_values;

  if (temp_vals.empty())
    error->one(FLERR, "The table {} is empty", keyword);

  int lo = 0;
  int hi = static_cast<int>(temp_vals.size()) - 1;

  if (temp <= temp_vals[lo])
    return y_vals[lo];
  if (temp >= temp_vals[hi])
    return y_vals[hi];

  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (temp < temp_vals[mid])
      hi = mid;
    else
      lo = mid;
  }

  return y_vals[lo] + (temp - temp_vals[lo]) * (dy_vals[lo] / dtemp_vals[lo]);
}

/* ----------------------------------------------------------------------
   performs the numerical integration of the specific heat
------------------------------------------------------------------------- */

double FixTTMCascade::integrated_ce(double Te) {
  if (!cetable_active) {
    return electronic_specific_heat * Te;
  }

  int lo = 0;
  int hi = static_cast<int>(ce_values.size()) - 1;
  int mid;

  if (Te <= temp_ce_values[lo])
    return ce_values[lo] * Te;
  if (Te >= temp_ce_values[hi])
    return ce_integral_values[hi] + (ce_values[hi] * (Te - temp_ce_values[hi]));

  while (hi - lo > 1) {
    mid = (lo + hi) / 2;
    if (Te < temp_ce_values[mid])
      hi = mid;
    else
      lo = mid;
  }

  double Ce_interpolated =
      ce_values[lo] +
      (Te - temp_ce_values[lo]) * (dce_values[lo] / dtemp_ce_values[lo]);

  return ce_integral_values[lo] +
         0.5 * (ce_values[lo] + Ce_interpolated) * (Te - temp_ce_values[lo]);
}

/* ----------------------------------------------------------------------
   computes the heat flux gradient according to the type of thermal conductivity.
   constant ke: k \nabla^2 T
   variable ke: k \nabla^2 T + \nabla k \nabla T
------------------------------------------------------------------------- */

double FixTTMCascade::heat_flux_gradient(int ix, int iy, int iz, double dxinv, double dyinv,
                       double dzinv) {

  double heat_flux_gradient;

  if (!ketable_active) {

    heat_flux_gradient = electronic_thermal_conductivity * ((T_electron_old[iz][iy][ix-1] + T_electron_old[iz][iy][ix+1] -
               2.0*T_electron_old[iz][iy][ix])*dxinv*dxinv +
              (T_electron_old[iz][iy-1][ix] + T_electron_old[iz][iy+1][ix] -
               2.0*T_electron_old[iz][iy][ix])*dyinv*dyinv +
              (T_electron_old[iz-1][iy][ix] + T_electron_old[iz+1][iy][ix] -
               2.0*T_electron_old[iz][iy][ix])*dzinv*dzinv);

    return heat_flux_gradient;

  }


  else {

    heat_flux_gradient = thermal_conductivity_grid[iz][iy][ix] * ((T_electron_old[iz][iy][ix-1] + T_electron_old[iz][iy][ix+1] -
            2.0*T_electron_old[iz][iy][ix])*dxinv*dxinv +
          (T_electron_old[iz][iy-1][ix] + T_electron_old[iz][iy+1][ix] -
            2.0*T_electron_old[iz][iy][ix])*dyinv*dyinv +
          (T_electron_old[iz-1][iy][ix] + T_electron_old[iz+1][iy][ix] -
            2.0*T_electron_old[iz][iy][ix])*dzinv*dzinv);

    heat_flux_gradient +=

    (thermal_conductivity_grid[iz][iy][ix+1] -
    thermal_conductivity_grid[iz][iy][ix-1])/2.0*dxinv *
    (T_electron_old[iz][iy][ix+1] - T_electron_old[iz][iy][ix-1])/2.0*dxinv +
    (thermal_conductivity_grid[iz][iy+1][ix] -
    thermal_conductivity_grid[iz][iy-1][ix])/2.0*dyinv *
    (T_electron_old[iz][iy+1][ix] - T_electron_old[iz][iy-1][ix])/2.0*dyinv +
    (thermal_conductivity_grid[iz+1][iy][ix] -
    thermal_conductivity_grid[iz-1][iy][ix])/2.0*dzinv *
    (T_electron_old[iz+1][iy][ix] - T_electron_old[iz-1][iy][ix])/2.0*dzinv;


    return heat_flux_gradient;
  }
}
