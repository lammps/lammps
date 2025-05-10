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

#include "compute_esp_grid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "update.h"
#include "domain.h"
#include "pair.h"
#include "utils.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeESPGridKokkos<DeviceType>::ComputeESPGridKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeESPGrid(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | Q_MASK | TYPE_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeESPGridKokkos<DeviceType>::~ComputeESPGridKokkos()
{
  if (copymode) return;
  
  // No need to free device memory - DualView destructor will handle it
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeESPGridKokkos<DeviceType>::init()
{
  // Call parent init to set up ReaxFF parameters if needed
  ComputeESPGrid::init();
  
  // Create DualView for bcut_acks2
  if (reaxflag) {
    k_bcut_acks2 = tdual_1d_double("esp/grid/kk:bcut_acks2", atom->ntypes+1);
    k_bcut_acks2.modify_host();
    for (int i = 1; i <= atom->ntypes; i++) {
      k_bcut_acks2.h_view(i) = bcut_acks2[i]; 
    }
    k_bcut_acks2.sync_device();
  } else {
    k_bcut_acks2 = tdual_1d_double("esp/grid/kk:bcut_acks2_dummy", 1);
    k_bcut_acks2.modify_host();
    k_bcut_acks2.h_view(0) = bcut_acks2[0];
    k_bcut_acks2.sync_device();
  }
  
  // Create DualViews for the grid data
  int nz = izhi-izlo+1;
  int ny = iyhi-iylo+1;
  int nx = ixhi-ixlo+1;
  
  k_esp = tdual_3d_double("esp/grid/kk:esp", nz, ny, nx);
  k_weight = tdual_3d_double("esp/grid/kk:weight", nz, ny, nx);
  k_reference = tdual_3d_double("esp/grid/kk:reference", nz, ny, nx);
  
  // Copy reference data to device
  k_reference.modify_host();
  for (int iz = izlo; iz <= izhi; iz++) {
    for (int iy = iylo; iy <= iyhi; iy++) {
      for (int ix = ixlo; ix <= ixhi; ix++) {
        k_reference.h_view(iz-izlo, iy-iylo, ix-ixlo) = reference[iz][iy][ix];
      }
    }
  }
  k_reference.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeESPGridKokkos<DeviceType>::compute_pergrid()
{
  if (invoked_pergrid == update->ntimestep) return;
  invoked_pergrid = update->ntimestep;
  
  // Get atom data
  atomKK->sync(execution_space, datamask_read);
  
  // Get Kokkos views of atom data
  d_x = atomKK->k_x.view<DeviceType>();
  d_q = atomKK->k_q.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  
  // Calculate grid spacing for Kokkos kernel
  double dx = (oxhi - oxlo) / std::max(1, ixhi - ixlo);
  double dy = (oyhi - oylo) / std::max(1, iyhi - iylo);
  double dz = (ozhi - ozlo) / std::max(1, izhi - izlo);
  
  // Debug grid spacing (on host)
  utils::logmesg(lmp, "*** Grid spacing: dx={:.6f}, dy={:.6f}, dz={:.6f}\n", dx, dy, dz);
  
  // Get local copies for lambda capture
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;
  
  int ixlo_d = ixlo;
  int iylo_d = iylo;
  int izlo_d = izlo;
  int ixhi_d = ixhi;
  int iyhi_d = iyhi;
  int izhi_d = izhi;
  
  double oxlo_d = oxlo;
  double oylo_d = oylo;
  double ozlo_d = ozlo;
  
  bool reaxflag_d = reaxflag;
  int ntypes = atom->ntypes;
  
  // Initialize device views
  k_esp.modify_device();
  k_weight.modify_device();
  auto d_esp = k_esp.d_view;
  auto d_weight = k_weight.d_view;
  auto d_reference = k_reference.d_view;
  auto d_bcut_acks2 = k_bcut_acks2.d_view;
  
  Kokkos::deep_copy(d_esp, 0.0);
  Kokkos::deep_copy(d_weight, 0.0);
  
  // Calculate total grid points and optimal team size
  int nz = izhi_d - izlo_d + 1;
  int ny = iyhi_d - iylo_d + 1;
  int nx = ixhi_d - ixlo_d + 1;
  int num_grid_points = nz * ny * nx;
  
  // Use hierarchical parallelism with TeamPolicy
  // Each team handles a chunk of grid points, vector lanes handle atom loops
  int team_size = 0; // Let Kokkos choose automatically
  int vector_length = 32; // Often 32 for NVIDIA GPUs, 16 for AMD

  // Define number of teams - enough to cover all grid points
  int points_per_team = 4;  // Process multiple grid points per team for better efficiency
  int num_teams = (num_grid_points + points_per_team - 1) / points_per_team;

  using team_policy_type = Kokkos::TeamPolicy<DeviceType>;
  using member_type = typename team_policy_type::member_type;
  team_policy_type team_policy(num_teams, team_size, vector_length);

  // Launch compute_pergrid kernel with TeamPolicy
  Kokkos::parallel_for("ComputeESPGridKokkos::compute_pergrid", 
    team_policy_type(num_teams, team_size, vector_length),
    KOKKOS_LAMBDA(const member_type& team_member) {
      // Each team processes multiple grid points
      const int team_idx = team_member.league_rank();
      const int start_idx = team_idx * points_per_team;
      const int end_idx = Kokkos::min(start_idx + points_per_team, num_grid_points);
      
      // Loop over assigned grid points
      for (int point_idx = start_idx; point_idx < end_idx; ++point_idx) {
        // Convert flat index to 3D grid coordinates
        int ix_idx = point_idx % nx;
        int iy_idx = (point_idx / nx) % ny;
        int iz_idx = point_idx / (nx * ny);
        
        // Compute actual grid indices
        int ix = ix_idx + ixlo_d - ixlo_d;
        int iy = iy_idx + iylo_d - iylo_d;
        int iz = iz_idx + izlo_d - izlo_d;
        
        // Calculate grid point coordinates
        double gx = oxlo_d + (ix + 0.5) * dx;
        double gy = oylo_d + (iy + 0.5) * dy;
        double gz = ozlo_d + (iz + 0.5) * dz;
        
        // Initialize grid values for this point
        d_esp(iz, iy, ix) = 0.0;
        d_weight(iz, iy, ix) = 0.0;
        
        // Find nearest atom including ghost atoms
        double rmin2 = 1e60;
        int tnear = -1;
        
        // Use vector parallelism for the atom loop to find nearest atom
        double thread_rmin2 = 1e60;
        int thread_tnear = -1;
        
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team_member, ntotal),
          [&] (const int& i, double& min_r2, int& nearest_type) {
            double dx_i = gx - d_x(i, 0);
            double dy_i = gy - d_x(i, 1);
            double dz_i = gz - d_x(i, 2);
            
            // Apply minimum image for periodic boundaries
            if (domain->xperiodic) {
              double xprd = domain->xprd;
              if (dx_i > 0.5*xprd) dx_i -= xprd;
              else if (dx_i < -0.5*xprd) dx_i += xprd;
            }
            if (domain->yperiodic) {
              double yprd = domain->yprd;
              if (dy_i > 0.5*yprd) dy_i -= yprd;
              else if (dy_i < -0.5*yprd) dy_i += yprd;
            }
            if (domain->zperiodic) {
              double zprd = domain->zprd;
              if (dz_i > 0.5*zprd) dz_i -= zprd;
              else if (dz_i < -0.5*zprd) dz_i += zprd;
            }
            
            double r2 = dx_i*dx_i + dy_i*dy_i + dz_i*dz_i;
            if (r2 < min_r2) {
              min_r2 = r2;
              nearest_type = d_type(i);
            }
          }, Kokkos::Min<double>(thread_rmin2), Kokkos::Min<int>(thread_tnear));
        
        // Update rmin2 and tnear with thread results
        if (thread_rmin2 < rmin2) {
          rmin2 = thread_rmin2;
          tnear = thread_tnear;
        }
        
        // Early return if no valid atom type found
        if (tnear < 1 || tnear > ntypes) continue;
        
        // Get cutoff radius
        double rcut;
        if (reaxflag_d) {
          rcut = d_bcut_acks2(tnear);
        } else {
          rcut = d_bcut_acks2(0);
        }
        
        double r = sqrt(rmin2);
        
        // Skip if distance is invalid
        if (r < 1.4 || r > rcut) {
          d_weight(iz, iy, ix) = 0.0;
          continue;
        }
        
        // Calculate weight
        d_weight(iz, iy, ix) = compute_weight(r, rcut);
        
        // Calculate ESP contributions using vector parallelism
        double esp_sum = 0.0;
        
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team_member, ntotal),
          [&] (const int& i, double& esp_local) {
            double dx_i = gx - d_x(i, 0);
            double dy_i = gy - d_x(i, 1);
            double dz_i = gz - d_x(i, 2);
            
            // Apply minimum image for periodic boundaries
            if (domain->xperiodic) {
              double xprd = domain->xprd;
              if (dx_i > 0.5*xprd) dx_i -= xprd;
              else if (dx_i < -0.5*xprd) dx_i += xprd;
            }
            if (domain->yperiodic) {
              double yprd = domain->yprd;
              if (dy_i > 0.5*yprd) dy_i -= yprd;
              else if (dy_i < -0.5*yprd) dy_i += yprd;
            }
            if (domain->zperiodic) {
              double zprd = domain->zprd;
              if (dz_i > 0.5*zprd) dz_i -= zprd;
              else if (dz_i < -0.5*zprd) dz_i += zprd;
            }
            
            double r2 = dx_i*dx_i + dy_i*dy_i + dz_i*dz_i + 1e-12;
            esp_local += d_q(i) / sqrt(r2);
          }, esp_sum);
        
        // Store result
        d_esp(iz, iy, ix) = esp_sum;
      }
    });

  
  // Sync DualViews
  k_esp.modify_device();
  k_weight.modify_device();
  k_esp.sync_host();
  k_weight.sync_host();
  
  // Update host arrays
  for (int iz = izlo; iz <= izhi; iz++) {
    for (int iy = iylo; iy <= iyhi; iy++) {
      for (int ix = ixlo; ix <= ixhi; ix++) {
        esp[iz][iy][ix] = k_esp.h_view(iz-izlo, iy-iylo, ix-ixlo);
        weight[iz][iy][ix] = k_weight.h_view(iz-izlo, iy-iylo, ix-ixlo);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeESPGridKokkos<DeviceType>::compute_scalar()
{
  if (invoked_pergrid != update->ntimestep) compute_pergrid();
  invoked_scalar = update->ntimestep;
  
  // Ensure device views are synced
  k_esp.sync_device();
  k_weight.sync_device();
  k_reference.sync_device();
  
  auto d_esp = k_esp.d_view;
  auto d_weight = k_weight.d_view;
  auto d_reference = k_reference.d_view;
  
  // Create TeamPolicy for reduction
  int team_size = 0; // Let Kokkos choose automatically
  int vector_length = 32; // Often 32 for NVIDIA GPUs, 16 for AMD

  // Calculate number of teams for reduction
  int points_per_team = 16;  // Process multiple grid points per team
  int num_grid_points = (izhi-izlo+1) * (iyhi-iylo+1) * (ixhi-ixlo+1);
  int num_teams = (num_grid_points + points_per_team - 1) / points_per_team;

  using team_policy_type = Kokkos::TeamPolicy<DeviceType>;
  using member_type = typename team_policy_type::member_type;

  double local_loss = 0.0;
  double local_weight = 0.0;

  Kokkos::parallel_reduce("ComputeESPGridKokkos::compute_scalar",
    team_policy_type(num_teams, team_size, vector_length),
    KOKKOS_LAMBDA(const member_type& team_member, double& team_loss, double& team_weight) {
      const int team_idx = team_member.league_rank();
      const int start_idx = team_idx * points_per_team;
      const int end_idx = Kokkos::min(start_idx + points_per_team, num_grid_points);

      double thread_loss = 0.0;
      double thread_weight = 0.0;

      for (int point_idx = start_idx; point_idx < end_idx; ++point_idx) {
        int ix_len = ixhi - ixlo + 1;
        int iy_len = iyhi - iylo + 1;

        int ix = point_idx % ix_len;
        int iy = (point_idx / ix_len) % iy_len;
        int iz = point_idx / (ix_len * iy_len);

        double w = d_weight(iz, iy, ix);
        if (w == 0.0) continue;

        double diff = d_esp(iz, iy, ix) - d_reference(iz, iy, ix);
        thread_loss += w * diff * diff;
        thread_weight += w;
      }

      team_loss += thread_loss;
      team_weight += thread_weight;
    },
    local_loss, local_weight
  );

  double local_result[2] = {local_loss, local_weight};
  double global_result[2];
  MPI_Allreduce(local_result, global_result, 2, MPI_DOUBLE, MPI_SUM, world);

  scalar = (global_result[1] > 0.0) ? global_result[0] / global_result[1] : 0.0;
  if (comm->me == 0) utils::logmesg(lmp, "*** scalar {}\n", scalar);
  
  return scalar;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double ComputeESPGridKokkos<DeviceType>::compute_weight(double r, double rcut) const
{
  double w = 1.0/(r*r);
  w *= 1.0-exp(-pow((r-1.4)/0.3, 6));
  w *= exp(-pow(r/rcut, 6));
  return w;
}

/* ---------------------------------------------------------------------- */

// Explicit instantiation for device types
namespace LAMMPS_NS {
template class ComputeESPGridKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeESPGridKokkos<LMPHostType>;
#endif
}
