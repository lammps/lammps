/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.

   Package-local PPPM source/sensor potential projection for xTB QM/MM.
------------------------------------------------------------------------- */

#include "pppm_tip4p_xtb.h"
#include "pppm_xtb.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "grid3d.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr FFT_SCALAR ZEROF = 0.0;

namespace LAMMPS_NS {

// The concrete KSpace style headers only depend on their core PPPM parents.
// This friend helper contains the shared package implementation while still
// allowing both derived styles to access PPPM's protected mesh state.
template <class Solver> class QMMMXTBPPPMHelper {
 public:
  static void compute_group_potential(Solver &, double *, int, int, bool);

 private:
  static void make_rho_group(Solver &, int, bool);
  static void project_group_potential(Solver &, double *, int);
};

}    // namespace LAMMPS_NS

/* ----------------------------------------------------------------------
   compute the mesh potential from one atom group at another atom group

   Sensors may have zero charge, unlike a potential reconstructed from
   per-atom energy.  A subsequent regular PPPM compute rebuilds the mesh data
   overwritten by this package-private source/sensor operation.
------------------------------------------------------------------------- */

template <class Solver>
void QMMMXTBPPPMHelper<Solver>::compute_group_potential(Solver &solver, double *potential,
                                                        int sensor_groupbit, int source_groupbit,
                                                        bool invert_source)
{
  if (solver.triclinic)
    solver.error->all(FLERR, "KSpace style pppm/xtb requires an orthorhombic box");

  solver.boxlo = solver.domain->boxlo;
  if (solver.atom->nmax > solver.nmax) {
    solver.memory->destroy(solver.part2grid);
    solver.nmax = solver.atom->nmax;
    solver.memory->create(solver.part2grid, solver.nmax, 3, "pppm/xtb:part2grid");
  }

  // u_brick is normally absent for analytic differentiation unless per-atom
  // energy is requested.  The scalar projection needs it temporarily.
  const bool temporary_u_brick = (solver.u_brick == nullptr);
  if (temporary_u_brick)
    solver.memory->create3d_offset(solver.u_brick, solver.nzlo_out, solver.nzhi_out,
                                   solver.nylo_out, solver.nyhi_out, solver.nxlo_out,
                                   solver.nxhi_out, "pppm/xtb:u_brick");

  solver.particle_map();
  make_rho_group(solver, source_groupbit, invert_source);
  solver.gc->reverse_comm(Grid3d::KSPACE, &solver, KSpace::REVERSE_RHO, 1, sizeof(FFT_SCALAR),
                          solver.gc_buf1, solver.gc_buf2, MPI_FFT_SCALAR);
  solver.brick2fft();

  for (int i = 0, n = 0; i < solver.nfft; ++i) {
    solver.work1[n++] = solver.density_fft[i];
    solver.work1[n++] = ZEROF;
  }
  solver.fft1->compute(solver.work1, solver.work1, FFT3d::FORWARD);

  for (int i = 0, n = 0; i < solver.nfft; ++i) {
    solver.work2[n] = solver.work1[n] * solver.greensfn[i];
    ++n;
    solver.work2[n] = solver.work1[n] * solver.greensfn[i];
    ++n;
  }
  solver.fft2->compute(solver.work2, solver.work2, FFT3d::BACKWARD);
  for (int k = solver.nzlo_in, n = 0; k <= solver.nzhi_in; ++k)
    for (int j = solver.nylo_in; j <= solver.nyhi_in; ++j)
      for (int i = solver.nxlo_in; i <= solver.nxhi_in; ++i) {
        solver.u_brick[k][j][i] = solver.work2[n];
        n += 2;
      }

  solver.gc->forward_comm(Grid3d::KSPACE, &solver, KSpace::FORWARD_AD, 1, sizeof(FFT_SCALAR),
                          solver.gc_buf1, solver.gc_buf2, MPI_FFT_SCALAR);
  project_group_potential(solver, potential, sensor_groupbit);

  if (temporary_u_brick)
    solver.memory->destroy3d_offset(solver.u_brick, solver.nzlo_out, solver.nylo_out,
                                    solver.nxlo_out);
}

/* ----------------------------------------------------------------------
   create a discretized density from a selected source group
------------------------------------------------------------------------- */

template <class Solver>
void QMMMXTBPPPMHelper<Solver>::make_rho_group(Solver &solver, int source_groupbit,
                                               bool invert_source)
{
  std::memset(&(solver.density_brick[solver.nzlo_out][solver.nylo_out][solver.nxlo_out]), 0,
              solver.ngrid * sizeof(FFT_SCALAR));

  for (int i = 0; i < solver.atom->nlocal; ++i) {
    const bool in_source = !!(solver.atom->mask[i] & source_groupbit) != invert_source;
    if (!in_source) continue;

    double site[3], weights[3];
    int indices[3];
    solver.get_charge_site(i, site, indices, weights);

    const int nx = solver.part2grid[i][0];
    const int ny = solver.part2grid[i][1];
    const int nz = solver.part2grid[i][2];
    const FFT_SCALAR dx = nx + solver.shiftone - (site[0] - solver.boxlo[0]) * solver.delxinv;
    const FFT_SCALAR dy = ny + solver.shiftone - (site[1] - solver.boxlo[1]) * solver.delyinv;
    const FFT_SCALAR dz = nz + solver.shiftone - (site[2] - solver.boxlo[2]) * solver.delzinv;

    solver.compute_rho1d(dx, dy, dz);
    const FFT_SCALAR z0 = solver.delvolinv * solver.atom->q[i];
    for (int n = solver.nlower; n <= solver.nupper; ++n) {
      const int mz = n + nz;
      const FFT_SCALAR y0 = z0 * solver.rho1d[2][n];
      for (int m = solver.nlower; m <= solver.nupper; ++m) {
        const int my = m + ny;
        const FFT_SCALAR x0 = y0 * solver.rho1d[1][m];
        for (int l = solver.nlower; l <= solver.nupper; ++l) {
          const int mx = l + nx;
          solver.density_brick[mz][my][mx] += x0 * solver.rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate the scalar mesh potential to selected sensor atoms
------------------------------------------------------------------------- */

template <class Solver>
void QMMMXTBPPPMHelper<Solver>::project_group_potential(Solver &solver, double *potential,
                                                        int sensor_groupbit)
{
  const bigint ngridtotal = static_cast<bigint>(solver.nx_pppm) * solver.ny_pppm * solver.nz_pppm;
  const double scaleinv = 1.0 / ngridtotal;

  for (int i = 0; i < solver.atom->nlocal; ++i) {
    if (!(solver.atom->mask[i] & sensor_groupbit)) continue;

    double site[3], weights[3];
    int indices[3];
    solver.get_charge_site(i, site, indices, weights);

    const int nix = solver.part2grid[i][0];
    const int niy = solver.part2grid[i][1];
    const int niz = solver.part2grid[i][2];
    const FFT_SCALAR dix = nix + solver.shiftone - (site[0] - solver.boxlo[0]) * solver.delxinv;
    const FFT_SCALAR diy = niy + solver.shiftone - (site[1] - solver.boxlo[1]) * solver.delyinv;
    const FFT_SCALAR diz = niz + solver.shiftone - (site[2] - solver.boxlo[2]) * solver.delzinv;
    solver.compute_rho1d(dix, diy, diz);

    double value = 0.0;
    for (int ni = solver.nlower; ni <= solver.nupper; ++ni) {
      const double iz0 = solver.rho1d[2][ni];
      const int miz = ni + niz;
      for (int mi = solver.nlower; mi <= solver.nupper; ++mi) {
        const double iy0 = iz0 * solver.rho1d[1][mi];
        const int miy = mi + niy;
        for (int li = solver.nlower; li <= solver.nupper; ++li) {
          const int mix = li + nix;
          value += iy0 * solver.rho1d[0][li] * solver.u_brick[miz][miy][mix];
        }
      }
    }
    potential[i] = value * scaleinv;
  }
}

PPPMXTB::PPPMXTB(LAMMPS *lmp) : PPPM(lmp)
{ xtbflag = 1; }

void PPPMXTB::compute_group_potential(double *potential, int sensor_groupbit, int source_groupbit,
                                      bool invert_source)
{
  QMMMXTBPPPMHelper<PPPMXTB>::compute_group_potential(*this, potential, sensor_groupbit,
                                                      source_groupbit, invert_source);
}

int PPPMXTB::get_charge_site(int i, double *site, int *indices, double *weights)
{
  for (int dim = 0; dim < 3; ++dim) site[dim] = atom->x[i][dim];
  indices[0] = i;
  weights[0] = 1.0;
  return 1;
}

PPPMTIP4PXTB::PPPMTIP4PXTB(LAMMPS *lmp) : PPPMTIP4P(lmp)
{ xtbflag = 1; }

void PPPMTIP4PXTB::compute_group_potential(double *potential, int sensor_groupbit,
                                           int source_groupbit, bool invert_source)
{
  QMMMXTBPPPMHelper<PPPMTIP4PXTB>::compute_group_potential(*this, potential, sensor_groupbit,
                                                           source_groupbit, invert_source);
}

/* ----------------------------------------------------------------------
   TIP4P M-site position and linear O/H/H force redistribution
------------------------------------------------------------------------- */

int PPPMTIP4PXTB::get_charge_site(int i, double *site, int *indices, double *weights)
{
  if (atom->type[i] != typeO) {
    for (int dim = 0; dim < 3; ++dim) site[dim] = atom->x[i][dim];
    indices[0] = i;
    weights[0] = 1.0;
    return 1;
  }

  int iH1 = atom->map(atom->tag[i] + 1);
  int iH2 = atom->map(atom->tag[i] + 2);
  if (iH1 == -1 || iH2 == -1) error->one(FLERR, "TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

  iH1 = domain->closest_image(i, iH1);
  iH2 = domain->closest_image(i, iH2);
  for (int dim = 0; dim < 3; ++dim) {
    const double midpoint_delta =
        0.5 * (atom->x[iH1][dim] + atom->x[iH2][dim] - 2.0 * atom->x[i][dim]);
    site[dim] = atom->x[i][dim] + alpha * midpoint_delta;
  }

  indices[0] = i;
  indices[1] = iH1;
  indices[2] = iH2;
  weights[0] = 1.0 - alpha;
  weights[1] = 0.5 * alpha;
  weights[2] = 0.5 * alpha;
  return 3;
}
