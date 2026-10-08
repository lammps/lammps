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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#include "pppm_electrode_tip4p.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "boundary_correction.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "remap_wrap.h"
#include "slab_dipole.h"
#include "update.h"
#include "wire_dipole.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr int OFFSET = 16384;
static constexpr FFT_SCALAR ZEROF = 0.0;

/* ---------------------------------------------------------------------- */

PPPMElectrodeTIP4P::PPPMElectrodeTIP4P(LAMMPS *lmp) :
    PPPMElectrode(lmp)
{

  group_group_enable = 0;
  tip4pflag = 1;

}
/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   TIP4P-aware slab boundary correction
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::compute_boundary_corr(double qsum_in, int eflag_atom_in,
                                               int eflag_global_in, double &energy_in,
                                               double *eatom_in)
{
  if (slabflag != 1) {
    PPPMElectrode::compute_boundary_corr(qsum_in, eflag_atom_in, eflag_global_in,
                                          energy_in, eatom_in);
    return;
  }

  const double zprd_slab = domain->zprd * slab_volfactor;
  const int nlocal = atom->nlocal;
  int *type = atom->type;
  double *q = atom->q;
  double **x = atom->x;
  double xM[3];
  double *xi;
  int iH1, iH2;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else {
      xi = x[i];
    }
    dipole += q[i] * xi[2];
  }

  double dipole_all;
  MPI_Allreduce(&dipole, &dipole_all, 1, MPI_DOUBLE, MPI_SUM, world);

  double dipole_r2 = 0.0;
  constexpr double SMALL = 0.00001;
  if (eflag_atom_in || fabs(qsum_in) > SMALL) {
    for (int i = 0; i < nlocal; i++) {
      if (type[i] == typeO) {
        find_M(i, iH1, iH2, xM);
        xi = xM;
      } else {
        xi = x[i];
      }
      dipole_r2 += q[i] * xi[2] * xi[2];
    }

    double tmp;
    MPI_Allreduce(&dipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    dipole_r2 = tmp;
  }

  const double e_slabcorr =
      MY_2PI * (dipole_all * dipole_all - qsum_in * dipole_r2 -
                qsum_in * qsum_in * zprd_slab * zprd_slab / 12.0) /
      volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global_in) energy_in += qscale * e_slabcorr;

  if (eflag_atom_in) {
    const double efact = qscale * MY_2PI / volume;
    for (int i = 0; i < nlocal; i++) {
      if (type[i] == typeO) {
        find_M(i, iH1, iH2, xM);
        const double e_pa =
            efact * q[i] *
            (xM[2] * dipole_all -
             0.5 * (dipole_r2 + qsum_in * xM[2] * xM[2]) -
             qsum_in * zprd_slab * zprd_slab / 12.0);
        eatom_in[i] += e_pa * (1.0 - alpha);
        eatom_in[iH1] += e_pa * alpha * 0.5;
        eatom_in[iH2] += e_pa * alpha * 0.5;
      } else {
        eatom_in[i] +=
            efact * q[i] *
            (x[i][2] * dipole_all -
             0.5 * (dipole_r2 + qsum_in * x[i][2] * x[i][2]) -
             qsum_in * zprd_slab * zprd_slab / 12.0);
      }
    }
  }

  const double ffact = qscale * (-4.0 * MY_PI / volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      const double fzi_corr = ffact * q[i] * (dipole_all - qsum_in * xM[2]);
      f[i][2] += fzi_corr * (1.0 - alpha);
      f[iH1][2] += 0.5 * alpha * fzi_corr;
      f[iH2][2] += 0.5 * alpha * fzi_corr;
    } else {
      f[i][2] += ffact * q[i] * (dipole_all - qsum_in * x[i][2]);
    }
  }
}

/* ----------------------------------------------------------------------
   TIP4P-aware boundary vector correction
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::compute_vector_boundary_corr(double *vec, int sensor_grpbit,
                                                      int source_grpbit, bool invert_source)
{
  if (slabflag != 1) {
    PPPMElectrode::compute_vector_boundary_corr(vec, sensor_grpbit, source_grpbit,
                                                invert_source);
    return;
  }

  const int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double *q = atom->q;
  double **x = atom->x;
  double xM[3];
  int iH1, iH2;
  double dipole = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!!(mask[i] & source_grpbit) != invert_source) {
      if (type[i] == typeO) {
        find_M(i, iH1, iH2, xM);
        dipole += q[i] * xM[2];
      } else {
        dipole += q[i] * x[i][2];
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &dipole, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole *= 4.0 * MY_PI / volume;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & sensor_grpbit) vec[i] += x[i][2] * dipole;
}


void PPPMElectrodeTIP4P::init_tip4p()
{
  int itmp = 0;

  // if kspace is TIP4P, extract TIP4P params from pair style
  // bond/angle are not yet init(), so ensure equilibrium request is valid

  qdist = 0.0;

  if (tip4pflag && wireflag)
    error->all(FLERR, "Kspace style pppm/electrode/tip4p does not support kspace_modify wire");

  if (tip4pflag) {

    auto *p_qdist = (double *) force->pair->extract("qdist", itmp);
    int *p_typeO = (int *) force->pair->extract("typeO", itmp);
    int *p_typeH = (int *) force->pair->extract("typeH", itmp);
    int *p_typeA = (int *) force->pair->extract("typeA", itmp);
    int *p_typeB = (int *) force->pair->extract("typeB", itmp);
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all(FLERR, "Pair style is incompatible with TIP4P KSpace style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == nullptr || force->bond == nullptr || force->angle->setflag == nullptr ||
        force->bond->setflag == nullptr)
      error->all(FLERR, "Bond and angle potentials must be defined for TIP4P");
    if (typeA < 1 || typeA > atom->nangletypes || force->angle->setflag[typeA] == 0)
      error->all(FLERR, "Bad TIP4P angle type for PPPM/TIP4P");
    if (typeB < 1 || typeB > atom->nbondtypes || force->bond->setflag[typeB] == 0)
      error->all(FLERR, "Bad TIP4P bond type for PPPM/TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5 * theta) * blen);
  }
}
/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the fictitious TIP4P charge site (M)
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::find_M(int i, int &iH1, int &iH2, double *xM)
{
  double **x = atom->x;

  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1)
    error->one(FLERR, "TIP4P hydrogen is missing");

  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

  iH1 = domain->closest_image(i,iH1);
  iH2 = domain->closest_image(i,iH2);

  // compute the fictitious M-site position from the O-H geometry

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];

  xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
}

/* ----------------------------------------------------------------------
   TIP4P-aware particle mapping
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::particle_map()
{

  int nx, ny, nz;
  int iH1, iH2;
  double *xi;
  double xM[3];

  int *type = atom->type;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable" + utils::errorurl(6));

  int flag = 0;

  for (int i = 0; i < nlocal; i++) {

    // use the fictitious M-site position for TIP4P oxygen atoms
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);
      xi = xM;
    } else {
      xi = x[i];
    }

    // determine the lower-left PPPM grid point for the interpolation stencil
    nx = static_cast<int>((xi[0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int>((xi[1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int>((xi[2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // verify that the complete PPPM interpolation stencil fits inside the local grid
    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag)
    error->one(FLERR, Error::NOLASTLINE,
               "Out of range atoms - cannot compute PPPM" + utils::errorurl(4));
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
-------------------------------------------------------------------------
*/


/* ----------------------------------------------------------------------
   compute estimated kspace force error
-------------------------------------------------------------------------
*/

/* ----------------------------------------------------------------------
   compute qopt
-------------------------------------------------------------------------
*/

/* ----------------------------------------------------------------------
   set local subset of PPPM/FFT grid that I own
   n xyz lo/hi in = 3d brick that I own (inclusive)
   n xyz lo/hi out = 3d brick + ghost cells in 6 directions (inclusive)
   n xyz lo/hi fft = FFT columns that I own (all of x dim, 2d decomp in yz)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   create discretized "density" of a set of particles (c.f. make_rho())
   in a specified scratch space.
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including
ghosts) in global grid
-------------------------------------------------------------------------
*/

/* ----------------------------------------------------------------------
   TIP4P-aware charge assignment
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::make_rho()
{

  int i, l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;

  FFT_SCALAR *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = ZEROF;

  int *type = atom->type;
  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int iH1, iH2;
  double *xi;
  double xM[3];

  for (i = 0; i < nlocal; i++) {

    // assign the TIP4P oxygen charge at the fictitious M-site position
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else {
      xi = x[i];
    }

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];

    dx = nx + shiftone - (xi[0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (xi[1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (xi[2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);

    // distribute the particle charge using the PPPM interpolation stencil
    z0 = delvolinv * q[i];

    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      y0 = z0 * rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        x0 = y0 * rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          density_brick[mz][my][mx] += x0 * rho1d[0][l];
        }
      }
    }
  }
}
/* ----------------------------------------------------------------------
   TIP4P-aware field interpolation (ik differentiation)
------------------------------------------------------------------------- */

void PPPMElectrodeTIP4P::fieldforce_ik()
{
  int i, l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;
  FFT_SCALAR ekx, eky, ekz;
  double *xi;
  int iH1, iH2;
  double xM[3];
  double fx, fy, fz;

  // interpolate electric field from nearby grid points

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    // evaluate per-atom energy/virial at the TIP4P M-site position
    // and redistribute the contribution to O and H atoms below
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else {
      xi = x[i];
    }

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];

    dx = nx + shiftone - (xi[0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (xi[1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (xi[2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);

    // interpolate the PPPM potential and virial quantities from the grid
    ekx = eky = ekz = ZEROF;

    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      z0 = rho1d[2][n];

      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        y0 = z0 * rho1d[1][m];

        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          x0 = y0 * rho1d[0][l];

          ekx -= x0 * vdx_brick[mz][my][mx];
          eky -= x0 * vdy_brick[mz][my][mx];
          ekz -= x0 * vdz_brick[mz][my][mx];
        }
      }
    }

    // convert the interpolated electric field into real-space forces
    const double qfactor = qqrd2e * scale * q[i];

    if (type[i] != typeO) {

      f[i][0] += qfactor * ekx;
      f[i][1] += qfactor * eky;

      if (slabflag != 2) f[i][2] += qfactor * ekz;

      // distribute the M-site contribution according to the TIP4P geometry
    } else {

      // redistribute the M-site force to the oxygen and hydrogen atoms
      fx = qfactor * ekx;
      fy = qfactor * eky;
      fz = qfactor * ekz;

      f[i][0] += fx * (1.0 - alpha);
      f[i][1] += fy * (1.0 - alpha);

      if (slabflag != 2) f[i][2] += fz * (1.0 - alpha);

      f[iH1][0] += 0.5 * alpha * fx;
      f[iH1][1] += 0.5 * alpha * fy;

      if (slabflag != 2) f[iH1][2] += 0.5 * alpha * fz;

      f[iH2][0] += 0.5 * alpha * fx;
      f[iH2][1] += 0.5 * alpha * fy;

      if (slabflag != 2) f[iH2][2] += 0.5 * alpha * fz;
    }
  }
}
void PPPMElectrodeTIP4P::fieldforce_ad()
{
  int i, l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz;
  FFT_SCALAR ekx, eky, ekz;
  double *xi;
  int iH1, iH2;
  double xM[3];
  double s1, s2, s3;
  double sf;
  double fx, fy, fz;

  double *prd = domain->prd;
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  double hx_inv = nx_pppm / xprd;
  double hy_inv = ny_pppm / yprd;
  double hz_inv = nz_pppm / zprd;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    // use the fictitious M-site position for TIP4P oxygen atoms
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else
      xi = x[i];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx + shiftone - (xi[0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (xi[1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (xi[2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);
    compute_drho1d(dx, dy, dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          ekx += drho1d[0][l] * rho1d[1][m] * rho1d[2][n] * u_brick[mz][my][mx];
          eky += rho1d[0][l] * drho1d[1][m] * rho1d[2][n] * u_brick[mz][my][mx];
          ekz += rho1d[0][l] * rho1d[1][m] * drho1d[2][n] * u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force and subtract self forces

    const double qfactor = qqrd2e * scale;

    s1 = xi[0] * hx_inv;
    s2 = xi[1] * hy_inv;
    s3 = xi[2] * hz_inv;
    sf = sf_coeff[0] * sin(2 * MY_PI * s1);
    sf += sf_coeff[1] * sin(4 * MY_PI * s1);
    sf *= 2.0 * q[i] * q[i];
    fx = qfactor * (ekx * q[i] - sf);

    sf = sf_coeff[2] * sin(2 * MY_PI * s2);
    sf += sf_coeff[3] * sin(4 * MY_PI * s2);
    sf *= 2.0 * q[i] * q[i];
    fy = qfactor * (eky * q[i] - sf);

    sf = sf_coeff[4] * sin(2 * MY_PI * s3);
    sf += sf_coeff[5] * sin(4 * MY_PI * s3);
    sf *= 2.0 * q[i] * q[i];
    fz = qfactor * (ekz * q[i] - sf);

    if (type[i] != typeO) {
      f[i][0] += fx;
      f[i][1] += fy;
      if (slabflag != 2) f[i][2] += fz;

    } else {

      // redistribute the M-site force to the oxygen and hydrogen atoms
      f[i][0] += fx * (1 - alpha);
      f[i][1] += fy * (1 - alpha);
      if (slabflag != 2) f[i][2] += fz * (1 - alpha);

      f[iH1][0] += 0.5 * alpha * fx;
      f[iH1][1] += 0.5 * alpha * fy;
      if (slabflag != 2) f[iH1][2] += 0.5 * alpha * fz;

      f[iH2][0] += 0.5 * alpha * fx;
      f[iH2][1] += 0.5 * alpha * fy;
      if (slabflag != 2) f[iH2][2] += 0.5 * alpha * fz;
    }
  }
}

void PPPMElectrodeTIP4P::fieldforce_peratom()
{
  int i, l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;
  double *xi;
  int iH1, iH2;
  double xM[3];
  FFT_SCALAR u_pa, v0, v1, v2, v3, v4, v5;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;

  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    // use the fictitious M-site position for TIP4P oxygen atoms
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else
      xi = x[i];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx + shiftone - (xi[0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (xi[1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (xi[2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);

    // interpolate per-atom energy and virial contributions from the PPPM grid
    u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        y0 = z0 * rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          x0 = y0 * rho1d[0][l];
          if (eflag_atom) u_pa += x0 * u_brick[mz][my][mx];
          if (vflag_atom) {
            v0 += x0 * v0_brick[mz][my][mx];
            v1 += x0 * v1_brick[mz][my][mx];
            v2 += x0 * v2_brick[mz][my][mx];
            v3 += x0 * v3_brick[mz][my][mx];
            v4 += x0 * v4_brick[mz][my][mx];
            v5 += x0 * v5_brick[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) {
      if (type[i] != typeO) {
        eatom[i] += q[i] * u_pa;
      } else {
        // redistribute the M-site contribution to the oxygen and hydrogen atoms
        eatom[i] += q[i] * u_pa * (1 - alpha);
        eatom[iH1] += q[i] * u_pa * alpha * 0.5;
        eatom[iH2] += q[i] * u_pa * alpha * 0.5;
      }
    }
    if (vflag_atom) {
      if (type[i] != typeO) {
        vatom[i][0] += v0 * q[i];
        vatom[i][1] += v1 * q[i];
        vatom[i][2] += v2 * q[i];
        vatom[i][3] += v3 * q[i];
        vatom[i][4] += v4 * q[i];
        vatom[i][5] += v5 * q[i];
      } else {
        vatom[i][0] += v0 * (1 - alpha) * q[i];
        vatom[i][1] += v1 * (1 - alpha) * q[i];
        vatom[i][2] += v2 * (1 - alpha) * q[i];
        vatom[i][3] += v3 * (1 - alpha) * q[i];
        vatom[i][4] += v4 * (1 - alpha) * q[i];
        vatom[i][5] += v5 * (1 - alpha) * q[i];
        vatom[iH1][0] += v0 * alpha * 0.5 * q[i];
        vatom[iH1][1] += v1 * alpha * 0.5 * q[i];
        vatom[iH1][2] += v2 * alpha * 0.5 * q[i];
        vatom[iH1][3] += v3 * alpha * 0.5 * q[i];
        vatom[iH1][4] += v4 * alpha * 0.5 * q[i];
        vatom[iH1][5] += v5 * alpha * 0.5 * q[i];
        vatom[iH2][0] += v0 * alpha * 0.5 * q[i];
        vatom[iH2][1] += v1 * alpha * 0.5 * q[i];
        vatom[iH2][2] += v2 * alpha * 0.5 * q[i];
        vatom[iH2][3] += v3 * alpha * 0.5 * q[i];
        vatom[iH2][4] += v4 * alpha * 0.5 * q[i];
        vatom[iH2][5] += v5 * alpha * 0.5 * q[i];
      }
    }
  }
}

void PPPMElectrodeTIP4P::make_rho_in_brick(int source_grpbit, FFT_SCALAR ***scratch_brick,
                                      bool invert_source)
{
  int l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;


  // clear 3d density array
  memset(&(scratch_brick[nzlo_out][nylo_out][nxlo_out]), 0, ngrid * sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int iH1, iH2;
  double *xi;
  double xM[3];

  for (int i = 0; i < nlocal; i++) {
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (!i_in_source) continue;

    // use the fictitious M-site position for TIP4P oxygen atoms
    if (type[i] == typeO) {
      find_M(i, iH1, iH2, xM);
      xi = xM;
    } else {
      xi = x[i];
    }

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];

    dx = nx + shiftone - (xi[0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (xi[1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (xi[2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);

    // spread the particle charge onto the surrounding PPPM interpolation stencil
    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      y0 = z0 * rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        x0 = y0 * rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          scratch_brick[mz][my][mx] += x0 * rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   group-group interactions
 -------------------------------------------------------------------------
*/
/* ----------------------------------------------------------------------
   compute b-vector EW3DC correction of constant potential approach
 -------------------------------------------------------------------------
*/
