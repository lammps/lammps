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

#include "fix_oxdna_lrf.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

FixOxdnaLRF::FixOxdnaLRF(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), nxyz(nullptr)
{
  comm_forward = 9;

  peratom_flag = 1;
  size_peratom_cols = 9;
  peratom_freq = 1;

  FixOxdnaLRF::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < size_peratom_cols; j++) nxyz[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixOxdnaLRF::~FixOxdnaLRF()
{
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(nxyz);
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::setmask()
{
  int mask = 0;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::min_setup_pre_force(int vflag)
{
  min_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::min_pre_force(int /*vflag*/)
{
  compute_lrf();
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::pre_force(int /*vflag*/)
{
  compute_lrf();
}

/* ---------------------------------------------------------------------- */

double FixOxdnaLRF::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * 9 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::grow_arrays(int nmax)
{
  memory->grow(nxyz, nmax, size_peratom_cols, "fix_oxdna/lrf:nxyz");
  array_atom = nxyz;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::copy_arrays(int i, int j, int /*delflag*/)
{
  memcpy(nxyz[j], nxyz[i], sizeof(double) * size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::set_arrays(int i)
{
  memset(nxyz[i], 0, sizeof(double) * size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = nxyz[i][0];
  buf[m++] = nxyz[i][1];
  buf[m++] = nxyz[i][2];
  buf[m++] = nxyz[i][3];
  buf[m++] = nxyz[i][4];
  buf[m++] = nxyz[i][5];
  buf[m++] = nxyz[i][6];
  buf[m++] = nxyz[i][7];
  buf[m++] = nxyz[i][8];

  return m;
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  nxyz[nlocal][0] = buf[m++];
  nxyz[nlocal][1] = buf[m++];
  nxyz[nlocal][2] = buf[m++];
  nxyz[nlocal][3] = buf[m++];
  nxyz[nlocal][4] = buf[m++];
  nxyz[nlocal][5] = buf[m++];
  nxyz[nlocal][6] = buf[m++];
  nxyz[nlocal][7] = buf[m++];
  nxyz[nlocal][8] = buf[m++];

  return m;
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = nxyz[j][0];
    buf[m++] = nxyz[j][1];
    buf[m++] = nxyz[j][2];
    buf[m++] = nxyz[j][3];
    buf[m++] = nxyz[j][4];
    buf[m++] = nxyz[j][5];
    buf[m++] = nxyz[j][6];
    buf[m++] = nxyz[j][7];
    buf[m++] = nxyz[j][8];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    nxyz[i][0] = buf[m++];
    nxyz[i][1] = buf[m++];
    nxyz[i][2] = buf[m++];
    nxyz[i][3] = buf[m++];
    nxyz[i][4] = buf[m++];
    nxyz[i][5] = buf[m++];
    nxyz[i][6] = buf[m++];
    nxyz[i][7] = buf[m++];
    nxyz[i][8] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::compute_lrf()
{
  using namespace MathExtra;    // q_to_exyz

  int nlocal = atom->nlocal;
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  // loop over all local atoms, calculation of local reference frame
  for (int i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit) {

      int n = ellipsoid[i];
      if (n < 0) continue;    // skip non-ellipsoid atoms

      // quaternion and Cartesian unit vectors in lab frame
      double *qn, nx_temp[3], ny_temp[3], nz_temp[3];

      qn = bonus[n].quat;

      q_to_exyz(qn, nx_temp, ny_temp, nz_temp);

      nxyz[i][0] = nx_temp[0];
      nxyz[i][1] = nx_temp[1];
      nxyz[i][2] = nx_temp[2];
      nxyz[i][3] = ny_temp[0];
      nxyz[i][4] = ny_temp[1];
      nxyz[i][5] = ny_temp[2];
      nxyz[i][6] = nz_temp[0];
      nxyz[i][7] = nz_temp[1];
      nxyz[i][8] = nz_temp[2];
    }
  }

  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */
