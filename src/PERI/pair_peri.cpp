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

#include "pair_peri.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_peri_neigh.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPeri::PairPeri(LAMMPS *_lmp) :
    Pair(_lmp), fix_peri_neigh(nullptr), bulkmodulus(nullptr), shearmodulus(nullptr),
    m_lambdai(nullptr), m_taubi(nullptr), m_yieldstress(nullptr), s00(nullptr), alpha(nullptr),
    cut(nullptr), kspring(nullptr), s0_new(nullptr), theta(nullptr), elastic_energy(nullptr)
{
  for (int i = 0; i < 6; i++) virial[i] = 0.0;
  no_virial_fdotr_compute = 1;
  single_enable = 0;
  nmax = -1;
}

/* ---------------------------------------------------------------------- */

PairPeri::~PairPeri()
{
  if (fix_peri_neigh) modify->delete_fix(fix_peri_neigh->id);

  if (allocated) {
    memory->destroy(bulkmodulus);
    memory->destroy(shearmodulus);
    memory->destroy(m_lambdai);
    memory->destroy(m_taubi);
    memory->destroy(m_yieldstress);
    memory->destroy(s00);
    memory->destroy(alpha);
    memory->destroy(cut);
    memory->destroy(cutsq);
    memory->destroy(setflag);
    memory->destroy(kspring);

    memory->destroy(s0_new);
    memory->destroy(theta);
    memory->destroy(elastic_energy);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPeri::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(bulkmodulus, n, n, "pair:bulkmodulus");
  memory->create(shearmodulus, n, n, "pair:shearmodulus");
  memory->create(s00, n, n, "pair:s00");
  memory->create(alpha, n, n, "pair:alpha");
  memory->create(cut, n, n, "pair:cut");
  memory->create(m_yieldstress, n, n, "pair:m_yieldstress");
  memory->create(m_lambdai, n, n, "pair:m_lambdai");
  memory->create(m_taubi, n, n, "pair:m_taubi");
  memory->create(kspring, n, n, "pair:m_taubi");
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PairPeri::memory_usage()
{
  double bytes = 2.0 * (double) nmax * sizeof(double);
  bytes += 10.0 * (double) atom->ntypes * atom->ntypes * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPeri::settings(int narg, char ** /*arg*/)
{
  if (narg) error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   init common to all peridynamics pair styles
------------------------------------------------------------------------- */

void PairPeri::init_style()
{
  // error checks

  if (!atom->peri_flag) error->all(FLERR, "Pair style peri requires atom style peri");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Pair peri requires an atom map, see atom_modify");

  if (domain->lattice == nullptr) error->all(FLERR, "Pair peri requires a lattice be defined");
  if (domain->lattice->xlattice != domain->lattice->ylattice ||
      domain->lattice->xlattice != domain->lattice->zlattice ||
      domain->lattice->ylattice != domain->lattice->zlattice)
    error->all(FLERR, "Pair peri lattice is not identical in x, y, and z");

  // if first init, create Fix needed for storing fixed neighbors

  if (!fix_peri_neigh)
    fix_peri_neigh = dynamic_cast<FixPeriNeigh *>(modify->add_fix("PERI_NEIGH all PERI_NEIGH"));

  neighbor->add_request(this);
}

/* ---------------------------------------------------------------------- */

void PairPeri::compute_dilatation(int ifrom, int ito)
{
  int i, j, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double xtmp0, ytmp0, ztmp0, delx0, dely0, delz0;
  double rsq, r, dr;
  double delta;

  double **x = atom->x;
  int *type = atom->type;
  double **x0 = atom->x0;
  double *vfrac = atom->vfrac;
  double vfrac_scale = 1.0;

  double lc = domain->lattice->xlattice;
  double half_lc = 0.5 * lc;

  double **r0 = fix_peri_neigh->r0;
  tagint **partner = fix_peri_neigh->partner;
  int *npartner = fix_peri_neigh->npartner;
  double *wvolume = fix_peri_neigh->wvolume;

  int periodic = domain->xperiodic || domain->yperiodic || domain->zperiodic;

  // compute the dilatation theta

  for (i = ifrom; i < ito; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    jnum = npartner[i];
    theta[i] = 0.0;
    itype = type[i];

    for (jj = 0; jj < jnum; jj++) {

      // if bond already broken, skip this partner
      if (partner[i][jj] == 0) continue;

      // look up local index of this partner particle
      j = atom->map(partner[i][jj]);

      // skip if particle is "lost"
      if (j < 0) continue;

      // compute force density and add to PD equation of motion
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (periodic) domain->minimum_image(delx, dely, delz);
      rsq = delx * delx + dely * dely + delz * delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0, dely0, delz0);

      r = sqrt(rsq);
      dr = r - r0[i][jj];
      if (fabs(dr) < NEAR_ZERO) dr = 0.0;

      jtype = type[j];
      delta = cut[itype][jtype];

      // scale vfrac[j] if particle j near the horizon

      if ((fabs(r0[i][jj] - delta)) <= half_lc)
        vfrac_scale =
            (-1.0 / (2 * half_lc)) * (r0[i][jj]) + (1.0 + ((delta - half_lc) / (2 * half_lc)));
      else
        vfrac_scale = 1.0;

      theta[i] += influence_function(delx0, dely0, delz0) * r0[i][jj] * dr * vfrac[j] * vfrac_scale;
    }

    // if wvolume[i] is zero, then particle i has no bonds
    // therefore, the dilatation is set to zero

    if (wvolume[i] != 0.0)
      theta[i] = (3.0 / wvolume[i]) * theta[i];
    else
      theta[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   communication routines
---------------------------------------------------------------------- */

int PairPeri::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = theta[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairPeri::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) { theta[i] = buf[m++]; }
}

/* ---------------------------------------------------------------------- */

void *PairPeri::extract(const char *name, int &dim)
{
  dim = 1;
  if (strcmp(name, "theta") == 0) return (void *) theta;
  if (strcmp(name, "elastic_energy") == 0) return (void *) elastic_energy;
  return nullptr;
}
