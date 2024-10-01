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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "fix_rheo_oxidation.h"

#include "atom.h"
#include "atom_vec.h"
#include "citeme.h"
#include "comm.h"
#include "compute_rheo_surface.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
enum { NONE, CONSTANT };

static const char cite_rheo_oxide[] =
    "@article{ApplMathModel.130.310,\n"
    " title = {A hybrid smoothed-particle hydrodynamics model of oxide skins on molten aluminum},\n"
    " journal = {Applied Mathematical Modelling},\n"
    " volume = {130},\n"
    " pages = {310-326},\n"
    " year = {2024},\n"
    " issn = {0307-904X},\n"
    " doi = {https://doi.org/10.1016/j.apm.2024.02.027},\n"
    " author = {Joel T. Clemmer and Flint Pierce and Thomas C. O'Connor and Thomas D. Nevins and "
    "Elizabeth M.C. Jones and Jeremy B. Lechman and John Tencer},\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixRHEOOxidation::FixRHEOOxidation(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), compute_surface(nullptr), fix_rheo(nullptr)
{
  if (narg != 6) error->all(FLERR, "Illegal fix rheo/oxidation command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  comm_forward = 3;

  cut = utils::numeric(FLERR, arg[3], false, lmp);
  if (cut <= 0.0) error->all(FLERR, "Illegal bond cutoff {} in fix rheo/oxidation", cut);

  btype = utils::inumeric(FLERR, arg[4], false, lmp);
  if (btype < 1 || btype > atom->nbondtypes)
    error->all(FLERR, "Illegal value {} for bond type in fix rheo/oxidation", btype);

  rsurf = utils::numeric(FLERR, arg[5], false, lmp);
  if (rsurf <= 0.0) error->all(FLERR, "Illegal surface distance {} in fix rheo/oxidation", rsurf);

  cutsq = cut * cut;

  if (lmp->citeme) lmp->citeme->add(cite_rheo_oxide);
}

/* ---------------------------------------------------------------------- */

FixRHEOOxidation::~FixRHEOOxidation() {}

/* ---------------------------------------------------------------------- */

int FixRHEOOxidation::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/oxidation");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (cut > fix_rheo->cut) error->all(FLERR, "Bonding length exceeds kernel cutoff");

  if (rsurf >= fix_rheo->cut) error->all(FLERR, "Surface distance must be less than kernel cutoff");

  if (!force->bond) error->all(FLERR, "Must define a bond style with fix rheo/oxidation");
  if (!atom->avec->bonds_allow) error->all(FLERR, "Fix rheo/oxidation requires atom bonds");

  int tmp1, tmp2;
  index_nb = atom->find_custom("shell_nbond", tmp1, tmp2);
  if (index_nb == -1) error->all(FLERR, "Must use bond style rheo/shell to use fix rheo/oxidation");
  nbond = atom->ivector[index_nb];

  // need a half neighbor list
  auto req = neighbor->add_request(this, NeighConst::REQ_FULL);
  req->set_cutoff(cut);
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::setup_pre_force(int /*vflag*/)
{
  // Not strictly required that this fix be after FixRHEO,
  // but enforce to be consistent with other RHEO fixes
  fix_rheo->oxidation_fix_defined = 1;

  if (!fix_rheo->surface_flag)
    error->all(FLERR, "fix rheo/oxidation requires surface calculation in fix rheo");
  compute_surface = fix_rheo->compute_surface;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::pre_force(int /*vflag*/) {}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::post_integrate()
{
  int i, j, n, ii, jj, inum, jnum, bflag, fluidi, fluidj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double delx, dely, delz, rsq;
  tagint tagi, tagj;

  int newton_bond = force->newton_bond;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  int *mask = atom->mask;
  int *status = atom->rheo_status;
  double *rsurface = compute_surface->rsurface;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Forward positions (after inititial integrate, before comm)
  // Note: surface designation lags one timestep, acceptable error
  comm->forward_comm(this);

  int added_bonds = 0;
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    // Exclude particles that aren't solid or surface
    fluidi = !(status[i] & PHASECHECK);
    if (fluidi && (rsurface[i] > rsurf)) continue;

    tagi = tag[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;

      fluidj = !(status[j] & PHASECHECK);
      if (fluidj && (rsurface[j] > rsurf)) continue;

      // Skip solid-solid, leaves surface-surface or surface-solid
      if ((!fluidi) && (!fluidj)) continue;

      tagj = tag[j];

      // Ensure pair is always ordered to ensure numerical operations
      // are identical to minimize the possibility that a bond straddling
      // an mpi grid (newton off) isn't created on one proc but not the other
      if (tagi < tagj) {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
      } else {
        delx = x[j][0] - x[i][0];
        dely = x[j][1] - x[i][1];
        delz = x[j][2] - x[i][2];
      }
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq) continue;

      // Check if already have an oxide bond
      bflag = 0;
      for (n = 0; n < num_bond[i]; n++) {
        if (bond_type[i][n] == btype && bond_atom[i][n] == tagj) {
          bflag = 1;
          break;
        }
      }
      if (bflag) continue;

      added_bonds += 1;

      // Add bonds to owned atoms
      // If newton bond off, add to both, otherwise add to whichever has a smaller tag

      if (!newton_bond || (tagi < tagj)) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom in fix rheo/oxidation for atom {}",
                     tagi);
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tagj;
        num_bond[i]++;
      }
    }
  }

  int added_bonds_all;
  MPI_Allreduce(&added_bonds, &added_bonds_all, 1, MPI_INT, MPI_SUM, world);

  if (added_bonds_all > 0) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::post_force(int /*vflag*/)
{
  int *status = atom->rheo_status;
  int *num_bond = atom->num_bond;
  for (int i = 0; i < atom->nlocal; i++)
    if (num_bond[i] != 0) status[i] |= STATUS_NO_SHIFT;
}

/* ---------------------------------------------------------------------- */

int FixRHEOOxidation::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i, j, m;
  double **x = atom->x;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = x[j][0];
    buf[m++] = x[j][1];
    buf[m++] = x[j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  double **x = atom->x;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}
