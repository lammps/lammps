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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  Support for bonds added by : Evangelos Voyiatzis (NovaMechanics)
  --------------------------------------------------------------------------*/

#include "compute_stress_mop_profile.h"

#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { X, Y, Z };
enum { LOWER, CENTER, UPPER, COORD };
enum { TOTAL, CONF, KIN, PAIR, BOND };

/* ---------------------------------------------------------------------- */

ComputeStressMopProfile::ComputeStressMopProfile(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute stress/mop/profile", error);

  bondflag = 0;

  // set compute mode and direction of plane(s) for pressure calculation

  if (strcmp(arg[3], "x") == 0) {
    dir = X;
  } else if (strcmp(arg[3], "y") == 0) {
    dir = Y;
  } else if (strcmp(arg[3], "z") == 0) {
    dir = Z;
  } else
    error->all(FLERR, "Illegal compute stress/mop/profile command");

  // bin parameters

  if (strcmp(arg[4], "lower") == 0)
    originflag = LOWER;
  else if (strcmp(arg[4], "center") == 0)
    originflag = CENTER;
  else if (strcmp(arg[4], "upper") == 0)
    originflag = UPPER;
  else
    originflag = COORD;
  if (originflag == COORD) origin = utils::numeric(FLERR, arg[4], false, lmp);
  delta = utils::numeric(FLERR, arg[5], false, lmp);
  invdelta = 1.0 / delta;

  // parse values until one isn't recognized

  which = new int[3 * (narg - 6)];
  nvalues = 0;
  int i;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "conf") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = CONF;
        nvalues++;
      }
    } else if (strcmp(arg[iarg], "kin") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = KIN;
        nvalues++;
      }
    } else if (strcmp(arg[iarg], "total") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = TOTAL;
        nvalues++;
      }
    } else if (strcmp(arg[iarg], "pair") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = PAIR;
        nvalues++;
      }
    } else if (strcmp(arg[iarg], "bond") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = BOND;
        nvalues++;
      }
    } else
      error->all(FLERR, "Illegal compute stress/mop/profile command");    //break;

    iarg++;
  }

  // check domain related errors

  // 3D only

  if (domain->dimension < 3)
    error->all(FLERR, "Compute stress/mop/profile incompatible with simulation dimension");

  // orthogonal simulation box

  if (domain->triclinic != 0)
    error->all(FLERR, "Compute stress/mop/profile incompatible with triclinic simulation box");

  // initialize some variables

  nbins = 0;
  coord = coordp = nullptr;
  values_local = values_global = array = nullptr;
  bond_local = nullptr;
  bond_global = nullptr;
  local_contribution = nullptr;

  // bin setup

  setup_bins();

  // this fix produces a global array

  memory->create(array, nbins, 1 + nvalues, "stress/mop/profile:array");
  size_array_rows = nbins;
  size_array_cols = 1 + nvalues;

  array_flag = 1;
  extarray = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStressMopProfile::~ComputeStressMopProfile()
{
  delete[] which;

  memory->destroy(coord);
  memory->destroy(coordp);
  memory->destroy(values_local);
  memory->destroy(values_global);
  memory->destroy(bond_local);
  memory->destroy(bond_global);
  memory->destroy(local_contribution);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeStressMopProfile::init()
{
  // conversion constants

  nktv2p = force->nktv2p;
  ftm2v = force->ftm2v;

  // plane area

  area = 1;
  int i;
  for (i = 0; i < 3; i++) {
    if (i != dir) area = area * domain->prd[i];
  }

  // timestep Value

  dt = update->dt;

  // Error checks:

  // Compute stress/mop/profile requires fixed simulation box

  if (domain->box_change_size || domain->box_change_shape || domain->deform_flag)
    error->all(FLERR, "Compute stress/mop/profile requires a fixed simulation box");

  // This compute requires a pair style with pair_single method implemented

  if (!force->pair) error->all(FLERR, "No pair style is defined for compute stress/mop/profile");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/mop/profile");

  // Errors

  if (comm->me == 0) {

    // Compute stress/mop/profile only accounts for pair interactions.
    // issue an error if any intramolecular potential or Kspace is defined.

    if (force->bond) bondflag = 1;

    if (force->angle)
      if ((strcmp(force->angle_style, "zero") != 0) && (strcmp(force->angle_style, "none") != 0))
        error->all(FLERR, "compute stress/mop/profile does not account for angle potentials");
    if (force->dihedral)
      if ((strcmp(force->dihedral_style, "zero") != 0) &&
          (strcmp(force->dihedral_style, "none") != 0))
        error->all(FLERR, "compute stress/mop/profile does not account for dihedral potentials");
    if (force->improper)
      if ((strcmp(force->improper_style, "zero") != 0) &&
          (strcmp(force->improper_style, "none") != 0))
        error->all(FLERR, "compute stress/mop/profile does not account for improper potentials");
    if (force->kspace)
      error->warning(FLERR, "compute stress/mop/profile does not account for kspace contributions");
  }

  // need an occasional half neighbor list

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeStressMopProfile::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   compute output array
   ------------------------------------------------------------------------- */

void ComputeStressMopProfile::compute_array()
{
  invoked_array = update->ntimestep;

  // Compute pressures on separate procs

  compute_pairs();

  // sum pressure contributions over all procs

  MPI_Allreduce(&values_local[0][0], &values_global[0][0], nbins * nvalues, MPI_DOUBLE, MPI_SUM,
                world);

  // Compute bond contribution on separate procs

  if (bondflag) {
    compute_bonds();
  } else {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++) { bond_local[m][i] = 0.0; }
    }
  }

  // sum bond contribution over all procs

  MPI_Allreduce(&bond_local[0][0], &bond_global[0][0], nbins * nvalues, MPI_DOUBLE, MPI_SUM, world);

  for (int ibin = 0; ibin < nbins; ibin++) {
    array[ibin][0] = coord[ibin][0];

    int mo = 1;
    int m = 0;
    while (m < nvalues) {
      array[ibin][m + mo] = values_global[ibin][m] + bond_global[ibin][m];
      m++;
    }
  }
}

/*------------------------------------------------------------------------
  compute pressure contribution of local proc
  -------------------------------------------------------------------------*/

void ComputeStressMopProfile::compute_pairs()

{
  int i, j, m, ii, jj, inum, jnum, itype, jtype, ibin;
  double delx, dely, delz;
  double rsq, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double pos, pos1;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // zero out arrays for one sample

  for (m = 0; m < nbins; m++) {
    for (i = 0; i < nvalues; i++) values_local[m][i] = 0.0;
  }

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  // parse values

  double xi[3];
  double vi[3];
  double fi[3];
  double xj[3];

  m = 0;
  while (m < nvalues) {
    if ((which[m] == CONF) || (which[m] == TOTAL) || (which[m] == PAIR)) {

      // Compute configurational contribution to pressure

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        xi[0] = atom->x[i][0];
        xi[1] = atom->x[i][1];
        xi[2] = atom->x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          factor_lj = special_lj[sbmask(j)];
          factor_coul = special_coul[sbmask(j)];
          j &= NEIGHMASK;

          // skip if neither I nor J are in group

          if (!(mask[i] & groupbit || mask[j] & groupbit)) continue;

          xj[0] = atom->x[j][0];
          xj[1] = atom->x[j][1];
          xj[2] = atom->x[j][2];
          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx * delx + dely * dely + delz * delz;
          jtype = type[j];
          if (rsq >= cutsq[itype][jtype]) continue;

          if (newton_pair || j < nlocal) {

            for (ibin = 0; ibin < nbins; ibin++) {
              pos = coord[ibin][0];
              pos1 = coordp[ibin][0];

              // check if ij pair is across plane, add contribution to pressure

              if (((xi[dir] > pos) && (xj[dir] < pos)) || ((xi[dir] > pos1) && (xj[dir] < pos1))) {

                pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                values_local[ibin][m] += fpair * (xi[0] - xj[0]) / area * nktv2p;
                values_local[ibin][m + 1] += fpair * (xi[1] - xj[1]) / area * nktv2p;
                values_local[ibin][m + 2] += fpair * (xi[2] - xj[2]) / area * nktv2p;

              } else if (((xi[dir] < pos) && (xj[dir] > pos)) ||
                         ((xi[dir] < pos1) && (xj[dir] > pos1))) {

                pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                values_local[ibin][m] -= fpair * (xi[0] - xj[0]) / area * nktv2p;
                values_local[ibin][m + 1] -= fpair * (xi[1] - xj[1]) / area * nktv2p;
                values_local[ibin][m + 2] -= fpair * (xi[2] - xj[2]) / area * nktv2p;
              }
            }
          } else {

            for (ibin = 0; ibin < nbins; ibin++) {
              pos = coord[ibin][0];
              pos1 = coordp[ibin][0];

              //check if ij pair is across plane, add contribution to pressure

              if (((xi[dir] > pos) && (xj[dir] < pos)) || ((xi[dir] > pos1) && (xj[dir] < pos1))) {

                pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                values_local[ibin][m] += fpair * (xi[0] - xj[0]) / area * nktv2p;
                values_local[ibin][m + 1] += fpair * (xi[1] - xj[1]) / area * nktv2p;
                values_local[ibin][m + 2] += fpair * (xi[2] - xj[2]) / area * nktv2p;
              }
            }
          }
        }
      }
    }

    // compute kinetic contribution to pressure
    // counts local particles transfers across the plane

    if ((which[m] == KIN) || (which[m] == TOTAL)) {

      double sgn;

      for (int i = 0; i < nlocal; i++) {

        // skip if I is not in group

        if (mask[i] & groupbit) {

          itype = type[i];

          // coordinates at t

          xi[0] = atom->x[i][0];
          xi[1] = atom->x[i][1];
          xi[2] = atom->x[i][2];

          // velocities at t

          vi[0] = atom->v[i][0];
          vi[1] = atom->v[i][1];
          vi[2] = atom->v[i][2];

          // forces at t

          fi[0] = atom->f[i][0];
          fi[1] = atom->f[i][1];
          fi[2] = atom->f[i][2];

          const double imass = (rmass) ? rmass[i] : mass[itype];
          const double iterm = 0.5 / imass * dt * ftm2v;

          // coordinates at t-dt (based on Velocity-Verlet alg.)

          xj[0] = xi[0] - vi[0] * dt + fi[0] * iterm * dt;
          xj[1] = xi[1] - vi[1] * dt + fi[1] * iterm * dt;
          xj[2] = xi[2] - vi[2] * dt + fi[2] * iterm * dt;

          for (ibin = 0; ibin < nbins; ibin++) {
            pos = coord[ibin][0];
            pos1 = coordp[ibin][0];

            if (((xi[dir] - pos) * (xj[dir] - pos) * (xi[dir] - pos1) * (xj[dir] - pos1) < 0)) {

              sgn = copysign(1.0, vi[dir]);

              // approximate crossing velocity by v(t-dt/2) (based on Velocity-Verlet alg.)

              double vcross[3];
              vcross[0] = vi[0] - fi[0] * iterm;
              vcross[1] = vi[1] - fi[1] * iterm;
              vcross[2] = vi[2] - fi[2] * iterm;

              values_local[ibin][m] += imass * vcross[0] * sgn / dt / area * nktv2p / ftm2v;
              values_local[ibin][m + 1] += imass * vcross[1] * sgn / dt / area * nktv2p / ftm2v;
              values_local[ibin][m + 2] += imass * vcross[2] * sgn / dt / area * nktv2p / ftm2v;
            }
          }
        }
      }
    }
    m += 3;
  }
}

/*------------------------------------------------------------------------
  compute bond contribution to pressure of local proc
  -------------------------------------------------------------------------*/

void ComputeStressMopProfile::compute_bonds()
{
  int i, nb, atom1, atom2, imol, iatom, btype;
  tagint tagprev;
  double rsq, fpair;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  int molecular = atom->molecular;

  Bond *bond = force->bond;

  double dx[3] = {0.0, 0.0, 0.0};
  double x_bond_1[3] = {0.0, 0.0, 0.0};
  double x_bond_2[3] = {0.0, 0.0, 0.0};

  // initialization

  for (int m = 0; m < nbins; m++) {
    for (int i = 0; i < nvalues; i++) { bond_local[m][i] = 0.0; }
    local_contribution[m][0] = 0.0;
    local_contribution[m][1] = 0.0;
    local_contribution[m][2] = 0.0;
  }

  // loop over all bonded atoms in the current proc

  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;

    if (molecular == 1)
      nb = num_bond[atom1];
    else {
      if (molindex[atom1] < 0) continue;
      imol = molindex[atom1];
      iatom = molatom[atom1];
      nb = onemols[imol]->num_bond[iatom];
    }

    for (i = 0; i < nb; i++) {
      if (molecular == 1) {
        btype = bond_type[atom1][i];
        atom2 = atom->map(bond_atom[atom1][i]);
      } else {
        tagprev = tag[atom1] - iatom - 1;
        btype = onemols[imol]->bond_type[iatom][i];
        atom2 = atom->map(onemols[imol]->bond_atom[iatom][i] + tagprev);
      }

      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype <= 0) continue;

      for (int ibin = 0; ibin < nbins; ibin++) {
        double pos = coord[ibin][0];

        // minimum image of atom1 with respect to the plane of interest

        dx[0] = x[atom1][0];
        dx[1] = x[atom1][1];
        dx[2] = x[atom1][2];
        dx[dir] -= pos;
        domain->minimum_image(dx[0], dx[1], dx[2]);
        x_bond_1[0] = dx[0];
        x_bond_1[1] = dx[1];
        x_bond_1[2] = dx[2];
        x_bond_1[dir] += pos;

        // minimum image of atom2 with respect to atom1

        dx[0] = x[atom2][0] - x_bond_1[0];
        dx[1] = x[atom2][1] - x_bond_1[1];
        dx[2] = x[atom2][2] - x_bond_1[2];
        domain->minimum_image(dx[0], dx[1], dx[2]);
        x_bond_2[0] = x_bond_1[0] + dx[0];
        x_bond_2[1] = x_bond_1[1] + dx[1];
        x_bond_2[2] = x_bond_1[2] + dx[2];

        // check if the bond vector crosses the plane of interest

        double tau = (x_bond_1[dir] - pos) / (x_bond_1[dir] - x_bond_2[dir]);
        if ((tau <= 1) && (tau >= 0)) {
          dx[0] = x_bond_1[0] - x_bond_2[0];
          dx[1] = x_bond_1[1] - x_bond_2[1];
          dx[2] = x_bond_1[2] - x_bond_2[2];
          rsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
          bond->single(btype, rsq, atom1, atom2, fpair);

          double sgn = copysign(1.0, x_bond_1[dir] - pos);
          local_contribution[ibin][0] += sgn * fpair * dx[0] / area * nktv2p;
          local_contribution[ibin][1] += sgn * fpair * dx[1] / area * nktv2p;
          local_contribution[ibin][2] += sgn * fpair * dx[2] / area * nktv2p;
        }
      }
    }
  }

  // loop over the keywords and if necessary add the bond contribution

  int m = 0;
  while (m < nvalues) {
    if ((which[m] == CONF) || (which[m] == TOTAL) || (which[m] == BOND)) {
      for (int ibin = 0; ibin < nbins; ibin++) {
        bond_local[ibin][m] = local_contribution[ibin][0];
        bond_local[ibin][m + 1] = local_contribution[ibin][1];
        bond_local[ibin][m + 2] = local_contribution[ibin][2];
      }
    }
    m += 3;
  }
}

/* ----------------------------------------------------------------------
   setup 1d bins and their extent and coordinates
   called at init()
   ------------------------------------------------------------------------- */

void ComputeStressMopProfile::setup_bins()
{
  int i, n;
  double lo = 0.0, hi = 0.0;

  double *boxlo, *boxhi;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  if (originflag == LOWER)
    origin = boxlo[dir];
  else if (originflag == UPPER)
    origin = boxhi[dir];
  else if (originflag == CENTER)
    origin = 0.5 * (boxlo[dir] + boxhi[dir]);

  if (origin < boxlo[dir]) {
    error->all(FLERR, "Origin of bins for compute stress/mop/profile is out of bounds");
  } else {
    n = static_cast<int>((origin - boxlo[dir]) * invdelta);
    lo = origin - n * delta;
  }
  if (origin < boxhi[dir]) {
    n = static_cast<int>((boxhi[dir] - origin) * invdelta);
    hi = origin + n * delta;
  } else {
    error->all(FLERR, "Origin of bins for compute stress/mop/profile is out of bounds");
  }

  offset = lo;
  nbins = static_cast<int>((hi - lo) * invdelta + 1.5);

  // allocate bin arrays

  memory->create(coord, nbins, 1, "stress/mop/profile:coord");
  memory->create(coordp, nbins, 1, "stress/mop/profile:coordp");
  memory->create(values_local, nbins, nvalues, "stress/mop/profile:values_local");
  memory->create(values_global, nbins, nvalues, "stress/mop/profile:values_global");
  memory->create(bond_local, nbins, nvalues, "stress/mop/profile:bond_local");
  memory->create(bond_global, nbins, nvalues, "stress/mop/profile:bond_global");
  memory->create(local_contribution, nbins, 3, "stress/mop/profile:local_contribution");

  // set bin coordinates

  for (i = 0; i < nbins; i++) {
    coord[i][0] = offset + i * delta;
    if (coord[i][0] < (domain->boxlo[dir] + domain->prd_half[dir])) {
      coordp[i][0] = coord[i][0] + domain->prd[dir];
    } else {
      coordp[i][0] = coord[i][0] - domain->prd[dir];
    }
  }
}
