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
  Support for bonds, angles and dihedrals added by : Evangelos Voyiatzis (NovaMechanics)
  --------------------------------------------------------------------------*/

#include "compute_stress_mop_profile.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
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

static constexpr double SMALL =     0.001;

enum { X, Y, Z };
enum { TOTAL, CONF, KIN, PAIR, BOND, ANGLE, DIHEDRAL };

/* ---------------------------------------------------------------------- */

ComputeStressMopProfile::ComputeStressMopProfile(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute stress/mop/profile", error);

  bondflag = 0;
  angleflag = 0;
  dihedralflag = 0;

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

  if (strcmp(arg[4], "lower") == 0) {
    origin = domain->boxlo[dir];
  } else if (strcmp(arg[4], "center") == 0) {
    origin = 0.5 * (domain->boxlo[dir] + domain->boxhi[dir]);
  } else if (strcmp(arg[4], "upper") == 0) {
    origin = domain->boxhi[dir];
  } else {
    origin = utils::numeric(FLERR, arg[4], false, lmp);
  }
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
    } else if (strcmp(arg[iarg], "angle") == 0) {
      for (i = 0; i < 3; i++) {
        which[nvalues] = ANGLE;
        nvalues++;
      }
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      for (i=0; i<3; i++) {
        which[nvalues] = DIHEDRAL;
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
  angle_local = nullptr;
  angle_global = nullptr;
  dihedral_local = nullptr;
  dihedral_global = nullptr;
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
  memory->destroy(angle_local);
  memory->destroy(angle_global);
  memory->destroy(dihedral_local);
  memory->destroy(dihedral_global);
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

    if (force->angle) {
      if (force->angle->born_matrix_enable == 0) {
        if ((strcmp(force->angle_style, "zero") != 0) && (strcmp(force->angle_style, "none") != 0))
          error->all(FLERR,"compute stress/mop/profile does not account for angle potentials");
      } else {
        angleflag = 1;
      }
    }

    if (force->dihedral) {
      if (force->dihedral->born_matrix_enable == 0) {
        if ((strcmp(force->dihedral_style, "zero") != 0) &&
            (strcmp(force->dihedral_style, "none") != 0))
          error->all(FLERR, "compute stress/mop/profile does not account for dihedral potentials");
      } else {
        dihedralflag = 1;
      }
    }

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

  if (angleflag) {
    //Compute angle contribution on separate procs
    compute_angles();
  } else {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++) {
        angle_local[m][i] = 0.0;
      }
    }
  }

  // sum angle contribution over all procs
  MPI_Allreduce(&angle_local[0][0],&angle_global[0][0],nbins*nvalues,MPI_DOUBLE,MPI_SUM,world);

  if (dihedralflag) {
    //Compute dihedral contribution on separate procs
    compute_dihedrals();
  } else {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++) {
        dihedral_local[m][i] = 0.0;
      }
    }
  }

  // sum dihedral contribution over all procs
  MPI_Allreduce(&dihedral_local[0][0],&dihedral_global[0][0],nbins*nvalues,MPI_DOUBLE,MPI_SUM,world);

  for (int ibin = 0; ibin < nbins; ibin++) {
    array[ibin][0] = coord[ibin];

    int mo = 1;
    int m = 0;
    while (m < nvalues) {
      array[ibin][m + mo] = values_global[ibin][m] + bond_global[ibin][m] + angle_global[ibin][m] + dihedral_global[ibin][m];
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
              pos = coord[ibin];
              pos1 = coordp[ibin];

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
              pos = coord[ibin];
              pos1 = coordp[ibin];

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
            pos = coord[ibin];
            pos1 = coordp[ibin];

            // minimum image of xi with respect to the plane
            xi[dir] -= pos;
            domain->minimum_image(xi[0], xi[1], xi[2]);
            xi[dir] += pos;

            // minimum image of xj with respect to xi
            xj[0] -= xi[0];
            xj[1] -= xi[1];
            xj[2] -= xi[2];
            domain->minimum_image(xi[0], xi[1], xi[2]);
            xj[0] += xi[0];
            xj[1] += xi[1];
            xj[2] += xi[2];

            double tau = (xi[dir] - pos) / (xi[dir] - xj[dir]);
            if ((tau <= 1) && (tau >= 0)) {

              sgn = copysign(1.0, vi[dir]);

              //approximate crossing velocity by v(t-dt/2) (based on Velocity-Verlet alg.)
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
        double pos = coord[ibin];

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

/*------------------------------------------------------------------------
  compute angle contribution to pressure of local proc
  -------------------------------------------------------------------------*/

void ComputeStressMopProfile::compute_angles()
{
  int na, atom1, atom2, atom3, imol, iatom, atype;
  tagint tagprev;
  double r1, r2, cos_theta;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int molecular = atom->molecular;

  // loop over all atoms and their angles
  Angle *angle = force->angle;

  double duang, du2ang;
  double dx[3] = {0.0, 0.0, 0.0};
  double dx_left[3] = {0.0, 0.0, 0.0};
  double dx_right[3] = {0.0, 0.0, 0.0};
  double x_angle_left[3] = {0.0, 0.0, 0.0};
  double x_angle_middle[3] = {0.0, 0.0, 0.0};
  double x_angle_right[3] = {0.0, 0.0, 0.0};
  double dcos_theta[3] = {0.0, 0.0, 0.0};

  // initialization
  for (int m = 0; m < nbins; m++) {
    for (int i = 0; i < nvalues; i++) {
      angle_local[m][i] = 0.0;
    }
    local_contribution[m][0] = 0.0;
    local_contribution[m][1] = 0.0;
    local_contribution[m][2] = 0.0;
  }


  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == 1)
      na = num_angle[atom2];
    else {
      if (molindex[atom2] < 0) continue;
      imol = molindex[atom2];
      iatom = molatom[atom2];
      na = onemols[imol]->num_angle[iatom];
    }

    for (int i = 0; i < na; i++) {
      if (molecular == 1) {
        if (tag[atom2] != angle_atom2[atom2][i]) continue;
        atype = angle_type[atom2][i];
        atom1 = atom->map(angle_atom1[atom2][i]);
        atom3 = atom->map(angle_atom3[atom2][i]);
      } else {
        if (tag[atom2] != onemols[imol]->angle_atom2[atom2][i]) continue;
        atype = onemols[imol]->angle_type[atom2][i];
        tagprev = tag[atom2] - iatom - 1;
        atom1 = atom->map(onemols[imol]->angle_atom1[atom2][i] + tagprev);
        atom3 = atom->map(onemols[imol]->angle_atom3[atom2][i] + tagprev);
      }

      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      if (atype <= 0) continue;

      for (int ibin = 0; ibin<nbins; ibin++) {
        double pos = coord[ibin];

        // minimum image of atom1 with respect to the plane of interest
        dx[0] = x[atom1][0];
        dx[1] = x[atom1][1];
        dx[2] = x[atom1][2];
        dx[dir] -= pos;
        domain->minimum_image(dx[0], dx[1], dx[2]);
        x_angle_left[0] = dx[0];
        x_angle_left[1] = dx[1];
        x_angle_left[2] = dx[2];
        x_angle_left[dir] += pos;

        // minimum image of atom2 with respect to atom1
        dx_left[0] = x[atom2][0] - x_angle_left[0];
        dx_left[1] = x[atom2][1] - x_angle_left[1];
        dx_left[2] = x[atom2][2] - x_angle_left[2];
        domain->minimum_image(dx_left[0], dx_left[1], dx_left[2]);
        x_angle_middle[0] = x_angle_left[0] + dx_left[0];
        x_angle_middle[1] = x_angle_left[1] + dx_left[1];
        x_angle_middle[2] = x_angle_left[2] + dx_left[2];

        // minimum image of atom3 with respect to atom2
        dx_right[0] = x[atom3][0] - x_angle_middle[0];
        dx_right[1] = x[atom3][1] - x_angle_middle[1];
        dx_right[2] = x[atom3][2] - x_angle_middle[2];
        domain->minimum_image(dx_right[0], dx_right[1], dx_right[2]);
        x_angle_right[0] = x_angle_middle[0] + dx_right[0];
        x_angle_right[1] = x_angle_middle[1] + dx_right[1];
        x_angle_right[2] = x_angle_middle[2] + dx_right[2];

        // check if any bond vector crosses the plane of interest
        double tau_right = (x_angle_right[dir] - pos) / (x_angle_right[dir] - x_angle_middle[dir]);
        double tau_left = (x_angle_middle[dir] - pos) / (x_angle_middle[dir] - x_angle_left[dir]);
        bool right_cross = ((tau_right >=0) && (tau_right  <= 1));
        bool left_cross = ((tau_left >=0) && (tau_left <= 1));

        // no bonds crossing the plane
        if (!right_cross && !left_cross) continue;

        // compute the cos(theta) of the angle
        r1 = sqrt(dx_left[0]*dx_left[0] + dx_left[1]*dx_left[1] + dx_left[2]*dx_left[2]);
        r2 = sqrt(dx_right[0]*dx_right[0] + dx_right[1]*dx_right[1] + dx_right[2]*dx_right[2]);
        cos_theta = -(dx_right[0]*dx_left[0] + dx_right[1]*dx_left[1] + dx_right[2]*dx_left[2])/(r1*r2);

        if (cos_theta >  1.0) cos_theta = 1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;

        // The method returns derivative with regards to cos(theta)
        angle->born_matrix(atype, atom1, atom2, atom3, duang, du2ang);
        // only right bond crossing the plane
        if (right_cross && !left_cross)
        {
          double sgn = copysign(1.0, x_angle_right[dir] - pos);
          dcos_theta[0] = sgn*(dx_right[0]*cos_theta/r2 + dx_left[0]/r1)/r2;
          dcos_theta[1] = sgn*(dx_right[1]*cos_theta/r2 + dx_left[1]/r1)/r2;
          dcos_theta[2] = sgn*(dx_right[2]*cos_theta/r2 + dx_left[2]/r1)/r2;
        }

        // only left bond crossing the plane
        if (!right_cross && left_cross)
        {
          double sgn = copysign(1.0, x_angle_left[dir] - pos);
          dcos_theta[0] = -sgn*(dx_left[0]*cos_theta/r1 + dx_right[0]/r2)/r1;
          dcos_theta[1] = -sgn*(dx_left[1]*cos_theta/r1 + dx_right[1]/r2)/r1;
          dcos_theta[2] = -sgn*(dx_left[2]*cos_theta/r1 + dx_right[2]/r2)/r1;
        }

        // both bonds crossing the plane
        if (right_cross && left_cross)
        {
          // due to right bond
          double sgn = copysign(1.0, x_angle_middle[dir] - pos);
          dcos_theta[0] = -sgn*(dx_right[0]*cos_theta/r2 + dx_left[0]/r1)/r2;
          dcos_theta[1] = -sgn*(dx_right[1]*cos_theta/r2 + dx_left[1]/r1)/r2;
          dcos_theta[2] = -sgn*(dx_right[2]*cos_theta/r2 + dx_left[2]/r1)/r2;

          // due to left bond
          dcos_theta[0] += sgn*(dx_left[0]*cos_theta/r1 + dx_right[0]/r2)/r1;
          dcos_theta[1] += sgn*(dx_left[1]*cos_theta/r1 + dx_right[1]/r2)/r1;
          dcos_theta[2] += sgn*(dx_left[2]*cos_theta/r1 + dx_right[2]/r2)/r1;
        }

        // final contribution of the given angle term
        local_contribution[ibin][0] += duang*dcos_theta[0]/area*nktv2p;
        local_contribution[ibin][1] += duang*dcos_theta[1]/area*nktv2p;
        local_contribution[ibin][2] += duang*dcos_theta[2]/area*nktv2p;
      }
    }
  }

  // loop over the keywords and if necessary add the angle contribution
  int m = 0;
  while (m < nvalues) {
    if (which[m] == CONF || which[m] == TOTAL || which[m] == ANGLE) {
      for (int ibin = 0; ibin < nbins; ibin++) {
        angle_local[ibin][m] = local_contribution[ibin][0];
        angle_local[ibin][m+1] = local_contribution[ibin][1];
        angle_local[ibin][m+2] = local_contribution[ibin][2];
      }
    }
    m += 3;
  }
}

/*------------------------------------------------------------------------
  compute dihedral contribution to pressure of local proc
  -------------------------------------------------------------------------*/

void ComputeStressMopProfile::compute_dihedrals()
{
  int i, nd, atom1, atom2, atom3, atom4, imol, iatom;
  tagint tagprev;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
  double vb2xm, vb2ym, vb2zm;
  double sb1, sb2, sb3, rb1, rb3, c0, b1mag2, b1mag, b2mag2;
  double b2mag, b3mag2, b3mag, c2mag, ctmp, r12c1, c1mag, r12c2;
  double s1, s2, s12, sc1, sc2, a11, a22, a33, a12, a13, a23;
  double df[3], f1[3], f2[3], f3[3], f4[3];
  double c, sx2, sy2, sz2, sin2;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int molecular = atom->molecular;

  // loop over all atoms and their dihedrals

  Dihedral *dihedral = force->dihedral;

  double dudih, du2dih;

  double diffx[3] = {0.0, 0.0, 0.0};
  double x_atom_1[3] = {0.0, 0.0, 0.0};
  double x_atom_2[3] = {0.0, 0.0, 0.0};
  double x_atom_3[3] = {0.0, 0.0, 0.0};
  double x_atom_4[3] = {0.0, 0.0, 0.0};

  // initialization
  for (int m = 0; m < nbins; m++) {
    for (int i = 0; i < nvalues; i++) {
      dihedral_local[m][i] = 0.0;
    }
    local_contribution[m][0] = 0.0;
    local_contribution[m][1] = 0.0;
    local_contribution[m][2] = 0.0;
  }

  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == Atom::MOLECULAR)
      nd = num_dihedral[atom2];
    else {
      if (molindex[atom2] < 0) continue;
      imol = molindex[atom2];
      iatom = molatom[atom2];
      nd = onemols[imol]->num_dihedral[iatom];
    }

    for (i = 0; i < nd; i++) {
      if (molecular == 1) {
        if (tag[atom2] != dihedral_atom2[atom2][i]) continue;
          atom1 = atom->map(dihedral_atom1[atom2][i]);
          atom3 = atom->map(dihedral_atom3[atom2][i]);
          atom4 = atom->map(dihedral_atom4[atom2][i]);
      } else {
        if (tag[atom2] != onemols[imol]->dihedral_atom2[atom2][i]) continue;
        tagprev = tag[atom2] - iatom - 1;
        atom1 = atom->map(onemols[imol]->dihedral_atom1[atom2][i] + tagprev);
        atom3 = atom->map(onemols[imol]->dihedral_atom3[atom2][i] + tagprev);
        atom4 = atom->map(onemols[imol]->dihedral_atom4[atom2][i] + tagprev);
      }

      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      for (int ibin = 0; ibin<nbins; ibin++) {
        double pos = coord[ibin];

        // minimum image of atom1 with respect to the plane of interest
        x_atom_1[0] = x[atom1][0];
        x_atom_1[1] = x[atom1][1];
        x_atom_1[2] = x[atom1][2];
        x_atom_1[dir] -= pos;
        domain->minimum_image(x_atom_1[0], x_atom_1[1], x_atom_1[2]);
        x_atom_1[dir] += pos;

        // minimum image of atom2 with respect to atom1
        diffx[0] = x[atom2][0] - x_atom_1[0];
        diffx[1] = x[atom2][1] - x_atom_1[1];
        diffx[2] = x[atom2][2] - x_atom_1[2];
        domain->minimum_image(diffx[0], diffx[1], diffx[2]);
        x_atom_2[0] = x_atom_1[0] + diffx[0];
        x_atom_2[1] = x_atom_1[1] + diffx[1];
        x_atom_2[2] = x_atom_1[2] + diffx[2];

        // minimum image of atom3 with respect to atom2
        diffx[0] = x[atom3][0] - x_atom_2[0];
        diffx[1] = x[atom3][1] - x_atom_2[1];
        diffx[2] = x[atom3][2] - x_atom_2[2];
        domain->minimum_image(diffx[0], diffx[1], diffx[2]);
        x_atom_3[0] = x_atom_2[0] + diffx[0];
        x_atom_3[1] = x_atom_2[1] + diffx[1];
        x_atom_3[2] = x_atom_2[2] + diffx[2];

        // minimum image of atom3 with respect to atom2
        diffx[0] = x[atom4][0] - x_atom_3[0];
        diffx[1] = x[atom4][1] - x_atom_3[1];
        diffx[2] = x[atom4][2] - x_atom_3[2];
        domain->minimum_image(diffx[0], diffx[1], diffx[2]);
        x_atom_4[0] = x_atom_3[0] + diffx[0];
        x_atom_4[1] = x_atom_3[1] + diffx[1];
        x_atom_4[2] = x_atom_3[2] + diffx[2];

        // check if any bond vector crosses the plane of interest
        double tau_right = (x_atom_2[dir] - pos) / (x_atom_2[dir] - x_atom_1[dir]);
        double tau_middle = (x_atom_3[dir] - pos) / (x_atom_3[dir] - x_atom_2[dir]);
        double tau_left = (x_atom_4[dir] - pos) / (x_atom_4[dir] - x_atom_3[dir]);
        bool right_cross = ((tau_right >=0) && (tau_right  <= 1));
        bool middle_cross = ((tau_middle >= 0) && (tau_middle <= 1));
        bool left_cross = ((tau_left >=0) && (tau_left <= 1));

        // no bonds crossing the plane
        if (!right_cross && !middle_cross && !left_cross) continue;

        dihedral->born_matrix(i, atom1, atom2, atom3, atom4, dudih, du2dih);

        // first bond
        vb1x = x_atom_1[0] - x_atom_2[0];
        vb1y = x_atom_1[1] - x_atom_2[1];
        vb1z = x_atom_1[2] - x_atom_2[2];

        // second bond
        vb2x = x_atom_3[0] - x_atom_2[0];
        vb2y = x_atom_3[1] - x_atom_2[1];
        vb2z = x_atom_3[2] - x_atom_2[2];

        vb2xm = -vb2x;
        vb2ym = -vb2y;
        vb2zm = -vb2z;

        // third bond
        vb3x = x_atom_4[0] - x_atom_3[0];
        vb3y = x_atom_4[1] - x_atom_3[1];
        vb3z = x_atom_4[2] - x_atom_3[2];

        // c0 calculation
        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

        rb1 = sqrt(sb1);
        rb3 = sqrt(sb3);

        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;
        // 1st and 2nd angle
        b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
        b1mag = sqrt(b1mag2);
        b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
        b2mag = sqrt(b2mag2);
        b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
        b3mag = sqrt(b3mag2);

        ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
        r12c1 = 1.0 / (b1mag*b2mag);
        c1mag = ctmp * r12c1;

        ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
        r12c2 = 1.0 / (b2mag*b3mag);
        c2mag = ctmp * r12c2;

        // cos and sin of 2 angles and final c
        sin2 = MAX(1.0 - c1mag*c1mag,0.0);
        sc1 = sqrt(sin2);
        if (sc1 < SMALL) sc1 = SMALL;
        sc1 = 1.0/sc1;

        sin2 = MAX(1.0 - c2mag*c2mag,0.0);
        sc2 = sqrt(sin2);
        if (sc2 < SMALL) sc2 = SMALL;
        sc2 = 1.0/sc2;

        s1 = sc1 * sc1;
        s2 = sc2 * sc2;
        s12 = sc1 * sc2;
        c = (c0 + c1mag*c2mag) * s12;

        // error check
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        // forces on each particle
        double a = dudih;
        c = c * a;
        s12 = s12 * a;
        a11 = c*sb1*s1;
        a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2));
        a33 = c*sb3*s2;
        a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
        a13 = -rb1*rb3*s12;
        a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);

        sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z;

        f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
        f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
        f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;

        f2[0] = -sx2 - f1[0];
        f2[1] = -sy2 - f1[1];
        f2[2] = -sz2 - f1[2];

        f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
        f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
        f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;

        f3[0] = sx2 - f4[0];
        f3[1] = sy2 - f4[1];
        f3[2] = sz2 - f4[2];

        // only right bond crossing the plane
        if (right_cross && !middle_cross && !left_cross)
        {
          double sgn = copysign(1.0, x_atom_1[dir] - pos);
          df[0] = sgn * f1[0];
          df[1] = sgn * f1[1];
          df[2] = sgn * f1[2];
        }

        // only middle bond crossing the plane
        if (!right_cross && middle_cross && !left_cross)
        {
          double sgn = copysign(1.0, x_atom_2[dir] - pos);
          df[0] = sgn * (f2[0] + f1[0]);
          df[1] = sgn * (f2[1] + f1[1]);
          df[2] = sgn * (f2[2] + f1[2]);
        }

        // only left bond crossing the plane
        if (!right_cross && !middle_cross && left_cross)
        {
          double sgn = copysign(1.0, x_atom_4[dir] - pos);
          df[0] = sgn * f4[0];
          df[1] = sgn * f4[1];
          df[2] = sgn * f4[2];
        }

        // only right & middle bonds crossing the plane
        if (right_cross && middle_cross && !left_cross)
        {
          double sgn = copysign(1.0, x_atom_2[dir] - pos);
          df[0] = sgn * f2[0];
          df[1] = sgn * f2[1];
          df[2] = sgn * f2[2];
        }

        // only right & left bonds crossing the plane
        if (right_cross && !middle_cross && left_cross)
        {
          double sgn = copysign(1.0, x_atom_1[dir] - pos);
          df[0] = sgn * (f1[0] + f4[0]);
          df[1] = sgn * (f1[1] + f4[1]);
          df[2] = sgn * (f1[2] + f4[2]);
        }

        // only middle & left bonds crossing the plane
        if (!right_cross && middle_cross && left_cross)
        {
          double sgn = copysign(1.0, x_atom_3[dir] - pos);
          df[0] = sgn * f3[0];
          df[1] = sgn * f3[1];
          df[2] = sgn * f3[2];
        }

        // all three bonds crossing the plane
        if (right_cross && middle_cross && left_cross)
        {
          double sgn = copysign(1.0, x_atom_1[dir] - pos);
          df[0] = sgn * (f1[0] + f3[0]);
          df[1] = sgn * (f1[1] + f3[1]);
          df[2] = sgn * (f1[2] + f3[2]);
        }

        local_contribution[ibin][0] += df[0]/area*nktv2p;
        local_contribution[ibin][1] += df[1]/area*nktv2p;
        local_contribution[ibin][2] += df[2]/area*nktv2p;
      }
    }
  }

  // loop over the keywords and if necessary add the dihedral contribution
  int m = 0;
  while (m < nvalues) {
    if ((which[m] == CONF) || (which[m] == TOTAL) || (which[m] == DIHEDRAL)) {
      for (int ibin = 0; ibin < nbins; ibin++) {
        dihedral_local[ibin][m] = local_contribution[ibin][0];
        dihedral_local[ibin][m+1] = local_contribution[ibin][1];
        dihedral_local[ibin][m+2] = local_contribution[ibin][2];
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

  if ((origin > domain->boxhi[dir]) || (origin < domain->boxlo[dir]))
    error->all(FLERR, "Origin of bins for compute stress/mop/profile is out of bounds");

  n = static_cast<int> ((origin - boxlo[dir]) * invdelta);
  lo = origin - n*delta;

  n = static_cast<int> ((boxhi[dir] - origin) * invdelta);
  hi = origin + n*delta;

  offset = lo;
  nbins = static_cast<int>((hi - lo) * invdelta + 1.5);

  //allocate bin arrays
  memory->create(coord, nbins, "stress/mop/profile:coord");
  memory->create(coordp, nbins, "stress/mop/profile:coordp");
  memory->create(values_local, nbins, nvalues, "stress/mop/profile:values_local");
  memory->create(values_global, nbins, nvalues, "stress/mop/profile:values_global");
  memory->create(bond_local, nbins, nvalues, "stress/mop/profile:bond_local");
  memory->create(bond_global, nbins, nvalues, "stress/mop/profile:bond_global");
  memory->create(angle_local, nbins, nvalues, "stress/mop/profile:angle_local");
  memory->create(angle_global, nbins, nvalues, "stress/mop/profile:angle_global");
  memory->create(dihedral_local,nbins,nvalues,"stress/mop/profile:dihedral_local");
  memory->create(dihedral_global,nbins,nvalues,"stress/mop/profile:dihedral_global");
  memory->create(local_contribution, nbins, 3, "stress/mop/profile:local_contribution");

  // set bin coordinates

  for (i = 0; i < nbins; i++) {
    coord[i] = offset + i * delta;
    if (coord[i] < (domain->boxlo[dir] + domain->prd_half[dir])) {
      coordp[i] = coord[i] + domain->prd[dir];
    } else {
      coordp[i] = coord[i] - domain->prd[dir];
    }
  }
}
