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

/* ----------------------------------------------------------------------
  Contributing Authors : Germain Clavier (TUe), Aidan Thompson (Sandia)
------------------------------------------------------------------------- */

#include "compute_born_matrix.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1000000000
#define SMALL 1e-16

// this table is used to pick the 3d rij vector indices used to
// compute the 6 indices long Voigt stress vector

static int constexpr sigma_albe[6][2] = {
    {0, 0},    // s11
    {1, 1},    // s22
    {2, 2},    // s33
    {1, 2},    // s44
    {0, 2},    // s55
    {0, 1},    // s66
};

// this table is used to pick the correct indices from the Voigt
// stress vector to compute the Cij matrix (21 terms, see doc) contribution

static int constexpr C_albe[21][2] = {
    {0, 0},    // C11
    {1, 1},    // C22
    {2, 2},    // C33
    {3, 3},    // C44
    {4, 4},    // C55
    {5, 5},    // C66
    {0, 1},    // C12
    {0, 2},    // C13
    {0, 3},    // C14
    {0, 4},    // C15
    {0, 5},    // C16
    {1, 2},    // C23
    {1, 3},    // C24
    {1, 4},    // C25
    {1, 5},    // C26
    {2, 3},    // C34
    {2, 4},    // C35
    {2, 5},    // C36
    {3, 4},    // C45
    {3, 5},    // C46
    {4, 5}     // C56
};

// this table is used to pick the 3d rij vector indices used to
// compute the 21 indices long Cij matrix

static int constexpr albemunu[21][4] = {
    {0, 0, 0, 0},    // C11
    {1, 1, 1, 1},    // C22
    {2, 2, 2, 2},    // C33
    {1, 2, 1, 2},    // C44
    {0, 2, 0, 2},    // C55
    {0, 1, 0, 1},    // C66
    {0, 0, 1, 1},    // C12
    {0, 0, 2, 2},    // C13
    {0, 0, 1, 2},    // C14
    {0, 0, 0, 2},    // C15
    {0, 0, 0, 1},    // C16
    {1, 1, 2, 2},    // C23
    {1, 1, 1, 2},    // C24
    {1, 1, 0, 2},    // C25
    {1, 1, 0, 1},    // C26
    {2, 2, 1, 2},    // C34
    {2, 2, 0, 2},    // C35
    {2, 2, 0, 1},    // C36
    {1, 2, 0, 2},    // C45
    {1, 2, 0, 1},    // C46
    {0, 1, 0, 2}     // C56
};

/* ---------------------------------------------------------------------- */

ComputeBornMatrix::ComputeBornMatrix(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), id_virial(nullptr), temp_x(nullptr), temp_f(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute born/matrix command");

  nvalues = 21;

  numflag = 0;
  numdelta = 0.0;

  pairflag = bondflag = angleflag = dihedflag = impflag = 0;
  if (narg == 3) {
    pairflag = bondflag = angleflag = dihedflag = impflag = 1;
  } else {
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "numdiff") == 0) {
        if (iarg + 3 > narg) error->all(FLERR, "Illegal compute born/matrix command");
        numflag = 1;
        numdelta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        if (numdelta <= 0.0) error->all(FLERR, "Illegal compute born/matrix command");
        id_virial = utils::strdup(arg[iarg + 2]);
        int icompute = modify->find_compute(id_virial);
        if (icompute < 0) error->all(FLERR, "Could not find compute born/matrix pressure ID");
        compute_virial = modify->compute[icompute];
        if (compute_virial->pressflag == 0)
          error->all(FLERR, "Compute born/matrix pressure ID does not compute pressure");
        iarg += 3;
      } else if (strcmp(arg[iarg], "pair") == 0) {
        pairflag = 1;
      } else if (strcmp(arg[iarg], "bond") == 0) {
        bondflag = 1;
      } else if (strcmp(arg[iarg], "angle") == 0) {
        angleflag = 1;
      } else if (strcmp(arg[iarg], "dihedral") == 0) {
        dihedflag = 1;
      } else if (strcmp(arg[iarg], "improper") == 0) {
        impflag = 1;
      } else {
        error->all(FLERR, "Illegal compute born/matrix command");
      }
      ++iarg;
    }
  }

  if (pairflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->pair) {
      if (force->pair->born_matrix_enable == 0)
        error->all(FLERR, "Pair style {} does not support compute born/matrix", force->pair_style);
    } else {
      pairflag = 0;
    }
  }

  if (bondflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->bond) {
      if (force->bond->born_matrix_enable == 0)
        error->all(FLERR, "Bond style {} does not support compute born/matrix", force->bond_style);
    } else {
      bondflag = 0;
    }
  }

  if (angleflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->angle) {
      if (force->angle->born_matrix_enable == 0)
        error->all(FLERR, "Angle style {} does not support compute born/matrix",
                   force->angle_style);
    } else {
      angleflag = 0;
    }
  }

  if (dihedflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->dihedral) {
      if (force->dihedral->born_matrix_enable == 0)
        error->all(FLERR, "Dihedral style {} does not support compute born/matrix",
                   force->dihedral_style);
    } else {
      dihedflag = 0;
    }
  }

  if (impflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->improper) {
      if (force->improper->born_matrix_enable == 0)
        error->all(FLERR, "Improper style {} does not support compute born/matrix",
                   force->improper_style);
    } else {
      impflag = 0;
    }
  }

  if (force->kspace) {
    if (!numflag && (comm->me == 0))
      error->all(FLERR, "KSpace contribution not supported by compute born/matrix");
  }

  // Initialize some variables

  values_local = values_global = vector = nullptr;

  // this fix produces a global vector

  memory->create(vector, nvalues, "born_matrix:vector");
  memory->create(values_global, nvalues, "born_matrix:values_global");
  size_vector = nvalues;

  vector_flag = 1;
  extvector = 0;
  maxatom = 0;

  if (!numflag) {
    memory->create(values_local, nvalues, "born_matrix:values_local");
  } else {

    reallocate();

    // set fixed-point to default = center of cell

    fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
    fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
    fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);

    // define the cartesian indices for each strain (Voigt order)

    dirlist[0][0] = 0;
    dirlist[0][1] = 0;
    dirlist[1][0] = 1;
    dirlist[1][1] = 1;
    dirlist[2][0] = 2;
    dirlist[2][1] = 2;

    dirlist[3][0] = 1;
    dirlist[3][1] = 2;
    dirlist[4][0] = 0;
    dirlist[4][1] = 2;
    dirlist[5][0] = 0;
    dirlist[5][1] = 1;
  }
}

/* ---------------------------------------------------------------------- */

ComputeBornMatrix::~ComputeBornMatrix()
{
  memory->destroy(values_global);
  memory->destroy(vector);
  if (!numflag) {
    memory->destroy(values_local);
  } else {
    memory->destroy(temp_x);
    memory->destroy(temp_f);
    delete[] id_virial;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBornMatrix::init()
{
  if (!numflag) {

    // need an occasional half neighbor list

    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  } else {

    // check for virial compute

    int icompute = modify->find_compute(id_virial);
    if (icompute < 0) error->all(FLERR, "Virial compute ID for compute born/matrix does not exist");
    compute_virial = modify->compute[icompute];

    // set up reverse index lookup
    // This table is used for consistency between numdiff and analytical
    // ordering of the terms.

    for (int m = 0; m < nvalues; m++) {
      int a = C_albe[m][0];
      int b = C_albe[m][1];
      revalbe[a][b] = m;
      revalbe[b][a] = m;
    }

    // reorder LAMMPS virial vector to Voigt order

    virialVtoV[0] = 0;
    virialVtoV[1] = 1;
    virialVtoV[2] = 2;
    virialVtoV[3] = 5;
    virialVtoV[4] = 4;
    virialVtoV[5] = 3;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBornMatrix::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   compute output vector
------------------------------------------------------------------------- */

void ComputeBornMatrix::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (!numflag) {

    // zero out arrays for one sample

    for (int m = 0; m < nvalues; m++) values_local[m] = 0.0;

    // compute Born contribution

    if (pairflag) compute_pairs();
    if (bondflag) compute_bonds();
    if (angleflag) compute_angles();
    if (dihedflag) compute_dihedrals();

    // sum Born contributions over all procs

    MPI_Allreduce(values_local, values_global, nvalues, MPI_DOUBLE, MPI_SUM, world);

  } else {

    // calculate Born matrix using stress finite differences

    compute_numdiff();

    // convert from pressure to energy units

    double inv_nktv2p = 1.0 / force->nktv2p;
    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (int m = 0; m < nvalues; m++) { values_global[m] *= inv_nktv2p * volume; }
  }

  for (int m = 0; m < nvalues; m++) vector[m] = values_global[m];
}

/* ----------------------------------------------------------------------
  compute Born contribution of local proc
------------------------------------------------------------------------- */

void ComputeBornMatrix::compute_pairs()

{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double rsq, factor_coul, factor_lj;
  double dupair, du2pair, rinv;
  int *ilist, *jlist, *numneigh, **firstneigh;

  int *type = atom->type;
  int *mask = atom->mask;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  // declares born values

  int a, b, c, d;
  double xi1, xi2, xi3;
  double xj1, xj2, xj3;
  double rij[3];
  double pair_pref;
  double r2inv;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xi1 = atom->x[i][0];
    xi2 = atom->x[i][1];
    xi3 = atom->x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      xj1 = atom->x[j][0];
      xj2 = atom->x[j][1];
      xj3 = atom->x[j][2];
      rij[0] = xj1 - xi1;
      rij[1] = xj2 - xi2;
      rij[2] = xj3 - xi3;
      rsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      r2inv = 1.0 / rsq;
      rinv = sqrt(r2inv);
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      // Add contribution to Born tensor

      pair_pref = dupair = du2pair = 0.0;
      pair->born_matrix(i, j, itype, jtype, rsq, factor_coul, factor_lj, dupair, du2pair);
      pair_pref = 0.5 * du2pair - dupair * rinv;

      // See albemunu in compute_born_matrix.h for indices order.

      a = b = c = d = 0;
      for (int m = 0; m < nvalues; m++) {
        a = albemunu[m][0];
        b = albemunu[m][1];
        c = albemunu[m][2];
        d = albemunu[m][3];
        values_local[m] += pair_pref * rij[a] * rij[b] * rij[c] * rij[d] * r2inv;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  compute Born matrix using virial stress finite differences
------------------------------------------------------------------------- */

void ComputeBornMatrix::compute_numdiff()
{
  int vec_index;

  // grow arrays if necessary

  int nall = atom->nlocal + atom->nghost;
  if (nall > maxatom) reallocate();

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) {
      temp_x[i][k] = x[i][k];
      temp_f[i][k] = f[i][k];
    }

  // loop over 6 strain directions
  // compute stress finite difference in each direction

  for (int idir = 0; idir < NDIR_VIRIAL; idir++) {

    // forward

    displace_atoms(nall, idir, 1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_index = revalbe[idir][jdir];
      values_global[vec_index] = compute_virial->vector[virialVtoV[jdir]];
    }
    restore_atoms(nall, idir);

    // backward

    displace_atoms(nall, idir, -1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_index = revalbe[idir][jdir];
      values_global[vec_index] -= compute_virial->vector[virialVtoV[jdir]];
    }
    restore_atoms(nall, idir);
  }

  // apply derivative factor

  double denominator = -0.5 / numdelta;
  for (int m = 0; m < nvalues; m++) values_global[m] *= denominator;

  // recompute virial so all virial and energy contributions are as before
  // also needed for virial stress addon contributions to Born matrix

  update_virial();

  // add on virial terms

  virial_addon();

  // restore original forces for owned and ghost atoms

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) f[i][k] = temp_f[i][k];
}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrix::displace_atoms(int nall, int idir, double magnitude)
{
  double **x = atom->x;

  // NOTE: virial_addon() expressions predicated on
  // shear strain fields (l != k) being symmetric here
  int k = dirlist[idir][0];
  int l = dirlist[idir][1];

  // axial strain

  if (l == k)
    for (int i = 0; i < nall; i++)
      x[i][k] = temp_x[i][k] + numdelta * magnitude * (temp_x[i][l] - fixedpoint[l]);

  // symmetric shear strain

  else
    for (int i = 0; i < nall; i++) {
      x[i][k] = temp_x[i][k] + 0.5 * numdelta * magnitude * (temp_x[i][l] - fixedpoint[l]);
      x[i][l] = temp_x[i][l] + 0.5 * numdelta * magnitude * (temp_x[i][k] - fixedpoint[k]);
    }
}

/* ----------------------------------------------------------------------
   restore position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrix::restore_atoms(int nall, int idir)
{

  // reset only idir coord

  int k = dirlist[idir][0];
  int l = dirlist[idir][1];
  double **x = atom->x;
  if (l == k)
    for (int i = 0; i < nall; i++) x[i][k] = temp_x[i][k];
  else
    for (int i = 0; i < nall; i++) {
      x[i][l] = temp_x[i][l];
      x[i][k] = temp_x[i][k];
    }
}

/* ----------------------------------------------------------------------
   evaluate potential forces and virial
   same logic as in Verlet
------------------------------------------------------------------------- */

void ComputeBornMatrix::update_virial()
{
  int eflag = 0;

  // this may not be completely general
  // but it works for lj/cut and sw pair styles
  // and compute stress/atom output is unaffected

  int vflag = VIRIAL_PAIR;

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  compute_virial->compute_vector();
}

/* ----------------------------------------------------------------------
   calculate virial stress addon terms to the Born matrix
   this is based on original code of Dr. Yubao Zhen
   described here: Comp. Phys. Comm. 183 (2012) 261-265
   as well as Yoshimoto et al., PRB, 71 (2005) 184108, Eq 15.and eq A3.
------------------------------------------------------------------------- */

void ComputeBornMatrix::virial_addon()
{
  double *sigv = compute_virial->vector;

  // This way of doing is not very elegant but is correct.
  // The complete Cijkl terms are the sum of symmetric terms
  // computed in compute_numdiff and virial stress terms.
  // The viral terms are not symmetric in the tensor computation.
  // For example:
  // C2323 = s33+D2323; C3232 = s22+D3232 etc...
  // However there are two symmetry breaking when reducing
  // the 4-rank tensor to a 2-rank tensor
  // Cijkl = (Bijkl+Bjikl+Bijlk+Bjilk)/4. = (Bijkl+Bjilk)/2.
  // and when computing only the 21 independant term.

  // these expressions have been verified
  // correct to about 1e-7 compared
  // to the analytic expressions for lj/cut,
  // predicated on shear strain fields being
  // symmetric in displace_atoms()

  values_global[0] += 2.0 * sigv[0];
  values_global[1] += 2.0 * sigv[1];
  values_global[2] += 2.0 * sigv[2];

  values_global[3] += 0.5 * (sigv[1] + sigv[2]);
  values_global[4] += 0.5 * (sigv[0] + sigv[2]);
  values_global[5] += 0.5 * (sigv[0] + sigv[1]);
  values_global[6] += 0.0;
  values_global[7] += 0.0;
  values_global[8] += 0.0;
  values_global[9] += sigv[4];
  values_global[10] += sigv[3];
  values_global[11] += 0.0;
  values_global[12] += sigv[5];
  values_global[13] += 0.0;
  values_global[14] += sigv[3];
  values_global[15] += sigv[5];
  values_global[16] += sigv[4];
  values_global[17] += 0.0;
  values_global[18] += 0.5 * sigv[3];
  values_global[19] += 0.5 * sigv[4];
  values_global[20] += 0.5 * sigv[5];
}

/* ----------------------------------------------------------------------
   clear forces needed
------------------------------------------------------------------------- */

void ComputeBornMatrix::force_clear(int nall)
{
  double **forces = atom->f;
  size_t nbytes = 3 * sizeof(double) * nall;
  if (nbytes) memset(&forces[0][0], 0, nbytes);
}

/* ----------------------------------------------------------------------
   reallocated local per-atoms arrays
------------------------------------------------------------------------- */

void ComputeBornMatrix::reallocate()
{
  memory->destroy(temp_x);
  memory->destroy(temp_f);
  maxatom = atom->nmax;
  memory->create(temp_x, maxatom, 3, "born/matrix:temp_x");
  memory->create(temp_f, maxatom, 3, "born/matrix:temp_f");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double ComputeBornMatrix::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) 2 * maxatom * 3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   Count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted or turned off (type <= 0)
   do not count or count contribution
---------------------------------------------------------------------- */

void ComputeBornMatrix::compute_bonds()
{
  int i, m, nb, atom1, atom2, imol, iatom, btype;
  tagint tagprev;
  double dx, dy, dz, rsq;

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

  int a, b, c, d;
  double rij[3];
  double rinv, r2inv;
  double pair_pref, dupair, du2pair;

  // loop over all atoms and their bonds

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

      dx = x[atom2][0] - x[atom1][0];
      dy = x[atom2][1] - x[atom1][1];
      dz = x[atom2][2] - x[atom1][2];
      domain->minimum_image(dx, dy, dz);
      rsq = dx * dx + dy * dy + dz * dz;
      rij[0] = dx;
      rij[1] = dy;
      rij[2] = dz;
      r2inv = 1.0 / rsq;
      rinv = sqrt(r2inv);

      pair_pref = dupair = du2pair = 0.0;
      bond->born_matrix(btype, rsq, atom1, atom2, dupair, du2pair);
      pair_pref = du2pair - dupair * rinv;

      a = b = c = d = 0;
      for (m = 0; m < 21; m++) {
        a = albemunu[m][0];
        b = albemunu[m][1];
        c = albemunu[m][2];
        d = albemunu[m][3];
        values_local[m] += pair_pref * rij[a] * rij[b] * rij[c] * rij[d] * r2inv;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   count angles and compute angle info on this proc
   only count if 2nd atom is the one storing the angle
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted or turned off (type <= 0)
   do not count or count contribution
---------------------------------------------------------------------- */

void ComputeBornMatrix::compute_angles()
{
  int i, m, na, atom1, atom2, atom3, imol, iatom, atype;
  tagint tagprev;
  double delx1, dely1, delz1, delx2, dely2, delz2;
  double rsq1, rsq2, r1, r2, cost;
  double r1r2, r1r2inv;
  double rsq1inv, rsq2inv, cinv;

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

  int a, b, c, d, e, f;
  double duang, du2ang;
  double del1[3], del2[3];
  double dcos[6];
  double d2cos;
  double d2lncos;

  // Initializing array for intermediate cos derivatives
  // w regard to strain
  for (i = 0; i < 6; i++) dcos[i] = 0;

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

    for (i = 0; i < na; i++) {
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

      delx1 = x[atom1][0] - x[atom2][0];
      dely1 = x[atom1][1] - x[atom2][1];
      delz1 = x[atom1][2] - x[atom2][2];
      domain->minimum_image(delx1, dely1, delz1);
      del1[0] = delx1;
      del1[1] = dely1;
      del1[2] = delz1;

      rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
      rsq1inv = 1.0 / rsq1;
      r1 = sqrt(rsq1);

      delx2 = x[atom3][0] - x[atom2][0];
      dely2 = x[atom3][1] - x[atom2][1];
      delz2 = x[atom3][2] - x[atom2][2];
      domain->minimum_image(delx2, dely2, delz2);
      del2[0] = delx2;
      del2[1] = dely2;
      del2[2] = delz2;

      rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
      rsq2inv = 1.0 / rsq2;
      r2 = sqrt(rsq2);

      r1r2 = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
      r1r2inv = 0;
      if (r1r2 == 0.) {
        // Worst case scenario. A 0 cosine has an undefined logarithm.
        // We then add a small amount of the third bond to the first one
        // so they are not orthogonal anymore and recompute.
        del1[0] += SMALL * del2[0];
        del1[1] += SMALL * del2[1];
        del1[2] += SMALL * del2[2];
        r1r2 = del1[0] * del2[0] + del1[1] * del2[1] + del1[2] * del2[2];
        r1r2inv = 1 / r1r2;

        // cost = cosine of angle
        cost = r1r2 / (r1 * r2);

      } else {
        r1r2inv = 1 / r1r2;
        cost = r1r2 / (r1 * r2);
      }

      if (cost > 1.0) cost = 1.0;
      if (cost < -1.0) cost = -1.0;
      cinv = 1.0 / cost;

      // The method must return derivative with regards
      // to cos(theta)!
      // Use the chain rule if needed:
      // dU(t)/de = dt/dcos(t)*dU(t)/dt*dcos(t)/de
      // with dt/dcos(t) = -1/sin(t)
      angle->born_matrix(atype, atom1, atom2, atom3, duang, du2ang);

      // Voigt notation
      // 1 = 11, 2 = 22, 3 = 33
      // 4 = 23, 5 = 13, 6 = 12
      a = b = c = d = 0;
      // clang-format off
      for (m = 0; m<6; m++) {
        a = sigma_albe[m][0];
        b = sigma_albe[m][1];
        dcos[m] = cost*((del1[a]*del2[b]+del1[b]*del2[a])*r1r2inv -
                        del1[a]*del1[b]*rsq1inv - del2[a]*del2[b]*rsq2inv);
      }
      for (m = 0; m<21; m++) {
        a = albemunu[m][0];
        b = albemunu[m][1];
        c = albemunu[m][2];
        d = albemunu[m][3];
        e = C_albe[m][0];
        f = C_albe[m][1];
        d2lncos = 2*(del1[a]*del1[b]*del1[c]*del1[d]*rsq1inv*rsq1inv +
                        del2[a]*del2[b]*del2[c]*del2[d]*rsq2inv*rsq2inv) -
                       (del1[a]*del2[b]+del1[b]*del2[a]) *
                       (del1[c]*del2[d]+del1[d]*del2[c]) *
                       r1r2inv*r1r2inv;
        d2cos =  cost*d2lncos + dcos[e]*dcos[f]*cinv;
        values_local[m] += duang*d2cos + du2ang*dcos[e]*dcos[f];
      }
      // clang-format on
    }
  }
}

/* ----------------------------------------------------------------------
   count dihedrals on this proc
   only count if 2nd atom is the one storing the dihedral
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if flag is set, compute requested info about dihedral
------------------------------------------------------------------------- */

void ComputeBornMatrix::compute_dihedrals()
{
  int i, m, nd, atom1, atom2, atom3, atom4, imol, iatom;
  tagint tagprev;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
  double ax, ay, az, bx, by, bz, rasq, rbsq, dotab, ra2inv, rb2inv, dotabinv, rabinv;
  double co;

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

  int al, be, mu, nu, e, f;
  double dudih, du2dih, b1sq, b2sq, b3sq, b1b2, b1b3, b2b3;
  double b1[3], b2[3], b3[3];

  // 1st and 2nd order derivatives of the dot products
  double dab[6], daa[6], dbb[6];
  double d2ab, d2aa, d2bb;

  double dcos[6], d2cos;

  for (i = 0; i < 6; i++) dab[i] = daa[i] = dbb[i] = dcos[i] = 0;

  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    // if (molecular == 1)
    //   nd = num_dihedral[atom2];
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

      // phi calculation from dihedral style harmonic

      // The method must return derivative with regards
      // to cos(phi)!
      // Use the chain rule if needed:
      // dU(t)/de = dt/dcos(t)*dU(t)/dt*dcos(t)/de
      // with dt/dcos(t) = -1/sin(t)

      dihedral->born_matrix(i, atom1, atom2, atom3, atom4, dudih, du2dih);

      vb1x = x[atom2][0] - x[atom1][0];
      vb1y = x[atom2][1] - x[atom1][1];
      vb1z = x[atom2][2] - x[atom1][2];
      domain->minimum_image(vb1x, vb1y, vb1z);
      b1[0] = vb1x;
      b1[1] = vb1y;
      b1[2] = vb1z;
      b1sq = b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2];

      vb2x = x[atom3][0] - x[atom2][0];
      vb2y = x[atom3][1] - x[atom2][1];
      vb2z = x[atom3][2] - x[atom2][2];
      domain->minimum_image(vb2x, vb2y, vb2z);
      b2[0] = vb2x;
      b2[1] = vb2y;
      b2[2] = vb2z;
      b2sq = b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2];

      vb3x = x[atom4][0] - x[atom3][0];
      vb3y = x[atom4][1] - x[atom3][1];
      vb3z = x[atom4][2] - x[atom3][2];
      domain->minimum_image(vb3x, vb3y, vb3z);
      b3[0] = vb3x;
      b3[1] = vb3y;
      b3[2] = vb3z;
      b3sq = b3[0] * b3[0] + b3[1] * b3[1] + b3[2] * b3[2];

      b1b2 = b1[0] * b2[0] + b1[1] * b2[1] + b1[2] * b2[2];
      b1b3 = b1[0] * b3[0] + b1[1] * b3[1] + b1[2] * b3[2];
      b2b3 = b2[0] * b3[0] + b2[1] * b3[1] + b2[2] * b3[2];

      // a = b1xb2; b = b2xb3
      ax = vb1y * vb2z - vb1z * vb2y;
      ay = vb1z * vb2x - vb1x * vb2z;
      az = vb1x * vb2y - vb1y * vb2x;
      bx = vb2y * vb3z - vb2z * vb3y;
      by = vb2z * vb3x - vb2x * vb3z;
      bz = vb2x * vb3y - vb2y * vb3x;

      rasq = ax * ax + ay * ay + az * az;
      rbsq = bx * bx + by * by + bz * bz;
      dotab = ax * bx + ay * by + az * bz;

      ra2inv = rb2inv = rabinv = dotabinv = 0.0;
      if (rasq > 0) ra2inv = 1.0 / rasq;
      if (rbsq > 0) rb2inv = 1.0 / rbsq;
      if (dotab != 0.) dotabinv = 1.0 / dotab;
      rabinv = sqrt(ra2inv * rb2inv);

      co = (ax * bx + ay * by + az * bz) * rabinv;

      if (co == 0.0) {
        // Worst case scenario. A 0 cosine has an undefined logarithm.
        // We then add a small amount of the third bond to the first one
        // so they are not orthogonal anymore and recompute.
        b1[0] += SMALL * b3[0];
        b1[1] += SMALL * b3[1];
        b1[2] += SMALL * b3[2];
        b1sq = b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2];
        b3sq = b3[0] * b3[0] + b3[1] * b3[1] + b3[2] * b3[2];
        b1b2 = b1[0] * b2[0] + b1[1] * b2[1] + b1[2] * b2[2];
        b1b3 = b1[0] * b3[0] + b1[1] * b3[1] + b1[2] * b3[2];
        b2b3 = b2[0] * b3[0] + b2[1] * b3[1] + b2[2] * b3[2];
        ax = b1[1] * b2[2] - b1[2] * b2[1];
        ay = b1[2] * b2[0] - b1[0] * b2[2];
        az = b1[0] * b2[1] - b1[1] * b2[0];
        bx = b2[1] * b3[2] - b2[2] * b3[1];
        by = b2[2] * b3[0] - b2[0] * b3[2];
        bz = b2[0] * b3[1] - b2[1] * b3[0];
        rasq = ax * ax + ay * ay + az * az;
        rbsq = bx * bx + by * by + bz * bz;
        dotab = ax * bx + ay * by + az * bz;
        ra2inv = rb2inv = rabinv = dotabinv = 0.0;
        if (rasq > 0) ra2inv = 1.0 / rasq;
        if (rbsq > 0) rb2inv = 1.0 / rbsq;
        if (dotab != 0.) dotabinv = 1.0 / dotab;
        rabinv = sqrt(ra2inv * rb2inv);
        co = (ax * bx + ay * by + az * bz) * rabinv;
      }
      if (co > 1.0) co = 1.0;
      if (co < -1.0) co = -1.0;

      al = be = mu = nu = e = f = 0;
      // clang-format off
      for (m = 0; m<6; m++) {
        al = sigma_albe[m][0];
        be = sigma_albe[m][1];

        daa[m] = 2*(b2sq*b1[al]*b1[be]
                  + b1sq*b2[al]*b2[be]
                  - b1b2*(b1[al]*b2[be]+b2[al]*b1[be]));

        dbb[m] = 2*(b2sq*b3[al]*b3[be]
                  + b3sq*b2[al]*b2[be]
                  - b2b3*(b3[al]*b2[be]+b2[al]*b3[be]));

        dab[m] = b1b2*(b2[al]*b3[be]+b3[al]*b2[be])
               + b2b3*(b1[al]*b2[be]+b2[al]*b1[be])
               - b1b3*(b2[al]*b2[be]+b2[al]*b2[be])
               - b2sq*(b1[al]*b3[be]+b3[al]*b1[be]);


        dcos[m] = 0.5*co*(2*dotabinv*dab[m] - ra2inv*daa[m] - rb2inv*dbb[m]);
      }

      for (m = 0; m<21; m++) {
        al = albemunu[m][0];
        be = albemunu[m][1];
        mu = albemunu[m][2];
        nu = albemunu[m][3];
        e = C_albe[m][0];
        f = C_albe[m][1];

        d2aa = 4*(b1[al]*b1[be]*b2[mu]*b2[nu] + b1[mu]*b1[nu]*b2[al]*b2[be])
             - 2*(b1[al]*b2[be]+b1[be]*b2[al])*(b1[mu]*b2[nu]+b1[nu]*b2[mu]);

        d2bb = 4*(b3[al]*b3[be]*b2[mu]*b2[nu] + b3[mu]*b3[nu]*b2[al]*b2[be])
             - 2*(b3[al]*b2[be]+b3[be]*b2[al])*(b3[mu]*b2[nu]+b3[nu]*b2[mu]);

        d2ab = (b1[al]*b2[be]+b2[al]*b1[be])*(b3[mu]*b2[nu]+b2[mu]*b3[nu])
             + (b3[al]*b2[be]+b2[al]*b3[be])*(b1[mu]*b2[nu]+b2[mu]*b1[nu])
             - (b1[al]*b3[be]+b3[al]*b1[be])*(b2[mu]*b2[nu]+b2[mu]*b2[nu])
             - (b2[al]*b2[be]+b2[al]*b2[be])*(b1[mu]*b3[nu]+b3[mu]*b1[nu]);

        d2cos = -2*dotabinv*dotabinv*dab[e]*dab[f] + 2*dotabinv*d2ab
                + ra2inv*ra2inv*daa[e]*daa[f] - ra2inv*d2aa
                + rb2inv*rb2inv*dbb[e]*dbb[f] - rb2inv*d2bb
                + 2*dcos[e]*dcos[f]/co/co;
        d2cos *= 0.5*co;

        values_local[m] += dudih*d2cos;
        values_local[m] += du2dih*dcos[e]*dcos[f];
      }
      // clang-format on
    }
  }
}
