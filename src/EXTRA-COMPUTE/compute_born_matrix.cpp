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

  MPI_Comm_rank(world, &me);

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
  dt = update->dt;

  if (!numflag) {

    // need an occasional half neighbor list

    neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    // int irequest = neighbor->request((void *) this);
    // neighbor->requests[irequest]->pair = 0;
    // neighbor->requests[irequest]->compute = 1;
    // neighbor->requests[irequest]->occasional = 1;

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
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

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
  double fi1, fi2, fi3;
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

      if (newton_pair || j < nlocal) {

        // Add contribution to Born tensor

        pair->born_matrix(i, j, itype, jtype, rsq, factor_coul, factor_lj, dupair, du2pair);
        pair_pref = du2pair - dupair * rinv;

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
}

/* ----------------------------------------------------------------------
  compute Born matrix using virial stress finite differences
------------------------------------------------------------------------- */

void ComputeBornMatrix::compute_numdiff()
{
  double energy;
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

  int flag, allflag;

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
  for (int m = 0; m < nvalues; m++)
    values_global[m] *= denominator;

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

  int kd, nd, id, jd;
  int m;

  double* sigv = compute_virial->vector;

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
  int i, j, m, n, nb, atom1, atom2, imol, iatom, btype, ivar;
  tagint tagprev;
  double dx, dy, dz, rsq;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
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

  m = 0;
  while (m < nvalues) {

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
        for (j = 0; j < 21; j++) {
          a = albemunu[j][0];
          b = albemunu[j][1];
          c = albemunu[j][2];
          d = albemunu[j][3];
          values_local[m + j] += pair_pref * rij[a] * rij[b] * rij[c] * rij[d] * r2inv;
        }
      }
    }
    m += 21;
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
  int i, j, m, n, na, atom1, atom2, atom3, imol, iatom, atype, ivar;
  tagint tagprev;
  double delx1, dely1, delz1, delx2, dely2, delz2;
  double rsq1, rsq2, r1, r2, cost;
  double r1r2, r1r2inv;
  double rsq1inv, rsq2inv, r1inv, r2inv, cinv;
  double *ptr;

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
      r1inv = 1.0 / r1;

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
      r2inv = 1.0 / r2;

      r1r2 = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
      r1r2inv = 1 / r1r2;
      // cost = cosine of angle

      cost = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
      cost /= r1 * r2;
      if (cost > 1.0) cost = 1.0;
      if (cost < -1.0) cost = -1.0;
      if (cost == 0.) {
          cinv = 1.0 / cost;
      } else {
          cinv = 0.;
      }

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
      for (j = 0; j<6; j++) {
          if (cost == 0) {
              dcos[j] = 0.;
          } else {
            a = sigma_albe[j][0];
            b = sigma_albe[j][1];
            dcos[j] = cost*((del1[a]*del2[b]+del1[b]*del2[a])*r1r2inv -
                            del1[a]*del1[b]*rsq1inv - del2[a]*del2[b]*rsq2inv);
          }
      }
      for (j = 0; j<21; j++) {
        a = albemunu[j][0];
        b = albemunu[j][1];
        c = albemunu[j][2];
        d = albemunu[j][3];
        e = C_albe[j][0];
        f = C_albe[j][1];
        d2lncos = 2*(del1[a]*del1[b]*del1[c]*del1[d]*rsq1inv*rsq1inv +
                        del2[a]*del2[b]*del2[c]*del2[d]*rsq2inv*rsq2inv) -
                       (del1[a]*del2[b]+del1[b]*del2[a]) *
                       (del1[c]*del2[d]+del1[d]*del2[c]) *
                       r1r2inv*r1r2inv;
        d2cos =  cost*d2lncos + dcos[e]*dcos[f]*cinv;
        values_local[m+j] += duang*d2cos + du2ang*dcos[e]*dcos[f];
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
  int i, j, m, n, nd, atom1, atom2, atom3, atom4, imol, iatom, dtype, ivar;
  tagint tagprev;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z, vb2xm, vb2ym, vb2zm;
  double ax, ay, az, bx, by, bz, rasq, rbsq, dotab, rgsq, rg, ra2inv, rb2inv, dotabinv, rabinv;
  double si, co, phi;
  double *ptr;

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
  int al, be, mu, nu, e, f;
  double b1sq;
  double b2sq;
  double b3sq;
  double b1b2;
  double b1b3;
  double b2b3;
  double b1[3];
  double b2[3];
  double b3[3];

  // Actually derivatives of the square of the
  // vectors dot product.
  double dab[6];
  double daa[6];
  double dbb[6];
  double d2ab;
  double d2aa;
  double d2bb;

  double dcos[6];
  double d2cos;

  for (i = 0; i < 6; i++) dab[i] = daa[i] = dbb[i] = dcos[i] = 0;

  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == 1)
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

      dihedral->born_matrix(nd, atom1, atom2, atom3, atom4, dudih, du2dih);
      printf("Energy: %f, %f\n", dudih, du2dih);

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
      printf("b1 :");
      for (i = 0; i<3; i++) printf(" %f ", b1[i]);
      printf("\n");
      printf("b2 :");
      for (i = 0; i<3; i++) printf(" %f ", b2[i]);
      printf("\n");
      printf("b3 :");
      for (i = 0; i<3; i++) printf(" %f ", b3[i]);
      printf("\n");

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
      rgsq = vb2x * vb2x + vb2y * vb2y + vb2z * vb2z;
      dotab = ax*bx + ay*by + az*bz;
      rg = sqrt(rgsq);

      ra2inv = rb2inv = rabinv = dotabinv = 0.0;
      if (rasq > 0) ra2inv = 1.0 / rasq;
      if (rbsq > 0) rb2inv = 1.0 / rbsq;
      dotabinv = 1.0/dotab;
      rabinv = sqrt(ra2inv * rb2inv);
      printf("a :");
      printf(" %f %f %f %f ", ax, ay, az, ra2inv);
      printf("b :");
      printf(" %f %f %f %f ", bx, by, bz, rb2inv);
      printf("rabinv :");
      printf(" %f", dotabinv);
      printf("\n");
      printf("TESTa1: %f, TESTa2: %f\n", rasq, b1sq*b2sq-b1b2*b1b2);
      printf("TESTb1: %f, TESTb2: %f\n", rbsq, b3sq*b2sq-b2b3*b2b3);
      printf("TESTab1: %f, TESTab2: %f\n", 1/dotabinv, b1b3*b2sq-b1b2*b2b3);

      co = (ax * bx + ay * by + az * bz) * rabinv;
      si = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);

      if (co > 1.0) co = 1.0;
      if (co < -1.0) co = -1.0;
      phi = atan2(si, co);
      printf("Cos: %18.15g, Sin: %18.15g, Phi: %18.15g\n", co, si, phi);

      al = be = mu = nu = e = f = 0;
      double d2lncos;
      int test1 = 19;
      int test2 = -1;
      int test3 = -1;
      // clang-format off
      for (j = 0; j<6; j++) {
        al = sigma_albe[j][0];
        be = sigma_albe[j][1];

        daa[j] = 2*(b2sq*b1[al]*b1[be]
                  + b1sq*b2[al]*b2[be]
                  - b1b2*(b1[al]*b2[be]+b2[al]*b1[be]));

        dbb[j] = 2*(b2sq*b3[al]*b3[be]
                  + b3sq*b2[al]*b2[be]
                  - b2b3*(b3[al]*b2[be]+b2[al]*b3[be]));

        dab[j] = b1b2*(b2[al]*b3[be]+b3[al]*b2[be])
               + b2b3*(b1[al]*b2[be]+b2[al]*b1[be])
               - b1b3*(b2[al]*b2[be]+b2[al]*b2[be])
               - b2sq*(b1[al]*b3[be]+b3[al]*b1[be]);

        if (j == test1) {
            printf("b3sq = %f, b2[al] = %f, b2[be] = %f\n", b3sq, b2[al], b3[be]);
            printf("daa1 = %18.15g, daa2 = %18.15g, daa3 = %18.15g\n", b2sq*b1[al]*b1[be], b1sq*b2[al]*b2[be], b1b2*(b1[al]*b2[be]+b2[al]*b1[be]));
            printf("dbb1 = %18.15g, dbb2 = %18.15g, dbb3 = %18.15g\n", b2sq*b3[al]*b3[be], b3sq*b2[al]*b2[be], b2b3*(b3[al]*b2[be]+b2[al]*b3[be]));
            printf("dab1 = %18.15g, dab2 = %18.15g, dab3 = %18.15g, dab4 = %18.15g\n", b1b2*(b2[al]*b3[be]+b3[al]*b2[be]), b2b3*(b1[al]*b2[be]+b2[al]*b1[be]), -b1b3*(b2[al]*b2[be]+b2[al]*b2[be]), -b2sq*(b1[al]*b3[be]+b3[al]*b1[be]));
        }

        dcos[j] = 0.5*co*(2*dotabinv*dab[j] - ra2inv*daa[j] - rb2inv*dbb[j]);
        if (j == test1 || j == test2) {
            printf("order 1 %d al: %d, be: %d\n", j+1, al, be);
            printf("daa = %18.15g, dbb = %18.15g, dab = %18.15g\n", daa[j], dbb[j], dab[j]);
            printf("daa/raa = %18.15g, dbb/rbb = %18.15g, dab/rab = %18.15g\n", ra2inv*daa[j], rb2inv*dbb[j], dotabinv*dab[j]);
            printf("dcos = %e\n", dcos[j]);
        }
      }
      printf("dcos:\n");
      printf("%e %e %e\n", dcos[0], dcos[1], dcos[2]);
      printf("%e %e %e\n", dcos[3], dcos[4], dcos[5]);
      printf("\n");
      printf("daa:\n");
      printf("%e %e %e\n", daa[0], daa[1], daa[2]);
      printf("%e %e %e\n", daa[3], daa[4], daa[5]);
      printf("\n");
      printf("dbb:\n");
      printf("%e %e %e\n", dbb[0], dbb[1], dbb[2]);
      printf("%e %e %e\n", dbb[3], dbb[4], dbb[5]);
      printf("\n");
      printf("dab:\n");
      printf("%e %e %e\n", dab[0], dab[1], dab[2]);
      printf("%e %e %e\n", dab[3], dab[4], dab[5]);
      printf("\n");

      for (j = 0; j<21; j++) {
        al = albemunu[j][0];
        be = albemunu[j][1];
        mu = albemunu[j][2];
        nu = albemunu[j][3];
        e = C_albe[j][0];
        f = C_albe[j][1];

        d2aa = 4*(b1[al]*b1[be]*b2[mu]*b2[nu] + b1[mu]*b1[nu]*b2[al]*b2[be])
             - 2*(b1[al]*b2[be]+b1[be]*b2[al])*(b1[mu]*b2[nu]+b1[nu]*b2[mu]);

        d2bb = 4*(b3[al]*b3[be]*b2[mu]*b2[nu] + b3[mu]*b3[nu]*b2[al]*b2[be])
             - 2*(b3[al]*b2[be]+b3[be]*b2[al])*(b3[mu]*b2[nu]+b3[nu]*b2[mu]);

        d2ab = (b1[al]*b2[be]+b2[al]*b1[be])*(b3[mu]*b2[nu]+b2[mu]*b3[nu])
             + (b3[al]*b2[be]+b2[al]*b3[be])*(b1[mu]*b2[nu]+b2[mu]*b1[nu])
             - (b1[al]*b3[be]+b3[al]*b1[be])*(b2[mu]*b2[nu]+b2[mu]*b2[nu])
             - (b2[al]*b2[be]+b2[al]*b2[be])*(b1[mu]*b3[nu]+b3[mu]*b1[nu]);

        if (j == test1 || j == test2 || j == test3) {
            printf("b1al = %g ", b1[al]);
            printf("b1be = %g ", b1[be]);
            printf("b1mu = %g ", b1[mu]);
            printf("b1nu = %g ", b1[nu]);
            printf("b2al = %g ", b2[al]);
            printf("b2be = %g ", b2[be]);
            printf("b2mu = %g ", b2[mu]);
            printf("b2nu = %g ", b2[nu]);
            printf("b3al = %g ", b3[al]);
            printf("b3be = %g ", b3[be]);
            printf("b3mu = %g ", b3[mu]);
            printf("b3nu = %g \n", b3[nu]);
            printf("d2aa details:\n");
            printf("1: %e\n", 4*(b1[al]*b1[be]*b2[mu]*b2[nu]+b1[mu]*b1[nu]*b2[al]*b2[be]));
            printf("2: %e\n", 2*(b1[al]*b2[be]*(b1[mu]*b2[nu] + b1[nu]*b2[mu])));
            printf("3: %e\n", 2*(b1[be]*b2[al]*(b1[mu]*b2[nu] + b1[nu]*b2[mu])));
            printf("d2bb details:\n");
            printf("1: %e\n", 4*(b3[al]*b3[be]*b2[mu]*b2[nu]+b3[mu]*b3[nu]*b2[al]*b2[be]));
            printf("2: %e\n", 2*(b3[al]*b2[be]*(b3[mu]*b2[nu] + b3[nu]*b2[mu])));
            printf("3: %e\n", 2*(b3[be]*b2[al]*(b3[mu]*b2[nu] + b3[nu]*b2[mu])));
            printf("d2ab details:\n");
            printf("1: %e\n", (b1[al]*b3[be]+b1[be]*b3[al])*(b2[mu]*b2[nu]+b2[nu]*b2[mu]));
            printf("2: %e\n", 2*b2[al]*b2[be]*(b1[mu]*b3[nu]+b1[nu]*b3[mu]));
            printf("3: %e\n", (b2[al]*b3[be]+b2[be]*b3[al])*(b1[mu]*b2[nu]+b1[nu]*b2[mu]));
            printf("4: %e\n", (b2[al]*b1[be]+b2[be]*b1[al])*(b2[mu]*b3[nu]+b2[nu]*b3[mu]));
        }

        // d2lncos = 2*(dotabinv*d2ab - dotabinv*dotabinv*dab[e]*dab[f]) + (ra2inv*ra2inv*daa[e]*daa[f] - ra2inv*d2aa) + (rb2inv*rb2inv*dbb[e]*dbb[f] - rb2inv*d2bb);
        // d2cos = 0.5*co*d2lncos + dcos[e]*dcos[f]/co/co;
        d2cos = 0.5*co*(-2*dotabinv*dotabinv*dab[e]*dab[f] + 2*dotabinv*d2ab + ra2inv*ra2inv*daa[e]*daa[f] - ra2inv*d2aa + rb2inv*rb2inv*dbb[e]*dbb[f] - rb2inv*d2bb + 2*dcos[e]*dcos[f]/co);

        values_local[m+j] += dudih*d2cos;
        values_local[m+j] += du2dih*dcos[e]*dcos[f];
        if (j == test1 || j == test2 || j == test3) {
            printf("order 2 %d al: %d, be: %d\n", j+1, al, be);
            printf("           mu: %d, nu: %d\n", mu, nu);
            printf("           e: %d, f: %d\n", e, f);
            printf("d2aa = %18.15e, d2bb = %18.15e, d2ab = %18.15e\n", d2aa, d2bb, d2ab);
            printf("(ab) t1 = %18.15e; t2 = %18.15e\n", dotabinv*d2ab, dotabinv*dotabinv*dab[e]*dab[f]);
            printf("(aa) t1 = %18.15e; t2 = %18.15e\n", ra2inv*d2aa, ra2inv*ra2inv*daa[e]*daa[f]);
            printf("(aa) 0.5*(t1-t2) = %18.15e\n", 0.5*(ra2inv*d2aa-ra2inv*ra2inv*daa[e]*daa[f]));
            printf("(bb) t1 = %18.15e; t2 = %18.15e\n", rb2inv*d2bb, rb2inv*rb2inv*dbb[e]*dbb[f]);
            printf("(bb) 0.5*(t1-t2) = %18.15e\n", 0.5*(rb2inv*d2bb - rb2inv*rb2inv*dbb[e]*dbb[f]));
            printf("D1 = %18.15e; D2 = %18.15e\n", dotabinv*d2ab+0.5*(ra2inv*d2aa+rb2inv*d2bb), dotabinv*dotabinv*dab[e]*dab[f]-0.5*(ra2inv*ra2inv*daa[e]*daa[f]+rb2inv*rb2inv*dbb[e]*dbb[f]));
            printf("d2lncos = %18.15e\n", d2lncos);
            printf("co*d2lncos = %18.15e; dcos*dcos/co = %18.15e\n", co*d2lncos, dcos[e]*dcos[f]/co);
            printf("d2cos = %e\n", d2cos);
            printf("dudih*d2cos = %e; du2dih*dcos*dcos = %e\n", dudih*d2cos, du2dih*dcos[e]*dcos[f]);
        }
      }
      printf("\n");
      // clang-format on
    }
  }
}
