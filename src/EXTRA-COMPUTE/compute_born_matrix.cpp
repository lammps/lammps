// clang-format off
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
#include "update.h"
#include "universe.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1000000000

// This table is used to pick the 3d rij vector indices used to
// compute the 6 indices long Voigt stress vector
static int constexpr sigma_albe[6][2] = {
    {0, 0},    // s11
    {1, 1},    // s22
    {2, 2},    // s33
    {1, 2},    // s44
    {0, 2},    // s55
    {0, 1},    // s66
};

// This table is used to pick the correct indices from the Voigt
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

// This table is used to pick the 3d rij vector indices used to
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
    Compute(lmp, narg, arg), id_virial(nullptr), temp_x(nullptr), 
    temp_f(nullptr)
{
  if (narg < 3) error->all(FLERR,"Illegal compute born/matrix command");

  MPI_Comm_rank(world, &me);
  // For now the matrix can be computed as a 21 element vector

  nvalues = 21;

  // Error check

  numflag = 0;
  numdelta = 0.0;

  pairflag = 0;
  bondflag = 0;
  angleflag = 0;
  dihedflag = 0;
  impflag = 0;
  if (narg == 3) {
    pairflag = 1;
    bondflag = 1;
    angleflag = 1;
    dihedflag = 1;
    impflag = 1;
  } else {
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"numdiff") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute born/matrix command");
        numflag = 1;
        numdelta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (numdelta <= 0.0) error->all(FLERR, "Illegal compute born/matrix command");
        id_virial = utils::strdup(arg[iarg+2]);
        int icompute = modify->find_compute(id_virial);
        if (icompute < 0) error->all(FLERR,"Could not find compute born/matrix pressure ID");
        compute_virial = modify->compute[icompute];
        if (compute_virial->pressflag == 0)
        error->all(FLERR,"Compute born/matrix pressure ID does not compute pressure");
        iarg += 3;
      } else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) impflag = 1;
      else error->all(FLERR,"Illegal compute born/matrix command");
      ++iarg;
    }
  }

  if (pairflag) {
    if (numflag) error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->pair) {
      if (force->pair->born_matrix_enable == 0) {
        if (comm->me == 0) error->warning(FLERR, "Pair style does not support compute born/matrix");
      }
    } else {
      pairflag = 0;
    }
  }
  if (bondflag) {
    if (numflag) error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->bond) {
      if (force->bond->born_matrix_enable == 0) {
        if (comm->me == 0) error->warning(FLERR, "Bond style does not support compute born/matrix");
      }
    } else {
      bondflag = 0;
    }
  }
  if (angleflag) {
    if (numflag) error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->angle) {
      if (force->angle->born_matrix_enable == 0) {
        if (comm->me == 0) error->warning(FLERR, "Angle style does not support compute born/matrix");
      }
    } else {
      angleflag = 0;
    }
  }
  if (dihedflag) {
    if (numflag) error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->dihedral) {
     if (force->dihedral->born_matrix_enable == 0) {
       if (comm->me == 0) error->warning(FLERR, "Dihedral style does not support compute born/matrix");
     }
    } else {
      dihedflag = 0;
    }
  }
  if (impflag) {
    if (numflag) error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->improper) {
      if (force->improper->born_matrix_enable == 0) {
        if (comm->me == 0) error->warning(FLERR, "Improper style does not support compute born/matrix");
      }
    } else {
      impflag = 0;
    }
  }
  if (force->kspace) {
    if (!numflag) error->warning(FLERR, "KSPACE contribution not supported by compute born/matrix");
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
  //Timestep value
  dt = update->dt;

  if (!numflag) {

    // need an occasional half neighbor list
    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;

  } else {

    // check for virial compute

    int icompute = modify->find_compute(id_virial);
    if (icompute < 0) error->all(FLERR, "Virial compute ID for compute born/matrix does not exist");
    compute_virial = modify->compute[icompute];

    // set up reverse index lookup

    for (int m = 0; m < nvalues; m++) {
      int a = C_albe[m][0];
      int b = C_albe[m][1];
      revalbe[a][b] = m;
      revalbe[b][a] = m;
    }

    // for (int a = 0; a < NDIR_VIRIAL; a++) {
    //   for (int b = 0; b < NDIR_VIRIAL; b++) {
    //     printf("%d ",revalbe[a][b]);
    //   }
    //   printf("\n");
    // }

    // voigt3VtoM notation in normal physics sense,
    // 3x3 matrix and vector indexing
    //  i-j:      (1-1), (2-2), (3-3), (2-3), (1-3), (1-2)
    //  voigt3VtoM: 1      2      3      4      5      6

    voigt3VtoM[0][0]=0;  //  for 1
    voigt3VtoM[0][1]=0;
    voigt3VtoM[1][0]=1;  //  for 2
    voigt3VtoM[1][1]=1;
    voigt3VtoM[2][0]=2;  //  for 3
    voigt3VtoM[2][1]=2;
    voigt3VtoM[3][0]=1;  //  for 4 
    voigt3VtoM[3][1]=2;
    voigt3VtoM[4][0]=0;  //  for 5 
    voigt3VtoM[4][1]=2;
    voigt3VtoM[5][0]=0;  //  for 6 
    voigt3VtoM[5][1]=1;

    // to convert to vector indexing:
    // matrix index to vector index, double -> single index
    // this is not used at all

    voigt3MtoV[0][0]=0;    voigt3MtoV[0][1]=5;    voigt3MtoV[0][2]=4;
    voigt3MtoV[1][0]=5;    voigt3MtoV[1][1]=1;    voigt3MtoV[1][2]=3;
    voigt3MtoV[2][0]=4;    voigt3MtoV[2][1]=3;    voigt3MtoV[2][2]=2;

    // this is just for the virial.
    // since they use the xx,yy,zz,xy,xz,yz
    // order not the ordinary voigt

    virialMtoV[0][0]=0;    virialMtoV[0][1]=3;    virialMtoV[0][2]=4;
    virialMtoV[1][0]=3;    virialMtoV[1][1]=1;    virialMtoV[1][2]=5;
    virialMtoV[2][0]=4;    virialMtoV[2][1]=5;    virialMtoV[2][2]=2;

    // reorder LAMMPS virial vector to Voigt order

    virialVtoV[0] = 0;
    virialVtoV[1] = 1;
    virialVtoV[2] = 2;
    virialVtoV[3] = 5;
    virialVtoV[4] = 4;
    virialVtoV[5] = 3;

    // the following is for 6x6 matrix and vector indexing converter
    // this is clearly different order form albe[][] and revalbe[]
    // should not be used

    int indcounter = 0;
    for(int row = 0; row < NDIR_VIRIAL; row++)
      for(int col = row; col< NDIR_VIRIAL; col++) {
        voigt6MtoV[row][col] = voigt6MtoV[col][row] = indcounter;
        indcounter++;
      }
    // printf("Voigt6MtoV:\n");
    // for (int a = 0; a < NDIR_VIRIAL; a++) {
    //   for (int b = 0; b < NDIR_VIRIAL; b++) {
    //     printf("%d ", voigt6MtoV[a][b]);
    //   }
    //   printf("\n");
    // }

    // set up 3x3 kronecker deltas

    for(int row = 0; row < NXYZ_VIRIAL; row++)
      for(int col = 0; col < NXYZ_VIRIAL; col++)
        kronecker[row][col] = 0;
    for(int row = 0; row < NXYZ_VIRIAL; row++)
      kronecker[row][row] = 1;
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

    // Compute Born contribution

    if (pairflag) compute_pairs();
    if (bondflag) compute_bonds();
    if (angleflag) compute_angles();
    if (dihedflag) compute_dihedrals();

    // sum Born contributions over all procs

    MPI_Allreduce(values_local, values_global, nvalues, MPI_DOUBLE, MPI_SUM, world);

    // // convert to pressure units
    // // As discussed, it might be better to keep it as energy units.
    // // but this is to be defined

    // double nktv2p = force->nktv2p;
    // double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    // for (int m = 0; m < nvalues; m++) {
    //   values_global[m] *= (nktv2p * inv_volume);
    // }
  } else {

    // calculate Born matrix using stress finite differences
    compute_numdiff();

    // for consistency this is returned in energy units
    double inv_nktv2p = 1.0/force->nktv2p;
    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (int m = 0; m < nvalues; m++) {
      values_global[m] *= inv_nktv2p * volume;
    }

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

  // Declares born values

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

        a = 0;
        b = 0;
        c = 0;
        d = 0;
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
  int vec_indice;

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
  // It must be noted that, as stated in Yoshimoto's eq. 15, eq 16.
  // and eq. A3, this tensor is NOT the true Cijkl tensor.
  // We have the relationship
  // C_ijkl=1./4.(\hat{C}_ijkl+\hat{C}_jikl+\hat{C}_ijlk+\hat{C}_jilk)

  int flag, allflag;

  for (int idir = 0; idir < NDIR_VIRIAL; idir++) {
    displace_atoms(nall, idir, 1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_indice = revalbe[idir][jdir];
      values_global[vec_indice] = compute_virial->vector[virialVtoV[jdir]];
    }
    restore_atoms(nall, idir);

    displace_atoms(nall, idir, -1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_indice = revalbe[idir][jdir];
      values_global[vec_indice] -= compute_virial->vector[virialVtoV[jdir]];
    }

    // End of the strain
    restore_atoms(nall, idir);
  }

  // apply derivative factor

  double denominator = -0.5 / numdelta;
  for (int m = 0; m < nvalues; m++) values_global[m] *= denominator;

  // recompute virial so all virial and energy contributions are as before
  // also needed for virial stress addon contributions to Born matrix
  // this will possibly break compute stress/atom, need to test

  update_virial();

  virial_addon();

  // restore original forces for owned and ghost atoms

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++)
      f[i][k] = temp_f[i][k];

}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrix::displace_atoms(int nall, int idir, double magnitude)
{
  double **x = atom->x;

  // A.T.
  // this works for vector indices  7,  8,  9, 12, 14, 18 and 15, 16, 17
  // corresponding i,j indices     12, 13, 14, 23, 25, 36 and 26, 34, 35
  // int k = dirlist[idir][1];
  // int l = dirlist[idir][0];

  // A.T.
  // this works for vector indices  7,  8,  9, 12, 14, 18 and 10, 11, 13
  // corresponding i,j indices     12, 13, 14, 23, 25, 36 and 15, 16, 24
  // G.C.:
  // I see no difference with a 0 step simulation between both
  // methods.
  int k = dirlist[idir][0];
  int l = dirlist[idir][1];
  for (int i = 0; i < nall; i++)
      x[i][k] = temp_x[i][k] + numdelta * magnitude * (temp_x[i][l] - fixedpoint[l]);
}

/* ----------------------------------------------------------------------
   restore position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrix::restore_atoms(int nall, int idir)
{

  // reset all coords, just to be safe, ignore idir

  double **x = atom->x;
  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++)
      x[i][k] = temp_x[i][k];
}

/* ----------------------------------------------------------------------
   evaluate potential forces and virial
   same logic as in Verlet
------------------------------------------------------------------------- */

void ComputeBornMatrix::update_virial()
{
  int eflag = 0;
  // int vflag = VIRIAL_FDOTR; // Need to generalize this 
  int vflag = 1;

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
------------------------------------------------------------------------- */

void ComputeBornMatrix::virial_addon()
{
  // compute the contribution due to perturbation
  // here the addon parts are put into born
  // delta_il sigv_jk +  delta_ik sigv_jl +
  // delta_jl sigv_ik +  delta_jk sigv_il
  // Note: in calculation kl is all there from 0 to 6, and ij=(id,jd)
  //       each time there are six numbers passed for (Dijkl+Djikl)
  //       and the term I need should be div by 2.
  //       Job is to arrange the 6 numbers with ij indexing to the 21
  //       element data structure.
  // the sigv is the virial stress at current time. It is never touched.
  // Note the symmetry of (i-j), (k-n), and (ij, kn)
  // so we only need to evaluate 6x6 matrix with symmetry

  int kd, nd, id, jd;
  int m;

  double* sigv = compute_virial->vector;
  double modefactor[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};

  // Back to the ugly way
  // You can compute these factor by looking at
  // every Dijkl terms and adding the proper virials
  // Take into account the symmetries. For example:
  // B2323 = s33+D2323; B3232= s22+D3232;
  // but D3232=D2323 (computed in compute_numdiff)
  // and Cijkl = (Bijkl+Bjikl+Bijlk+Bjilk)/4. = (Bijkl+Bjilk)/2.
  // see Yoshimoto eq 15.and eq A3.
  values_global[0]  += 2.0*sigv[0];
  values_global[1]  += 2.0*sigv[1];
  values_global[2]  += 2.0*sigv[2];
  // values_global[3]  += 0.5*(sigv[1]+sigv[2]);
  // values_global[4]  += 0.5*(sigv[0]+sigv[2]);
  // values_global[5]  += 0.5*(sigv[0]+sigv[1]);
  values_global[3]  += sigv[2];
  values_global[4]  += sigv[2];
  values_global[5]  += sigv[1];
  values_global[6]  += 0.0;
  values_global[7]  += 0.0;
  values_global[8]  += 0.0;
  //  values_global[9]  += sigv[4];
  values_global[9]  += 2.0*sigv[4];
  //  values_global[10] += sigv[3];
  values_global[10] += 2.0*sigv[3];
  values_global[11] += 0.0;
  //  values_global[12] += sigv[5];
  values_global[12] += 2.0*sigv[5];
  values_global[13] += 0.0;
  //  values_global[14] += sigv[3];
  values_global[14] += 0.0;
  //  values_global[15] += sigv[5];
  values_global[15] += 0.0;
  //  values_global[16] += sigv[4];
  values_global[16] += 0.0;
  values_global[17] += 0.0;
  values_global[18] += 0.0;
  //  values_global[19] += sigv[4];
  values_global[19] += 0.0;
  values_global[20] += sigv[5];

  // This loop is actually bogus.
  //
  // for (int idir = 0; idir < NDIR_VIRIAL; idir++) {
  //   int ijvgt = idir; // this is it.
  //   double addon;

  //   // extract the two indices composing the voigt representation

  //   id  = voigt3VtoM[ijvgt][0];
  //   jd  = voigt3VtoM[ijvgt][1];

  //   for (int knvgt=ijvgt; knvgt < NDIR_VIRIAL; knvgt++) {
  //     kd = voigt3VtoM[knvgt][0];
  //     nd = voigt3VtoM[knvgt][1];
  //     addon = kronecker[id][nd]*sigv[virialMtoV[jd][kd]] +
  //             kronecker[id][kd]*sigv[virialMtoV[jd][nd]];
  //     if(id != jd)
  //     addon += kronecker[jd][nd]*sigv[virialMtoV[id][kd]] +
  //              kronecker[jd][kd]*sigv[virialMtoV[id][nd]];

  //     m = revalbe[ijvgt][knvgt];

  //     values_global[revalbe[ijvgt][knvgt]] += 0.5*modefactor[idir]*addon;
  //   }
  // }
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
  int i,m,n,nb,atom1,atom2,imol,iatom,btype,ivar;
  tagint tagprev;
  double dx,dy,dz,rsq;

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

  int a,b,c,d;
  double rij[3];
  double rinv, r2inv;
  double pair_pref, dupair, du2pair;

  // loop over all atoms and their bonds

  m = 0;
  while (m<nvalues) {

    for (atom1 = 0; atom1 < nlocal; atom1++) {
      if (!(mask[atom1] & groupbit)) continue;

      if (molecular == 1) nb = num_bond[atom1];
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
          atom2 = atom->map(onemols[imol]->bond_atom[iatom][i]+tagprev);
        }

        if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
        if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
        if (btype <= 0) continue;

        dx = x[atom2][0] - x[atom1][0];
        dy = x[atom2][1] - x[atom1][1];
        dz = x[atom2][2] - x[atom1][2];
        domain->minimum_image(dx,dy,dz);
        rsq = dx*dx + dy*dy + dz*dz;
        rij[0] = dx;
        rij[1] = dy;
        rij[2] = dz;
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        pair_pref = 0.0;
        dupair = 0.0;
        du2pair = 0.0;
        bond->born_matrix(btype,rsq,atom1,atom2,dupair,du2pair);
        pair_pref = du2pair - dupair*rinv;

        a = 0;
        b = 0;
        c = 0;
        d = 0;
        for (i = 0; i<21; i++) {
          a = albemunu[i][0];
          b = albemunu[i][1];
          c = albemunu[i][2];
          d = albemunu[i][3];
          values_local[m+i] += pair_pref*rij[a]*rij[b]*rij[c]*rij[d]*r2inv;
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
  int i,m,n,na,atom1,atom2,atom3,imol,iatom,atype,ivar;
  tagint tagprev;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,cost;
  double r1r2, r1r2inv;
  double rsq1inv,rsq2inv,r1inv,r2inv,cinv;
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

  int a,b,c,d,e,f;
  double duang, du2ang;
  double del1[3], del2[3];
  double dcos[6];
  double d2cos[21];
  double d2lncos[21];

  // Initializing array for intermediate cos derivatives
  // w regard to strain
  for (i = 0; i < 6; i++) {
    dcos[i] = 0;
  }
  for (i = 0; i < 21; i++) {
    d2cos[i] = 0;
    d2lncos[i] = 0;
  }

  m = 0;
  while (m < nvalues) {
    for (atom2 = 0; atom2 < nlocal; atom2++) {
      if (!(mask[atom2] & groupbit)) continue;

      if (molecular == 1) na = num_angle[atom2];
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
          atom1 = atom->map(onemols[imol]->angle_atom1[atom2][i]+tagprev);
          atom3 = atom->map(onemols[imol]->angle_atom3[atom2][i]+tagprev);
        }

        if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
        if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
        if (atype <= 0) continue;

        delx1 = x[atom1][0] - x[atom2][0];
        dely1 = x[atom1][1] - x[atom2][1];
        delz1 = x[atom1][2] - x[atom2][2];
        domain->minimum_image(delx1,dely1,delz1);
        del1[0] = delx1;
        del1[1] = dely1;
        del1[2] = delz1;

        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
        rsq1inv = 1.0/rsq1;
        r1 = sqrt(rsq1);
        r1inv = 1.0/r1;

        delx2 = x[atom3][0] - x[atom2][0];
        dely2 = x[atom3][1] - x[atom2][1];
        delz2 = x[atom3][2] - x[atom2][2];
        domain->minimum_image(delx2,dely2,delz2);
        del2[0] = delx2;
        del2[1] = dely2;
        del2[2] = delz2;

        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
        rsq2inv = 1.0/rsq2;
        r2 = sqrt(rsq2);
        r2inv = 1.0/r2;

        r1r2 = delx1*delx2 + dely1*dely2 + delz1*delz2;
        r1r2inv = 1/r1r2;
        // cost = cosine of angle

        cost = delx1*delx2 + dely1*dely2 + delz1*delz2;
        cost /= r1*r2;
        if (cost > 1.0) cost = 1.0;
        if (cost < -1.0) cost = -1.0;
        cinv = 1.0/cost;

        // The method must return derivative with regards
        // to cos(theta)!
        // Use the chain rule if needed:
        // dU(t)/de = dt/dcos(t)*dU(t)/dt*dcos(t)/de
        // with dt/dcos(t) = -1/sin(t)
        angle->born_matrix(atype,atom1,atom2,atom3,duang,du2ang);

        // Voigt notation
        // 1 = 11, 2 = 22, 3 = 33
        // 4 = 23, 5 = 13, 6 = 12
        a = 0;
        b = 0;
        c = 0;
        d = 0;
        for (i = 0; i<6; i++) {
          a = sigma_albe[i][0];
          b = sigma_albe[i][1];
          dcos[i] = cost*(del1[a]*del2[b]+del1[b]*del2[a]*r1r2inv -
                          del1[a]*del1[b]*rsq1inv - del2[a]*del2[b]*rsq2inv);
        }
        for (i = 0; i<21; i++) {
          a = albemunu[i][0];
          b = albemunu[i][1];
          c = albemunu[i][2];
          d = albemunu[i][3];
          e = C_albe[i][0];
          f = C_albe[i][1];
          d2lncos[i] = 2*(del1[a]*del1[b]*del1[c]*del1[d]*rsq1inv*rsq1inv +
                          del2[a]*del2[b]*del2[c]*del2[d]*rsq2inv*rsq2inv) -
                         (del1[a]*del2[b]+del1[b]*del2[a]) *
                         (del1[c]*del2[d]+del1[d]*del2[c]) *
                         r1r2inv*r1r2inv;
          d2cos[i] =  cost*d2lncos[i] + dcos[e]*dcos[f]*cinv;
          values_local[m+i] += duang*d2cos[i] + du2ang*dcos[e]*dcos[f];
        }
      }
    }
  m+=21;
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
  int i,m,n,nd,atom1,atom2,atom3,atom4,imol,iatom,dtype,ivar;
  tagint tagprev;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,ra2inv,rb2inv,rabinv;
  double si,co,phi;
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

  double dudih,du2dih;
  int a,b,c,d,e,f;
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
  double dmn[6];
  double dmm[6];
  double dnn[6];
  double d2mn[21];
  double d2mm[21];
  double d2nn[21];

  double dcos[6];
  double d2cos[21];

  for (i = 0; i < 6; i++) {
    dmn[i] =0;
    dmm[i] = 0;
    dnn[i] = 0;
    dcos[i] = 0;
  }
  for (i = 0; i < 21; i++) {
    d2mn[i] = 0;
    d2mm[i] = 0;
    d2nn[i] = 0;
    d2cos[i] = 0;
  }

  m = 0;
  while (m < nvalues) {
    for (atom2 = 0; atom2 < nlocal; atom2++) {
      if (!(mask[atom2] & groupbit)) continue;

      if (molecular == 1) nd = num_dihedral[atom2];
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
          atom1 = atom->map(onemols[imol]->dihedral_atom1[atom2][i]+tagprev);
          atom3 = atom->map(onemols[imol]->dihedral_atom3[atom2][i]+tagprev);
          atom4 = atom->map(onemols[imol]->dihedral_atom4[atom2][i]+tagprev);
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

        dihedral->born_matrix(nd,atom1,atom2,atom3,atom4,dudih,du2dih);

        vb1x = x[atom1][0] - x[atom2][0];
        vb1y = x[atom1][1] - x[atom2][1];
        vb1z = x[atom1][2] - x[atom2][2];
        domain->minimum_image(vb1x,vb1y,vb1z);
        b1[0] = vb1x;
        b1[1] = vb1y;
        b1[2] = vb1z;
        b1sq = b1[0]*b1[0]+b1[1]*b1[1]+b1[2]*b1[2];

        vb2x = x[atom3][0] - x[atom2][0];
        vb2y = x[atom3][1] - x[atom2][1];
        vb2z = x[atom3][2] - x[atom2][2];
        domain->minimum_image(vb2x,vb2y,vb2z);
        b2[0] = vb2x;
        b2[1] = vb2y;
        b2[2] = vb2z;
        b2sq = b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2];

        vb2xm = -vb2x;
        vb2ym = -vb2y;
        vb2zm = -vb2z;
        domain->minimum_image(vb2xm,vb2ym,vb2zm);

        vb3x = x[atom4][0] - x[atom3][0];
        vb3y = x[atom4][1] - x[atom3][1];
        vb3z = x[atom4][2] - x[atom3][2];
        domain->minimum_image(vb3x,vb3y,vb3z);
        b3[0] = vb3x;
        b3[1] = vb3y;
        b3[2] = vb3z;
        b3sq = b3[0]*b3[0]+b3[1]*b3[1]+b3[2]*b3[2];

        b1b2 = b1[0]*b2[0]+b1[1]*b2[1]+b1[2]*b2[2];
        b1b3 = b1[0]*b3[0]+b1[1]*b3[1]+b1[2]*b3[2];
        b2b3 = b2[0]*b3[0]+b2[1]*b3[1]+b2[2]*b3[2];

        ax = vb1y*vb2zm - vb1z*vb2ym;
        ay = vb1z*vb2xm - vb1x*vb2zm;
        az = vb1x*vb2ym - vb1y*vb2xm;
        bx = vb3y*vb2zm - vb3z*vb2ym;
        by = vb3z*vb2xm - vb3x*vb2zm;
        bz = vb3x*vb2ym - vb3y*vb2xm;

        rasq = ax*ax + ay*ay + az*az;
        rbsq = bx*bx + by*by + bz*bz;
        rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
        rg = sqrt(rgsq);

        ra2inv = rb2inv = 0.0;
        if (rasq > 0) ra2inv = 1.0/rasq;
        if (rbsq > 0) rb2inv = 1.0/rbsq;
        rabinv = sqrt(ra2inv*rb2inv);

        co = (ax*bx + ay*by + az*bz)*rabinv;
        si = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

        if (co > 1.0) co = 1.0;
        if (co < -1.0) co = -1.0;
        phi = atan2(si,co);

        // above a and b are m and n vectors
        // here they are integers indices
        a = 0;
        b = 0;
        c = 0;
        d = 0;
        e = 0;
        f = 0;
        for (i = 0; i<6; i++) {
          a = sigma_albe[i][0];
          b = sigma_albe[i][1];
          dmm[i] = 2*(b2sq*b1[a]*b1[b]+b1sq*b2[a]*b2[b] - b1b2*(b1[a]*b2[b]+b1[b]*b2[a]));
          dnn[i] = 2*(b3sq*b2[a]*b2[b]+b2sq*b3[a]*b3[b] - b2b3*(b2[a]*b3[b]+b2[b]*b3[a]));
          dmn[i] = b1b2*(b2[a]*b3[b]+b2[b]*b3[a]) + b2b3*(b1[a]*b2[b]+b1[b]*b2[a])
                 - 2*(b1b3*b2[a]*b2[b]) - b2sq*(b1[a]*b3[b]+b1[b]*b3[a]);
          dcos[i] = co*(rabinv*rabinv*dmn[i] - ra2inv*dmm[i] - rb2inv*dnn[i])/2.;
        }
        for (i = 0; i<21; i++) {
          a = albemunu[i][0];
          b = albemunu[i][1];
          c = albemunu[i][2];
          d = albemunu[i][3];
          e = C_albe[i][0];
          f = C_albe[i][1];
          d2mm[i] = 4*(b1[a]*b1[b]*b2[c]*b2[d] + b1[c]*b1[d]*b2[a]*b2[b])
                  - 8*(b1[a]*b2[b]+b1[b]*b2[a])*(b1[c]*b2[d]+b1[d]*b2[c]);
          d2nn[i] = 4*(b2[a]*b2[b]*b3[c]*b3[d] + b2[c]*b2[d]*b3[a]*b3[b])
                  - 8*(b2[a]*b3[b]+b2[b]*b3[a])*(b2[c]*b3[d]+b2[d]*b3[c]);
          d2mn[i] = (b1[a]*b2[b]+b1[b]*b2[a])*(b2[c]*b3[d]+b2[d]*b3[c])
                  + (b2[a]*b3[b]+b2[b]*b3[a])*(b1[c]*b2[d]+b1[d]*b2[d])
                  - (b1[a]*b3[b]+b1[b]*b3[a])*(b2[c]*b2[d]+b2[c]*b2[d])
                  - (b1[c]*b3[d]+b1[d]*b3[c])*(b2[a]*b2[b]+b2[a]*b2[b]);
          d2cos[i] = co/2.*(
                  rabinv*rabinv*d2mn[i]
                  - rabinv*rabinv*rabinv*rabinv*dmn[e]*dmn[f]
                  + ra2inv*ra2inv*dmm[e]*dmm[f]
                  - ra2inv*d2mm[i]
                  + rb2inv*rb2inv*dnn[e]*dnn[f]
                  - rb2inv*d2nn[i]);
          values_local[m+i] += dudih*d2cos[i] + du2dih*dcos[e]*dcos[f];
        }
      }
    }
  m+=21;
  }
}
