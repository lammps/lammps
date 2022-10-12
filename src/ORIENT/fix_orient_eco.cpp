// clang-format off
/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratdir_veces
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Adrian A. Schratt and Volker Mohles (ICAMS)
------------------------------------------------------------------------- */

#include "fix_orient_eco.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static const char cite_fix_orient_eco[] =
  "fix orient/eco command: doi:j.commatsci.2020.109774\n\n"
  "@Article{Schratt20,\n"
  " author = {A. A. Schratt and V. Mohles},\n"
  " title = {Efficient Calculation of the {ECO} Driving Force for Atomistic Simulations of Grain Boundary Motion},\n"
  " journal = {Computational Materials Science},\n"
  " volume = {182},\n"
  " year = {2020},\n"
  " pages = {109774},\n"
  " doi = {j.commatsci.2020.109774},\n"
  " url = {https://doi.org/10.1016/j.commatsci.2020.109774}\n"
  "}\n\n";

struct FixOrientECO::Nbr {
  double duchi;            // potential derivative
  double real_phi[2][3];   // real part of wave function
  double imag_phi[2][3];   // imaginary part of wave function
};

/* ---------------------------------------------------------------------- */

FixOrientECO::FixOrientECO(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    dir_filename(nullptr), order(nullptr), nbr(nullptr), list(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_orient_eco);

  MPI_Comm_rank(world, &me);

  if (narg != 7) error->all(FLERR, "Illegal fix orient/eco command");

  scalar_flag = 1;          // computes scalar
  global_freq = 1;          // values can be computed at every timestep
  extscalar = 1;            // scalar scales with # of atoms
  peratom_flag = 1;         // quantities are per atom quantities
  size_peratom_cols = 2;    // # of per atom quantities
  peratom_freq = 1;
  energy_global_flag = 1;

  // parse input parameters

  u_0 = utils::numeric(FLERR, arg[3],false,lmp);
  sign = (u_0 >= 0.0 ? 1 : -1);
  eta = utils::numeric(FLERR, arg[4],false,lmp);
  r_cut = utils::numeric(FLERR, arg[5],false,lmp);

  // read reference orientations from file
  // work on rank 0 only

  dir_filename = utils::strdup(arg[6]);
  if (me == 0) {
    char line[IMGMAX];
    char *result;
    int count;

    FILE *infile = utils::open_potential(dir_filename,lmp,nullptr);
    if (infile == nullptr)
      error->one(FLERR,"Cannot open fix orient/eco file {}: {}",
                                   dir_filename, utils::getsyserror());
    for (int i = 0; i < 6; ++i) {
      result = fgets(line, IMGMAX, infile);
      if (!result) error->one(FLERR, "Fix orient/eco file read failed");
      count = sscanf(line, "%lg %lg %lg", &dir_vec[i][0], &dir_vec[i][1], &dir_vec[i][2]);
      if (count != 3) error->one(FLERR, "Fix orient/eco file read failed");
    }
    fclose(infile);

    // calculate reciprocal lattice vectors
    get_reciprocal();

    squared_cutoff = r_cut * r_cut;
    inv_squared_cutoff = 1.0 / squared_cutoff;
    half_u = 0.5 * u_0;
    inv_eta = 1.0 / eta;
  }

  // initializations
  MPI_Bcast(&dir_vec[0][0], 18, MPI_DOUBLE, 0, world);                  // communicate direct lattice vectors
  MPI_Bcast(&reciprocal_vectors[0][0][0], 18, MPI_DOUBLE, 0, world);    // communicate reciprocal vectors
  MPI_Bcast(&squared_cutoff, 1, MPI_DOUBLE, 0, world);                  // communicate squared cutoff radius
  MPI_Bcast(&inv_squared_cutoff, 1, MPI_DOUBLE, 0, world);              // communicate inverse squared cutoff radius
  MPI_Bcast(&half_u, 1, MPI_DOUBLE, 0, world);                          // communicate half potential energy
  MPI_Bcast(&inv_eta, 1, MPI_DOUBLE, 0, world);                         // communicate inverse threshold

  // set comm size needed by this Fix
  if (u_0 != 0) comm_forward = sizeof(Nbr) / sizeof(double);

  added_energy = 0.0;  // initial energy

  nmax = atom->nmax;
  nbr = (Nbr *) memory->smalloc(nmax * sizeof(Nbr), "orient/eco:nbr");
  memory->create(order, nmax, 2, "orient/eco:order");
  array_atom = order;

  // zero the array since a variable may access it before first run
  for (int i = 0; i < atom->nlocal; ++i) order[i][0] = order[i][1] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixOrientECO::~FixOrientECO() {
  memory->destroy(order);
  memory->sfree(nbr);
  delete[] dir_filename;
}

/* ---------------------------------------------------------------------- */

int FixOrientECO::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::init() {
  // get this processors rank
  MPI_Comm_rank(world, &me);

  // compute normalization factor
  int neigh = get_norm();
  if (me == 0)
    utils::logmesg(lmp,"  fix orient/eco: cutoff={} norm_fac={} neighbors={}\n",
                   r_cut, norm_fac, neigh);

  inv_norm_fac = 1.0 / norm_fac;

  // check parameters
  if (r_cut > force->pair->cutforce) {
    error->all(FLERR, "Cutoff radius used by fix orient/eco must be smaller than force cutoff");
  }

  // communicate normalization factor
  MPI_Bcast(&norm_fac, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&inv_norm_fac, 1, MPI_DOUBLE, 0, world);

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }

  // need a full perpetual neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::init_list(int /* id */, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::setup(int vflag) {
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::post_force(int /* vflag */) {
  // set local variables
  int ii, i, jj, j;                                             // loop variables and atom IDs
  int k;                                                        // variable to loop over 3 reciprocal directions
  int lambda;                                                   // variable to loop over 2 crystals
  int dim;                                                      // variable to loop over 3 spatial components
  double dx, dy, dz;                                            // stores current interatomic vector
  double squared_distance;                                      // stores current squared distance
  double chi;                                                   // stores current order parameter
  double weight;                                                // stores current weight function
  double scalar_product;                                        // stores current scalar product
  double omega;                                                 // phase of sine transition
  double omega_pre = MY_PI2 * inv_eta;                          // prefactor for omega
  double duchi_pre = half_u * MY_PI * inv_eta * inv_norm_fac;   // prefactor for duchi
  double sin_om;                                                // stores the value of the sine transition

  // initialize added energy at this step
  added_energy = 0.0;

  // set local pointers and neighbor lists
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;

  int inum = list->inum;
  int jnum;
  int *ilist = list->ilist;
  int *jlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // insure nbr and order data structures are adequate size
  if (nall > nmax) {
    nmax = nall;
    memory->destroy(nbr);
    memory->destroy(order);
    nbr = (Nbr *) memory->smalloc(nmax * sizeof(Nbr), "orient/eco:nbr");
    memory->create(order, nmax, 2, "orient/eco:order");
    array_atom = order;
  }

  // loop over owned atoms and compute order parameter
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // initializations
    chi = 0.0;
    for (k = 0; k < 3; ++k) {
      nbr[i].real_phi[0][k] = nbr[i].real_phi[1][k] = 0.0;
      nbr[i].imag_phi[0][k] = nbr[i].imag_phi[1][k] = 0.0;
    }

    // loop over all neighbors of atom i
    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // vector difference
      dx = x[i][0] - x[j][0];
      dy = x[i][1] - x[j][1];
      dz = x[i][2] - x[j][2];
      squared_distance = dx * dx + dy * dy + dz * dz;

      if (squared_distance < squared_cutoff) {
        squared_distance *= inv_squared_cutoff;
        weight = squared_distance * (squared_distance - 2.0) + 1.0;

        for (lambda = 0; lambda < 2; ++lambda) {
          for (k = 0; k < 3; ++k) {
            scalar_product = reciprocal_vectors[lambda][k][0] * dx + reciprocal_vectors[lambda][k][1] * dy + reciprocal_vectors[lambda][k][2] * dz;
            nbr[i].real_phi[lambda][k] += weight * cos(scalar_product);
            nbr[i].imag_phi[lambda][k] += weight * sin(scalar_product);
          }
        }
      }
    }

    // collect contributions
    for (k = 0; k < 3; ++k) {
      chi += (nbr[i].real_phi[0][k] * nbr[i].real_phi[0][k] + nbr[i].imag_phi[0][k] * nbr[i].imag_phi[0][k] -
                nbr[i].real_phi[1][k] * nbr[i].real_phi[1][k] - nbr[i].imag_phi[1][k] * nbr[i].imag_phi[1][k]);
    }
    chi *= inv_norm_fac;
    order[i][0] = chi;

    // compute normalized order parameter
    // and potential energy
    if (chi > eta) {
      added_energy += half_u;
      nbr[i].duchi = 0.0;
      order[i][1] = sign;
    } else if (chi < -eta) {
      added_energy -= half_u;
      nbr[i].duchi = 0.0;
      order[i][1] = -sign;
    } else {
      omega = omega_pre * chi;
      sin_om = sin(omega);

      added_energy += half_u * sin_om;
      nbr[i].duchi = duchi_pre * cos(omega);
      order[i][1] = sign * sin_om;
    }

    // compute product with potential derivative
    for (k = 0; k < 3; ++k) {
      for (lambda = 0; lambda < 2; ++lambda) {
        nbr[i].real_phi[lambda][k] *= nbr[i].duchi;
        nbr[i].imag_phi[lambda][k] *= nbr[i].duchi;
      }
    }
  }

  // set additional local variables
  double gradient_ii_cos[2][3][3];      // gradient ii cosine term
  double gradient_ii_sin[2][3][3];      // gradient ii sine term
  double gradient_ij_vec[2][3][3];      // gradient ij vector term
  double gradient_ij_sca[2][3];         // gradient ij scalar term
  double weight_gradient_prefactor;     // gradient prefactor
  double weight_gradient[3];            // gradient of weight
  double cos_scalar_product;            // cosine of scalar product
  double sin_scalar_product;            // sine of scalar product
  double gcos_scalar_product;           // gradient weight function * cosine of scalar product
  double gsin_scalar_product;           // gradient weight function * sine of scalar product

  // compute force only if synthetic
  // potential is not zero
  if (u_0 != 0.0) {
    // communicate to acquire nbr data for ghost atoms
    comm->forward_comm(this);

    // loop over all atoms
    for (ii = 0; ii < inum; ++ii) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      const bool no_boundary_atom = (nbr[i].duchi == 0.0);

      // skip atoms not in group
      if (!(mask[i] & groupbit)) continue;

      // initializations
      for (k = 0; k < 3; ++k) {
        for (lambda = 0; lambda < 2; ++lambda) {
          for (dim = 0; dim < 3; ++dim) {
            gradient_ii_cos[lambda][k][dim] = 0.0;
            gradient_ii_sin[lambda][k][dim] = 0.0;
            gradient_ij_vec[lambda][k][dim] = 0.0;
          }
          gradient_ij_sca[lambda][k] = 0.0;
        }
      }
      // loop over all neighbors of atom i
      // for those within squared_cutoff, compute force
      for (jj = 0; jj < jnum; ++jj) {
        j = jlist[jj];
        j &= NEIGHMASK;

        // do not compute force on atom i if it is far from boundary
        if ((nbr[j].duchi == 0.0) && no_boundary_atom) continue;

        // vector difference
        dx = x[i][0] - x[j][0];
        dy = x[i][1] - x[j][1];
        dz = x[i][2] - x[j][2];
        squared_distance = dx * dx + dy * dy + dz * dz;

        if (squared_distance < squared_cutoff) {
          // compute force on atom i
          // need weight and its gradient
          squared_distance *= inv_squared_cutoff;
          weight = squared_distance * (squared_distance - 2.0) + 1.0;
          weight_gradient_prefactor = 4.0 * (squared_distance - 1.0) * inv_squared_cutoff;
          weight_gradient[0] = weight_gradient_prefactor * dx;
          weight_gradient[1] = weight_gradient_prefactor * dy;
          weight_gradient[2] = weight_gradient_prefactor * dz;

          // (1) compute scalar product and sine and cosine of it
          // (2) compute product of sine and cosine with gradient of weight function
          // (3) compute gradient_ii_cos and gradient_ii_sin by summing up result of (2)
          // (4) compute gradient_ij_vec and gradient_ij_sca
          for (lambda = 0; lambda < 2; ++lambda) {
            for (k = 0; k < 3; ++k) {
              scalar_product = reciprocal_vectors[lambda][k][0] * dx + reciprocal_vectors[lambda][k][1] * dy + reciprocal_vectors[lambda][k][2] * dz;
              cos_scalar_product = cos(scalar_product);
              sin_scalar_product = sin(scalar_product);
              for (dim = 0; dim < 3; ++dim) {
                gradient_ii_cos[lambda][k][dim] += (gcos_scalar_product = weight_gradient[dim] * cos_scalar_product);
                gradient_ii_sin[lambda][k][dim] += (gsin_scalar_product = weight_gradient[dim] * sin_scalar_product);
                gradient_ij_vec[lambda][k][dim] += (nbr[j].real_phi[lambda][k] * gcos_scalar_product - nbr[j].imag_phi[lambda][k] * gsin_scalar_product);
              }
              gradient_ij_sca[lambda][k] += weight * (nbr[j].real_phi[lambda][k] * sin_scalar_product + nbr[j].imag_phi[lambda][k] * cos_scalar_product);
            }
          }
        }
      }

      // sum contributions
      for (k = 0; k < 3; ++k) {
        for (dim = 0; dim < 3; ++dim) {
          f[i][dim] -= (nbr[i].real_phi[0][k] * gradient_ii_cos[0][k][dim] + nbr[i].imag_phi[0][k] * gradient_ii_sin[0][k][dim] + gradient_ij_vec[0][k][dim] + reciprocal_vectors[1][k][dim] * gradient_ij_sca[1][k]);
          f[i][dim] += (nbr[i].real_phi[1][k] * gradient_ii_cos[1][k][dim] + nbr[i].imag_phi[1][k] * gradient_ii_sin[1][k][dim] + gradient_ij_vec[1][k][dim] + reciprocal_vectors[0][k][dim] * gradient_ij_sca[0][k]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::post_force_respa(int vflag, int ilevel, int /* iloop */) {
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixOrientECO::compute_scalar() {
  double added_energy_total;
  MPI_Allreduce(&added_energy, &added_energy_total, 1, MPI_DOUBLE, MPI_SUM, world);
  return added_energy_total;
}

/* ---------------------------------------------------------------------- */

int FixOrientECO::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int */*pbc*/) {
  int i, j;
  int m = 0;

  for (i = 0; i < n; ++i) {
    j = list[i];
    *((Nbr*)&buf[m]) = nbr[j];
    m += sizeof(Nbr)/sizeof(double);
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixOrientECO::unpack_forward_comm(int n, int first, double *buf) {
  int i;
  int last = first + n;
  int m = 0;

  for (i = first; i < last; ++i) {
    nbr[i] = *((Nbr*) &buf[m]);
    m += sizeof(Nbr) / sizeof(double);
  }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double FixOrientECO::memory_usage() {
  double bytes = (double)nmax * sizeof(Nbr);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
 reciprocal lattice vectors from file input
 ------------------------------------------------------------------------- */

void FixOrientECO::get_reciprocal() {
  // volume of unit cell A
  double vol = 0.5 * (dir_vec[0][0] * (dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2]) + dir_vec[1][0] * (dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2]) + dir_vec[2][0] * (dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2])) / MY_PI;
  double i_vol = 1.0 / vol;

  // grain A: reciprocal_vectors 0
  reciprocal_vectors[0][0][0] = (dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2]) * i_vol;
  reciprocal_vectors[0][0][1] = (dir_vec[1][2] * dir_vec[2][0] - dir_vec[2][2] * dir_vec[1][0]) * i_vol;
  reciprocal_vectors[0][0][2] = (dir_vec[1][0] * dir_vec[2][1] - dir_vec[2][0] * dir_vec[1][1]) * i_vol;

  // grain A: reciprocal_vectors 1
  reciprocal_vectors[0][1][0] = (dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2]) * i_vol;
  reciprocal_vectors[0][1][1] = (dir_vec[2][2] * dir_vec[0][0] - dir_vec[0][2] * dir_vec[2][0]) * i_vol;
  reciprocal_vectors[0][1][2] = (dir_vec[2][0] * dir_vec[0][1] - dir_vec[0][0] * dir_vec[2][1]) * i_vol;

  // grain A: reciprocal_vectors 2
  reciprocal_vectors[0][2][0] = (dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2]) * i_vol;
  reciprocal_vectors[0][2][1] = (dir_vec[0][2] * dir_vec[1][0] - dir_vec[1][2] * dir_vec[0][0]) * i_vol;
  reciprocal_vectors[0][2][2] = (dir_vec[0][0] * dir_vec[1][1] - dir_vec[1][0] * dir_vec[0][1]) * i_vol;

  // volume of unit cell B
  vol = 0.5 * (dir_vec[3][0] * (dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2]) + dir_vec[4][0] * (dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2]) + dir_vec[5][0] * (dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2])) / MY_PI;
  i_vol = 1.0 / vol;

  // grain B: reciprocal_vectors 0
  reciprocal_vectors[1][0][0] = (dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2]) * i_vol;
  reciprocal_vectors[1][0][1] = (dir_vec[4][2] * dir_vec[5][0] - dir_vec[5][2] * dir_vec[4][0]) * i_vol;
  reciprocal_vectors[1][0][2] = (dir_vec[4][0] * dir_vec[5][1] - dir_vec[5][0] * dir_vec[4][1]) * i_vol;

  // grain B: reciprocal_vectors 1
  reciprocal_vectors[1][1][0] = (dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2]) * i_vol;
  reciprocal_vectors[1][1][1] = (dir_vec[5][2] * dir_vec[3][0] - dir_vec[3][2] * dir_vec[5][0]) * i_vol;
  reciprocal_vectors[1][1][2] = (dir_vec[5][0] * dir_vec[3][1] - dir_vec[3][0] * dir_vec[5][1]) * i_vol;

  // grain B: reciprocal_vectors 2
  reciprocal_vectors[1][2][0] = (dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2]) * i_vol;
  reciprocal_vectors[1][2][1] = (dir_vec[3][2] * dir_vec[4][0] - dir_vec[4][2] * dir_vec[3][0]) * i_vol;
  reciprocal_vectors[1][2][2] = (dir_vec[3][0] * dir_vec[4][1] - dir_vec[4][0] * dir_vec[3][1]) * i_vol;
}

/* ----------------------------------------------------------------------
 normalization factor
 ------------------------------------------------------------------------- */

int FixOrientECO::get_norm() {
  // set up local variables
  double delta[3];                        // relative position
  double squared_distance;                           // squared distance of atoms i and j
  double weight;                           // weight function for atoms i and j
  double wsum = 0.0;                    // sum of all weight functions
  double scalar_product;                          // scalar product
  double reesum[3] = {0.0, 0.0, 0.0};   // sum of real part
  double imesum[3] = {0.0, 0.0, 0.0};   // sum of imaginary part

  int max_co = 4;                       // will produce wrong results for rcut > 3 * lattice constant

  int neigh = 0;                            // count number of neighbors used

  // loop over ideal lattice positions
  int i, k, idx[3];
  for (idx[0] = -max_co; idx[0] <= max_co; ++idx[0]) {
    for (idx[1] = -max_co; idx[1] <= max_co; ++idx[1]) {
      for (idx[2] = -max_co; idx[2] <= max_co; ++idx[2]) {
        // distance of atoms
        for (i = 0; i < 3; ++i) {
          delta[i] = dir_vec[0][i] * idx[0] + dir_vec[1][i] * idx[1] + dir_vec[2][i] * idx[2];
        }
        squared_distance = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

        // check if atom is within cutoff region
        if ((squared_distance != 0.0) && (squared_distance < squared_cutoff)) {
          ++neigh;
          squared_distance *= inv_squared_cutoff;

          // weight
          weight = squared_distance * (squared_distance - 2.0) + 1.0;
          wsum += weight;

          // three reciprocal directions
          for (k = 0; k < 3; ++k) {
            scalar_product = reciprocal_vectors[1][k][0] * delta[0] + reciprocal_vectors[1][k][1] * delta[1] + reciprocal_vectors[1][k][2] * delta[2];
            reesum[k] += weight * cos(scalar_product);
            imesum[k] -= weight * sin(scalar_product);
          }
        }
      }
    }
  }

  // compute normalization
  norm_fac = 3.0 * wsum * wsum;
  for (k = 0; k < 3; ++k) {
    norm_fac -= reesum[k] * reesum[k] + imesum[k] * imesum[k];
  }
  return neigh;
}



