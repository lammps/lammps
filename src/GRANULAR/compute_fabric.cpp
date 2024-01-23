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
   Contributing authors: Joel Clemmer (SNL), Ishan Srivastava (LBNL)
------------------------------------------------------------------------- */

#include "compute_fabric.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { OTHER, GRANULAR };
enum { TYPE, RADIUS };
enum { CN, BR, FN, FT };

/* ---------------------------------------------------------------------- */

ComputeFabric::ComputeFabric(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), tensor_style(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute fabric command");

  if (strcmp(arg[3], "type") == 0)
    cutstyle = TYPE;
  else if (strcmp(arg[3], "radius") == 0)
    cutstyle = RADIUS;
  else
    error->all(FLERR, "Illegal compute fabric command");

  if (cutstyle == RADIUS && !atom->radius_flag)
    error->all(FLERR, "Compute fabric radius style requires atom attribute radius");

  // If optional arguments included, this will be oversized
  ntensors = narg - 4;
  tensor_style = new int[ntensors];

  cn_flag = 0;
  br_flag = 0;
  fn_flag = 0;
  ft_flag = 0;
  type_filter = nullptr;

  ntensors = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "contact") == 0) {
      cn_flag = 1;
      tensor_style[ntensors++] = CN;
    } else if (strcmp(arg[iarg], "branch") == 0) {
      br_flag = 1;
      tensor_style[ntensors++] = BR;
    } else if (strcmp(arg[iarg], "force/normal") == 0) {
      fn_flag = 1;
      tensor_style[ntensors++] = FN;
    } else if (strcmp(arg[iarg], "force/tangential") == 0) {
      ft_flag = 1;
      tensor_style[ntensors++] = FT;
    } else if (strcmp(arg[iarg], "type/include") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in compute fabric command");
      int ntypes = atom->ntypes;

      int i, j, itype, jtype;
      int inlo, inhi, jnlo, jnhi;
      if (!type_filter) {
        memory->create(type_filter, ntypes + 1, ntypes + 1, "compute/fabric:type_filter");

        for (i = 0; i <= ntypes; i++) {
          for (j = 0; j <= ntypes; j++) { type_filter[i][j] = 0; }
        }
      }

      std::vector<std::string> iwords = Tokenizer(arg[iarg + 1], ",").as_vector();
      std::vector<std::string> jwords = Tokenizer(arg[iarg + 2], ",").as_vector();
      for (const auto &ifield : iwords) {
        utils::bounds(FLERR, ifield, 1, ntypes, inlo, inhi, error);

        for (const auto &jfield : jwords) {
          utils::bounds(FLERR, jfield, 1, ntypes, jnlo, jnhi, error);

          for (itype = inlo; itype <= inhi; itype++) {
            for (jtype = jnlo; jtype <= jnhi; jtype++) {
              type_filter[itype][jtype] = 1;
              type_filter[jtype][itype] = 1;
            }
          }
        }
      }
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute fabric command");
    iarg++;
  }

  vector_flag = 1;
  size_vector = ntensors * 6;
  extvector = 0;

  scalar_flag = 1;
  extscalar = 1;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeFabric::~ComputeFabric()
{
  delete[] vector;
  delete[] tensor_style;
  memory->destroy(type_filter);
}

/* ---------------------------------------------------------------------- */

void ComputeFabric::init()
{
  if (force->pair == nullptr) error->all(FLERR, "No pair style is defined for compute fabric");
  if (force->pair->single_enable == 0 && (fn_flag || ft_flag))
    error->all(FLERR, "Pair style does not support compute fabric normal or tangential force");

  // Find if granular or gran
  pstyle = OTHER;
  if (force->pair_match("^granular", 0) || force->pair_match("^gran/", 0)) pstyle = GRANULAR;

  if (pstyle != GRANULAR && ft_flag)
    error->all(FLERR, "Pair style does not calculate tangential forces for compute fabric");

  if (force->pair->beyond_contact)
    error->all(FLERR, "Compute fabric does not support pair styles that extend beyond contact");

  // need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list (e.g. for granular pstyle)

  auto pairrequest = neighbor->find_request(force->pair);
  if (pairrequest && pairrequest->get_size())
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_OCCASIONAL);
  else
    neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeFabric::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeFabric::compute_vector()
{
  invoked_vector = update->ntimestep;

  int i, j, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r, rinv, rsq, radsum, fpair;

  double nx, ny, nz;
  double ncinv, denom, fn, ft, prefactor;
  double br_tensor[6], ft_tensor[6], fn_tensor[6];
  double trace_phi, trace_D, trace_Xfn, trace_Xft;
  double phi_ij[6] = {0.0};
  double Ac_ij[6] = {0.0};
  double D_ij[6] = {0.0};
  double Xfn_ij[6] = {0.0};
  double Xft_ij[6] = {0.0};
  double temp_dbl[6];

  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  // invoke compute_scalar() to update the number of contacts, if needed
  nc = compute_scalar();

  // If no contacts, everything will be zero
  if (nc == 0) {
    for (i = 0; i < size_vector; i++) vector[i] = 0.0;
    return;
  }
  ncinv = 1.0 / nc;

  // First loop through and calculate contact fabric tensor
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      jtype = type[j];

      if (type_filter)
        if (type_filter[itype][jtype] == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (cutstyle == TYPE) {
        if (rsq >= cutsq[itype][jtype]) continue;
      } else {
        radsum = radius[i] + radius[j];
        if (rsq >= radsum * radsum) continue;
      }

      r = sqrt(rsq);
      rinv = 1.0 / r;
      nx = delx * rinv;
      ny = dely * rinv;
      nz = delz * rinv;

      phi_ij[0] += nx * nx;
      phi_ij[1] += ny * ny;
      phi_ij[2] += nz * nz;
      phi_ij[3] += nx * ny;
      phi_ij[4] += nx * nz;
      phi_ij[5] += ny * nz;
    }
  }

  //Sum phi across processors
  MPI_Allreduce(phi_ij, temp_dbl, 6, MPI_DOUBLE, MPI_SUM, world);
  for (i = 0; i < 6; i++) phi_ij[i] = temp_dbl[i] * ncinv;

  trace_phi = (1.0 / 3.0) * (phi_ij[0] + phi_ij[1] + phi_ij[2]);

  Ac_ij[0] = (15.0 / 2.0) * (phi_ij[0] - trace_phi);
  Ac_ij[1] = (15.0 / 2.0) * (phi_ij[1] - trace_phi);
  Ac_ij[2] = (15.0 / 2.0) * (phi_ij[2] - trace_phi);
  Ac_ij[3] = (15.0 / 2.0) * (phi_ij[3]);
  Ac_ij[4] = (15.0 / 2.0) * (phi_ij[4]);
  Ac_ij[5] = (15.0 / 2.0) * (phi_ij[5]);

  // If needed, loop through and calculate other fabric tensors
  if (br_flag || fn_flag || ft_flag) {

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itag = tag[i];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        if (!(mask[j] & groupbit)) continue;

        // itag = jtag is possible for long cutoffs that include images of self

        if (newton_pair == 0 && j >= nlocal) {
          jtag = tag[j];
          if (itag > jtag) {
            if ((itag + jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag + jtag) % 2 == 1) continue;
          } else {
            if (x[j][2] < ztmp) continue;
            if (x[j][2] == ztmp) {
              if (x[j][1] < ytmp) continue;
              if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
            }
          }
        }

        jtype = type[j];

        if (type_filter)
          if (type_filter[itype][jtype] == 0) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (cutstyle == TYPE) {
          if (rsq >= cutsq[itype][jtype]) continue;
        } else {
          radsum = radius[i] + radius[j];
          if (rsq >= radsum * radsum) continue;
        }

        if (fn_flag || ft_flag) pair->single(i, j, itype, jtype, rsq, 1.0, 1.0, fpair);

        r = sqrt(rsq);
        rinv = 1.0 / r;
        nx = delx * rinv;
        ny = dely * rinv;
        nz = delz * rinv;

        denom = 1 + Ac_ij[0] * nx * nx + Ac_ij[1] * ny * ny + Ac_ij[2] * nz * nz;
        denom += 2 * Ac_ij[3] * nx * ny + 2 * Ac_ij[4] * nx * nz + 2 * Ac_ij[5] * ny * nz;
        prefactor = ncinv / denom;

        if (br_flag) {
          D_ij[0] += prefactor * nx * nx * r;
          D_ij[1] += prefactor * ny * ny * r;
          D_ij[2] += prefactor * nz * nz * r;
          D_ij[3] += prefactor * nx * ny * r;
          D_ij[4] += prefactor * nx * nz * r;
          D_ij[5] += prefactor * ny * nz * r;
        }

        if (fn_flag || ft_flag) {
          fn = r * fpair;

          Xfn_ij[0] += prefactor * nx * nx * fn;
          Xfn_ij[1] += prefactor * ny * ny * fn;
          Xfn_ij[2] += prefactor * nz * nz * fn;
          Xfn_ij[3] += prefactor * nx * ny * fn;
          Xfn_ij[4] += prefactor * nx * nz * fn;
          Xfn_ij[5] += prefactor * ny * nz * fn;

          if (ft_flag) {
            ft = force->pair->svector[3];

            Xft_ij[0] += prefactor * nx * nx * ft;
            Xft_ij[1] += prefactor * ny * ny * ft;
            Xft_ij[2] += prefactor * nz * nz * ft;
            Xft_ij[3] += prefactor * nx * ny * ft;
            Xft_ij[4] += prefactor * nx * nz * ft;
            Xft_ij[5] += prefactor * ny * nz * ft;
          }
        }
      }
    }
  }

  // Output results

  if (cn_flag) {
    for (i = 0; i < ntensors; i++) {
      if (tensor_style[i] == CN) {
        for (j = 0; j < 6; j++) vector[6 * i + j] = Ac_ij[j];
      }
    }
  }

  if (br_flag) {
    MPI_Allreduce(D_ij, temp_dbl, 6, MPI_DOUBLE, MPI_SUM, world);
    for (i = 0; i < 6; i++) D_ij[i] = temp_dbl[i];

    trace_D = (1.0 / 3.0) * (D_ij[0] + D_ij[1] + D_ij[2]);

    br_tensor[0] = (15.0 / (6.0 * trace_D)) * (D_ij[0] - trace_D);
    br_tensor[1] = (15.0 / (6.0 * trace_D)) * (D_ij[1] - trace_D);
    br_tensor[2] = (15.0 / (6.0 * trace_D)) * (D_ij[2] - trace_D);
    br_tensor[3] = (15.0 / (6.0 * trace_D)) * (D_ij[3]);
    br_tensor[4] = (15.0 / (6.0 * trace_D)) * (D_ij[4]);
    br_tensor[5] = (15.0 / (6.0 * trace_D)) * (D_ij[5]);

    for (i = 0; i < ntensors; i++) {
      if (tensor_style[i] == BR) {
        for (j = 0; j < 6; j++) vector[6 * i + j] = br_tensor[j];
      }
    }
  }

  if (fn_flag || ft_flag) {
    MPI_Allreduce(Xfn_ij, temp_dbl, 6, MPI_DOUBLE, MPI_SUM, world);
    for (i = 0; i < 6; i++) Xfn_ij[i] = temp_dbl[i];

    trace_Xfn = (1.0 / 3.0) * (Xfn_ij[0] + Xfn_ij[1] + Xfn_ij[2]);
  }

  if (fn_flag) {

    fn_tensor[0] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[0] - trace_Xfn);
    fn_tensor[1] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[1] - trace_Xfn);
    fn_tensor[2] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[2] - trace_Xfn);
    fn_tensor[3] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[3]);
    fn_tensor[4] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[4]);
    fn_tensor[5] = (15.0 / (6.0 * trace_Xfn)) * (Xfn_ij[5]);

    for (i = 0; i < ntensors; i++) {
      if (tensor_style[i] == FN) {
        for (j = 0; j < 6; j++) vector[6 * i + j] = fn_tensor[j];
      }
    }
  }

  if (ft_flag) {
    MPI_Allreduce(Xft_ij, temp_dbl, 6, MPI_DOUBLE, MPI_SUM, world);
    for (i = 0; i < 6; i++) Xft_ij[i] = temp_dbl[i];

    trace_Xft = (1.0 / 3.0) * (Xft_ij[0] + Xft_ij[1] + Xft_ij[2]);

    ft_tensor[0] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[0] - trace_Xft);
    ft_tensor[1] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[1] - trace_Xft);
    ft_tensor[2] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[2] - trace_Xft);
    ft_tensor[3] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[3]);
    ft_tensor[4] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[4]);
    ft_tensor[5] = (15.0 / (9.0 * trace_Xfn)) * (Xft_ij[5]);

    for (i = 0; i < ntensors; i++) {
      if (tensor_style[i] == FT) {
        for (j = 0; j < 6; j++) vector[6 * i + j] = ft_tensor[j];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeFabric::compute_scalar()
{
  // Skip if already calculated on this timestep
  if (invoked_scalar == update->ntimestep) return nc;

  invoked_scalar = update->ntimestep;

  int i, j, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, radsum, temp_dbl;

  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double **cutsq = force->pair->cutsq;

  // First loop through and calculate contact fabric tensor
  nc = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      jtype = type[j];

      if (type_filter)
        if (type_filter[itype][jtype] == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (cutstyle == TYPE) {
        if (rsq >= cutsq[itype][jtype]) continue;
      } else {
        radsum = radius[i] + radius[j];
        if (rsq >= radsum * radsum) continue;
      }

      nc += 1.0;
    }
  }
  //Count total contacts across processors
  MPI_Allreduce(&nc, &temp_dbl, 1, MPI_DOUBLE, MPI_SUM, world);
  nc = temp_dbl;

  scalar = nc;
  return nc;
}
