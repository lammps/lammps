// clang-format off
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
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "compute_continuum_chunk.h"

#include "arg_info.h"
#include "atom.h"
#include "citeme.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace NeighConst;

enum { OTHER, GRANULAR };
enum {
  NATOMS,
  DENSITY,
  VOLFRAC,
  MOMENTUM,
  VELOCITY,
  MGRAD,
  VGRAD,
  STRAINRATE,
  STRESS,
  STRESSKE,
  STRESSCON,
  IFD,
  FABRIC,
  TEMPERATURE
};
enum { BOUNDARY_NONE, BOUNDARY_FIX, BOUNDARY_ATOM, BOUNDARY_BOTH };

static const char cite_continuum[] =
    "Coarse-graining procedure: doi:10.1007/s10035-010-0181-z\n\n"
    "@Article{Goldhirsch2010,\n"
    " author = {Goldhirsch, Isaac},\n"
    " title = {{Stress, stress asymmetry and couple stress: From discrete particles to continuous fields}},\n"
    " journal = {Granular Matter},\n"
    " year =    2010,\n"
    " volume =  12,\n"
    " number =  3,\n"
    " pages =   {239--252}\n"
    "}\n\n";

static const char cite_boundary[] =
    "Boundary corrections: doi:10.1007/s10035-012-0317-4\n\n"
    "@Article{Weinhart2012,\n"
    " author = {Weinhart, Thomas and Thornton, Anthony R. and Luding, Stefan and Bokhove, Onno},\n"
    " title = {{From discrete particles to continuum fields near a boundary}},\n"
    " journal = {Granular Matter},\n"
    " year =    2012,\n"
    " volume =  14,\n"
    " number =  2,\n"
    " pages =   {289--294}\n"
    "}\n\n";

inline double ComputeContinuumChunk::calc_w(double r) const
{
  if (r > w_cut) {
    return 0.0;
  } else {
    return w_scale * exp(-(r * r) / (2.0 * w_sd_sq)) - w_offset;
  }
}

inline double ComputeContinuumChunk::calc_w_int(double *dr, double *rij) const
{
  double dr_sq = MathExtra::lensq3(dr);
  double dr_dot_rij = MathExtra::dot3(dr, rij);
  double rij_sq = MathExtra::lensq3(rij);

  double tmp = dr_dot_rij * dr_dot_rij - dr_sq * rij_sq + rij_sq * w_cut_sq;
  if (tmp < 0.0) return 0.0;
  tmp = sqrt(tmp);
  double smin = MAX(0.0, (dr_dot_rij - tmp) / rij_sq);
  double smax = MIN(1.0, (dr_dot_rij + tmp) / rij_sq);

  if ((smin >= 1.0) || (smax <= 0.0)) return 0.0;

  double rij_mag = sqrt(rij_sq);
  tmp = MY_SQRT2 * rij_mag * w_sd;
  double w_int = erf((rij_sq * smax - dr_dot_rij) / tmp) - erf((rij_sq * smin - dr_dot_rij) / tmp);
  w_int *= exp((dr_dot_rij * dr_dot_rij - dr_sq * rij_sq) / (tmp * tmp));
  w_int *= sqrt(0.5 * MY_PI) * w_sd / rij_mag;

  return w_scale * w_int - w_offset * (smax - smin);
}

/* ---------------------------------------------------------------------- */

ComputeContinuumChunk::ComputeContinuumChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), list(nullptr), delta(nullptr), ncoord(0), reducedflag(0),
    nlayers(nullptr), chunk_dim(nullptr), values_local(nullptr), values_global(nullptr),
    density_local(nullptr), density_global(nullptr), momentum_local(nullptr),
    momentum_global(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute continuum/chunk", error);

  w_cut = utils::numeric(FLERR, arg[4], false, lmp);
  w_sd = utils::numeric(FLERR, arg[5], false, lmp);
  if (w_cut <= 0.0) error->all(FLERR, "Illegal compute continuum/chunk cutoff value: {}", w_cut);
  if (w_sd <= 0.0) error->all(FLERR, "Illegal compute continuum/chunk width value: {}", w_sd);

  dim = domain->dimension;
  calculate_pair = 0;
  calculate_2_loops = 0;
  boundaryflag = BOUNDARY_NONE;
  boundary_group_flag = 0;
  boundary_groupbit = 0;
  radius_required = 0;
  index_density = -1;
  for (int a = 0; a < 3; a++) {
    index_momentum[a] = -1;
    index_velocity[a] = -1;
    for (int b = 0; b < 3; b++) index_vgrad[a][b] = -1;
  }

  int need_momentum = 0;
  int need_density = 0;
  int need_velocity = 0;
  int need_vgrad = 0;
  int need_radius = 0;

  int iarg = 6;
  values.clear();
  labels.clear();
  while (iarg < narg) {
    if (strcmp(arg[iarg], "natoms") == 0) {
      values.push_back(std::make_pair(NATOMS, -1));
      labels.push_back("natoms");
    } else if (strcmp(arg[iarg], "density") == 0) {
      values.push_back(std::make_pair(DENSITY, -1));
      labels.push_back("density");
      index_density = static_cast<int>(values.size()) - 1;
    } else if (strcmp(arg[iarg], "volume/fraction") == 0) {
      values.push_back(std::make_pair(VOLFRAC, -1));
      labels.push_back("volume/fraction");
      need_radius = 1;
    } else if (utils::strmatch(arg[iarg], "^momentum/.$")) {
      add_vector_component(arg[iarg], MOMENTUM);
    } else if (utils::strmatch(arg[iarg], "^velocity/.$")) {
      add_vector_component(arg[iarg], VELOCITY);
      need_density = 1;
      need_momentum = 1;
    } else if (utils::strmatch(arg[iarg], "^momentum/grad/")) {
      add_tensor_component(arg[iarg], MGRAD);
      need_momentum = 1;
      calculate_2_loops = 1;
    } else if (utils::strmatch(arg[iarg], "^velocity/grad/")) {
      add_tensor_component(arg[iarg], VGRAD);
      need_density = 1;
      need_momentum = 1;
      need_velocity = 1;
      calculate_2_loops = 1;
    } else if (utils::strmatch(arg[iarg], "^strain/rate/")) {
      add_tensor_component(arg[iarg], STRAINRATE);
      need_density = 1;
      need_momentum = 1;
      need_velocity = 1;
      need_vgrad = 1;
      calculate_2_loops = 1;
    } else if (utils::strmatch(arg[iarg], "^stress/.$") ||
               utils::strmatch(arg[iarg], "^stress/..$")) {
      add_tensor_component(arg[iarg], STRESS);
      calculate_pair = 1;
      need_radius = 1;
      need_density = 1;
      need_momentum = 1;
      calculate_2_loops = 1;
    } else if (utils::strmatch(arg[iarg], "^stress/ke/")) {
      add_tensor_component(arg[iarg], STRESSKE);
      need_density = 1;
      need_momentum = 1;
      calculate_2_loops = 1;
    } else if (utils::strmatch(arg[iarg], "^stress/contacts/")) {
      add_tensor_component(arg[iarg], STRESSCON);
      calculate_pair = 1;
      need_radius = 1;
    } else if (utils::strmatch(arg[iarg], "^boundary/force")) {
      add_vector_component(arg[iarg], IFD);
      calculate_pair = 1;
      need_radius = 1;
    } else if (utils::strmatch(arg[iarg], "^fabric/")) {
      add_tensor_component(arg[iarg], FABRIC);
      calculate_pair = 1;
      need_radius = 1;
    } else if (strcmp(arg[iarg], "temperature") == 0) {
      values.push_back(std::make_pair(TEMPERATURE, -1));
      labels.push_back("temperature");
      need_density = 1;
      need_momentum = 1;
      calculate_2_loops = 1;
    } else {
      break;
    }
    iarg++;
  }

  // Add any necessary intermediate values, won't be saved as output
  nskip = 0;
  if ((need_density) && (index_density == -1)) {
    values.push_back(std::make_pair(DENSITY, -1));
    index_density = static_cast<int>(values.size()) - 1;
    labels.push_back("density/internal");
    nskip += 1;
  }

  if (need_momentum) {
    for (int a = 0; a < dim; a++) {
      if (index_momentum[a] == -1) {
        values.push_back(std::make_pair(MOMENTUM, a));
        index_momentum[a] = static_cast<int>(values.size()) - 1;
        labels.push_back("momentum/internal");
        nskip += 1;
      }
    }
  }

  if (need_velocity) {
    for (int a = 0; a < dim; a++) {
      if (index_velocity[a] == -1) {
        values.push_back(std::make_pair(VELOCITY, a));
        index_velocity[a] = static_cast<int>(values.size()) - 1;
        labels.push_back("velocity/internal");
        nskip += 1;
      }
    }
  }

  if (need_vgrad) {
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if ((dim == 2) && ((b == 2) || (a == 2))) continue;
        if (index_vgrad[a][b] == -1) {
          values.push_back(std::make_pair(VGRAD, a * 3 + b));
          index_vgrad[a][b] = static_cast<int>(values.size()) - 1;
          labels.push_back("vgrad/internal");
          nskip += 1;
        }
      }
    }
  }

  nvalues = static_cast<int>(values.size());
  if (nvalues == 0) error->all(FLERR, "No values in compute continuum/chunk command");

  while (iarg < narg) {
    if (strcmp(arg[iarg], "boundary/fix") == 0) {
      if (boundaryflag == BOUNDARY_ATOM)
        boundaryflag = BOUNDARY_BOTH;
      else
        boundaryflag = BOUNDARY_FIX;
      iarg += 1;
    } else if (strcmp(arg[iarg], "boundary/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("compute continuum/chunk ") + arg[iarg], error);
      if (boundaryflag == BOUNDARY_FIX)
        boundaryflag = BOUNDARY_BOTH;
      else
        boundaryflag = BOUNDARY_ATOM;
      boundary_group_flag = 1;
      boundary_groupbit = group->get_bitmask_by_id(FLERR, arg[iarg + 1], "compute continuum/chunk");
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown compute continuum/chunk keyword: {}", arg[iarg]);
    }
  }

  if ((boundaryflag == BOUNDARY_NONE)) {
    for (auto &val : values)
      if (val.first == IFD)
        error->all(FLERR, "Must specify boundary/atom and/or boundary/fix to compute boundary/force");
  }

  int which = cchunk->get_which();
  if (which == ArgInfo::BIN1D) {
    bin_dim = 1;
  } else if (which == ArgInfo::BIN2D) {
    bin_dim = 2;
  } else if (which == ArgInfo::BIN3D) {
    bin_dim = 3;
  } else {
    error->all(FLERR, "Can only use bin chunk/atom styles with compute continuum/chunk");
  }

  w_cut_sq = w_cut * w_cut;
  w_sd_sq = w_sd * w_sd;

  // Normalization factor for truncated Gaussian
  double exp_cut = exp(-w_cut_sq / (2.0 * w_sd_sq));
  if (bin_dim == 1) {
    w_scale = sqrt(2.0 * MY_PI) * w_sd * erf(w_cut / (MY_SQRT2 * w_sd));
    w_scale -= 2.0 * w_cut * exp_cut;
    w_scale = 1.0 / w_scale;
  } else if (bin_dim == 2) {
    w_scale = -0.5 * w_cut_sq * exp_cut + w_sd_sq * (1.0 - exp_cut);
    w_scale = 1.0 / (2.0 * MY_PI * w_scale);
  } else {
    w_scale = -THIRD * w_cut * exp_cut * (w_cut_sq + 3.0 * w_sd_sq);
    w_scale += sqrt(0.5 * MY_PI) * w_sd_sq * w_sd * erf(w_cut / (MY_SQRT2 * w_sd));
    w_scale = 1.0 / (4.0 * MY_PI * w_scale);
  }
  w_offset = w_scale * exp(-w_cut_sq / (2.0 * w_sd_sq));

  if (need_radius && !atom->radius_flag)
    error->all(FLERR, "Compute continuum/chunk requires atom attribute radius");
  radius_required = need_radius;

  array_flag = 1;
  size_array_cols = nvalues - nskip;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;
  thermo_modify_colname = 1;

  allocate();

  if (lmp->citeme) {
    lmp->citeme->add(cite_continuum);
    if (boundaryflag) lmp->citeme->add(cite_boundary);
  }
}

/* ---------------------------------------------------------------------- */

ComputeContinuumChunk::~ComputeContinuumChunk()
{
  memory->destroy(array);
  memory->destroy(values_local);
  memory->destroy(values_global);
  memory->destroy(density_local);
  memory->destroy(density_global);
  memory->destroy(momentum_local);
  memory->destroy(momentum_global);
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::init()
{
  ComputeChunk::init();

  if ((boundaryflag == BOUNDARY_FIX) || (boundaryflag == BOUNDARY_BOTH)) {
    auto wall_fixes = modify->get_fix_by_style("wall/gran");
    if (wall_fixes.size() == 0)
      error->all(FLERR, "Could not find any instances of fix wall/gran for boundary corrections");
    for (auto fix : wall_fixes)
      if (!fix->peratom_flag)
        error->all(FLERR, "Must use contacts keyword in fix wall/gran {} for boundary corrections",
                   fix->id);
  }

  if (calculate_pair) {
    if (force->pair == nullptr)
      error->all(FLERR, "No pair style is defined for compute continuum/chunk stress calculation");
    if (force->pair->single_enable == 0)
      error->all(FLERR, "Pair style does not support compute continuum/chunk stress calculation");

    // Find if granular or gran, need to include tangential forces
    pstyle = OTHER;
    if (force->pair_match("^granular", 0) || force->pair_match("^gran/", 0)) pstyle = GRANULAR;

    auto *pairrequest = neighbor->find_request(force->pair);
    if (pairrequest && pairrequest->get_size())
      neighbor->add_request(this, REQ_SIZE | REQ_OCCASIONAL | REQ_FULL);
    else
      neighbor->add_request(this, REQ_OCCASIONAL | REQ_FULL);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::init_list(int /*id*/, NeighList *ptr)
{
  if (calculate_pair) list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::compute_array()
{
  int i, j, m, mtmp;

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  build_stencil();
  int *cdim = cchunk->get_dim();
  double *chunk_delta = cchunk->get_delta();
  double **coord = cchunk->coord;
  int chunk_ncoord = cchunk->ncoord;
  int chunk_reducedflag = cchunk->get_reducedflag();

  for (m = 0; m < nchunk; m++) {
    for (i = 0; i < nvalues; i++) values_local[m][i] = 0.0;
    for (i = 0; i < size_array_cols; i++) array[m][i] = 0.0;
  }

  int a, b, itype, style, component, field_index, jboundary;
  double w, wc, mi, voli, rsq_atom_bin, rsq_cont_bin, rsq_pair, r_pair;
  double f_norm, w_int_tmp;
  double coordx[3], xbin0[3], xbin[3], xbin2[3], xcont[3], f_pair[3], f_wall[3];
  double dx_pair[3], dx_pair_filtered[3], dx_atom_bin[3], dx_bin_cont[3], dx_atom_cont[3];

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int jj, jnum;
  int *jlist, *numneigh, **firstneigh;
  double **array_atom_fix;

  if (calculate_pair) {
    neighbor->build_one(list);
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }

  Pair *pair = force->pair;
  auto wall_fixes = modify->get_fix_by_style("wall/gran");

  for (i = 0; i < nlocal; i++) {
    if ((mask[i] & groupbit) && (ichunk[i] > 0)) {
      m = ichunk[i] - 1;

      if (boundary_group_flag && (mask[i] & boundary_groupbit)) continue;

      MathExtra::copy3(x[i], xbin0);
      for (a = 0; a < chunk_ncoord; a++) {
        if (chunk_reducedflag) {
          domain->lamda2x(coord[m], coordx);
          xbin0[cdim[a]] = coordx[a];
        } else {
          xbin0[cdim[a]] = coord[m][a];
        }
      }

      itype = type[i];
      if (rmass)
        mi = rmass[i];
      else
        mi = mass[itype];
      voli = 0.0;
      if (radius_required) {
        voli = MY_PI * radius[i] * radius[i];
        if (dim == 3) voli *= 4.0 * THIRD * radius[i];
      }

      for (auto &stencil_offset : stencil) {
        xbin[0] = xbin0[0] + stencil_offset.dx[0];
        xbin[1] = xbin0[1] + stencil_offset.dx[1];
        xbin[2] = xbin0[2] + stencil_offset.dx[2];

        mtmp = shifted_bin(m, stencil_offset.dn);
        if (mtmp == -1) continue;

        MathExtra::sub3(x[i], xbin, dx_atom_bin);
        rsq_atom_bin = MathExtra::lensq3(dx_atom_bin);
        w = calc_w(sqrt(rsq_atom_bin));

        field_index = 0;
        for (auto &val : values) {
          style = val.first;
          component = val.second;

          a = component % 3;
          b = (component - a) / 3;

          if (style == NATOMS) {
            values_local[mtmp][field_index] += 1.0;
          } else if (style == DENSITY) {
            values_local[mtmp][field_index] += mi * w;
          } else if (style == VOLFRAC) {
            values_local[mtmp][field_index] += voli * w;
          } else if (style == MOMENTUM) {
            values_local[mtmp][field_index] += mi * v[i][component] * w;
          }

          if (boundaryflag && ((style == STRESS) || (style == STRESSCON))) {
            for (auto wall_fix : wall_fixes) {
              array_atom_fix = wall_fix->array_atom;

              if (array_atom_fix[i][0] < 0.5) continue;
              f_wall[0] = array_atom_fix[i][1];
              f_wall[1] = array_atom_fix[i][2];
              f_wall[2] = array_atom_fix[i][3];
              xcont[0] = array_atom_fix[i][4];
              xcont[1] = array_atom_fix[i][5];
              xcont[2] = array_atom_fix[i][6];

              MathExtra::sub3(x[i], xcont, dx_atom_cont);
              w_int_tmp = calc_w_int(dx_atom_bin, dx_atom_cont);

              values_local[mtmp][field_index] -= f_wall[a] * dx_atom_cont[b] * w_int_tmp;
            }
          }

          field_index++;
        }

        if (calculate_pair) {
          jlist = firstneigh[i];
          jnum = numneigh[i];
          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;

            if (!(mask[j] & groupbit)) continue;

            if (boundary_group_flag && (mask[j] & boundary_groupbit))
              jboundary = 1;
            else
              jboundary = 0;

            MathExtra::sub3(x[i], x[j], dx_pair);
            rsq_pair = MathExtra::lensq3(dx_pair);
            r_pair = sqrt(rsq_pair);
            pair->single(i, j, itype, type[j], rsq_pair, 1.0, 1.0, f_norm);

            MathExtra::scale3(f_norm / r_pair, dx_pair, f_pair);
            if (pstyle == GRANULAR) {
              f_pair[0] += force->pair->svector[0];
              f_pair[1] += force->pair->svector[1];
              f_pair[2] += force->pair->svector[2];
            }

            if (MathExtra::lensq3(f_pair) == 0.0) continue;

            if (jboundary) {
              MathExtra::add3(x[i], x[j], xcont);
              MathExtra::scaleadd3((radius[j] - radius[i]) / r_pair, dx_pair, xcont, xcont);
              MathExtra::scale3(0.5, xcont);

              MathExtra::copy3(xcont, xbin2);
              for (a = 0; a < chunk_ncoord; a++) xbin2[cdim[a]] = xbin[cdim[a]];
              MathExtra::sub3(xbin2, xcont, dx_bin_cont);

              rsq_cont_bin = MathExtra::lensq3(dx_bin_cont);
              wc = calc_w(sqrt(rsq_cont_bin));

              MathExtra::sub3(x[i], xcont, dx_atom_cont);
              w_int_tmp = calc_w_int(dx_atom_bin, dx_atom_cont);
            } else {
              MathExtra::zero3(dx_pair_filtered);
              for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                dx_pair_filtered[cdim[coord_index]] = dx_pair[cdim[coord_index]];
              w_int_tmp = calc_w_int(dx_atom_bin, dx_pair_filtered);
            }

            field_index = 0;
            for (auto &val : values) {
              style = val.first;
              component = val.second;

              a = component % 3;
              b = (component - a) / 3;

              // Factors of 1/2 b/c this is a full nlist EXCEPT for boundary which filters if one atom is a boundary
              if ((style == STRESS) || (style == STRESSCON)) {
                if (jboundary) {
                  values_local[mtmp][field_index] -= 0.5 * f_pair[a] * dx_atom_cont[b] * w_int_tmp;
                } else {
                  values_local[mtmp][field_index] -= 0.5 * f_pair[a] * dx_pair[b] * w_int_tmp;
                }
              } else if (style == IFD) {
                if (!jboundary) continue;
                values_local[mtmp][field_index] -= f_pair[a] * wc;
              } else if (style == FABRIC) {
                if (jboundary) continue;
                values_local[mtmp][field_index] +=
                    0.5 * voli * dx_pair[a] * dx_pair[b] * w_int_tmp / rsq_pair;
              }

              field_index++;
            }
          }
        }
      }
    }
  }

  if (calculate_2_loops) {
    for (m = 0; m < nchunk; m++) {
      density_local[m] = values_local[m][index_density];
      for (a = 0; a < 3; a++) {
        if (a < dim)
          momentum_local[m][a] = values_local[m][index_momentum[a]];
        else
          momentum_local[m][a] = 0.0;
      }
    }

    MPI_Allreduce(&density_local[0], &density_global[0], nchunk, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&momentum_local[0][0], &momentum_global[0][0], nchunk * 3, MPI_DOUBLE, MPI_SUM,
                  world);

    double dtemp, vtemp[3];
    for (i = 0; i < nlocal; i++) {
      if ((mask[i] & groupbit) && (ichunk[i] > 0)) {
        m = ichunk[i] - 1;

        MathExtra::copy3(x[i], xbin0);
        for (a = 0; a < chunk_ncoord; a++) {
          if (chunk_reducedflag) {
            domain->lamda2x(coord[m], coordx);
            xbin0[cdim[a]] = coordx[a];
          } else {
            xbin0[cdim[a]] = coord[m][a];
          }
        }

        for (auto &stencil_offset : stencil) {
          xbin[0] = xbin0[0] + stencil_offset.dx[0];
          xbin[1] = xbin0[1] + stencil_offset.dx[1];
          xbin[2] = xbin0[2] + stencil_offset.dx[2];

          mtmp = shifted_bin(m, stencil_offset.dn);
          if (mtmp == -1) continue;

          MathExtra::sub3(x[i], xbin, dx_atom_bin);
          rsq_atom_bin = MathExtra::lensq3(dx_atom_bin);
          w = calc_w(sqrt(rsq_atom_bin));

          dtemp = density_global[m];
          MathExtra::copy3(momentum_global[m], vtemp);
          if (dtemp != 0.0) MathExtra::scale3(1.0 / dtemp, vtemp);

          MathExtra::sub3(vtemp, v[i], vtemp);
          itype = type[i];
          if (rmass)
            mi = rmass[i];
          else
            mi = mass[itype];

          field_index = 0;
          for (auto &val : values) {
            style = val.first;
            component = val.second;

            a = component % 3;
            b = (component - a) / 3;

            if (style == TEMPERATURE) {
              values_local[mtmp][field_index] += 0.5 * mi * MathExtra::lensq3(vtemp) * w;
            } else if ((style == STRESS) || (style == STRESSKE)) {
              values_local[mtmp][field_index] -= mi * vtemp[a] * vtemp[b] * w;
            }

            field_index++;
          }
        }
      }
    }
  }

  MPI_Allreduce(&values_local[0][0], &values_global[0][0], nchunk * nvalues, MPI_DOUBLE, MPI_SUM,
                world);

  // Calculate trivially derived values, in the order used
  // velocity
  double dtemp;
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = val.first;
      component = val.second;

      if (style == VELOCITY) {
        dtemp = values_global[m][index_density];
        if (dtemp != 0.0) values_global[m][field_index] = values_global[m][index_momentum[component]] / dtemp;
      }

      field_index++;
    }
  }

  // gradients
  int shift[3], mp, mm, ac;
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = val.first;
      component = val.second;

      if ((style != MGRAD) && (style != VGRAD)) {
        field_index++;
        continue;
      }

      a = component % 3;
      b = (component - a) / 3;

      ac = -1;
      for (int c = 0; c < ncoord; c++)
        if (chunk_dim[c] == a) ac = c;

      if (ac == -1) {
        values_global[m][field_index] = 0.0;
        field_index++;
        continue;
      }

      shift[0] = shift[1] = shift[2] = 0;
      shift[ac] = 1;
      mp = shifted_bin(m, shift);
      shift[ac] = -1;
      mm = shifted_bin(m, shift);

      if ((mp == -1) || (mm == -1)) {
        values_global[m][field_index] = 0.0;
        field_index++;
        continue;
      }

      if (style == MGRAD) {
        values_global[m][field_index] =
            (values_global[mp][index_momentum[b]] - values_global[mm][index_momentum[b]]) /
            (2.0 * chunk_delta[a]);
      } else if (style == VGRAD) {
        values_global[m][field_index] =
            (values_global[mp][index_velocity[b]] - values_global[mm][index_velocity[b]]) /
            (2.0 * chunk_delta[a]);
      }

      field_index++;
    }
  }

  // strain rate
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = val.first;
      component = val.second;
      a = component % 3;
      b = (component - a) / 3;

      if (style == STRAINRATE)
        values_global[m][field_index] =
            0.5 * (values_global[m][index_vgrad[a][b]] + values_global[m][index_vgrad[b][a]]);

      field_index++;
    }
  }

  // Normalize by any unused dimensions
  if (bin_dim != dim) {
    int unused_dim[3] = {1, 1, 1};
    for (a = 0; a < ncoord; a++) unused_dim[cdim[a]] = 0;

    for (a = 0; a < dim; a++)
      if (unused_dim[a])
        for (m = 0; m < nchunk; m++)
          for (int n = 0; n < nvalues; n++)
            values_global[m][n] /= domain->prd[a];
  }

  for (m = 0; m < nchunk; m++)
    for (i = 0; i < size_array_cols; i++) array[m][i] = values_global[m][i];
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(array);
  memory->destroy(values_local);
  memory->destroy(values_global);
  memory->destroy(density_local);
  memory->destroy(density_global);
  memory->destroy(momentum_local);
  memory->destroy(momentum_global);

  maxchunk = nchunk;
  memory->create(array, maxchunk, size_array_cols, "continuum/chunk:array");
  memory->create(values_local, maxchunk, nvalues, "continuum/chunk:values_local");
  memory->create(values_global, maxchunk, nvalues, "continuum/chunk:values_global");
  if (calculate_2_loops) {
    memory->create(density_local, maxchunk, "continuum/chunk:density_local");
    memory->create(density_global, maxchunk, "continuum/chunk:density_global");
    memory->create(momentum_local, maxchunk, 3, "continuum/chunk:momentum_local");
    memory->create(momentum_global, maxchunk, 3, "continuum/chunk:momentum_global");
  }
}

/* ---------------------------------------------------------------------- */

double ComputeContinuumChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += 2.0 * maxchunk * nvalues * sizeof(double);
  bytes += (double) maxchunk * size_array_cols * sizeof(double);
  if (calculate_2_loops) {
    bytes += 2.0 * maxchunk * sizeof(double);
    bytes += 2.0 * maxchunk * 3 * sizeof(double);
  }
  return bytes;
}

/* ---------------------------------------------------------------------- */

std::string ComputeContinuumChunk::get_thermo_colname(int m)
{
  if ((m >= 0) && (m < size_array_cols)) return labels[m];
  return {};
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::add_tensor_component(char *option, int variable)
{
  if (((std::string) option).back() == '*') {
    std::vector<std::string> suffices = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if ((dim == 2) && ((b == 2) || (a == 2))) continue;
        values.push_back(std::make_pair(variable, a * 3 + b));
        labels.push_back(trimmed_option + suffices[a * 3 + b]);
        if (variable == VGRAD) index_vgrad[a][b] = static_cast<int>(values.size()) - 1;
      }
    }
  } else {
    int index = -1;
    int dim_error = 0;

    if (utils::strmatch(option, "xx$")) {
      index = 0;
    } else if (utils::strmatch(option, "xy$")) {
      index = 1;
    } else if (utils::strmatch(option, "xz$")) {
      index = 2;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "yx$")) {
      index = 3;
    } else if (utils::strmatch(option, "yy$")) {
      index = 4;
    } else if (utils::strmatch(option, "yz$")) {
      index = 5;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zx$")) {
      index = 6;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zy$")) {
      index = 7;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zz$")) {
      index = 8;
      if (dim == 2) dim_error = 1;
    } else {
      error->all(FLERR, "Invalid compute continuum/chunk property {}", option);
    }

    if (dim_error) error->all(FLERR, "Invalid compute continuum/chunk property {} in 2D", option);

    values.push_back(std::make_pair(variable, index));
    labels.push_back(option);
    if (variable == VGRAD) {
      int a = index % 3;
      int b = (index - a) / 3;
      index_vgrad[a][b] = static_cast<int>(values.size()) - 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::add_vector_component(char *option, int variable)
{
  if (((std::string) option).back() == '*') {
    std::vector<std::string> suffices = {"x", "y", "z"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < dim; a++) {
      values.push_back(std::make_pair(variable, a));
      labels.push_back(trimmed_option + suffices[a]);
      if (variable == MOMENTUM) index_momentum[a] = static_cast<int>(values.size()) - 1;
      if (variable == VELOCITY) index_velocity[a] = static_cast<int>(values.size()) - 1;
    }
  } else {
    int index = -1;
    if (utils::strmatch(option, "x$")) {
      index = 0;
    } else if (utils::strmatch(option, "y$")) {
      index = 1;
    } else if (utils::strmatch(option, "z$")) {
      if (dim == 2)
        error->all(FLERR, "Invalid compute continuum/chunk property {} in 2D", option);
      index = 2;
    } else {
      error->all(FLERR, "Invalid compute continuum/chunk property {}", option);
    }

    values.push_back(std::make_pair(variable, index));
    labels.push_back(option);
    if (variable == MOMENTUM) index_momentum[index] = static_cast<int>(values.size()) - 1;
    if (variable == VELOCITY) index_velocity[index] = static_cast<int>(values.size()) - 1;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeContinuumChunk::shifted_bin(int origin_bin, int *dn) const
{
  int x[3] = {0, 0, 0};

  if (ncoord == 1) {
    x[0] = origin_bin + dn[0];
  } else if (ncoord == 2) {
    x[0] = origin_bin / nlayers[1] + dn[0];
    x[1] = origin_bin % nlayers[1] + dn[1];
  } else if (ncoord == 3) {
    x[0] = origin_bin / (nlayers[1] * nlayers[2]) + dn[0];
    x[1] = (origin_bin / nlayers[2]) % nlayers[1] + dn[1];
    x[2] = origin_bin % nlayers[2] + dn[2];
  }

  for (int a = 0; a < ncoord; a++) {
    if (!domain->periodicity[chunk_dim[a]]) {
      if ((x[a] < 0) || (x[a] >= nlayers[a])) return -1;
      continue;
    }
    while (x[a] < 0) x[a] += nlayers[a];
    while (x[a] >= nlayers[a]) x[a] -= nlayers[a];
  }

  int new_bin = -1;
  if (ncoord == 1) {
    new_bin = x[0];
  } else if (ncoord == 2) {
    new_bin = x[0] * nlayers[1] + x[1];
  } else if (ncoord == 3) {
    new_bin = x[0] * nlayers[1] * nlayers[2] + x[1] * nlayers[2] + x[2];
  }

  if ((new_bin < 0) || (new_bin >= nchunk))
    error->one(FLERR, "Bad chunk index {} shifted by {} {} {}", origin_bin, dn[0], dn[1], dn[2]);

  return new_bin;
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::build_stencil()
{
  stencil.clear();
  int stencil_size[3] = {0, 0, 0};
  double width[3] = {0.0, 0.0, 0.0};

  nlayers = cchunk->get_nlayers();
  delta = cchunk->get_delta();
  chunk_dim = cchunk->get_dim();
  ncoord = cchunk->ncoord;
  reducedflag = cchunk->get_reducedflag();

  for (int a = 0; a < ncoord; a++) {
    width[a] = delta[a];
    if (reducedflag) width[a] *= domain->prd[chunk_dim[a]];
    stencil_size[a] = static_cast<int>(ceil(w_cut / width[a]));
  }

  for (int dn0 = -stencil_size[0]; dn0 <= stencil_size[0]; dn0++) {
    for (int dn1 = -stencil_size[1]; dn1 <= stencil_size[1]; dn1++) {
      for (int dn2 = -stencil_size[2]; dn2 <= stencil_size[2]; dn2++) {
        StencilOffset offset;
        offset.dx[0] = 0.0;
        offset.dx[1] = 0.0;
        offset.dx[2] = 0.0;
        if (ncoord >= 1) offset.dx[chunk_dim[0]] = dn0 * width[0];
        if (ncoord >= 2) offset.dx[chunk_dim[1]] = dn1 * width[1];
        if (ncoord >= 3) offset.dx[chunk_dim[2]] = dn2 * width[2];

        double r_sq = MathExtra::lensq3(offset.dx);
        if (r_sq <= w_cut_sq) {
          offset.dn[0] = dn0;
          offset.dn[1] = dn1;
          offset.dn[2] = dn2;
          stencil.push_back(offset);
        }
      }
    }
  }
}
