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

static constexpr double EPSILON = 1.0e-8;

static const char cite_continuum[] = "Coarse-graining procedure: doi:10.1007/s10035-010-0181-z\n\n"
                                     "@Article{Goldhirsch2010,\n"
                                     " author = {Goldhirsch, Isaac},\n"
                                     " title = {{Stress, stress asymmetry and couple stress: From "
                                     "discrete particles to continuous fields}},\n"
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

  // In case atoms are only separated along a dimension which is not being binned
  if (rij_sq < EPSILON * w_sd_sq) return calc_w(sqrt(dr_sq));

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
    ComputeChunk(lmp, narg, arg), list(nullptr), nlayers(nullptr), chunk_dim(nullptr),
    delta(nullptr), values_local(nullptr), values_global(nullptr), density_local(nullptr),
    density_global(nullptr), momentum_local(nullptr), momentum_global(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute continuum/chunk", error);

  w_cut = utils::numeric(FLERR, arg[4], false, lmp);
  w_sd = utils::numeric(FLERR, arg[5], false, lmp);
  if (w_cut <= 0.0) error->all(FLERR, 4, "Illegal compute continuum/chunk cutoff value: {}", w_cut);
  if (w_sd <= 0.0) error->all(FLERR, 5, "Illegal compute continuum/chunk width value: {}", w_sd);

  ncoord = 0;
  reducedflag = 0;
  dim = domain->dimension;
  calculate_pair = 0;
  calculate_2_loops = 0;
  boundaryflag = BOUNDARY_NONE;
  boundary_group_flag = 0;
  boundary_groupbit = 0;
  radius_required = 0;
  index_density = -1;
  pstyle = OTHER;
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
      no_norm.insert(values.size());
      values.emplace_back(std::make_tuple(NATOMS, 0, 0));
      labels.emplace_back("natoms");
    } else if (strcmp(arg[iarg], "density") == 0) {
      values.emplace_back(std::make_tuple(DENSITY, 0, 0));
      labels.emplace_back("density");
      index_density = static_cast<int>(values.size()) - 1;
    } else if (strcmp(arg[iarg], "volume/fraction") == 0) {
      values.emplace_back(std::make_tuple(VOLFRAC, 0, 0));
      labels.emplace_back("volume/fraction");
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
      values.emplace_back(std::make_tuple(TEMPERATURE, 0, 0));
      labels.emplace_back("temperature");
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
    values.emplace_back(std::make_tuple(DENSITY, 0, 0));
    index_density = static_cast<int>(values.size()) - 1;
    labels.emplace_back("density/internal");
    nskip += 1;
  }

  if (need_momentum) {
    for (int a = 0; a < dim; a++) {
      if (index_momentum[a] == -1) {
        values.emplace_back(std::make_tuple(MOMENTUM, 1, a));
        index_momentum[a] = static_cast<int>(values.size()) - 1;
        labels.emplace_back("momentum/internal");
        nskip += 1;
      }
    }
  }

  if (need_velocity) {
    for (int a = 0; a < dim; a++) {
      if (index_velocity[a] == -1) {
        values.emplace_back(std::make_tuple(VELOCITY, 1, a));
        index_velocity[a] = static_cast<int>(values.size()) - 1;
        labels.emplace_back("velocity/internal");
        nskip += 1;
      }
    }
  }

  if (need_vgrad) {
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if ((dim == 2) && ((b == 2) || (a == 2))) continue;
        if (index_vgrad[a][b] == -1) {
          values.emplace_back(std::make_tuple(VGRAD, 2, a * 3 + b));
          index_vgrad[a][b] = static_cast<int>(values.size()) - 1;
          labels.emplace_back("vgrad/internal");
          nskip += 1;
        }
      }
    }
  }

  nvalues = static_cast<int>(values.size());
  if (nvalues == 0) error->all(FLERR, 6, "No values in compute continuum/chunk command");

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
      error->all(FLERR, iarg, "Unknown compute continuum/chunk keyword: {}", arg[iarg]);
    }
  }

  if (boundaryflag == BOUNDARY_NONE) {
    for (auto &val : values)
      if (std::get<0>(val) == IFD)
        error->all(FLERR,
                   "Must specify how boundary/force is calculated: "
                   "using boundary/atom and/or boundary/fix");
  }

  int which = cchunk->get_which();
  if (which == ArgInfo::BIN1D) {
    bin_dim = 1;
  } else if (which == ArgInfo::BIN2D) {
    bin_dim = 2;
  } else if (which == ArgInfo::BIN3D) {
    bin_dim = 3;
  } else {
    error->all(FLERR, 3, "Can only use bin chunk/atom styles with compute continuum/chunk");
  }

  // Can't use discard if bound is set, so this is stricter than needed but assumed by position_to_bin()
  if (cchunk->get_discard() == 0)
    error->all(FLERR, "The compute chunk/atom discard no option is not supported");
  if (cchunk->compress)
    error->all(FLERR, "The compute chunk/atom compress option is not supported");
  if (cchunk->get_limit())
    error->all(FLERR, "The compute chunk/atom limit option is not supported");

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
    error->all(FLERR, Error::NOLASTLINE, "Compute continuum/chunk requires atom attribute radius");
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
    if (wall_fixes.empty())
      error->all(FLERR, Error::NOLASTLINE,
                 "Could not find any instances of fix wall/gran for boundary corrections");
    for (auto *fix : wall_fixes)
      if (!fix->peratom_flag)
        error->all(FLERR, Error::NOLASTLINE,
                   "Must use contacts keyword in fix wall/gran {} for boundary corrections",
                   fix->id);
  }

  if (calculate_pair) {
    if (force->pair == nullptr)
      error->all(FLERR, Error::NOLASTLINE,
                 "No pair style is defined for compute continuum/chunk stress calculation");
    if (force->pair->single_enable == 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Pair style does not support compute continuum/chunk stress calculation");

    // Find if granular or gran, need to include tangential forces
    if (force->pair_match("^granular", 0) || force->pair_match("^gran/", 0)) pstyle = GRANULAR;

    // As in pair/local, create occasional list vs using actual pair list (could be half/full)
    //   Note, will need to update if any new granular pair styles use history with a full list
    auto *pairrequest = neighbor->find_request(force->pair);
    if (pairrequest && pairrequest->get_size())
      neighbor->add_request(this, REQ_SIZE | REQ_OCCASIONAL);
    else
      neighbor->add_request(this, REQ_OCCASIONAL);
  }

  if (domain->triclinic)
    error->all(FLERR, "Compute continuum/chunk does not support triclinic simulation boxes");

  // compute chunk does not wrap bins across periodic boundaries
  double *chunk_delta = cchunk->get_delta();
  int chunk_ncoord = cchunk->ncoord;
  int chunk_reducedflag = cchunk->get_reducedflag();
  chunk_dim = cchunk->get_dim();

  double *prd = domain->prd;
  int *periodicity = domain->periodicity;
  if (chunk_reducedflag) {
    for (int a = 0; a < chunk_ncoord; a++)
      if (periodicity[chunk_dim[a]]) {
        double nbins = 1.0 / chunk_delta[a];
        if (std::fabs(nbins - std::round(nbins)) > EPSILON)
          error->warning(FLERR,
                         "Bins do not evenly divide the simulation box,"
                         " results on the boundary may be incorrect");
      }
  } else {
    for (int a = 0; a < chunk_ncoord; a++)
      if (periodicity[chunk_dim[a]]) {
        double nbins = prd[chunk_dim[a]] / chunk_delta[a];
        if (std::fabs(nbins - std::round(nbins)) > EPSILON)
          error->warning(FLERR,
                         "Bins do not evenly divide the simulation box,"
                         " results on the boundary may be incorrect");
      }
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
  int nlocal = atom->nlocal;

  ComputeChunk::compute_array();

  int *ichunk = cchunk->ichunk;

  build_stencil();
  int *cdim = cchunk->get_dim();
  double *chunk_delta = cchunk->get_delta();
  double **coord = cchunk->coord;
  int chunk_ncoord = cchunk->ncoord;
  int chunk_reducedflag = cchunk->get_reducedflag();

  // Check lower bound of bin aligns with domain in periodic directions
  const double *boxlo = reducedflag ? domain->boxlo_lamda : domain->boxlo;

  for (int a = 0; a < ncoord; a++) {
    const int idim = chunk_dim[a];

    if (!domain->periodicity[idim]) continue;

    // coord[0][a] is the center of the first bin along this coordinate.
    const double offset = cchunk->coord[0][a] - 0.5 * delta[a];
    if (std::fabs(boxlo[idim] - offset) > EPSILON)
      error->warning(FLERR,
                     "Bins do not start at lower edge of simulation box."
                     " Results on the boundary may be incorrect");
  }

  for (m = 0; m < nchunk; m++) {
    for (i = 0; i < nvalues; i++) values_local[m][i] = 0.0;
    for (i = 0; i < size_array_cols; i++) array[m][i] = 0.0;
  }

  int a = 0;
  int b = 0;
  int itype, style, vtype, component, field_index, iboundary, jboundary;
  double w, wc, massi, voli, volj, rsq_atom_bin, rsq_cont_bin, rsq_pair, r_pair, r_cont;
  double f_norm, w_int_tmp, factor_lj;
  double xbin0[3], xbinc[3], xbin[3], xbin2[3], xcont[3], f_pair[3], f_wall[3], dx_pair[3],
      xj_near[3];
  double dx_pair_filtered[3], dx_atom_bin[3], dx_bin_cont[3], dx_atom_cont[3];
  double dx_atom_cont_filtered[3];
  double **array_atom_fix;

  double pair_stencil_reach;
  int mc, stencil_size[3];

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  std::unordered_set<int> visited_bins;

  auto wall_fixes = modify->get_fix_by_style("wall/gran");

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    m = ichunk[i] - 1;
    if (m < 0) continue;

    if (boundary_group_flag && (mask[i] & boundary_groupbit)) continue;

    if (chunk_reducedflag) {
      double lamda[3];
      domain->x2lamda(x[i], lamda);
      for (a = 0; a < chunk_ncoord; a++) lamda[cdim[a]] = coord[m][a];
      domain->lamda2x(lamda, xbin0);
    } else {
      MathExtra::copy3(x[i], xbin0);
      for (a = 0; a < chunk_ncoord; a++) xbin0[cdim[a]] = coord[m][a];
    }

    itype = type[i];
    if (rmass)
      massi = rmass[i];
    else
      massi = mass[itype];
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
        style = std::get<0>(val);
        vtype = std::get<1>(val);
        component = std::get<2>(val);

        if (vtype == 1) {
          a = component;
        } else {
          a = component / 3;
          b = component % 3;
        }

        if (style == NATOMS) {
          if (rsq_atom_bin < w_cut_sq) values_local[mtmp][field_index] += 1.0;
        } else if (style == DENSITY) {
          values_local[mtmp][field_index] += massi * w;
        } else if (style == VOLFRAC) {
          values_local[mtmp][field_index] += voli * w;
        } else if (style == MOMENTUM) {
          values_local[mtmp][field_index] += massi * v[i][component] * w;
        }

        field_index++;
      }
    }

    if (boundaryflag == BOUNDARY_FIX || boundaryflag == BOUNDARY_BOTH) {

      // Use custom stencil because a bin may overlap with contact point but not atom i
      //   and the line integral needs to add that contribution

      for (auto *wall_fix : wall_fixes) {
        array_atom_fix = wall_fix->array_atom;

        if (array_atom_fix[i][0] < 0.5) continue;
        f_wall[0] = array_atom_fix[i][1];
        f_wall[1] = array_atom_fix[i][2];
        f_wall[2] = array_atom_fix[i][3];
        xcont[0] = array_atom_fix[i][4];
        xcont[1] = array_atom_fix[i][5];
        xcont[2] = array_atom_fix[i][6];

        mc = position_to_bin(xcont);
        if (mc < 0) continue;

        MathExtra::sub3(x[i], xcont, dx_atom_cont);
        r_cont = MathExtra::len3(dx_atom_cont);

        if (chunk_reducedflag) {
          double lamda[3];
          domain->x2lamda(xcont, lamda);
          for (a = 0; a < chunk_ncoord; a++) lamda[cdim[a]] = coord[mc][a];
          domain->lamda2x(lamda, xbinc);
        } else {
          MathExtra::copy3(xcont, xbinc);
          for (a = 0; a < chunk_ncoord; a++) xbinc[cdim[a]] = coord[mc][a];
        }

        pair_stencil_reach = w_cut + r_cont + 0.5 * bin_diagonal;

        stencil_size[0] = stencil_size[1] = stencil_size[2] = 0;
        for (a = 0; a < ncoord; a++)
          stencil_size[a] = static_cast<int>(std::ceil(pair_stencil_reach / bin_width[a]));

        visited_bins.clear();
        for (int dn0 = -stencil_size[0]; dn0 <= stencil_size[0]; dn0++) {
          for (int dn1 = -stencil_size[1]; dn1 <= stencil_size[1]; dn1++) {
            for (int dn2 = -stencil_size[2]; dn2 <= stencil_size[2]; dn2++) {

              MathExtra::copy3(xbinc, xbin);
              if (ncoord >= 1) xbin[chunk_dim[0]] += dn0 * bin_width[0];
              if (ncoord >= 2) xbin[chunk_dim[1]] += dn1 * bin_width[1];
              if (ncoord >= 3) xbin[chunk_dim[2]] += dn2 * bin_width[2];

              mtmp = position_to_bin(xbin);
              if (mtmp == -1) continue;

              // can loop around depending on value of r_cont, so ensure bins only visited once
              if (visited_bins.find(mtmp) != visited_bins.end()) continue;
              visited_bins.insert(mtmp);

              MathExtra::copy3(x[i], xbin2);
              for (int c = 0; c < chunk_ncoord; ++c) xbin2[cdim[c]] = xbin[cdim[c]];

              MathExtra::sub3(x[i], xbin2, dx_atom_bin);
              rsq_atom_bin = MathExtra::lensq3(dx_atom_bin);
              w = calc_w(sqrt(rsq_atom_bin));

              field_index = 0;
              for (auto &val : values) {
                style = std::get<0>(val);
                vtype = std::get<1>(val);
                component = std::get<2>(val);

                if (vtype == 1) {
                  a = component;
                } else {
                  a = component / 3;
                  b = component % 3;
                }

                MathExtra::zero3(dx_atom_cont_filtered);
                for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                  dx_atom_cont_filtered[cdim[coord_index]] = dx_atom_cont[cdim[coord_index]];

                if ((style == STRESS) || (style == STRESSCON)) {
                  w_int_tmp = calc_w_int(dx_atom_bin, dx_atom_cont_filtered);
                  values_local[mtmp][field_index] -= f_wall[a] * dx_atom_cont[b] * w_int_tmp;
                } else if (style == IFD) {
                  MathExtra::copy3(xcont, xbin2);
                  for (int c = 0; c < chunk_ncoord; c++) xbin2[cdim[c]] = xbin[cdim[c]];
                  MathExtra::sub3(xbin2, xcont, dx_bin_cont);
                  rsq_cont_bin = MathExtra::lensq3(dx_bin_cont);
                  wc = calc_w(sqrt(rsq_cont_bin));
                  values_local[mtmp][field_index] -= f_wall[a] * wc;
                }

                field_index++;
              }
            }
          }
        }
      }
    }
  }

  if (calculate_pair) {
    Pair *pair = force->pair;
    double **cutsq = force->pair->cutsq;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;

    tagint itag, jtag;
    tagint *tag = atom->tag;

    int ii, jj, jnum, *jlist;

    neighbor->build_one(list);

    int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];

      if (!(mask[i] & groupbit)) continue;

      voli = 0.0;
      if (radius_required) {
        voli = MY_PI * radius[i] * radius[i];
        if (dim == 3) voli *= 4.0 * THIRD * radius[i];
      }

      if (boundary_group_flag && (mask[i] & boundary_groupbit))
        iboundary = 1;
      else
        iboundary = 0;

      jlist = firstneigh[i];
      jnum = numneigh[i];
      itag = tag[i];
      itype = type[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        factor_lj = special_lj[sbmask(j)];
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
            if (x[j][2] < x[i][2]) continue;
            if (x[j][2] == x[i][2]) {
              if (x[j][1] < x[i][1]) continue;
              if (x[j][1] == x[i][1] && x[j][0] < x[i][0]) continue;
            }
          }
        }

        volj = 0.0;
        if (radius_required) {
          volj = MY_PI * radius[j] * radius[j];
          if (dim == 3) volj *= 4.0 * THIRD * radius[j];
        }

        if (boundary_group_flag && (mask[j] & boundary_groupbit))
          jboundary = 1;
        else
          jboundary = 0;

        // ensure no boundary-boundary interactions (ideally should be pruned from neighbor list)
        if (iboundary && jboundary) continue;

        MathExtra::sub3(x[i], x[j], dx_pair);
        domain->minimum_image(FLERR, dx_pair[0], dx_pair[1], dx_pair[2]);

        rsq_pair = MathExtra::lensq3(dx_pair);
        if (rsq_pair >= cutsq[itype][type[j]]) continue;

        r_pair = sqrt(rsq_pair);

        pair->single(i, j, itype, type[j], rsq_pair, 1.0, 1.0, f_norm);

        MathExtra::scale3(f_norm, dx_pair, f_pair);
        if (pstyle == GRANULAR) {
          f_pair[0] += force->pair->svector[0];
          f_pair[1] += force->pair->svector[1];
          f_pair[2] += force->pair->svector[2];
        }

        MathExtra::scale3(factor_lj, f_pair, f_pair);
        if (MathExtra::lensq3(f_pair) == 0.0) continue;

        MathExtra::sub3(x[i], dx_pair, xj_near);
        MathExtra::add3(x[i], xj_near, xcont);
        MathExtra::scaleadd3((radius[j] - radius[i]) / r_pair, dx_pair, xcont, xcont);
        MathExtra::scale3(0.5, xcont);

        mc = position_to_bin(xcont);
        if (mc < 0) continue;

        if (chunk_reducedflag) {
          double lamda[3];
          domain->x2lamda(xcont, lamda);
          for (a = 0; a < chunk_ncoord; a++) lamda[cdim[a]] = coord[mc][a];
          domain->lamda2x(lamda, xbin0);
        } else {
          MathExtra::copy3(xcont, xbin0);
          for (a = 0; a < chunk_ncoord; a++) xbin0[cdim[a]] = coord[mc][a];
        }

        // create custom stencil and loop over for this pair style
        //  cannot loop i and j's stencil in case some bin contains midpoint (say) but not i or j
        //  in future, could use a bounding box or save stencils based on discretely binned pair
        //    distances to improve performance

        pair_stencil_reach = w_cut + r_pair + 0.5 * bin_diagonal;

        stencil_size[0] = stencil_size[1] = stencil_size[2] = 0;
        for (a = 0; a < ncoord; a++)
          stencil_size[a] = static_cast<int>(std::ceil(pair_stencil_reach / bin_width[a]));

        visited_bins.clear();
        for (int dn0 = -stencil_size[0]; dn0 <= stencil_size[0]; dn0++) {
          for (int dn1 = -stencil_size[1]; dn1 <= stencil_size[1]; dn1++) {
            for (int dn2 = -stencil_size[2]; dn2 <= stencil_size[2]; dn2++) {

              MathExtra::copy3(xbin0, xbin);
              if (ncoord >= 1) xbin[chunk_dim[0]] += dn0 * bin_width[0];
              if (ncoord >= 2) xbin[chunk_dim[1]] += dn1 * bin_width[1];
              if (ncoord >= 3) xbin[chunk_dim[2]] += dn2 * bin_width[2];

              mtmp = position_to_bin(xbin);
              if (mtmp == -1) continue;

              // can loop around depending on value of r_pair, so ensure bins only visited once
              if (visited_bins.find(mtmp) != visited_bins.end()) continue;
              visited_bins.insert(mtmp);

              MathExtra::copy3(x[i], xbin2);
              for (a = 0; a < chunk_ncoord; a++) xbin2[cdim[a]] = xbin[cdim[a]];
              MathExtra::sub3(x[i], xbin2, dx_atom_bin);
              domain->minimum_image(FLERR, dx_atom_bin[0], dx_atom_bin[1], dx_atom_bin[2]);
              rsq_atom_bin = MathExtra::lensq3(dx_atom_bin);
              w = calc_w(sqrt(rsq_atom_bin));

              // contributions from i

              if (jboundary) {
                MathExtra::copy3(xcont, xbin2);
                for (a = 0; a < chunk_ncoord; a++) xbin2[cdim[a]] = xbin[cdim[a]];
                MathExtra::sub3(xbin2, xcont, dx_bin_cont);

                rsq_cont_bin = MathExtra::lensq3(dx_bin_cont);
                wc = calc_w(sqrt(rsq_cont_bin));

                MathExtra::sub3(x[i], xcont, dx_atom_cont);
                MathExtra::zero3(dx_atom_cont_filtered);
                for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                  dx_atom_cont_filtered[cdim[coord_index]] = dx_atom_cont[cdim[coord_index]];
                w_int_tmp = calc_w_int(dx_atom_bin, dx_atom_cont_filtered);
              } else {
                MathExtra::zero3(dx_pair_filtered);
                for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                  dx_pair_filtered[cdim[coord_index]] = dx_pair[cdim[coord_index]];
                w_int_tmp = calc_w_int(dx_atom_bin, dx_pair_filtered);
              }

              field_index = 0;
              for (auto &val : values) {
                style = std::get<0>(val);
                vtype = std::get<1>(val);
                component = std::get<2>(val);

                if (vtype == 1) {
                  a = component;
                } else {
                  a = component / 3;
                  b = component % 3;
                }

                if ((style == STRESS) || (style == STRESSCON)) {
                  if (jboundary) {
                    values_local[mtmp][field_index] -= f_pair[a] * dx_atom_cont[b] * w_int_tmp;
                  } else if (!iboundary) {
                    values_local[mtmp][field_index] -=
                        0.5 * f_pair[a] * dx_pair[b] * w_int_tmp;    // half from each
                  }
                } else if (style == IFD) {
                  if (jboundary) values_local[mtmp][field_index] -= f_pair[a] * wc;
                } else if (style == FABRIC) {
                  if (!iboundary && !jboundary)
                    values_local[mtmp][field_index] +=
                        0.5 * voli * dx_pair[a] * dx_pair[b] * w_int_tmp / rsq_pair;
                }

                field_index++;
              }

              // contributions from j

              MathExtra::copy3(x[j], xbin2);
              for (a = 0; a < chunk_ncoord; a++) xbin2[cdim[a]] = xbin[cdim[a]];
              MathExtra::sub3(x[j], xbin2, dx_atom_bin);
              domain->minimum_image(FLERR, dx_atom_bin[0], dx_atom_bin[1], dx_atom_bin[2]);
              rsq_atom_bin = MathExtra::lensq3(dx_atom_bin);
              w = calc_w(sqrt(rsq_atom_bin));

              if (iboundary) {
                MathExtra::copy3(xcont, xbin2);
                for (a = 0; a < chunk_ncoord; a++) xbin2[cdim[a]] = xbin[cdim[a]];
                MathExtra::sub3(xbin2, xcont, dx_bin_cont);

                rsq_cont_bin = MathExtra::lensq3(dx_bin_cont);
                wc = calc_w(sqrt(rsq_cont_bin));

                MathExtra::sub3(x[j], xcont, dx_atom_cont);
                MathExtra::zero3(dx_atom_cont_filtered);
                for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                  dx_atom_cont_filtered[cdim[coord_index]] = dx_atom_cont[cdim[coord_index]];
                w_int_tmp = calc_w_int(dx_atom_bin, dx_atom_cont_filtered);
              } else {
                MathExtra::zero3(dx_pair_filtered);
                for (int coord_index = 0; coord_index < chunk_ncoord; coord_index++)
                  dx_pair_filtered[cdim[coord_index]] = dx_pair[cdim[coord_index]];

                MathExtra::negate3(dx_pair_filtered);
                w_int_tmp = calc_w_int(dx_atom_bin, dx_pair_filtered);
              }

              field_index = 0;
              for (auto &val : values) {
                style = std::get<0>(val);
                vtype = std::get<1>(val);
                component = std::get<2>(val);

                if (vtype == 1) {
                  a = component;
                } else {
                  a = component / 3;
                  b = component % 3;
                }

                if ((style == STRESS) || (style == STRESSCON)) {
                  if (iboundary) {
                    values_local[mtmp][field_index] -= (-f_pair[a]) * dx_atom_cont[b] * w_int_tmp;
                  } else if (!jboundary) {
                    values_local[mtmp][field_index] -=
                        0.5 * (-f_pair[a]) * (-dx_pair[b]) * w_int_tmp;    // half from each
                  }
                } else if (style == IFD) {
                  if (iboundary) values_local[mtmp][field_index] -= (-f_pair[a]) * wc;
                } else if (style == FABRIC) {
                  if (!iboundary && !jboundary)
                    values_local[mtmp][field_index] +=
                        0.5 * volj * dx_pair[a] * dx_pair[b] * w_int_tmp / rsq_pair;
                }

                field_index++;
              }
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

      if (boundary_group_flag && (mask[i] & boundary_groupbit)) continue;

      if (mask[i] & groupbit) {
        m = ichunk[i] - 1;
        if (m < 0) continue;

        if (chunk_reducedflag) {
          double lamda[3];
          domain->x2lamda(x[i], lamda);
          for (a = 0; a < chunk_ncoord; a++) lamda[cdim[a]] = coord[m][a];
          domain->lamda2x(lamda, xbin0);
        } else {
          MathExtra::copy3(x[i], xbin0);
          for (a = 0; a < chunk_ncoord; a++) xbin0[cdim[a]] = coord[m][a];
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
            massi = rmass[i];
          else
            massi = mass[itype];

          field_index = 0;
          for (auto &val : values) {
            style = std::get<0>(val);
            vtype = std::get<1>(val);
            component = std::get<2>(val);

            if (vtype == 1) {
              a = component;
            } else {
              a = component / 3;
              b = component % 3;
            }

            if (style == TEMPERATURE) {
              values_local[mtmp][field_index] += 0.5 * massi * MathExtra::lensq3(vtemp) * w;
            } else if ((style == STRESS) || (style == STRESSKE)) {
              values_local[mtmp][field_index] -= massi * vtemp[a] * vtemp[b] * w;
            }

            field_index++;
          }
        }
      }
    }
  }

  MPI_Allreduce(&values_local[0][0], &values_global[0][0], nchunk * nvalues, MPI_DOUBLE, MPI_SUM,
                world);

  // Normalize by any unused dimensions as needed
  //   no_norm irrelevant for derived values defined below (like gradients)
  if (bin_dim != dim) {
    int unused_dim[3] = {1, 1, 1};
    for (a = 0; a < ncoord; a++) unused_dim[cdim[a]] = 0;

    for (a = 0; a < dim; a++)
      if (unused_dim[a])
        for (m = 0; m < nchunk; m++)
          for (int n = 0; n < nvalues; n++)
            if (no_norm.find(n) == no_norm.end()) values_global[m][n] /= domain->prd[a];
  }

  // Calculate trivially derived values, in the order used

  // velocity
  double dtemp;
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = std::get<0>(val);
      vtype = std::get<1>(val);
      component = std::get<2>(val);

      if (vtype == 1) {
        a = component;
      } else {
        a = component / 3;
        b = component % 3;
      }

      if (style == VELOCITY) {
        dtemp = values_global[m][index_density];
        if (dtemp != 0.0)
          values_global[m][field_index] = values_global[m][index_momentum[component]] / dtemp;
      }

      field_index++;
    }
  }

  // gradients

  double width[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < chunk_ncoord; a++) {
    width[a] = chunk_delta[a];
    if (chunk_reducedflag) width[a] *= domain->prd[chunk_dim[a]];
  }

  int shift[3], mp, mm, bc;
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = std::get<0>(val);
      vtype = std::get<1>(val);
      component = std::get<2>(val);

      if (vtype == 1) {
        a = component;
      } else {
        a = component / 3;
        b = component % 3;
      }

      if ((style != MGRAD) && (style != VGRAD)) {
        field_index++;
        continue;
      }

      bc = -1;
      for (int c = 0; c < ncoord; c++)
        if (chunk_dim[c] == b) bc = c;

      if (bc == -1) {
        values_global[m][field_index] = 0.0;
        field_index++;
        continue;
      }

      shift[0] = shift[1] = shift[2] = 0;
      shift[bc] = 1;
      mp = shifted_bin(m, shift);
      shift[bc] = -1;
      mm = shifted_bin(m, shift);

      if ((mp == -1) && (mm == -1)) {
        values_global[m][field_index] = 0.0;
        field_index++;
        continue;
      }

      if (mp == -1) {
        if (style == MGRAD) {
          values_global[m][field_index] =
              (values_global[m][index_momentum[a]] - values_global[mm][index_momentum[a]]) /
              (width[bc]);
        } else if (style == VGRAD) {
          values_global[m][field_index] =
              (values_global[m][index_velocity[a]] - values_global[mm][index_velocity[a]]) /
              (width[bc]);
        }
      } else if (mm == -1) {
        if (style == MGRAD) {
          values_global[m][field_index] =
              (values_global[mp][index_momentum[a]] - values_global[m][index_momentum[a]]) /
              (width[bc]);
        } else if (style == VGRAD) {
          values_global[m][field_index] =
              (values_global[mp][index_velocity[a]] - values_global[m][index_velocity[a]]) /
              (width[bc]);
        }
      } else {
        if (style == MGRAD) {
          values_global[m][field_index] =
              (values_global[mp][index_momentum[a]] - values_global[mm][index_momentum[a]]) /
              (2.0 * width[bc]);
        } else if (style == VGRAD) {
          values_global[m][field_index] =
              (values_global[mp][index_velocity[a]] - values_global[mm][index_velocity[a]]) /
              (2.0 * width[bc]);
        }
      }

      field_index++;
    }
  }

  // strain rate
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = std::get<0>(val);
      vtype = std::get<1>(val);
      component = std::get<2>(val);

      if (vtype == 1) {
        a = component;
      } else {
        a = component / 3;
        b = component % 3;
      }

      if (style == STRAINRATE)
        values_global[m][field_index] =
            0.5 * (values_global[m][index_vgrad[a][b]] + values_global[m][index_vgrad[b][a]]);

      field_index++;
    }
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
  if (std::string(option).back() == '*') {
    std::vector<std::string> suffices = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if ((dim == 2) && ((b == 2) || (a == 2))) continue;
        values.emplace_back(std::make_tuple(variable, 2, a * 3 + b));
        labels.emplace_back(trimmed_option + suffices[a * 3 + b]);
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

    values.emplace_back(std::make_tuple(variable, 2, index));
    labels.emplace_back(option);
    if (variable == VGRAD) {
      int a = index / 3;
      int b = index % 3;
      index_vgrad[a][b] = static_cast<int>(values.size()) - 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeContinuumChunk::add_vector_component(char *option, int variable)
{
  if (std::string(option).back() == '*') {
    std::vector<std::string> suffices = {"x", "y", "z"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < dim; a++) {
      values.emplace_back(std::make_tuple(variable, 1, a));
      labels.emplace_back(trimmed_option + suffices[a]);
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
      if (dim == 2) error->all(FLERR, "Invalid compute continuum/chunk property {} in 2D", option);
      index = 2;
    } else {
      error->all(FLERR, "Invalid compute continuum/chunk property {}", option);
    }

    values.emplace_back(std::make_tuple(variable, 1, index));
    labels.emplace_back(option);
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
  bin_width[0] = bin_width[1] = bin_width[2] = 0.0;

  nlayers = cchunk->get_nlayers();
  delta = cchunk->get_delta();
  ncoord = cchunk->ncoord;
  reducedflag = cchunk->get_reducedflag();

  // The distance here is center-to-center between bins, so add diagonal distance to kernel cutoff
  bin_diagonal = 0.0;
  for (int a = 0; a < ncoord; a++) {
    bin_width[a] = delta[a];
    if (reducedflag) bin_width[a] *= domain->prd[chunk_dim[a]];
    bin_diagonal += bin_width[a] * bin_width[a];
  }

  bin_diagonal = sqrt(bin_diagonal);
  double cut = w_cut + 0.5 * bin_diagonal;
  double cut_sq = cut * cut;
  for (int a = 0; a < ncoord; a++) stencil_size[a] = static_cast<int>(ceil(cut / bin_width[a]));

  for (int dn0 = -stencil_size[0]; dn0 <= stencil_size[0]; dn0++) {
    for (int dn1 = -stencil_size[1]; dn1 <= stencil_size[1]; dn1++) {
      for (int dn2 = -stencil_size[2]; dn2 <= stencil_size[2]; dn2++) {
        StencilOffset offset;
        offset.dx[0] = 0.0;
        offset.dx[1] = 0.0;
        offset.dx[2] = 0.0;
        if (ncoord >= 1) offset.dx[chunk_dim[0]] = dn0 * bin_width[0];
        if (ncoord >= 2) offset.dx[chunk_dim[1]] = dn1 * bin_width[1];
        if (ncoord >= 3) offset.dx[chunk_dim[2]] = dn2 * bin_width[2];

        double r_sq = MathExtra::lensq3(offset.dx);
        if (r_sq <= cut_sq) {
          offset.dn[0] = dn0;
          offset.dn[1] = dn1;
          offset.dn[2] = dn2;
          stencil.push_back(offset);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return zero-based Cartesian-bin index for a position, or -1 if outside
------------------------------------------------------------------------- */

int ComputeContinuumChunk::position_to_bin(double *position)
{
  double xcoord[3];
  if (reducedflag)
    domain->x2lamda(position, xcoord);
  else
    MathExtra::copy3(position, xcoord);

  const double *boxlo = reducedflag ? domain->boxlo_lamda : domain->boxlo;
  const double *boxhi = reducedflag ? domain->boxhi_lamda : domain->boxhi;
  const double *prd = reducedflag ? domain->prd_lamda : domain->prd;

  int layer[3] = {0, 0, 0};
  for (int a = 0; a < ncoord; a++) {
    const int idim = chunk_dim[a];
    double xremap = xcoord[idim];

    if (domain->periodicity[idim]) {
      while (xremap < boxlo[idim]) xremap += prd[idim];
      while (xremap >= boxhi[idim]) xremap -= prd[idim];
    }

    // coord[0][a] is the center of the first bin along this coordinate.
    const double offset = cchunk->coord[0][a] - 0.5 * delta[a];
    layer[a] = static_cast<int>((xremap - offset) / delta[a]);
    if (xremap < offset) layer[a]--;

    if (layer[a] < 0 || layer[a] >= nlayers[a]) return -1;
  }

  if (ncoord == 1)
    return layer[0];
  else if (ncoord == 2)
    return layer[0] * nlayers[1] + layer[1];
  else
    return layer[0] * nlayers[1] * nlayers[2] + layer[1] * nlayers[2] + layer[2];
}
