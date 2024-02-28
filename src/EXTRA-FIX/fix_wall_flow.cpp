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
   Contributing authors: Vladislav Galigerov (HSE),
                         Daniil Pavlov (MIPT)
------------------------------------------------------------------------- */

#include "fix_wall_flow.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

static const char cite_fix_wall_flow_c[] =
    "fix wall/flow command: doi:10.1177/10943420231213013\n\n"
    "@Article{Pavlov-etal-IJHPCA-2024,\n"
    " author = {Daniil Pavlov and Vladislav Galigerov and Daniil Kolotinskii and Vsevolod "
    "Nikolskiy and Vladimir Stegailov},\n"
    " title = {GPU-based molecular dynamics of fluid flows: Reaching for turbulence},\n"
    " journal = {The International Journal of High Performance Computing Applications},\n"
    " year =    2024,\n"
    " volume =  38,\n"
    " number =  1,\n"
    " pages =   34-49\n"
    "}\n\n";

FixWallFlow::FixWallFlow(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), flowax(FlowAxis::AX_X), flowvel(0.0), flowdir(0), rndseed(0),
    current_segment(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_wall_flow_c);
  if (narg < 9) utils::missing_cmd_args(FLERR, "fix wall/flow", error);

  if (domain->triclinic != 0)
    error->all(FLERR, "Fix wall/flow cannot be used with triclinic simulation box");

  dynamic_group_allow = 1;
  bool do_abort = false;

  int iarg = 3;
  // parsing axis
  if (strcmp(arg[iarg], "x") == 0)
    flowax = FlowAxis::AX_X;
  else if (strcmp(arg[iarg], "y") == 0)
    flowax = FlowAxis::AX_Y;
  else if (strcmp(arg[iarg], "z") == 0)
    flowax = FlowAxis::AX_Z;
  else
    error->all(FLERR, "Illegal fix wall/flow argument: axis must by x or y or z, but {} specified",
               arg[iarg]);

  if (domain->periodicity[flowax] != 1)
    error->all(FLERR,
               "Fix wall/flow cannot be used with a non-periodic boundary along the flow axis");

  ++iarg;
  // parsing velocity
  flowvel = utils::numeric(FLERR, arg[iarg], do_abort, lmp);
  if (flowvel == 0.0) error->all(FLERR, "Illegal fix wall/flow argument: velocity cannot be 0");
  if (flowvel > 0.0)
    flowdir = 1;
  else
    flowdir = -1;
  if (flowdir < 0)
    error->all(FLERR, "Illegal fix wall/flow argument: negative direction is not supported yet");

  ++iarg;
  // parsing temperature
  double flowtemp = utils::numeric(FLERR, arg[iarg], do_abort, lmp);
  kT = lmp->force->boltz * flowtemp / force->mvv2e;

  ++iarg;
  // parsing seed
  rndseed = utils::inumeric(FLERR, arg[iarg], do_abort, lmp);
  if (rndseed <= 0)
    error->all(FLERR, "Illegal fix wall/flow argument: random seed must be positive integer");

  ++iarg;
  // parsing wall count
  int wallcount = utils::inumeric(FLERR, arg[iarg], do_abort, lmp);
  if (wallcount <= 0)
    error->all(FLERR, "Illegal fix wall/flow argument: wall count must be positive integer");

  ++iarg;
  // parsing walls
  if (narg - iarg != wallcount && narg - iarg != wallcount + 2)
    error->all(FLERR, "Wrong fix wall/flow wall count");

  double scale = 0.0;
  if (flowax == FlowAxis::AX_X)
    scale = domain->lattice->xlattice;
  else if (flowax == FlowAxis::AX_Y)
    scale = domain->lattice->ylattice;
  else if (flowax == FlowAxis::AX_Z)
    scale = domain->lattice->zlattice;

  if (narg - iarg == wallcount + 2) {
    if (strcmp(arg[narg - 2], "units") != 0) error->all(FLERR, "Wrong fix wall/flow units command");
    if (strcmp(arg[narg - 1], "box") == 0)
      scale = 1.0;
    else if (strcmp(arg[narg - 1], "lattice") != 0)
      error->all(FLERR, "Wrong fix wall/flow units command");
  }

  walls.resize(wallcount + 2);
  walls.front() = domain->boxlo[flowax];
  for (int w = 1; w <= wallcount; ++w, ++iarg) {
    walls[w] = utils::numeric(FLERR, arg[iarg], do_abort, lmp) * scale;
  }
  walls.back() = domain->boxhi[flowax];
  if (!std::is_sorted(walls.begin(), walls.end(), std::less_equal<double>())) {
    error->all(FLERR,
               "Wrong fix wall/flow wall ordering or some walls are outside simulation domain");
  }

  if (std::adjacent_find(walls.begin(), walls.end()) != walls.end()) {
    error->all(FLERR,
               "Wrong fix wall/flow wall coordinates: some walls have the same coordinates or lie "
               "on the boundary");
  }

  memory->grow(current_segment, atom->nmax, "WallFlow::current_segment");
  atom->add_callback(Atom::GROW);
  if (restart_peratom) atom->add_callback(Atom::RESTART);

  maxexchange = 1;

  random = new RanMars(lmp, rndseed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixWallFlow::~FixWallFlow()
{
  if (copymode) return;
  atom->delete_callback(id, Atom::GROW);
  if (restart_peratom) atom->delete_callback(id, Atom::RESTART);
  memory->destroy(current_segment);

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixWallFlow::setmask()
{
  int mask = 0;

  mask |= END_OF_STEP;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallFlow::init()
{
  if (domain->triclinic != 0)
    error->all(FLERR, "Fix wall/flow cannot be used with triclinic simulation box");

  int nrigid = 0;
  int box_change_flowax = 0;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) nrigid++;
    switch (flowax) {
      case FlowAxis::AX_X:
        if (ifix->box_change & Fix::BOX_CHANGE_X) box_change_flowax++;
        break;
      case FlowAxis::AX_Y:
        if (ifix->box_change & Fix::BOX_CHANGE_Y) box_change_flowax++;
        break;
      case FlowAxis::AX_Z:
        if (ifix->box_change & Fix::BOX_CHANGE_Z) box_change_flowax++;
        break;
    }
  }

  if (nrigid) error->all(FLERR, "Fix wall/flow is not compatible with rigid bodies");
  if (box_change_flowax)
    error->all(
        FLERR,
        "Fix wall/flow is not compatible with simulation box size changing along flow direction");

  for (int i = 0; i < atom->nlocal; ++i) {
    double pos = atom->x[i][flowax];
    current_segment[i] = compute_current_segment(pos);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallFlow::end_of_step()
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      double pos = x[i][flowax];
      int prev_segment = current_segment[i];
      current_segment[i] = compute_current_segment(pos);

      if (prev_segment != current_segment[i]) generate_velocity(i);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallFlow::generate_velocity(int atom_i)
{
  const int newton_iteration_count = 10;
  double *vel = atom->v[atom_i];

  double *prmass = atom->rmass;
  double *pmass = atom->mass;
  double mass = 0.0;
  if (prmass)
    mass = prmass[atom_i];
  else
    mass = pmass[atom->type[atom_i]];

  const double gamma = 1.0 / std::sqrt(2.0 * kT / mass);
  double delta = gamma * flowvel;

  const double edd = std::exp(-delta * delta) / MathConst::MY_PIS + delta * std::erf(delta);
  const double probability_threshold = 0.5f * (1.f + delta / edd);

  double direction = 1.0;

  if (random->uniform() > probability_threshold) {
    delta = -delta;
    direction = -direction;
  }

  const double xi_0 = random->uniform();
  const double F_inf = edd + delta;
  const double xi = xi_0 * F_inf;
  const double x_0 = (std::sqrt(delta * delta + 2) - delta) * 0.5;
  double x = x_0;
  for (int i = 0; i < newton_iteration_count; ++i) {
    x -= (std::exp(x * x) * MathConst::MY_PIS * (xi - delta * std::erfc(x)) - 1.0) / (x + delta) *
        0.5;
  }

  const double nu = x + delta;
  const double v = nu / gamma;

  vel[flowax] = v * direction;
  vel[(flowax + 1) % 3] = random->gaussian() / (gamma * MathConst::MY_SQRT2);
  vel[(flowax + 2) % 3] = random->gaussian() / (gamma * MathConst::MY_SQRT2);
}

/* ---------------------------------------------------------------------- */

int FixWallFlow::compute_current_segment(double pos) const
{
  int result = 0;
  for (; result < (int)walls.size() - 1; ++result) {
    if (pos >= walls[result] && pos < walls[result + 1]) { return result; }
  }
  return -1;    // -1 is "out of box" region
}

/* ---------------------------------------------------------------------- */

void FixWallFlow::grow_arrays(int nmax)
{
  memory->grow(current_segment, nmax, "WallFlow::current_segment");
}

/* ---------------------------------------------------------------------- */

void FixWallFlow::copy_arrays(int i, int j, int)
{
  current_segment[j] = current_segment[i];
}

/* ---------------------------------------------------------------------- */

int FixWallFlow::pack_exchange(int i, double *buf)
{
  buf[0] = static_cast<double>(current_segment[i]);
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixWallFlow::unpack_exchange(int i, double *buf)
{
  current_segment[i] = static_cast<int>(buf[0]);
  return 1;
}
