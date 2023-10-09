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
   Contributing authors: Vladislav Galigerov (HSE),
                         Daniil Pavlov (MIPT)
------------------------------------------------------------------------- */

#include "fix_wall_flow.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "update.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "force.h"
#include "math_const.h"

#include <cstring>
#include <iostream>
#include <algorithm>
#include <functional>
#include <fstream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

static const char cite_fix_wall_flow_c[] =
  "fix wall/flow command: doi:{tba}\n\n"
  "@Article{Pavlov-etal-IJHPCA-2023,\n"
  " author = {Daniil Pavlov and Vladislav Galigerov and Daniil Kolotinskii and Vsevolod Nikolskiy and Vladimir Stegailov},\n"
  " title = {GPU-based Molecular Dynamics of Fluid Flows: Reaching for Turbulence},\n"
  " journal = {International Journal of High Performance Computing Applications},\n"
  " year =    2023,\n"
  " volume =  {tba},\n"
  " pages =   {tba}\n”
  "}\n\n";

FixWallFlow::FixWallFlow(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  flowax(FlowAxis::AX_X),
  flowvel(0.0),
  flowdir(0),
  rndseed(0),
  current_segment(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_wall_flow_c);
  if (narg < 9) utils::missing_cmd_args(FLERR, "fix wall/flow", error);

  dynamic_group_allow = 1;
  bool do_abort = false;

  int iarg = 3;
  // parsing axis
  if (strcmp(arg[iarg], "x") == 0) flowax = FlowAxis::AX_X;
  else if (strcmp(arg[iarg],"y") == 0) flowax = FlowAxis::AX_Y;
  else if (strcmp(arg[iarg],"z") == 0) flowax = FlowAxis::AX_Z;
  else error->all(FLERR,"Illegal fix wall/flow argument: axis must by x or y or z, but {} specified", arg[iarg]);

  ++iarg;
  // parsing velocity
  flowvel = utils::numeric(FLERR,arg[iarg],do_abort,lmp);
  if (flowvel == 0.0) error->all(FLERR,"Illegal fix wall/flow argument: velocity cannot be 0");
  if (flowvel > 0.0) flowdir = 1;
  else flowdir = -1;
  if(flowdir < 0) error->all(FLERR, "Negative direction is not supported yet");

  ++iarg;
  // parsing temperature
  double flowtemp = utils::numeric(FLERR,arg[iarg],do_abort,lmp);
  kT = lmp->force->boltz * flowtemp / force->mvv2e;

  ++iarg;
  // parsing seed
  rndseed = utils::inumeric(FLERR, arg[iarg],do_abort,lmp);
  if(rndseed <= 0) error->all(FLERR, "Random seed must be positive!");

  ++iarg;
  // parsing wall count
  int wallcount = utils::inumeric(FLERR,arg[iarg],do_abort,lmp);
  if(wallcount <= 0) error->all(FLERR,"Illegal fix wall/flow argument: wall count must be positive");
  
  ++iarg;
  // parsing walls
  if(narg - iarg != wallcount) error->all(FLERR, "Wrong fix wall/flow wall count {},"
                                                     " must be {}",
                                                     wallcount, narg - iarg);
  walls.resize(wallcount + 2);
  walls.front() = domain->boxlo[flowax];
  for (size_t w = 1; w <= wallcount; ++w, ++iarg)
  {
    walls[w] = utils::numeric(FLERR,arg[iarg],do_abort,lmp);
  }
  walls.back() = domain->boxhi[flowax];
  if (!std::is_sorted(walls.begin(), walls.end(), std::less_equal<double>()))
  {
    error->all(FLERR, "Wrong fix wall/flow wall ordering or some walls are outside simulation domain");
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
  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0)
    error->warning(FLERR,"FixWallFlow is not compatible with rigid bodies");

  for (int i = 0; i < atom->nlocal; ++i)
  {
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

  for (int i = 0; i < nlocal; ++i)
  {
    if (mask[i] & groupbit)
    {
       double pos = x[i][flowax];
       int prev_segment = current_segment[i];
       current_segment[i] = compute_current_segment(pos);

       if (prev_segment != current_segment[i])
       {
          generate_velocity(i);
       }
    }
  }
}

void FixWallFlow::generate_velocity(int atom_i)
{
  const int newton_iteration_count = 10;
  double *vel = atom->v[atom_i];
  double mass = atom->mass[atom->type[atom_i]];
  const double gamma = 1.0 / std::sqrt(2.0 * kT / mass);
  double delta = gamma * flowvel;

  const double edd = std::exp(-delta*delta) / MathConst::MY_PIS + delta * std::erf(delta);
  const double probability_threshold = 0.5f * (1.f + delta / edd);

  double direction = 1.0;

  if (random->uniform() > probability_threshold)
  {
      delta = -delta;
      direction = -direction;
  }

  const double xi_0 = random->uniform();
  const double F_inf = edd + delta;
  const double xi = xi_0 * F_inf;
  const double x_0 = (std::sqrt(delta*delta + 2) - delta) * 0.5;
  double x = x_0;
  for (int i = 0; i < newton_iteration_count; ++i)
  {
    x -= (std::exp(x*x) * MathConst::MY_PIS * (xi - delta * std::erfc(x)) - 1.0) / (x + delta) * 0.5;
  }

  const double nu = x + delta;
  const double v = nu / gamma;

  vel[flowax] = v * direction;
  vel[(flowax + 1) % 3] = random->gaussian() / (gamma * MathConst::MY_SQRT2);
  vel[(flowax + 2) % 3] = random->gaussian() / (gamma * MathConst::MY_SQRT2);
}

int FixWallFlow::compute_current_segment(double pos) const
{
  int result = 0;
  for (; result < walls.size()-1; ++result)
  {
    if (pos >= walls[result] && pos < walls[result + 1])
    {
      return result;
    }
  }
  return -1; // -1 is "out of box" region
}

void FixWallFlow::grow_arrays(int nmax)
{
  memory->grow(current_segment, nmax, "WallFlow::current_segment");
}

void FixWallFlow::copy_arrays(int i, int j, int)
{
  current_segment[j] = current_segment[i];
}

int FixWallFlow::pack_exchange(int i, double* buf)
{
  buf[0] = static_cast<double>(current_segment[i]);
  return 1;
}

int FixWallFlow::unpack_exchange(int i, double* buf)
{
  current_segment[i] = static_cast<int>(buf[0]);
  return 1;
}
