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
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)

   Contact models of the GRANULAR package for the contacts between rounded
   polyhedra and a flat wall: each rounded vertex is a sphere on the wall,
   with the normal force and the normal damping of the normal and damping
   sub-models, and the friction force of the wall on a body, from all its
   vertices in contact, is that of the tangential sub-model, as for pair
   style body/rounded/polyhedron/granular.  The wall has the same material
   properties as the body, as for fix wall/gran.
------------------------------------------------------------------------- */

#include "fix_wall_body_polyhedron_granular.h"

#include "atom.h"
#include "error.h"
#include "gran_sub_mod.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "math_extra.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace Granular_NS;

/* ----------------------------------------------------------------------
   fix ID group wall/body/polyhedron/granular normal_model coeffs keyword
     values ... wallstyle lo hi keyword values ...
------------------------------------------------------------------------- */

FixWallBodyPolyhedronGranular::FixWallBodyPolyhedronGranular(LAMMPS *lmp, int narg, char **arg) :
    FixWallBodyPolyhedron(lmp, narg, arg, 1), model(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix wall/body/polyhedron/granular", error);

  // contact model of a vertex with the wall, as for fix wall/gran granular:
  // the sub-models up to the wall style

  model = new GranularModel(lmp);
  model->contact_type = WALL;

  int iarg = model->add_sub_model(arg, 3, narg, NORMAL);
  while (iarg < narg) {
    if (strcmp(arg[iarg], "damping") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
    } else if (strcmp(arg[iarg], "tangential") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
    } else if ((strcmp(arg[iarg], "xplane") == 0) || (strcmp(arg[iarg], "yplane") == 0) ||
               (strcmp(arg[iarg], "zplane") == 0)) {
      break;
    } else {
      error->all(FLERR, iarg, "Unknown fix wall/body/polyhedron/granular keyword {}", arg[iarg]);
    }
  }

  // the damping of Wang et al. is proportional to the normal relative velocity

  if (!model->damping_model) model->construct_sub_model("velocity", DAMPING);
  model->init();

  const std::string &normal = model->normal_model->name;
  const std::string &damping = model->damping_model->name;
  const std::string &tangential = model->tangential_model->name;
  if ((normal != "hooke") && (normal != "hertz") && (normal != "hertz/material"))
    error->all(FLERR, "Fix wall/body/polyhedron/granular does not support the normal model {}",
               normal);
  if ((damping != "velocity") && (damping != "mass_velocity"))
    error->all(FLERR, "Fix wall/body/polyhedron/granular does not support the damping model {}",
               damping);
  if ((tangential != "linear_nohistory") && (tangential != "linear_history"))
    error->all(FLERR,
               "Fix wall/body/polyhedron/granular does not support the tangential model {}",
               tangential);
  if ((model->sub_models[ROLLING]->name != "none") ||
      (model->sub_models[TWISTING]->name != "none") ||
      (model->sub_models[HEAT]->name != "none") || model->limit_damping ||
      model->synchronized_verlet)
    error->all(FLERR, "Fix wall/body/polyhedron/granular supports only the normal, damping, "
               "and tangential sub-models");

  // the history of the sub-models of a body at a wall is stored in the
  // 3 values per wall of the tangential deformation, see FixWallBodyPolyhedron

  int next_index = 0;
  for (int k = 0; k < NSUBMODELS; k++) {
    model->sub_models[k]->history_index = next_index;
    next_index += model->sub_models[k]->size_history;
  }
  if (next_index > 3)
    error->all(FLERR, "Fix wall/body/polyhedron/granular supports at most 3 history values");

  // wall style and keywords, the friction force of a body at the wall always
  // comes from the tangential sub-model, with its history if it has one

  if (iarg >= narg) utils::missing_cmd_args(FLERR, "fix wall/body/polyhedron/granular", error);
  iarg = parse_wall(narg, arg, iarg);
  while (iarg < narg) {
    if (strcmp(arg[iarg], "wiggle") == 0) {
      iarg = parse_wiggle(narg, arg, iarg);
    } else {
      error->all(FLERR, iarg, "Unknown fix wall/body/polyhedron/granular keyword {}", arg[iarg]);
    }
  }

  history = 1;
  setup_wall();
}

/* ---------------------------------------------------------------------- */

FixWallBodyPolyhedronGranular::~FixWallBodyPolyhedronGranular()
{
  delete model;
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolyhedronGranular::init()
{
  FixWallBodyPolyhedron::init();
  if (!atom->rmass_flag)
    error->all(FLERR, Error::NOLASTLINE, "Fix wall/body/polyhedron/granular requires per-atom mass");
}

/* ----------------------------------------------------------------------
   geometry of the contact of a rounded vertex of body i with the wall at
   the surface separation R: a sphere with the rounded radius of the body
   on a flat wall, with the overlap delta = -R
------------------------------------------------------------------------- */

void FixWallBodyPolyhedronGranular::contact_geometry(int i, double R)
{
  model->i = i;
  model->itype = atom->type[i];
  model->radi = rounded_radius[i];
  model->radj = 0.0;
  model->radsum = model->radi;
  model->Reff = model->radi;
  model->delta = -R;
  model->dR = MAX(model->delta, 0.0) * model->Reff;
  model->meff = atom->rmass[i];
  model->dt = dt;
  model->history_update = update->setupflag ? 0 : 1;
  if (model->normal_model->contact_radius_flag)
    model->contact_radius = model->normal_model->calculate_contact_radius();
}

/* ----------------------------------------------------------------------
   elastic normal force and normal damping of a vertex that overlaps the
   wall by -R > 0, along the inward normal n of the wall, the tangential
   damping is part of the tangential sub-model, see tangential_spring()
   returns the elastic normal force
------------------------------------------------------------------------- */

double FixWallBodyPolyhedronGranular::vertex_force(int i, double R, double vnnr, const double *n,
                                                   const double * /*vn*/, const double * /*vt*/,
                                                   double *fw)
{
  contact_geometry(i, R);
  model->vnnr = vnnr;
  double fne = model->normal_model->calculate_forces();
  double fdamp = model->damping_model->calculate_forces();
  for (int k = 0; k < 3; k++) fw[k] = (fne + fdamp) * n[k];
  return fne;
}

/* ----------------------------------------------------------------------
   friction force ft of the wall on body i from the tangential sub-model,
   with the average tangential velocity vt of its vertices in contact,
   limited by mu times the sum of their elastic normal forces fne, with the
   history xi of the body at the wall, as in FixWallBodyPolyhedron::post_force()
------------------------------------------------------------------------- */

void FixWallBodyPolyhedronGranular::tangential_spring(int i, const double *n, const double *vt,
                                                      double fne, double *xi, double *ft)
{
  contact_geometry(i, 0.0);
  for (int k = 0; k < 3; k++) {
    model->nx[k] = model->nx_unrotated[k] = n[k];
    model->vtr[k] = vt[k];
  }
  model->vrel = MathExtra::len3(vt);
  model->vnnr = 0.0;
  model->history = xi;

  // the same order as in GranularModel::calculate_forces(): the damping
  // sub-model sets the prefactor of the tangential damping, and the normal
  // sub-model the critical force of the tangential sub-model

  model->Fntot = fne;
  model->damping_model->calculate_forces();
  model->normal_model->set_fncrit();
  model->tangential_model->calculate_forces();

  for (int k = 0; k < 3; k++) ft[k] = model->fs[k];
}
