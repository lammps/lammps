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
   polyhedra: the normal force and the normal damping of each contact, and
   the tangential force of the single friction force of a pair of bodies,
   see PairBodyRoundedPolyhedron::history_friction(), are those of the
   sub-models of pair style granular.  The detection of the contacts, their
   weights, the cohesion of Wang et al., and the work done by the forces that
   do not derive from the energy are those of body/rounded/polyhedron.

   PROTOTYPE: the sub-models are called one by one, in the same order as in
   GranularModel::calculate_forces(), and the layout of the history is the
   same as in PairGranular::init_style().  Keep both in sync.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polyhedron_granular.h"

#include "atom.h"
#include "error.h"
#include "gran_sub_mod.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "math_extra.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace Granular_NS;

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronGranular::PairBodyRoundedPolyhedronGranular(LAMMPS *lmp) :
    PairBodyRoundedPolyhedron(lmp), model_index(nullptr)
{
  // the friction force of a pair of bodies is always computed from all its
  // contacts, with the history of the tangential sub-model, see settings()

  c_n = c_t = mu = 0.0;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronGranular::~PairBodyRoundedPolyhedronGranular()
{
  for (auto *model : models) delete model;
  memory->destroy(model_index);
}

/* ----------------------------------------------------------------------
   global settings: A_ua cut_inner
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::settings(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR, "Pair style body/rounded/polyhedron/granular requires A_ua and cut_inner");

  A_ua = utils::numeric(FLERR, arg[0], false, lmp);
  cut_inner = utils::numeric(FLERR, arg[1], false, lmp);
  if (A_ua < 0) A_ua = 1;

  // the history is always stored, its size is set in init_style()

  if (!history) history_dummy_fix(1);
  history = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs:
     normal sub-model with its coefficients, followed by the keywords
     damping, tangential (as in pair style granular), and cohesion k_na
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ntypes = atom->ntypes;
  if (!model_index) {
    memory->create(model_index, ntypes + 1, ntypes + 1, "pair:model_index");
    for (int i = 0; i <= ntypes; i++)
      for (int j = 0; j <= ntypes; j++) model_index[i][j] = -1;
  }

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, ntypes, jlo, jhi, error);

  GranularModel *model = add_model();
  int m = models.size() - 1;
  double k_na_one = parse_model(model, narg, arg);
  model_args[m].assign(arg, arg + narg);
  check_model(m);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      model_index[i][j] = m;
      k_n[i][j] = k_t[i][j] = 0.0;
      k_na[i][j] = k_na_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ---------------------------------------------------------------------- */

GranularModel *PairBodyRoundedPolyhedronGranular::add_model()
{
  auto *model = new GranularModel(lmp);
  model->contact_type = PAIR;
  models.push_back(model);
  energy_factor.push_back(0.0);
  model_args.emplace_back();
  mix_itype.push_back(0);
  mix_jtype.push_back(0);
  return model;
}

/* ----------------------------------------------------------------------
   sub-models of a model from the pair_coeff arguments arg[2] to arg[narg-1],
   return the strength of the cohesion
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronGranular::parse_model(GranularModel *model, int narg, char **arg)
{
  double k_na_one = 0.0;

  int iarg = model->add_sub_model(arg, 2, narg, NORMAL);
  while (iarg < narg) {
    if (strcmp(arg[iarg], "damping") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
    } else if (strcmp(arg[iarg], "tangential") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
    } else if (strcmp(arg[iarg], "cohesion") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_coeff cohesion", error);
      k_na_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (k_na_one < 0.0) error->all(FLERR, iarg + 1, "Cohesion strength must not be negative");
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown pair_coeff keyword {} for pair style "
                 "body/rounded/polyhedron/granular", arg[iarg]);
    }
  }

  // the damping of Wang et al. is proportional to the normal relative velocity

  if (!model->damping_model) model->construct_sub_model("velocity", DAMPING);
  model->init();
  return k_na_one;
}

/* ----------------------------------------------------------------------
   a new copy of model m, made the same way as model m, from its pair_coeff
   arguments or by mixing the models of two types, with the same layout of
   the history, see init_style()
   the coefficients are not copied, since some sub-models also keep whether
   their coefficients were mixed
------------------------------------------------------------------------- */

GranularModel *PairBodyRoundedPolyhedronGranular::clone_model(int m)
{
  auto *model = new GranularModel(lmp);
  model->contact_type = PAIR;

  if (!model_args[m].empty()) {
    std::vector<char *> args;
    for (auto &a : model_args[m]) args.push_back(const_cast<char *>(a.c_str()));
    parse_model(model, (int) args.size(), args.data());
  } else {
    model->mix_coeffs(type_model(mix_itype[m], mix_itype[m]),
                      type_model(mix_jtype[m], mix_jtype[m]));
    model->init();
  }

  for (int k = 0; k < NSUBMODELS; k++)
    model->sub_models[k]->history_index = models[m]->sub_models[k]->history_index;
  return model;
}

/* ----------------------------------------------------------------------
   accept only the sub-models that are well defined for the contacts of
   rounded polyhedra and for the single friction force of a pair of bodies:
   the normal force depends on the overlap of each contact only, and the
   prefactors of the damping and of the tangential force do not depend on
   the geometry of the contact, since the friction force of a pair of bodies
   combines all its contacts
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::check_model(int m)
{
  GranularModel *model = models[m];
  const std::string &normal = model->normal_model->name;
  const std::string &damping = model->damping_model->name;
  const std::string &tangential = model->tangential_model->name;

  // elastic energy: Fne = k delta for hooke, Fne ~ delta^3/2 for hertz

  if (normal == "hooke") energy_factor[m] = 0.5;
  else if ((normal == "hertz") || (normal == "hertz/material")) energy_factor[m] = 0.4;
  else
    error->all(FLERR, "Pair style body/rounded/polyhedron/granular does not support "
               "the normal model {}", normal);

  if ((damping != "velocity") && (damping != "mass_velocity"))
    error->all(FLERR, "Pair style body/rounded/polyhedron/granular does not support "
               "the damping model {}", damping);

  if ((tangential != "linear_nohistory") && (tangential != "linear_history"))
    error->all(FLERR, "Pair style body/rounded/polyhedron/granular does not support "
               "the tangential model {}", tangential);

  if ((model->sub_models[ROLLING]->name != "none") ||
      (model->sub_models[TWISTING]->name != "none") ||
      (model->sub_models[HEAT]->name != "none") || model->limit_damping ||
      model->synchronized_verlet)
    error->all(FLERR, "Pair style body/rounded/polyhedron/granular supports only the "
               "normal, damping, and tangential sub-models");
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::init_style()
{
  if (!atom->rmass_flag)
    error->all(FLERR, "Pair body/rounded/polyhedron/granular requires per-atom mass");

  // layout of the history of a pair of bodies, as in PairGranular::init_style():
  // the sub-models of each type follow each other with the largest size among
  // the models, at least one value, which is required by fix NEIGH_HISTORY

  int size_max[NSUBMODELS] = {0};
  nondefault_history_transfer = 0;
  for (auto *model : models) {
    for (int k = 0; k < NSUBMODELS; k++)
      size_max[k] = MAX(size_max[k], model->sub_models[k]->size_history);
    if (model->nondefault_history_transfer) nondefault_history_transfer = 1;
  }

  size_history = 0;
  for (int k = 0; k < NSUBMODELS; k++) size_history += size_max[k];
  size_history = MAX(size_history, 1);

  for (auto *model : models) {
    int next_index = 0;
    for (int k = 0; k < NSUBMODELS; k++) {
      model->sub_models[k]->history_index = next_index;
      next_index += size_max[k];
    }
  }

  PairBodyRoundedPolyhedron::init_style();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   the models of unset pairs of types are mixed, as in pair style granular
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronGranular::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    if ((model_index[i][i] < 0) || (model_index[j][j] < 0))
      error->all(FLERR, "All pair coeffs are not set");

    GranularModel *model1 = type_model(i, i);
    GranularModel *model2 = type_model(j, j);

    // replace the model mixed at a previous init, since the models of
    // types i and j may have changed since then

    GranularModel *model;
    int m = model_index[i][j];
    if (m >= 0) {
      delete models[m];
      model = models[m] = new GranularModel(lmp);
      model->contact_type = PAIR;
    } else {
      model = add_model();
      m = models.size() - 1;
    }

    int error_code = model->mix_coeffs(model1, model2);
    if (error_code != -1)
      error->all(FLERR, "Granular sub-models {} and {} of types {} and {} cannot be mixed, "
                 "set the coefficients of this pair of types explicitly",
                 model1->sub_models[error_code]->name, model2->sub_models[error_code]->name,
                 i, j);
    model->init();
    check_model(m);
    for (int k = 0; k < NSUBMODELS; k++)
      model->sub_models[k]->history_index = model1->sub_models[k]->history_index;
    model_args[m].clear();
    mix_itype[m] = i;
    mix_jtype[m] = j;

    model_index[i][j] = m;
    k_n[i][j] = k_t[i][j] = 0.0;
    k_na[i][j] = sqrt(k_na[i][i] * k_na[j][j]);
  }

  model_index[j][i] = model_index[i][j];
  return PairBodyRoundedPolyhedron::init_one(i, j);
}

/* ----------------------------------------------------------------------
   history of a pair of bodies seen from the other body, see
   PairGranular::transfer_history(), only called if a sub-model has a
   history that is not just the negative of that of the other body
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::transfer_history(double *source, double *target,
                                                         int itype, int jtype)
{
  GranularModel *model = type_model(itype, jtype);
  for (int n = 0; n < size_history; n++) target[n] = -source[n];
  for (int k = 0; k < NSUBMODELS; k++) {
    GranSubMod *sub = model->sub_models[k];
    if (!sub->nondefault_history_transfer) continue;
    for (int n = 0; n < sub->size_history; n++)
      target[sub->history_index + n] =
          sub->transfer_history_factor[n] * source[sub->history_index + n];
  }
}

/* ----------------------------------------------------------------------
   effective radius of curvature of a contact of the given kind, see
   ContactKind, between the rounded surfaces of ibody with the rounded
   radius radi and jbody with radj, for the Hertzian normal models:
   1/Reff is the geometric mean of the principal relative curvatures, so
   that the area of the elliptic contact is that of the circular one of
   radius sqrt(Reff delta), see Johnson, Contact Mechanics (1985), Ch. 4
     a rounded vertex is a sphere, so that 1/Reff = 1/radi + 1/radj
     a rounded face is flat, so that Reff is the radius of the vertex
     a rounded edge is a cylinder: on a vertex, the relative curvatures are
       1/r_v and 1/r_v + 1/r_e
     two edges are taken as crossing at a right angle, with the relative
       curvatures 1/radi and 1/radj, since a contact of nearly parallel edges
       is that of their vertices with the faces instead, see
       PairBodyRoundedPolyhedron::vertex_face_weight()
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronGranular::effective_radius(int kind, double radi, double radj)
{
  if ((radi <= 0.0) || (radj <= 0.0)) return 0.0;

  switch (kind) {
    case VERTEX_FACE:
      return radi;
    case FACE_VERTEX:
      return radj;
    case VERTEX_EDGE:
      return radi * sqrt(radj / (radi + radj));
    case EDGE_VERTEX:
      return radj * sqrt(radi / (radi + radj));
    case EDGE_EDGE:
      return sqrt(radi * radj);
    case VERTEX_VERTEX:
    default:
      return radi * radj / (radi + radj);
  }
}

/* ----------------------------------------------------------------------
   geometry of a contact of the given kind at the surface separation R
   between the rounded surfaces of ibody and jbody for the sub-models of
   the GRANULAR package: the overlap delta = -R, and the effective radius
   of curvature, see effective_radius()
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::contact_geometry(GranularModel *model, int ibody,
                                                         int jbody, int kind, double R)
{
  double *rmass = atom->rmass;

  model->i = ibody;
  model->j = jbody;
  model->itype = atom->type[ibody];
  model->jtype = atom->type[jbody];
  model->radi = rounded_radius[ibody];
  model->radj = rounded_radius[jbody];
  model->radsum = model->radi + model->radj;
  model->Reff = effective_radius(kind, model->radi, model->radj);
  model->delta = -R;
  model->dR = MAX(model->delta, 0.0) * model->Reff;
  model->meff = rmass[ibody] * rmass[jbody] / (rmass[ibody] + rmass[jbody]);
  model->dt = dt;
  model->history_update = update->setupflag ? 0 : 1;
  if (model->normal_model->contact_radius_flag)
    model->contact_radius = model->normal_model->calculate_contact_radius();
}

/* ----------------------------------------------------------------------
   elastic force of the normal sub-model of a contact of the given kind at
   the surface separation R < 0, and the cohesion of Wang et al. with the
   strength k_na, see PairBodyRoundedPolyhedron::normal_force()
   return the energy
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronGranular::normal_force(int ibody, int jbody, int kind, double R,
                                                       int itype, int jtype, double &fe,
                                                       double &fc)
{
  fe = 0.0;
  double energy = cohesive_force(R, itype, jtype, fc);

  if (R < 0.0) {
    GranularModel *model = contact_model(itype, jtype);
    contact_geometry(model, ibody, jbody, kind, R);
    fe = model->normal_model->calculate_forces();
    energy += energy_factor[model_index[itype][jtype]] * fe * model->delta;
  }
  return energy;
}

/* ----------------------------------------------------------------------
   normal damping of the damping sub-model of a contact of the given kind at
   the surface separation R with the weight w, the tangential damping is part
   of the tangential sub-model of the pair of bodies, see tangential_spring()
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::contact_damping(int ibody, int jbody, int kind, double R,
                                                        double w, double vnnr, const double *n,
                                                        const double * /*vn*/,
                                                        const double * /*vt*/, double *fdiss)
{
  GranularModel *model = contact_model(atom->type[ibody], atom->type[jbody]);
  contact_geometry(model, ibody, jbody, kind, R);
  model->vnnr = vnnr;
  double fdamp = w * model->damping_model->calculate_forces();
  for (int k = 0; k < 3; k++) fdiss[k] = fdamp * n[k];
}

/* ----------------------------------------------------------------------
   tangential force of the tangential sub-model for the single friction
   force of a pair of bodies, with the average normal n and tangential
   relative velocity vt of its contacts, limited by mu times the sum of
   the elastic normal forces fne, see
   PairBodyRoundedPolyhedron::history_friction()
   the history refers to body s.shear_i, the force on ibody is added to fs
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranular::tangential_spring(int ibody, int jbody, const double *n,
                                                          const double *vt, double fne,
                                                          Scratch &s, double *fs)
{
  GranularModel *model = contact_model(atom->type[ibody], atom->type[jbody]);
  double sign = (ibody == s.shear_i) ? 1.0 : -1.0;

  // the force is that of all contacts of the pair, so there is no single kind
  // of contact, but the supported tangential models do not depend on it

  contact_geometry(model, ibody, jbody, VERTEX_VERTEX, 0.0);
  for (int k = 0; k < 3; k++) {
    model->nx[k] = model->nx_unrotated[k] = sign * n[k];
    model->vtr[k] = sign * vt[k];
  }
  model->vrel = MathExtra::len3(model->vtr);
  model->vnnr = 0.0;
  model->history = s.shear;

  // the same order as in GranularModel::calculate_forces(): the damping
  // sub-model sets the prefactor of the tangential damping, and the normal
  // sub-model the critical force of the tangential sub-model

  model->Fntot = fne;
  model->damping_model->calculate_forces();
  model->normal_model->set_fncrit();
  model->tangential_model->calculate_forces();

  for (int k = 0; k < 3; k++) fs[k] += sign * model->fs[k];
  s.touched = 1;
}
