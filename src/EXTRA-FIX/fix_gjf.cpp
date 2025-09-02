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
   Contributing authors: Tim Linke & Niels Gronbech-Jensen (UC Davis)
------------------------------------------------------------------------- */

#include "fix_gjf.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NOBIAS, BIAS };
enum { CONSTANT, EQUAL, ATOM };

static const char cite_gjf[] =
    "GJ methods: doi:10.1080/00268976.2019.1662506\n\n"
    "@Article{gronbech-jensen_complete_2020,\n"
    "title = {Complete set of stochastic Verlet-type thermostats for correct Langevin "
    "simulations},\n"
    "volume = {118},\n"
    "number = {8},\n"
    "url = {https://www.tandfonline.com/doi/full/10.1080/00268976.2019.1662506},\n"
    "doi = {10.1080/00268976.2019.1662506},\n"
    "journal = {Molecular Physics},\n"
    "author = {Grønbech-Jensen, Niels},\n"
    "year = {2020}\n"
    "}\n\n";

static const char cite_gjf_7[] = "GJ-VII method: doi:10.1063/5.0066008\n\n"
                                 "@Article{finkelstein_2021,\n"
                                 "title = {Bringing discrete-time Langevin splitting methods into "
                                 "agreement with thermodynamics},\n"
                                 "volume = {155},\n"
                                 "number = {18},\n"
                                 "url = {https://doi.org/10.1063/5.0066008},\n"
                                 "doi = {10.1063/5.0066008},\n"
                                 "journal = {J. Chem. Phys.},\n"
                                 "author = {Finkelstein, Joshua and Cheng, Chungho and Fiorin, "
                                 "Giacomo and Seibold, Benjamin and Grønbech-Jensen, Niels},\n"
                                 "year = {2021},\n"
                                 "pages = {184104}\n"
                                 "}\n\n";

static const char cite_gjf_8[] =
    "GJ-VIII method: doi:10.1007/s10955-024-03345-1\n\n"
    "@Article{gronbech_jensen_2024,\n"
    "title = {On the Definition of Velocity in Discrete-Time, Stochastic Langevin Simulations},\n"
    "volume = {191},\n"
    "number = {10},\n"
    "url = {https://doi.org/10.1007/s10955-024-03345-1},\n"
    "doi = {10.1007/s10955-024-03345-1},\n"
    "journal = {J. Stat. Phys.},\n"
    "author = {Gronbech-Jensen, Niels},\n"
    "year = {2024},\n"
    "pages = {137}\n"
    "}\n\n";

static const char cite_gjf_vhalf[] =
    "GJ-I vhalf method: doi:10.1080/00268976.2019.1570369\n\n"
    "@Article{jensen_accurate_2019,\n"
    "title = {Accurate configurational and kinetic statistics in discrete-time Langevin systems},\n"
    "volume = {117},\n"
    "url = {https://www.tandfonline.com/doi/full/10.1080/00268976.2019.1570369},\n"
    "doi = {10.1080/00268976.2019.1570369},\n"
    "number = {18},\n"
    "journal = {Molecular Physics},\n"
    "author = {Jensen, Lucas Frese Grønbech and Grønbech-Jensen, Niels},\n"
    "year = {2019}\n"
    "}\n\n";

static const char cite_gjf_vfull[] =
    "GJ-I vfull method: doi:10.1080/00268976.2012.760055\n\n"
    "@Article{gronbech-jensen_simple_2013,\n"
    "title = {A simple and effective Verlet-type algorithm for simulating Langevin dynamics},\n"
    "volume = {111},\n"
    "url = {http://www.tandfonline.com/doi/abs/10.1080/00268976.2012.760055},\n"
    "doi = {10.1080/00268976.2012.760055},\n"
    "pages = {983-991},\n"
    "number = {8},\n"
    "journal = {Molecular Physics},\n"
    "author = {Grønbech-Jensen, Niels and Farago, Oded},\n"
    "year = {2013}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixGJF::FixGJF(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), tstr(nullptr), tforce(nullptr), lv(nullptr), id_temp(nullptr),
    random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_gjf);
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix gjf", error);

  time_integrate = 1;
  global_freq = 1;
  nevery = 1;

  if (utils::strmatch(arg[3], "^v_")) {
    tstr = utils::strdup(arg[3] + 2);
  } else {
    t_start = utils::numeric(FLERR, arg[3], false, lmp);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = utils::numeric(FLERR, arg[4], false, lmp);
  t_period = utils::numeric(FLERR, arg[5], false, lmp);
  seed = utils::inumeric(FLERR, arg[6], false, lmp);

  if (t_period <= 0.0) error->all(FLERR, 5, "Fix gjf period must be > 0.0");
  if (seed <= 0) error->all(FLERR, 6, "Illegal fix gjf command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp, seed + comm->me);

  int GJmethods = 8;    // number of currently implemented GJ methods
  maxatom = 0;

  // optional args
  // per default, use half step and GJ-I

  osflag = 0;
  GJmethod = 1;
  lv_allocated = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "vel") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix gjf vel", error);
      if (strcmp(arg[iarg + 1], "vfull") == 0) {
        osflag = 1;
      } else if (strcmp(arg[iarg + 1], "vhalf") == 0) {
        osflag = 0;
      } else
        error->all(FLERR, iarg + 1, "Unknown fix gjf vel keyword {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "method") == 0) {
      GJmethod = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (GJmethod == 7) {
        if (iarg + 3 > narg) error->all(FLERR, "Illegal fix gjf command for GJ-VII");
        gjfc2 = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        if (gjfc2 < 0 || gjfc2 > 1) error->all(FLERR, "Choice of c2 in GJ-VII must be 0≤c2≤1");
        iarg += 3;
        if (lmp->citeme) lmp->citeme->add(cite_gjf_7);
      } else {
        if (iarg + 2 > narg) error->all(FLERR, "Illegal fix gjf command");
        if (GJmethod < 0 || GJmethod > GJmethods)
          error->all(FLERR, "Invalid GJ method choice in gjf command");
        if (GJmethod == 8)
          if (lmp->citeme) lmp->citeme->add(cite_gjf_8);
        iarg += 2;
      }
    } else
      error->all(FLERR, "Illegal fix gjf command");
  }
  if (GJmethod == 1 && osflag == 0)
    if (lmp->citeme) lmp->citeme->add(cite_gjf_vhalf);
  if (GJmethod == 1 && osflag == 1)
    if (lmp->citeme) lmp->citeme->add(cite_gjf_vfull);

  // set temperature = nullptr, user can override via fix_modify if wants bias
  id_temp = nullptr;
  temperature = nullptr;

  lv = nullptr;
  tforce = nullptr;

  // setup atom-based array for lv
  // register with Atom class
  // no need to set peratom_flag, b/c data is for internal use only

  FixGJF::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // initialize lv to onsite velocity
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    lv[i][0] = 0.0;
    lv[i][1] = 0.0;
    lv[i][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixGJF::~FixGJF()
{
  if (copymode) return;

  delete random;
  delete[] tstr;
  delete[] id_temp;
  memory->destroy(tforce);

  memory->destroy(lv);
  if (modify->get_fix_by_id(id)) atom->delete_callback(id, Atom::GROW);
}

/* ---------------------------------------------------------------------- */

int FixGJF::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  if (!osflag) mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGJF::init()
{
  if (id_temp) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature) {
      error->all(FLERR, Error::NOLASTLINE, "Temperature compute ID {} for fix {} does not exist",
                 id_temp, style);
    } else {
      if (temperature->tempflag == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Compute ID {} for fix {} does not compute temperature", id_temp, style);
    }
  }
  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix gjf does not exist", tstr);
    if (input->variable->equalstyle(tvar))
      tstyle = EQUAL;
    else if (input->variable->atomstyle(tvar))
      tstyle = ATOM;
    else
      error->all(FLERR, Error::NOLASTLINE, "Variable {} for fix gjf is invalid style", tstr);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    error->all(FLERR, Error::NOLASTLINE, "Fix gjf and run style respa are not compatible");
  }

  if (temperature && temperature->tempbias)
    tbiasflag = BIAS;
  else
    tbiasflag = NOBIAS;

  // Complete set of thermostats is given in Gronbech-Jensen, Molecular Physics, 118 (2020)
  switch (GJmethod) {
    case 1:
      gjfc2 = (1.0 - update->dt / 2.0 / t_period) / (1.0 + update->dt / 2.0 / t_period);
      break;
    case 2:
      gjfc2 = exp(-update->dt / t_period);
      break;
    case 3:
      gjfc2 = 1.0 - update->dt / t_period;
      break;
    case 4:
      gjfc2 = (sqrt(1.0 + 4.0 * (update->dt / t_period)) - 1.0) / (2.0 * update->dt / t_period);
      break;
    case 5:
      gjfc2 = 1.0 / (1.0 + update->dt / t_period);
      break;
    case 6:
      gjfc2 =
          (1.0 / (1.0 + update->dt / 2.0 / t_period)) * (1.0 / (1.0 + update->dt / 2.0 / t_period));
      break;
    case 7:    // provided in Finkelstein (2021)
      update->dt = (1.0 + gjfc2) / (1.0 - gjfc2) * log(gjfc2) * log(gjfc2) * 0.5 * t_period;
      break;
    case 8:    // provided in Gronbech-Jensen (2024)
      gjfc2 = sqrt((update->dt / t_period) * (update->dt / t_period) + 1.0) - update->dt / t_period;
      break;
    case 0:
      gjfc2 = 0.0;
      break;
    default:
      error->all(FLERR, "Fix gjf method not found");
      break;
  }
  gjfc1 = (1.0 + gjfc2) / 2.0;
  gjfc3 = (1.0 - gjfc2) * t_period / update->dt;
}

/* ----------------------------------------------------------------------
  integrate position and velocity according to the GJ methods
  in Grønbech-Jensen, J Stat Phys 191, 137 (2024). The general workflow is
    1. GJ Initial Integration
    2. Force Update
    3. GJ Final Integration
    4. Velocity Choice in end_of_step()
------------------------------------------------------------------------- */

void FixGJF::initial_integrate(int /* vflag */)
{
  // This function provides the integration of the GJ formulation 24 a-e
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fran[3];

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double dtf = 0.5 * dt * ftm2v;
  double dtfm;
  double c1sqrt = sqrt(gjfc1);
  double c3sqrt = sqrt(gjfc3);
  double csq = sqrt(gjfc3 / gjfc1);
  double m, beta;

  // If user elected vhalf, v needs to be reassigned to onsite velocity for integration
  if (!osflag && lv_allocated) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        // lv is Eq. 24f from previous time step
        v[i][0] = lv[i][0];
        v[i][1] = lv[i][1];
        v[i][2] = lv[i][2];
      }
  }

  compute_target();
  if (tbiasflag) temperature->compute_scalar();

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
        m = rmass[i];
        beta = tsqrt * sqrt(2.0 * dt * m * boltz / t_period / mvv2e) / ftm2v;

        fran[0] = beta * random->gaussian();
        fran[1] = beta * random->gaussian();
        fran[2] = beta * random->gaussian();

        // First integration delivers Eq. 24a and 24b:
        dtfm = dtf / m;
        v[i][0] += csq * dtfm * f[i][0];
        v[i][1] += csq * dtfm * f[i][1];
        v[i][2] += csq * dtfm * f[i][2];
        x[i][0] += 0.5 * csq * dt * v[i][0];
        x[i][1] += 0.5 * csq * dt * v[i][1];
        x[i][2] += 0.5 * csq * dt * v[i][2];

        if (tbiasflag) temperature->remove_bias(i, v[i]);

        // Calculate Eq. 24c:
        lv[i][0] = c1sqrt * v[i][0] + ftm2v * (c3sqrt / (2.0 * m)) * fran[0];
        lv[i][1] = c1sqrt * v[i][1] + ftm2v * (c3sqrt / (2.0 * m)) * fran[1];
        lv[i][2] = c1sqrt * v[i][2] + ftm2v * (c3sqrt / (2.0 * m)) * fran[2];

        // Calculate Eq. 24d
        v[i][0] = (gjfc2 / c1sqrt) * lv[i][0] + ftm2v * csq * (0.5 / m) * fran[0];
        v[i][1] = (gjfc2 / c1sqrt) * lv[i][1] + ftm2v * csq * (0.5 / m) * fran[1];
        v[i][2] = (gjfc2 / c1sqrt) * lv[i][2] + ftm2v * csq * (0.5 / m) * fran[2];

        if (tbiasflag) temperature->restore_bias(i, v[i]);
        if (tbiasflag) temperature->restore_bias(i, lv[i]);

        // Calculate Eq. 24e. Final integrator then calculates Eq. 24f after force update.
        x[i][0] += 0.5 * csq * dt * v[i][0];
        x[i][1] += 0.5 * csq * dt * v[i][1];
        x[i][2] += 0.5 * csq * dt * v[i][2];
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
        m = mass[type[i]];
        beta = tsqrt * sqrt(2.0 * dt * m * boltz / t_period / mvv2e) / ftm2v;

        fran[0] = beta * random->gaussian();
        fran[1] = beta * random->gaussian();
        fran[2] = beta * random->gaussian();

        // First integration delivers Eq. 24a and 24b:
        dtfm = dtf / m;
        v[i][0] += csq * dtfm * f[i][0];
        v[i][1] += csq * dtfm * f[i][1];
        v[i][2] += csq * dtfm * f[i][2];
        x[i][0] += 0.5 * csq * dt * v[i][0];
        x[i][1] += 0.5 * csq * dt * v[i][1];
        x[i][2] += 0.5 * csq * dt * v[i][2];

        if (tbiasflag) temperature->remove_bias(i, v[i]);

        // Calculate Eq. 24c:
        lv[i][0] = c1sqrt * v[i][0] + ftm2v * (c3sqrt / (2.0 * m)) * fran[0];
        lv[i][1] = c1sqrt * v[i][1] + ftm2v * (c3sqrt / (2.0 * m)) * fran[1];
        lv[i][2] = c1sqrt * v[i][2] + ftm2v * (c3sqrt / (2.0 * m)) * fran[2];

        // Calculate Eq. 24d
        v[i][0] = (gjfc2 / c1sqrt) * lv[i][0] + ftm2v * csq * (0.5 / m) * fran[0];
        v[i][1] = (gjfc2 / c1sqrt) * lv[i][1] + ftm2v * csq * (0.5 / m) * fran[1];
        v[i][2] = (gjfc2 / c1sqrt) * lv[i][2] + ftm2v * csq * (0.5 / m) * fran[2];

        if (tbiasflag) temperature->restore_bias(i, v[i]);
        if (tbiasflag) temperature->restore_bias(i, lv[i]);

        // Calculate Eq. 24e. Final integrator then calculates Eq. 24f after force update.
        x[i][0] += 0.5 * csq * dt * v[i][0];
        x[i][1] += 0.5 * csq * dt * v[i][1];
        x[i][2] += 0.5 * csq * dt * v[i][2];
      }
    }
  }
}

void FixGJF::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfm;
  double dtf = 0.5 * update->dt * force->ftm2v;
  double csq = sqrt(gjfc3 / gjfc1);

  // Calculate Eq. 24f.
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += csq * dtfm * f[i][0];
        v[i][1] += csq * dtfm * f[i][1];
        v[i][2] += csq * dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += csq * dtfm * f[i][0];
        v[i][1] += csq * dtfm * f[i][1];
        v[i][2] += csq * dtfm * f[i][2];
      }
  }

  lv_allocated = 1;
}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixGJF::compute_target()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // if variable temp, evaluate variable, wrap with clear/add
  // reallocate tforce array if necessary

  if (tstyle == CONSTANT) {
    t_target = t_start + delta * (t_stop - t_start);
    tsqrt = sqrt(t_target);
  } else {
    modify->clearstep_compute();
    if (tstyle == EQUAL) {
      t_target = input->variable->compute_equal(tvar);
      if (t_target < 0.0) error->one(FLERR, "Fix gjf variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (atom->nmax > maxatom) {
        maxatom = atom->nmax;
        memory->destroy(tforce);
        memory->create(tforce, maxatom, "gjf:tforce");
      }
      input->variable->compute_atom(tvar, igroup, tforce, 1, 0);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          if (tforce[i] < 0.0) error->one(FLERR, "Fix gjf variable returned negative temperature");
    }
    modify->addstep_compute(update->ntimestep + 1);
  }
}

/* ----------------------------------------------------------------------
   select velocity for GJ
------------------------------------------------------------------------- */

void FixGJF::end_of_step()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // After the final integrator delivers 24f, either the on-site or half-step
  // velocity is used in remaining simulation tasks, depending on user input
  double tmp[3];
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      // v is Eq. 24f
      tmp[0] = v[i][0];
      tmp[1] = v[i][1];
      tmp[2] = v[i][2];
      // Move on with half-step velocity
      v[i][0] = lv[i][0];
      v[i][1] = lv[i][1];
      v[i][2] = lv[i][2];
      // store Eq. 24f in lv for next timestep
      lv[i][0] = tmp[0];
      lv[i][1] = tmp[1];
      lv[i][2] = tmp[2];
    }
}

// clang-format on
/* ---------------------------------------------------------------------- */

void FixGJF::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

void FixGJF::reset_dt()
{
  // Complete set of thermostats is given in Gronbech-Jensen, Molecular Physics, 118 (2020)
  switch (GJmethod) {
    case 1:
      gjfc2 = (1.0 - update->dt / 2.0 / t_period) / (1.0 + update->dt / 2.0 / t_period);
      break;
    case 2:
      gjfc2 = exp(-update->dt / t_period);
      break;
    case 3:
      gjfc2 = 1.0 - update->dt / t_period;
      break;
    case 4:
      gjfc2 = (sqrt(1.0 + 4.0 * (update->dt / t_period)) - 1.0) / (2.0 * update->dt / t_period);
      break;
    case 5:
      gjfc2 = 1.0 / (1.0 + update->dt / t_period);
      break;
    case 6:
      gjfc2 =
          (1.0 / (1.0 + update->dt / 2.0 / t_period)) * (1.0 / (1.0 + update->dt / 2.0 / t_period));
      break;
    case 7:    // provided in Finkelstein (2021)
      update->dt = (1.0 + gjfc2) / (1.0 - gjfc2) * log(gjfc2) * log(gjfc2) * 0.5 * t_period;
      break;
    case 8:    // provided in Gronbech-Jensen (2024)
      gjfc2 = sqrt((update->dt / t_period) * (update->dt / t_period) + 1.0) - update->dt / t_period;
      break;
    case 0:
      gjfc2 = 0.0;
      break;
    default:
      error->all(FLERR, "Fix gjf method not found");
      break;
  }
  gjfc1 = (1.0 + gjfc2) / 2.0;
  gjfc3 = (1.0 - gjfc2) * t_period / update->dt;
}

/* ---------------------------------------------------------------------- */

int FixGJF::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR, "Could not find fix_modify temperature compute ID: {}", id_temp);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature compute {} does not compute temperature", id_temp);
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group: {} vs {}",
                     group->names[igroup], group->names[temperature->igroup]);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixGJF::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "t_target") == 0) { return &t_target; }
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of tally array
------------------------------------------------------------------------- */

double FixGJF::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) atom->nmax * 3 * sizeof(double);
  if (tforce) bytes += (double) atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for lv
------------------------------------------------------------------------- */

void FixGJF::grow_arrays(int nmax)
{
  memory->grow(lv, nmax, 3, "fix_gjf:lv");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixGJF::copy_arrays(int i, int j, int /*delflag*/)
{
  lv[j][0] = lv[i][0];
  lv[j][1] = lv[i][1];
  lv[j][2] = lv[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixGJF::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = lv[i][0];
  buf[n++] = lv[i][1];
  buf[n++] = lv[i][2];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixGJF::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  lv[nlocal][0] = buf[n++];
  lv[nlocal][1] = buf[n++];
  lv[nlocal][2] = buf[n++];
  return n;
}
