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

   Force scaling fix for gREM.
   Cite: https://doi.org/10.1063/1.3432176
   Cite: https://doi.org/10.1021/acs.jpcb.5b07614

------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Edyta Malolepsza (Broad Institute)
                         David Stelter (Boston University)
                         Tom Keyes (Boston University)
------------------------------------------------------------------------- */

#include "fix_grem.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixGrem::FixGrem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix grem command");

  scalar_flag = 1;
  extscalar = 0;
  global_freq = 1;

  scale_grem = 1.0;

  // tbath - temp of bath, the same as defined in thermostat

  lambda = utils::numeric(FLERR,arg[3],false,lmp);
  eta = utils::numeric(FLERR,arg[4],false,lmp);
  h0 = utils::numeric(FLERR,arg[5],false,lmp);
  id_nh = utils::strdup(arg[6]);

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all PRESSURE/GREM {} {}", id_press, id_temp, id));

  // create a new compute ke style
  // id = fix-ID + ke

  id_ke = utils::strdup(std::string(id) + "_ke");
  modify->add_compute(fmt::format("{} all ke",id_ke));

  // create a new compute pe style
  // id = fix-ID + pe

  id_pe = utils::strdup(std::string(id) + "_pe");
  modify->add_compute(fmt::format("{} all pe",id_pe));

  int ifix = modify->find_fix(id_nh);
  if (ifix < 0)
    error->all(FLERR,"Fix id for nvt or npt fix does not exist");
  Fix *nh = modify->fix[ifix];

  pressflag = 0;
  int *p_flag = (int *)nh->extract("p_flag",ifix);
  if ((p_flag == nullptr) || (ifix != 1) || (p_flag[0] == 0)
      || (p_flag[1] == 0) || (p_flag[2] == 0)) {
    pressflag = 0;
  } else if ((p_flag[0] == 1) && (p_flag[1] == 1)
             && (p_flag[2] == 1) && (ifix == 1)) {
    pressflag = 1;
    char *modargs[2];
    modargs[0] = (char *) "press";
    modargs[1] = id_press;
    nh->modify_param(2,modargs);
  }
}

/* ---------------------------------------------------------------------- */

FixGrem::~FixGrem()
{
  // delete temperature, pressure and energies if fix created them

  modify->delete_compute(id_temp);
  modify->delete_compute(id_press);
  modify->delete_compute(id_ke);
  modify->delete_compute(id_pe);
  delete [] id_temp;
  delete [] id_press;
  delete [] id_ke;
  delete [] id_pe;
  delete [] id_nh;
}

/* ---------------------------------------------------------------------- */

int FixGrem::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGrem::init()
{

  if (domain->triclinic)
    error->all(FLERR,"Triclinic cells are not supported");

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature compute ID for fix grem does not exist");
  temperature = modify->compute[icompute];

  icompute = modify->find_compute(id_ke);
  if (icompute < 0)
    error->all(FLERR,"KE compute ID for fix grem does not exist");
  ke = modify->compute[icompute];

  icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"PE compute ID for fix grem does not exist");
  pe = modify->compute[icompute];

  int ifix = modify->find_fix(id_nh);
  if (ifix < 0)
    error->all(FLERR,"Fix id for nvt or npt fix does not exist");
  Fix *nh = modify->fix[ifix];

  auto t_start = (double *)nh->extract("t_start",ifix);
  auto t_stop = (double *)nh->extract("t_stop",ifix);
  if ((t_start != nullptr) && (t_stop != nullptr) && (ifix == 0)) {
    tbath = *t_start;
    if (*t_start != *t_stop)
      error->all(FLERR,"Thermostat temperature ramp not allowed");
  } else
    error->all(FLERR,"Problem extracting target temperature from fix nvt or npt");

  pressref = 0.0;
  if (pressflag) {
    int *p_flag = (int *)nh->extract("p_flag",ifix);
    auto p_start = (double *) nh->extract("p_start",ifix);
    auto p_stop = (double *) nh->extract("p_stop",ifix);
    if ((p_flag != nullptr) && (p_start != nullptr) && (p_stop != nullptr)
        && (ifix == 1)) {
      ifix = 0;
      pressref = p_start[0];
      if ((p_start[0] != p_stop[0]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[1] != p_stop[1]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[2] != p_stop[2]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[0] != p_start[1]) || (p_start[1] != p_start[2])) ++ifix;
      if ((p_flag[3] != 0) || (p_flag[4] != 0) || (p_flag[5] != 0)) ++ifix;
      if (ifix > 0)
        error->all(FLERR,"Unsupported pressure settings in fix npt");
    } else
      error->all(FLERR,"Problem extracting target pressure from fix npt");
  }
}

/* ---------------------------------------------------------------------- */

void FixGrem::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Run style 'respa' is not supported");
}

/* ---------------------------------------------------------------------- */

void FixGrem::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGrem::post_force(int /*vflag*/)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double tmpvolume = domain->xprd * domain->yprd * domain->zprd;
  double tmppe = pe->compute_scalar();
  // potential energy
  double tmpenthalpy = tmppe+pressref*tmpvolume/(force->nktv2p);

  double teffective = lambda+eta*(tmpenthalpy-h0);
  scale_grem = tbath/teffective;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] *= scale_grem;
      f[i][1] *= scale_grem;
      f[i][2] *= scale_grem;
    }
  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

double FixGrem::compute_scalar()
{
  return tbath / scale_grem;
}

/* ----------------------------------------------------------------------
   extract scale factor
------------------------------------------------------------------------- */

void *FixGrem::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"scale_grem") == 0) {
    return &scale_grem;
  }
  return nullptr;
}
