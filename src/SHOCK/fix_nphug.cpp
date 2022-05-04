// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nphug.h"

#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ISO,ANISO,TRICLINIC}; // same as fix_nh.cpp

/* ---------------------------------------------------------------------- */

FixNPHug::FixNPHug(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg), pe(nullptr), id_pe(nullptr)
{
  // Prevent masses from being updated every timestep

  eta_mass_flag = 0;
  omega_mass_flag = 0;
  etap_mass_flag = 0;

  // extend vector of base-class computes

  size_vector += 3;

  // turn off deviatoric flag and remove strain energy from vector

  deviatoric_flag = 0;
  size_vector -= 1;

  // use initial state as reference state

  v0_set = p0_set = e0_set = 0;

  // check pressure settings

  if (p_start[0] != p_stop[0] ||
      p_start[1] != p_stop[1] ||
      p_start[2] != p_stop[2])
    error->all(FLERR,"Pstart and Pstop must have the same value");

  // uniaxial = 0 means hydrostatic compression
  // uniaxial = 1 means uniaxial compression
  //      in x, y, or z (idir = 0, 1, or 2)

  // isotropic hydrostatic compression

  if (pstyle == ISO) {
    uniaxial = 0;

    // anisotropic compression

  } else if (pstyle == ANISO) {

    // anisotropic hydrostatic compression

    if (p_start[0] == p_start[1] &&
        p_start[0] == p_start[2] )
      uniaxial = 0;

    // uniaxial compression

    else if (p_flag[0] == 1 && p_flag[1] == 0
        && p_flag[2] == 0) {
      uniaxial = 1;
      idir = 0;
    } else if (p_flag[0] == 0 && p_flag[1] == 1
           && p_flag[2] == 0) {
      uniaxial = 1;
      idir = 1;
    } else if (p_flag[0] == 0 && p_flag[1] == 0
               && p_flag[2] == 1) {
      uniaxial = 1;
      idir = 2;

    } else error->all(FLERR,"Specified target stress must be uniaxial or hydrostatic");

    // triclinic hydrostatic compression

  } else if (pstyle == TRICLINIC) {

    if (p_start[0] == p_start[1] &&
        p_start[0] == p_start[2] &&
        p_start[3] == 0.0 &&
        p_start[4] == 0.0 &&
        p_start[5] == 0.0 )
      uniaxial = 0;

    else error->all(FLERR,"For triclinic deformation, specified target stress must be hydrostatic");
  }

  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nphug");
  if (!pstat_flag)
    error->all(FLERR,"Pressure control must be used with fix nphug");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  // and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id)+"_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));
  tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id)+"_press");
  modify->add_compute(fmt::format("{} all pressure {}",id_press,id_temp));
  pcomputeflag = 1;

  // create a new compute potential energy compute

  id_pe = utils::strdup(std::string(id)+"_pe");
  modify->add_compute(fmt::format("{} all pe",id_pe));
  peflag = 1;
}

/* ---------------------------------------------------------------------- */

FixNPHug::~FixNPHug()
{
  // temp and press computes handled by base class
  // delete pe compute

  if (peflag) modify->delete_compute(id_pe);
  delete [] id_pe;
}

/* ---------------------------------------------------------------------- */

void FixNPHug::init()
{
  // Call base class init()

  FixNH::init();

  // set pe ptr

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Potential energy ID for fix nvt/nph/npt does not exist");
  pe = modify->compute[icompute];
}


/* ----------------------------------------------------------------------
   compute initial state before integrator starts
------------------------------------------------------------------------- */

void FixNPHug::setup(int vflag)
{
  FixNH::setup(vflag);

  if (v0_set == 0) {
    v0 = compute_vol();
    v0_set = 1;
  }

  if (p0_set == 0) {
    p0_set = 1;
    if (uniaxial == 1)
      p0 = p_current[idir];
    else
      p0 = (p_current[0]+p_current[1]+p_current[2])/3.0;
  }

  if (e0_set == 0) {
    e0 = compute_etotal();
    e0_set = 1;
  }

  double masstot = group->mass(igroup);
  rho0 = nktv2p*force->mvv2e*masstot/v0;

  t_target = 0.01;
  ke_target = tdof*boltz*t_target;

  pe->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   compute target temperature and kinetic energy
-----------------------------------------------------------------------*/

void FixNPHug::compute_temp_target()
{
  t_target = t_current + compute_hugoniot();
  ke_target = tdof * boltz * t_target;
  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

double FixNPHug::compute_etotal()
{
  if (!pe) return 0.0;

  double epot,ekin,etot;
  epot = pe->compute_scalar();
  ekin = temperature->compute_scalar();
  ekin *= 0.5 * tdof * force->boltz;
  etot = epot+ekin;
  return etot;
}

/* ---------------------------------------------------------------------- */

double FixNPHug::compute_vol()
{
  if (domain->dimension == 3)
    return domain->xprd * domain->yprd * domain->zprd;
  else
    return domain->xprd * domain->yprd;
}

/* ----------------------------------------------------------------------
   Computes the deviation of the current point
   from the Hugoniot in temperature units.
------------------------------------------------------------------------- */

double FixNPHug::compute_hugoniot()
{
  if (!temperature) return 0.0;

  double v,e,p;
  double dhugo;

  e = compute_etotal();

  temperature->compute_vector();


  if (uniaxial == 1) {
    pressure->compute_vector();
    p = pressure->vector[idir];
  } else
    p = pressure->compute_scalar();

  v = compute_vol();

  dhugo = (0.5 * (p + p0 ) * ( v0 - v)) /
    force->nktv2p + e0 - e;

  dhugo /= tdof * boltz;

  return dhugo;
}

/* ----------------------------------------------------------------------
   Compute shock velocity is distance/time units
------------------------------------------------------------------------- */

double FixNPHug::compute_us()
{
  if (!temperature) return 0.0;

  double v,p;
  double eps,us;

  temperature->compute_vector();

  if (uniaxial == 1) {
    pressure->compute_vector();
    p = pressure->vector[idir];
  } else
    p = pressure->compute_scalar();

  v = compute_vol();

  // Us^2 = (p-p0)/(rho0*eps)

  eps = 1.0 - v/v0;
  if (eps < 1.0e-10) us = 0.0;
  else if (p < p0) us = 0.0;
  else us = sqrt((p-p0)/(rho0*eps));

  return us;
}

/* ----------------------------------------------------------------------
   Compute particle velocity is distance/time units
------------------------------------------------------------------------- */

double FixNPHug::compute_up()
{
  double v;
  double eps,us,up;

  v = compute_vol();
  us = compute_us();

  // u = eps*Us

  eps = 1.0 - v/v0;
  up = us*eps;

  return up;
}

/* ----------------------------------------------------------------------- */

// look for index in local class
// if index not found, look in base class

double FixNPHug::compute_vector(int n)
{
  int ilen;

  // n = 0: Hugoniot energy difference (temperature units)

  ilen = 1;
  if (n < ilen) return compute_hugoniot();
  n -= ilen;

  // n = 1: Shock velocity

  ilen = 1;
  if (n < ilen) return compute_us();
  n -= ilen;

  // n = 2: Particle velocity

  ilen = 1;
  if (n < ilen) return compute_up();
  n -= ilen;

  // index not found, look in base class

  return FixNH::compute_vector(n);
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

int FixNPHug::pack_restart_data(double *list)
{
  int n = 0;

  list[n++] = e0;
  list[n++] = v0;
  list[n++] = p0;

  // call the base class function

  n += FixNH::pack_restart_data(list+n);

  return n;
}

/* ----------------------------------------------------------------------
   calculate the number of data to be packed
------------------------------------------------------------------------- */

int FixNPHug::size_restart_global()
{
  int nsize = 3;

  // call the base class function

  nsize += FixNH::size_restart_global();

  return nsize;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixNPHug::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;
  e0 = list[n++];
  v0 = list[n++];
  p0 = list[n++];

  e0_set = 1;
  v0_set = 1;
  p0_set = 1;

  // call the base class function

  buf += n*sizeof(double);
  FixNH::restart(buf);

}

/* ---------------------------------------------------------------------- */

int FixNPHug::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"e0") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix nphug command");
    e0 = utils::numeric(FLERR,arg[1],false,lmp);
    e0_set = 1;
    return 2;
  } else if (strcmp(arg[0],"v0") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix nphug command");
    v0 = utils::numeric(FLERR,arg[1],false,lmp);
    v0_set = 1;
    return 2;
  } else if (strcmp(arg[0],"p0") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix nphug command");
    p0 = utils::numeric(FLERR,arg[1],false,lmp);
    p0_set = 1;
    return 2;
  }

  return 0;
}
