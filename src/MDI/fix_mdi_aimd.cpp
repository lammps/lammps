/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "error.h"
#include "fix_mdi_aimd.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDIAimd::FixMDIAimd(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix mdi/aimd command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // connect to MDI engine

  MDI_Accept_communicator(&engine);
}

/* ---------------------------------------------------------------------- */

FixMDIAimd::~FixMDIAimd()
{
  // send exit command to engine

  MDI_Send_command("EXIT",engine);
}

/* ---------------------------------------------------------------------- */

int FixMDIAimd::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to request energy
------------------------------------------------------------------------- */

void FixMDIAimd::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::post_force(int vflag)
{
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // send current coords to MDI engine
  // NOTE: need to gather them to proc 0

  MDI_Send_command(">COORDS",engine);
  MDI_Send(&atom->x[0][0],3*atom->nlocal,MDI_DOUBLE,engine);

  // trigger engine to evaluate forces,energy,pressure for current system

  MDI_Send_command("EVAL",engine);

  // request energy from MDI engine

  if (eflag_global) {
    MDI_Send_command("<PE",engine);
    MDI_Recv(&engine_energy,1,MDI_DOUBLE,engine);
  }

  // request pressure tensor from MDI engine, convert to virial

  if (vflag_global) {
    double ptensor[6];
    MDI_Send_command("<PTENSOR",engine);
    MDI_Recv(ptensor,6,MDI_DOUBLE,engine);

    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (int i = 0; i < 6; i++)
      virial[i] = ptensor[i] * volume / force->nktv2p;
  }

  // request forces from MDI engine
  // NOTE: need to scatter to procs, then add to my local forces

  MDI_Send_command("<FORCES",engine);
  MDI_Recv(&atom->f[0][0],3*atom->nlocal,MDI_DOUBLE,engine);
}

/* ---------------------------------------------------------------------- */

void FixMDIAimd::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy from MDI engine
------------------------------------------------------------------------- */

double FixMDIAimd::compute_scalar()
{
  return engine_energy;
}
