/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com


LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov


Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.


See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "loadbalance.h"
#include "style_lb.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Loadbalance::Loadbalance(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Loadbalance::command(int narg, char **arg)
{
  if (narg < 1)
    error->all(FLERR,"Illegal loadbalance command, not enough arguments");

  if (strcmp(arg[0],"off") == 0) {
    if (domain->lbalance) delete domain->lbalance;
    domain->lbalance = NULL;
  } else {
    if (comm->nprocs == 1) {
      error->warning(FLERR,"Running in serial, loadbalance command has no effects");
      return;
    }


    if (0) return;

#define LB_CLASS
#define LBStyle(key,Class)						\
    else if (strcmp(arg[0],#key) == 0) domain->lbalance = new Class(lmp,narg,arg);
#include "style_lb.h"
#undef LB_CLASS

    else error->all(FLERR,"Illegal loadbalance command: Unknown style");
  }
}
