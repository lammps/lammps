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

/* ----------------------------------------------------------------------
   Contributing author: Vaibhav Palkar (Kuksenok Lab, Clemson University)
   based on the pair srp code by Timothy Sirk (ARL)

This pair style interfaces the pair style srp with the stochastic reaction
fixes bond/break and bond/create by updating pseudo beads corresponding to
bonds as bond breaking and formation takes place. This is useful in
simulation of reactions in polymers with soft potentials such as DPD.

See the doc page for pair_style srp/react command for usage instructions.

There is an example script for this package in examples/PACKAGES/srp_react/.
------------------------------------------------------------------------- */

#include "pair_srp_react.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "fix_srp_react.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static const char cite_srpreact[] =
  "pair srp/react style: doi:10.1021/acs.jpcb.1c09570\n\n"
  "@Article{palkar2022\n"
  " author = {Palkar, Vaibhav and Kuksenok, Olga},\n"
  " title = {Controlling Degradation and Erosion of Polymer Networks: Insights from Mesoscale Modeling},\n"
  " journal = {J.~Phys.\\ Chem.~B},\n"
  " year = 2022,\n"
  " volume = 126,\n"
  " number = 1,\n"
  " pages = {336--346}\n"
  "}\n\n";

static int srp_instance = 0;

/* ----------------------------------------------------------------------
 constructor
 ---------------------------------------------------------------------- */

PairSRPREACT::PairSRPREACT(LAMMPS *lmp) :
  PairSRP(lmp), idbreak(nullptr), idcreate(nullptr), bond_break(false), bond_create(false)
{

  if (lmp->citeme) lmp->citeme->add(cite_srpreact);

  // pair srp/react has its own fix, hence delete fix srp instance
  // created in the constructor of pair srp
  for (auto &ifix : modify->get_fix_by_style("SRP"))
    modify->delete_fix(ifix->id);

  // similar to fix SRP, create fix SRP REACT instance here with unique fix id
  f_srp = (FixSRPREACT *) modify->add_fix(fmt::format("{:02d}_FIX_SRP_REACT all SRPREACT",srp_instance));
  ++srp_instance;
}

PairSRPREACT::~PairSRPREACT()
{
  delete[] idbreak;
  delete[] idcreate;
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSRPREACT::settings(int narg, char **arg)
{
  if (narg < 3 || narg > 8)
    error->all(FLERR,"Illegal pair_style command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair_style srp/react requires atom IDs");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  // wildcard
  if (strcmp(arg[1],"*") == 0) {
    btype = 0;
  } else {
    btype = utils::inumeric(FLERR,arg[1],false,lmp);
    if ((btype > atom->nbondtypes) || (btype <= 0))
      error->all(FLERR,"Illegal pair_style command");
  }

  // settings
  midpoint = false;
  min = false;

  if (strcmp(arg[2],"min") == 0) min = true;
  else if (strcmp(arg[2],"mid") == 0) midpoint = true;
  else
    error->all(FLERR,"Illegal pair_style command");

  // default for bond/break and bond/create settings
  bond_create=false;
  bond_break=false;
  idbreak = nullptr;
  idcreate= nullptr;

  // find whether id is of bond/break or bond/create
  const char *reactid = arg[3];
  if (utils::strmatch(modify->get_fix_by_id(reactid)->style,"^bond/break")) {
    bond_break = true;
    idbreak = utils::strdup(reactid);
  } else if (utils::strmatch(modify->get_fix_by_id(reactid)->style,"^bond/create")) {
    bond_create = true;
    idcreate = utils::strdup(reactid);
  } else error->all(FLERR,"Illegal pair_style command");

  int iarg = 4;
  // default exclude 1-2
  // scaling for 1-2, etc not supported
  exclude = 1;

  // use last atom type by default for bond particles
  bptype = atom->ntypes;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair srp command");
      exclude = utils::logical(FLERR, arg[iarg+1], false, lmp);
      if (min && !exclude) error->all(FLERR,"Illegal exclude option in pair srp command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"bptype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair srp command");
      bptype = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if ((bptype < 1) || (bptype > atom->ntypes))
        error->all(FLERR,"Illegal bond particle type for srp");
      iarg += 2;
   } else error->all(FLERR,"Illegal pair srp command");
  }

  // reset cutoffs if explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= bptype; i++)
      for (j = i; j <= bptype; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSRPREACT::init_style()
{
  if (!force->newton_pair)
    error->all(FLERR,"Pair srp/react requires newton pair on");

  // verify that fix SRP is still defined and has not been changed.

  if (strcmp(f_srp->style,"SRPREACT") != 0)
    error->all(FLERR,"Fix SRPREACT has been changed unexpectedly");

  if (comm->me == 0)
    utils::logmesg(lmp,"Using type {} for bond particles\n",bptype);

  // set bond and bond particle types in fix srp
  // bonds of this type will be represented by bond particles
  // if bond type is 0, then all bonds have bond particles
  // btype = bond type

  char c0[20];
  char* arg0[2];
  sprintf(c0, "%d", btype);
  arg0[0] = (char *) "btype";
  arg0[1] = c0;
  f_srp->modify_params(2, arg0);

  // bptype = bond particle type
  sprintf(c0, "%d", bptype);
  arg0[0] = (char *) "bptype";
  arg0[1] = c0;
  f_srp->modify_params(2, arg0);

  // if using fix bond/break, set id of fix bond/break in fix srp
  // idbreak = id of fix bond break
  if (bond_break) {
    sprintf(c0, "%s", idbreak);
    arg0[0] = (char *) "bond/break";
    arg0[1] = c0;
    f_srp->modify_params(2, arg0);
  }

  // if using fix bond/create, set id of fix bond/create in fix srp
  // idcreate = id of fix bond break
  if (bond_create) {
    sprintf(c0, "%s", idcreate);
    arg0[0] = (char *) "bond/create";
    arg0[1] = c0;
    f_srp->modify_params(2, arg0);
  }

  // bond particles do not contribute to energy or virial
  // bond particles do not belong to group all
  // but thermo normalization is by nall
  // therefore should turn off normalization
  char *arg1[2];
  arg1[0] = (char *) "norm";
  arg1[1] = (char *) "no";
  output->thermo->modify_params(2, arg1);
  if (comm->me == 0) error->message(FLERR,"Thermo normalization turned off by pair srp/react");

  neighbor->request(this,instance_me);
}
