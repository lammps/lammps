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

#include "fix_balance.h"

#include "atom.h"
#include "balance.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "irregular.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "rcb.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SHIFT,BISECTION};

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), balance(nullptr), irregular(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal fix balance command");

  box_change = BOX_CHANGE_DOMAIN;
  pre_exchange_migrate = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  // parse required arguments

  int dimension = domain->dimension;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Illegal fix balance command");
  thresh = utils::numeric(FLERR,arg[4],false,lmp);

  if (strcmp(arg[5],"shift") == 0) lbstyle = SHIFT;
  else if (strcmp(arg[5],"rcb") == 0) lbstyle = BISECTION;
  else error->all(FLERR,"Illegal fix balance command");

  int iarg = 5;
  if (lbstyle == SHIFT) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[iarg+1]) > Balance::BSTR_SIZE)
      error->all(FLERR,"Illegal fix balance command");
    strncpy(bstr,arg[iarg+1], Balance::BSTR_SIZE+1);
    nitermax = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
    if (nitermax <= 0) error->all(FLERR,"Illegal fix balance command");
    stopthresh = utils::numeric(FLERR,arg[iarg+3],false,lmp);
    if (stopthresh < 1.0) error->all(FLERR,"Illegal fix balance command");
    iarg += 4;

  } else if (lbstyle == BISECTION) {
    iarg++;
  }

  // error checks

  if (lbstyle == SHIFT) {
    int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Fix balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Fix balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Fix balance shift string is invalid");
    }
  }

  if (lbstyle == BISECTION && comm->style == Comm::BRICK)
    error->all(FLERR,"Fix balance rcb cannot be used with comm_style brick");

  // create instance of Balance class
  // if SHIFT, initialize it with params
  // process remaining optional args via Balance

  balance = new Balance(lmp);
  if (lbstyle == SHIFT) balance->shift_setup(bstr,nitermax,thresh);
  balance->options(iarg,narg,arg,0);
  wtflag = balance->wtflag;
  sortflag = balance->sortflag;

  if (balance->varflag && nevery == 0)
    error->all(FLERR,"Fix balance nevery = 0 cannot be used with weight var");

  // create instance of Irregular class

  irregular = new Irregular(lmp);

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;
  lastbalance = -1;
  next_reneighbor = -1;

  // compute initial outputs

  itercount = 0;
  pending = 0;
  imbfinal = imbprev = maxloadperproc = 0.0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::post_constructor()
{
  if (wtflag) balance->weight_storage(id);
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  balance->init_imbalance(1);
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup(int /*vflag*/)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // ensure atoms are in current box & update box via shrink-wrap
  // has to be be done before rebalance() invokes Irregular::migrate_atoms()
  //   since it requires atoms be inside simulation box
  //   even though pbc() will be done again in Verlet::run()
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshold exceeded

  balance->set_weights();
  imbnow = balance->imbalance_factor(maxloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // ensure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshold exceeded
  // if weight variable is used, wrap weight setting in clear/add compute

  if (balance->varflag) modify->clearstep_compute();
  balance->set_weights();
  if (balance->varflag) modify->addstep_compute(update->ntimestep + nevery);

  imbnow = balance->imbalance_factor(maxloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occurred
------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_factor(maxloadperproc);
  pending = 0;

  // set disable = 1, so weights no longer migrate with atoms

  if (wtflag) balance->fixstore->disable = 1;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::rebalance()
{
  imbprev = imbnow;

  // invoke balancer and reset comm->uniform flag

  int *sendproc;
  if (lbstyle == SHIFT) {
    itercount = balance->shift();
    comm->layout = Comm::LAYOUT_NONUNIFORM;
  } else if (lbstyle == BISECTION) {
    sendproc = balance->bisection();
    comm->layout = Comm::LAYOUT_TILED;
  }

  // reset proc sub-domains

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();

  // check and warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in comm->exchange()

  domain->subbox_too_small_check(neighbor->skin);

  // output of new decomposition

  if (balance->outflag) balance->dumpout(update->ntimestep);

  // move atoms to new processors via irregular()
  // for non-RCB only needed if migrate_check() says an atom moves too far
  // else allow caller's comm->exchange() to do it
  // set disable = 0, so weights migrate with atoms
  //   important to delay disable = 1 until after pre_neighbor imbfinal calc
  //   b/c atoms may migrate again in comm->exchange()
  // sortflag determines whether irregular sorts its
  //   comm messages for reproducibility
  //   if not, message order is random, atom order is non-deterministic

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (wtflag) balance->fixstore->disable = 0;
  if (lbstyle == BISECTION) irregular->migrate_atoms(sortflag,1,sendproc);
  else if (irregular->migrate_check()) irregular->migrate_atoms(sortflag);
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // notify all classes that store distributed grids
  // so they can adjust to new proc sub-domains

  modify->reset_grid();
  if (force->pair) force->pair->reset_grid();
  if (force->kspace) force->kspace->reset_grid();

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in comm->exchange()

  pending = 1;
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return maxloadperproc;
  if (i == 1) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = irregular->memory_usage();
  if (balance->rcb) bytes += balance->rcb->memory_usage();
  return bytes;
}
