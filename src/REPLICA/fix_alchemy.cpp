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

#include "fix_alchemy.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAlchemy::FixAlchemy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), commbuf(nullptr)
{
  if (narg != 4) error->universe_all(FLERR, "Incorrect number of arguments for fix alchemy");
  if (universe->nworlds != 2) error->universe_all(FLERR, "Must use exactly two partitions");
  if (utils::strmatch(arg[3], "^v_"))
    id_lambda = arg[3] + 2;
  else
    error->universe_all(FLERR, "Must use variable as lambda argument to fix alchemy");

  lambda = epot[0] = epot[1] = epot[2] = 0.0;
  progress = 0;
  for (int i = 0; i < 6; ++i) pressure[i] = 0.0;

  no_change_box = 1;
  time_depend = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 1;
  nmax = 6;
  ivar = -1;
  sync_box = 0;

  // set up rank-to-rank communicator for inter-partition communication

  int color = comm->me;
  int key = universe->iworld;
  MPI_Comm_split(universe->uworld, color, key, &samerank);

  // check that we have the same domain decomposition on all ranks

  int my_nlocal[2] = {0, 0};
  int all_nlocal[2] = {0, 0};
  my_nlocal[universe->iworld] = atom->nlocal;
  MPI_Allreduce(my_nlocal, all_nlocal, 2, MPI_INT, MPI_SUM, samerank);
  int fail = (all_nlocal[0] == all_nlocal[1]) ? 0 : 1;
  int allfail = 0;
  MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (allfail)
    error->universe_all(FLERR,
                        "Number of atoms and domain decomposition must be the same "
                        "on all partitions");

  id_pe = std::string(id) + "_pe";
  pe = modify->add_compute(id_pe + " all pe");
  pe->addstep(update->ntimestep);
  id_temp = std::string(id) + "_temp";
  temp = modify->add_compute(id_temp + " all temp");
  temp->addstep(update->ntimestep);
  id_press = std::string(id) + "_press";
  press = modify->add_compute(id_press + " all pressure " + id_temp);
  press->addstep(update->ntimestep);
}

/* ---------------------------------------------------------------------- */

FixAlchemy::~FixAlchemy()
{
  MPI_Comm_free(&samerank);
  modify->delete_compute(id_pe);
  modify->delete_compute(id_temp);
  modify->delete_compute(id_press);
  memory->destroy(commbuf);
}

/* ---------------------------------------------------------------------- */

int FixAlchemy::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
   check consistency of owned atom count and ordering
   compare each pair of replica procs
   checked before each exchange of atom coords or forces
     to ensure the replicas have not become out of sync
---------------------------------------------------------------------- */

void FixAlchemy::check_consistency_atoms()
{
  // check that owned atom count is same for each pair of replica procs

  const int nlocal = atom->nlocal;
  int my_nlocal[2] = {0, 0};
  int all_nlocal[2] = {0, 0};
  my_nlocal[universe->iworld] = nlocal;
  MPI_Allreduce(my_nlocal, all_nlocal, 2, MPI_INT, MPI_SUM, samerank);

  int fail = (all_nlocal[0] == all_nlocal[1]) ? 0 : 1;
  int allfail = 0;
  MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (allfail) error->universe_all(FLERR, "Fix alchemy local atom count is inconsistent");

  // check that owned atom ordering is same for each pair of replica procs
  // re-use communication buffer for positions and forces

  tagint *tagbuf = (tagint *) commbuf;
  tagint *tag = atom->tag;
  if (universe->iworld == 0) {
    for (int i = 0; i < nlocal; ++i) tagbuf[i] = tag[i];
  }
  MPI_Bcast(tagbuf, nlocal, MPI_LMP_TAGINT, 0, samerank);

  fail = allfail = 0;
  if (universe->iworld > 0) {
    for (int i = 0; i < nlocal; ++i)
      if (tag[i] != tagbuf[i]) fail = 1;
  }
  MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (allfail) error->universe_all(FLERR, "Fix alchemy local atom ordering is inconsistent");
}

/* ----------------------------------------------------------------------
   force simulation box size and shape to be identical for 2 replicas
   invoked by post_integrate() after integration may have changed box
---------------------------------------------------------------------- */

static void synchronize_box(Domain *domain, MPI_Comm samerank)
{
  MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->yz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xy, 1, MPI_DOUBLE, 0, samerank);
  domain->set_global_box();
  domain->set_local_box();
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::init()
{
  int onenmax = MAX(nmax, 3 * atom->nmax);
  MPI_Allreduce(&onenmax, &nmax, 1, MPI_INT, MPI_MAX, universe->uworld);
  memory->destroy(commbuf);
  memory->create(commbuf, sizeof(double) * nmax, "alchemy:nmax");

  if (modify->get_fix_by_style("^balance").size() > 0)
    error->universe_all(FLERR, "Fix alchemy is not compatible with load balancing");

  if (modify->get_fix_by_style("^alchemy").size() > 1)
    error->universe_all(FLERR, "There may only one fix alchemy at a time");

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->universe_all(FLERR, "Must not use run style respa with fix alchemy");

  ivar = input->variable->find(id_lambda.c_str());
  if (ivar < 0)
    error->universe_one(FLERR, fmt::format("Fix alchemy variable {} does not exist", id_lambda));
  if (!input->variable->equalstyle(ivar))
    error->universe_one(FLERR, fmt::format("Fix alchemy variable {} is invalid style", id_lambda));
  lambda = input->variable->compute_equal(ivar);

  // synchronize box dimensions, determine if resync during run will be needed.

  synchronize_box(domain, samerank);

  sync_box = 0;
  for (auto ifix : modify->get_fix_list())
    if (ifix->box_change) sync_box = 1;
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::setup(int vflag)
{
  if (universe->me == 0) {
    progress = 0;
    auto msg = fmt::format("Starting alchemical run\n");
    if (universe->uscreen) fmt::print(universe->uscreen, msg);
    if (universe->ulogfile) fmt::print(universe->ulogfile, msg);
  }

  // recheck domain decomposition, atom ordering, and synchronize positions

  post_integrate();

  // mix initial forces

  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_integrate()
{
  // check owned atom count and ordering between replicas

  check_consistency_atoms();

  // synchronize atom positions

  const int nall = atom->nlocal;
  MPI_Bcast(&atom->x[0][0], 3 * nall, MPI_DOUBLE, 0, samerank);

  // synchronize box dimensions, if needed

  if (sync_box) synchronize_box(domain, samerank);
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_force(int /*vflag*/)
{
  // grow commbuf if necessary

  if (3 * atom->nmax > nmax) {
    nmax = 3 * atom->nmax;
    memory->grow(commbuf, sizeof(double) * atom->nmax, "alchemy:commbuf");
  }

  // check owned atom count and ordering between replicas

  check_consistency_atoms();

  // evaluate lambda variable

  lambda = input->variable->compute_equal(ivar);

  // sum forces multiplied by lambda across 2 replicas

  const int nall = 3 * atom->nlocal;
  double *f = &atom->f[0][0];
  for (int i = 0; i < nall; ++i) commbuf[i] = f[i] * lambda;
  MPI_Allreduce(commbuf, f, nall, MPI_DOUBLE, MPI_SUM, samerank);

  // sum up potential energy

  const double scalefac = 1.0 / comm->nprocs;
  commbuf[0] = commbuf[1] = commbuf[2] = 0.0;
  commbuf[universe->iworld] = scalefac * pe->compute_scalar();
  commbuf[2] = lambda * scalefac * pe->compute_scalar();
  MPI_Allreduce(commbuf, epot, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);
  pe->addstep(update->ntimestep + 1);

  // sum up pressure

  press->compute_vector();
  for (int i = 0; i < 6; ++i) commbuf[i] = lambda * scalefac * press->vector[i];
  MPI_Allreduce(commbuf, pressure, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
  press->addstep(update->ntimestep + 1);

  // print progress info to universe screen/logfile

  if (universe->me == 0) {
    double delta = update->ntimestep - update->beginstep;
    if ((delta != 0.0) && (update->beginstep != update->endstep))
      delta /= update->endstep - update->beginstep;
    int status = static_cast<int>(delta * 100.0);
    if ((status / 10) > (progress / 10)) {
      progress = status;
      auto msg = fmt::format("  Alchemical run progress: {:>3d}%\n", progress);
      if (universe->uscreen) fmt::print(universe->uscreen, msg);
      if (universe->ulogfile) fmt::print(universe->ulogfile, msg);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixAlchemy::compute_scalar()
{
  return lambda;
}

/* ---------------------------------------------------------------------- */

double FixAlchemy::compute_vector(int n)
{
  return epot[n];
}

/* ---------------------------------------------------------------------- */

void *FixAlchemy::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "lambda") == 0) { return &lambda; }
  if (strcmp(str, "pe") == 0) { return &epot[2]; }
  dim = 1;
  if (strcmp(str, "pressure") == 0) { return pressure; }
  return nullptr;
}
