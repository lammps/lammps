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

#include "fix_halt.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "input.h"
#include "modify.h"
#include "neighbor.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { BONDMAX, TLIMIT, DISKFREE, VARIABLE };
enum { LT, LE, GT, GE, EQ, NEQ, XOR };
enum { HARD, SOFT, CONTINUE };
enum { NOMSG = 0, YESMSG = 1 };
static constexpr int UTAG = 999;

/* ---------------------------------------------------------------------- */

FixHalt::FixHalt(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idvar(nullptr), dlimit_path(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix halt", error);
  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix halt command: nevery must be > 0");

  // comparison args

  int iarg = 4;

  if (strcmp(arg[iarg], "tlimit") == 0) {
    attribute = TLIMIT;
  } else if (strcmp(arg[iarg], "diskfree") == 0) {
    attribute = DISKFREE;
    dlimit_path = utils::strdup(".");
  } else if (strcmp(arg[iarg], "bondmax") == 0) {
    attribute = BONDMAX;
  } else if (utils::strmatch(arg[iarg], "^v_")) {
    ArgInfo argi(arg[iarg], ArgInfo::VARIABLE);

    if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_type() == ArgInfo::NONE) ||
        (argi.get_dim() != 0))
      error->all(FLERR, iarg, "Invalid fix halt attribute {}", arg[iarg]);

    attribute = VARIABLE;
    idvar = argi.copy_name();
    ivar = input->variable->find(idvar);

    if (ivar < 0) error->all(FLERR, iarg, "Could not find fix halt variable name {}", idvar);
    if (input->variable->equalstyle(ivar) == 0)
      error->all(FLERR, iarg, "Fix halt variable is not equal-style variable");
  } else {
    error->all(FLERR, iarg, "Unknown fix halt keyword {}", arg[iarg]);
  }

  // clang-format off
  ++iarg;
  if (strcmp(arg[iarg],"<") == 0) operation = LT;
  else if (strcmp(arg[iarg],"<=") == 0) operation = LE;
  else if (strcmp(arg[iarg],">") == 0) operation = GT;
  else if (strcmp(arg[iarg],">=") == 0) operation = GE;
  else if (strcmp(arg[iarg],"==") == 0) operation = EQ;
  else if (strcmp(arg[iarg],"!=") == 0) operation = NEQ;
  else if (strcmp(arg[iarg],"|^") == 0) operation = XOR;
  else error->all(FLERR, iarg, "Invalid fix halt operator {}", arg[iarg]);

  ++iarg;
  value = utils::numeric(FLERR, arg[iarg], false, lmp);

  // parse optional args

  eflag = SOFT;
  msgflag = YESMSG;
  uflag = NOMSG;
  ++iarg;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "error") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt error", error);
      if (strcmp(arg[iarg + 1], "hard") == 0) eflag = HARD;
      else if (strcmp(arg[iarg + 1], "soft") == 0) eflag = SOFT;
      else if (strcmp(arg[iarg + 1], "continue") == 0) eflag = CONTINUE;
      else error->all(FLERR, iarg + 1, "Unknown fix halt error condition {}", arg[iarg]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "message") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt message", error);
      msgflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "universe") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt universe", error);
      uflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "path") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt path", error);
      ++iarg;
      delete[] dlimit_path;
      // strip off outer quotes, if present
      int len = strlen(arg[iarg]) + 1;
      if (((arg[iarg][0] == '"') || (arg[iarg][0] == '\'')) &&
          (arg[iarg][0] == arg[iarg][len - 2])) {
        arg[iarg][len - 2] = '\0';
        dlimit_path = utils::strdup(arg[iarg] + 1);
      } else
        dlimit_path = utils::strdup(arg[iarg]);
      ++iarg;
    } else error->all(FLERR, "Unknown fix halt keyword {}", arg[iarg]);
  }
  // clang-format on

  // add nfirst to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  if (attribute == VARIABLE) {
    const bigint nfirst = (update->ntimestep / nevery) * nevery + nevery;
    modify->addstep_compute_all(nfirst);
  }
}

/* ---------------------------------------------------------------------- */

FixHalt::~FixHalt()
{
  delete[] idvar;
  delete[] dlimit_path;
}

/* ---------------------------------------------------------------------- */

int FixHalt::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHalt::init()
{
  // set ivar from current variable list

  if (attribute == VARIABLE) {
    ivar = input->variable->find(idvar);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Could not find fix halt variable {}", idvar);
    if (input->variable->equalstyle(ivar) == 0)
      error->all(FLERR, Error::NOLASTLINE, "Fix halt variable {} is not equal-style", idvar);
  }

  // settings used by TLIMIT

  nextstep = (update->ntimestep / nevery) * nevery + nevery;
  thisstep = -1;
  tratio = 0.5;

  // check if disk limit is supported

  if (attribute == DISKFREE) {
    if (!dlimit_path || platform::disk_free(dlimit_path) < 0.0)
      error->all(FLERR, Error::NOLASTLINE, "Disk limit not supported by OS or illegal path");
  }
}

/* ---------------------------------------------------------------------- */

void FixHalt::min_post_force(int /* vflag */)
{
  if (update->ntimestep == thisstep) return;
  if ((update->ntimestep % nevery) == 0) end_of_step();
  thisstep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixHalt::end_of_step()
{
  // check if another partition has exited and we need to exit, too.

  if (uflag) {
    MPI_Status status;
    int partition = -1;
    int flag = 0;
    if (comm->me == 0) {

      // probe if any stop request from another partition is pending

      MPI_Iprobe(MPI_ANY_SOURCE, UTAG, universe->uworld, &flag, &status);

      if (flag) {
        // determine which partition sent the stop request and receive the message
        for (int i = 0; i < universe->nworlds; ++i)
          if (universe->root_proc[i] == status.MPI_SOURCE) partition = i + 1;

        MPI_Recv(&flag, 1, MPI_INT, status.MPI_SOURCE, UTAG, universe->uworld, MPI_STATUS_IGNORE);
      }
    }

    // broadcast stop request partition to all processes in our partition

    MPI_Bcast(&partition, 1, MPI_INT, 0, world);

    // exit request pending handle the same as below

    if (partition > 0) {

      // hard halt -> exit LAMMPS
      // soft/continue halt -> trigger timer to break from run loop
      // print message with ID of fix halt in case multiple instances

      auto message = fmt::format("Received universe halt request from partition {} for fix-id {} on step {}",
                                 partition, id, update->ntimestep);
      if (eflag == HARD) {
        error->all(FLERR, message);
      } else if ((eflag == SOFT) || (eflag == CONTINUE)) {
        if ((comm->me == 0) && (msgflag == YESMSG)) error->message(FLERR, message);
        timer->force_timeout();
      }
    }
  }

  // variable evaluation may invoke computes so wrap with clear/add

  double attvalue;

  if (attribute == TLIMIT) {
    if (update->ntimestep != nextstep) return;
    attvalue = tlimit();
  } else if (attribute == DISKFREE) {
    attvalue = platform::disk_free(dlimit_path) / 1048576.0;    // MBytes
  } else if (attribute == BONDMAX) {
    attvalue = bondmax();
  } else {
    modify->clearstep_compute();
    attvalue = input->variable->compute_equal(ivar);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // ensure that the attribute is *exactly* the same on all ranks

  MPI_Bcast(&attvalue, 1, MPI_DOUBLE, 0, world);

  // check if halt is triggered, else just return

  if (operation == LT) {
    if (attvalue >= value) return;
  } else if (operation == LE) {
    if (attvalue > value) return;
  } else if (operation == GT) {
    if (attvalue <= value) return;
  } else if (operation == GE) {
    if (attvalue < value) return;
  } else if (operation == EQ) {
    if (attvalue != value) return;
  } else if (operation == NEQ) {
    if (attvalue == value) return;
  } else if (operation == XOR) {
    if ((attvalue == 0.0 && value == 0.0) || (attvalue != 0.0 && value != 0.0)) return;
  }

  // send message to all other root processes to trigger exit across universe, if requested

  if (uflag && (comm->me == 0)) {
    MPI_Request *req = new MPI_Request[universe->nworlds];
    for (int i = 0; i < universe->nworlds; ++i) {
      if (universe->me == universe->root_proc[i]) continue;
      MPI_Isend(&eflag, 1, MPI_INT, universe->root_proc[i], UTAG, universe->uworld, req + i);
    }

    // wait for all sends to complete, so MPI_Finalize() will be happy
    for (int i = 0; i < universe->nworlds; ++i) {
      if (universe->me == universe->root_proc[i]) continue;
      MPI_Wait(req + i, MPI_STATUS_IGNORE);
      MPI_Request_free(req + i);
    }
  }

  // hard halt -> exit LAMMPS
  // soft/continue halt -> trigger timer to break from run loop
  // print message with ID of fix halt in case multiple instances

  std::string message = fmt::format("Fix halt condition for fix-id {} met on step {} with value {}",
                                    id, update->ntimestep, attvalue);
  if (eflag == HARD) {
    error->all(FLERR, message);
  } else if ((eflag == SOFT) || (eflag == CONTINUE)) {
    if ((comm->me == 0) && (msgflag == YESMSG)) error->message(FLERR, message);
    timer->force_timeout();
  }
}

/* ----------------------------------------------------------------------
   reset expired timer setting to original value, if requested
------------------------------------------------------------------------- */

void FixHalt::post_run()
{
  // continue halt -> subsequent runs are allowed

  if (eflag == CONTINUE) timer->reset_timeout();
}

/* ----------------------------------------------------------------------
   compute max length of any bond using Neighbor bondlist for each proc
------------------------------------------------------------------------- */

double FixHalt::bondmax()
{
  double **x = atom->x;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int i1, i2;
  double delx, dely, delz, rsq;
  double maxone = 0.0;

  for (int n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    maxone = MAX(rsq, maxone);
  }

  double maxall;
  MPI_Allreduce(&maxone, &maxall, 1, MPI_DOUBLE, MPI_MAX, world);

  return sqrt(maxall);
}

/* ----------------------------------------------------------------------
   compute synced elapsed time
   reset nextstep = estimate of timestep when run will end
   first project to 1/2 the run time, thereafter to end of run
------------------------------------------------------------------------- */

double FixHalt::tlimit()
{
  double cpu = timer->elapsed(Timer::TOTAL);
  MPI_Bcast(&cpu, 1, MPI_DOUBLE, 0, world);

  if (cpu < value) {
    bigint elapsed = update->ntimestep - update->firststep;
    bigint final = update->firststep + static_cast<bigint>(tratio * value / cpu * elapsed);
    nextstep = (final / nevery) * nevery + nevery;
    if (nextstep == update->ntimestep) nextstep += nevery;
    tratio = 1.0;
  }

  return cpu;
}
