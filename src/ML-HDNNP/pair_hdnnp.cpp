/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   This file initially came from n2p2 (https://github.com/CompPhysVienna/n2p2)
   Copyright (2018) Andreas Singraber (University of Vienna)

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andreas Singraber
------------------------------------------------------------------------- */

#include "pair_hdnnp.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

#include "InterfaceLammps.h"    // n2p2 interface header

using namespace LAMMPS_NS;

static const char cite_user_hdnnp_package[] =
    "ML-HDNNP package: 10.1021/acs.jctc.8b00770\n\n"
    "@Article{Singraber19,\n"
    " author = {Singraber, Andreas and Behler, J{\"o}rg and Dellago, Christoph},\n"
    " title = {Library-{{Based LAMMPS Implementation}} of {{High}}-{{Dimensional Neural Network "
    "Potentials}}},\n"
    " year = {2019},\n"
    " month = mar,\n"
    " volume = {15},\n"
    " pages = {1827--1840},\n"
    " doi = {10.1021/acs.jctc.8b00770},\n"
    " journal = {J.~Chem.~Theory~Comput.},\n"
    " number = {3}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairHDNNP::PairHDNNP(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_hdnnp_package);

  single_enable = 0;    // 1 if single() routine exists
  restartinfo = 0;      // 1 if pair style writes restart info
  one_coeff = 1;        // 1 if allows only one coeff * * call
  manybody_flag = 1;    // 1 if a manybody potential
  unit_convert_flag =
      0;    // TODO: Check possible values. value != 0 indicates support for unit conversion.
  reinitflag = 0;    // 1 if compatible with fix adapt and alike

  interface = new nnp::InterfaceLammps();
}

/* ---------------------------------------------------------------------- */

PairHDNNP::~PairHDNNP()
{
  delete interface;
  memory->destroy(setflag);
  memory->destroy(cutsq);
}

/* ---------------------------------------------------------------------- */

void PairHDNNP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // Set number of local atoms and add element.
  interface->setLocalAtoms(atom->nlocal, atom->type);
  // Transfer tags separately. Interface::setLocalTags is overloaded internally
  // to work with both -DLAMMPS_SMALLBIG (tagint = int) and -DLAMMPS_BIGBIG
  // (tagint = int64_t)
  interface->setLocalTags(atom->tag);

  // Transfer local neighbor list to n2p2 interface.
  transferNeighborList();

  // Compute symmetry functions, atomic neural networks and add up energy.
  interface->process();

  // Do all stuff related to extrapolation warnings.
  if (showew == true || showewsum > 0 || maxew >= 0) { handleExtrapolationWarnings(); }

  // Calculate forces of local and ghost atoms.
  interface->getForces(atom->f);

  // Add energy contribution to total energy.
  if (eflag_global)
    ev_tally(0, 0, atom->nlocal, 1, interface->getEnergy(), 0.0, 0.0, 0.0, 0.0, 0.0);

  // Add atomic energy if requested (CAUTION: no physical meaning!).
  if (eflag_atom)
    for (int i = 0; i < atom->nlocal; ++i) eatom[i] = interface->getAtomicEnergy(i);

  // If virial needed calculate via F dot r.
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairHDNNP::settings(int narg, char **arg)
{
  int iarg = 0;

  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  maxCutoffRadius = utils::numeric(FLERR, arg[0], false, lmp);
  iarg++;

  // default settings
  directory = utils::strdup("hdnnp/");
  showew = true;
  showewsum = 0;
  maxew = 0;
  resetew = false;
  cflength = 1.0;
  cfenergy = 1.0;
  numExtrapolationWarningsTotal = 0;
  numExtrapolationWarningsSummary = 0;

  while (iarg < narg) {
    // set HDNNP directory
    if (strcmp(arg[iarg], "dir") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      delete[] directory;
      directory = utils::strdup(arg[iarg + 1]);
      iarg += 2;
      // show extrapolation warnings
    } else if (strcmp(arg[iarg], "showew") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        showew = true;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        showew = false;
      else
        error->all(FLERR, "Illegal pair_style command");
      iarg += 2;
      // show extrapolation warning summary
    } else if (strcmp(arg[iarg], "showewsum") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      showewsum = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      // maximum allowed extrapolation warnings
    } else if (strcmp(arg[iarg], "maxew") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      maxew = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      // reset extrapolation warning counter
    } else if (strcmp(arg[iarg], "resetew") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        resetew = true;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        resetew = false;
      else
        error->all(FLERR, "Illegal pair_style command");
      iarg += 2;
      // length unit conversion factor
    } else if (strcmp(arg[iarg], "cflength") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cflength = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      // energy unit conversion factor
    } else if (strcmp(arg[iarg], "cfenergy") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cfenergy = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHDNNP::coeff(int narg, char **arg)
{
  int n = atom->ntypes;

  if (!allocated) allocate();

  if (narg != 2 + n) error->all(FLERR, "Incorrect args for pair coefficients");

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  int *map = new int[n + 1];
  for (int i = 0; i < n; i++) map[i] = 0;

  emap = "";
  for (int i = 2; i < narg; i++) {
    if (strcmp(arg[i], "NULL") != 0) {
      if (i != 2) emap += ",";
      emap += std::to_string(i - 1) + ":" + arg[i];
      map[i - 1] = 1;
    }
  }

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] > 0 && map[j] > 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

  delete[] map;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHDNNP::init_style()
{
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // Return immediately if HDNNP setup is already completed.
  if (interface->isInitialized()) return;

  // Activate screen and logfile output only for rank 0.
  if (comm->me == 0) {
    if (lmp->screen != nullptr) interface->log.registerCFilePointer(&(lmp->screen));
    if (lmp->logfile != nullptr) interface->log.registerCFilePointer(&(lmp->logfile));
  }

  // Initialize interface on all processors.
  interface->initialize(directory, emap.c_str(), showew, resetew, showewsum, maxew, cflength,
                        cfenergy, maxCutoffRadius, atom->ntypes, comm->me);

  // LAMMPS cutoff radius (given via pair_coeff) should not be smaller than
  // maximum symmetry function cutoff radius.
  if (maxCutoffRadius < interface->getMaxCutoffRadius())
    error->all(FLERR, "Inconsistent cutoff radius");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHDNNP::init_one(int, int)
{
  return maxCutoffRadius;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairHDNNP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");    // TODO: Is this required?
}

void PairHDNNP::transferNeighborList()
{
  // Transfer neighbor list to n2p2.
  double rc2 = maxCutoffRadius * maxCutoffRadius;
  for (int ii = 0; ii < list->inum; ++ii) {
    int i = list->ilist[ii];
    for (int jj = 0; jj < list->numneigh[i]; ++jj) {
      int j = list->firstneigh[i][jj];
      j &= NEIGHMASK;
      double dx = atom->x[i][0] - atom->x[j][0];
      double dy = atom->x[i][1] - atom->x[j][1];
      double dz = atom->x[i][2] - atom->x[j][2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= rc2) { interface->addNeighbor(i, j, atom->tag[j], atom->type[j], dx, dy, dz, d2); }
    }
  }
}

void PairHDNNP::handleExtrapolationWarnings()
{
  // Get number of extrapolation warnings for local atoms.
  long numCurrentEW = (long) interface->getNumExtrapolationWarnings();

  // Update (or set, resetew == true) total warnings counter.
  if (resetew)
    numExtrapolationWarningsTotal = numCurrentEW;
  else
    numExtrapolationWarningsTotal += numCurrentEW;

  // Update warnings summary counter.
  if (showewsum > 0) { numExtrapolationWarningsSummary += numCurrentEW; }

  // If requested write extrapolation warnings.
  // Requires communication of all symmetry functions statistics entries to
  // rank 0.
  if (showew > 0) {
    // First collect an overview of extrapolation warnings per process.
    long *numEWPerProc = nullptr;
    if (comm->me == 0) numEWPerProc = new long[comm->nprocs];
    MPI_Gather(&numCurrentEW, 1, MPI_LONG, numEWPerProc, 1, MPI_LONG, 0, world);

    if (comm->me == 0) {
      for (int i = 1; i < comm->nprocs; i++) {
        if (numEWPerProc[i] > 0) {
          long bs = 0;
          MPI_Status ms;
          // Get buffer size.
          MPI_Recv(&bs, 1, MPI_LONG, i, 0, world, &ms);
          char *buf = new char[bs];
          // Receive buffer.
          MPI_Recv(buf, bs, MPI_BYTE, i, 0, world, &ms);
          interface->extractEWBuffer(buf, bs);
          delete[] buf;
        }
      }
      interface->writeExtrapolationWarnings();
    } else if (numCurrentEW > 0) {
      // Get desired buffer length for all extrapolation warning entries.
      long bs = interface->getEWBufferSize();
      // Allocate and fill buffer.
      char *buf = new char[bs];
      interface->fillEWBuffer(buf, bs);
      // Send buffer size and buffer.
      MPI_Send(&bs, 1, MPI_LONG, 0, 0, world);
      MPI_Send(buf, bs, MPI_BYTE, 0, 0, world);
      delete[] buf;
    }

    if (comm->me == 0) delete[] numEWPerProc;
  }

  // If requested gather number of warnings to display summary.
  if (showewsum > 0 && update->ntimestep % showewsum == 0) {
    long globalEW = 0;
    // Communicate the sum over all processors to proc 0.
    MPI_Reduce(&numExtrapolationWarningsSummary, &globalEW, 1, MPI_LONG, MPI_SUM, 0, world);
    // Write to screen or logfile.
    if (comm->me == 0)
      utils::logmesg(lmp, "### NNP EW SUMMARY ### TS: {:10d} EW {:10d} EWPERSTEP {:10.3e}\n",
                     update->ntimestep, globalEW, double(globalEW) / showewsum);
    // Reset summary counter.
    numExtrapolationWarningsSummary = 0;
  }

  // Stop if maximum number of extrapolation warnings is exceeded.
  if (numExtrapolationWarningsTotal > maxew) {
    error->one(FLERR, "Too many extrapolation warnings");
  }

  // Reset internal extrapolation warnings counters.
  interface->clearExtrapolationWarnings();
}
