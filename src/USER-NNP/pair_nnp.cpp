// n2p2 - A neural network potential package
// Copyright (C) 2018 Andreas Singraber (University of Vienna)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <mpi.h>
#include <string.h>
#include "pair_nnp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairNNP::PairNNP(LAMMPS *lmp) : Pair(lmp)
{
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairNNP::~PairNNP()
{
}

/* ---------------------------------------------------------------------- */

void PairNNP::compute(int eflag, int vflag)
{
  if(eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // Set number of local atoms and add index and element.
  interface.setLocalAtoms(atom->nlocal,atom->tag,atom->type);

  // Transfer local neighbor list to NNP interface.
  transferNeighborList();

  // Compute symmetry functions, atomic neural networks and add up energy.
  interface.process();

  // Do all stuff related to extrapolation warnings.
  if(showew == true || showewsum > 0 || maxew >= 0) {
    handleExtrapolationWarnings();
  }

  // Calculate forces of local and ghost atoms.
  interface.getForces(atom->f);

  // Add energy contribution to total energy.
  if (eflag_global)
     ev_tally(0,0,atom->nlocal,1,interface.getEnergy(),0.0,0.0,0.0,0.0,0.0);

  // Add atomic energy if requested (CAUTION: no physical meaning!).
  if (eflag_atom)
    for (int i = 0; i < atom->nlocal; ++i)
      eatom[i] = interface.getAtomicEnergy(i);

  // If virial needed calculate via F dot r.
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNNP::settings(int narg, char **arg)
{
  int iarg = 0;

  if (narg == 0) error->all(FLERR,"Illegal pair_style command");

  // default settings
  int len = strlen("nnp/") + 1;
  directory = new char[len];
  strcpy(directory,"nnp/");
  showew = true;
  showewsum = 0;
  maxew = 0;
  resetew = false;
  cflength = 1.0;
  cfenergy = 1.0;
  len = strlen("") + 1;
  emap = new char[len];
  strcpy(emap,"");
  numExtrapolationWarningsTotal = 0;
  numExtrapolationWarningsSummary = 0;

  while(iarg < narg) {
    // set NNP directory
    if (strcmp(arg[iarg],"dir") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      delete[] directory;
      len = strlen(arg[iarg+1]) + 2;
      directory = new char[len];
      sprintf(directory, "%s/", arg[iarg+1]);
      iarg += 2;
    // element mapping
    } else if (strcmp(arg[iarg],"emap") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      delete[] emap;
      len = strlen(arg[iarg+1]) + 1;
      emap = new char[len];
      sprintf(emap, "%s", arg[iarg+1]);
      iarg += 2;
    // show extrapolation warnings
    } else if (strcmp(arg[iarg],"showew") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      if (strcmp(arg[iarg+1],"yes") == 0)
        showew = true;
      else if (strcmp(arg[iarg+1],"no") == 0)
        showew = false;
      else
        error->all(FLERR,"Illegal pair_style command");
      iarg += 2;
    // show extrapolation warning summary
    } else if (strcmp(arg[iarg],"showewsum") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      showewsum = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    // maximum allowed extrapolation warnings
    } else if (strcmp(arg[iarg],"maxew") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      maxew = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    // reset extrapolation warning counter
    } else if (strcmp(arg[iarg],"resetew") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      if (strcmp(arg[iarg+1],"yes") == 0)
        resetew = true;
      else if (strcmp(arg[iarg+1],"no") == 0)
        resetew = false;
      else
        error->all(FLERR,"Illegal pair_style command");
      iarg += 2;
    // length unit conversion factor
    } else if (strcmp(arg[iarg],"cflength") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      cflength = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    // energy unit conversion factor
    } else if (strcmp(arg[iarg],"cfenergy") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command");
      cfenergy = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNNP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  maxCutoffRadius = utils::numeric(FLERR,arg[2],false,lmp);

  // TODO: Check how this flag is set.
  int count = 0;
  for(int i=ilo; i<=ihi; i++) {
    for(int j=MAX(jlo,i); j<=jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNNP::init_style()
{
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // Return immediately if NNP setup is already completed.
  if (interface.isInitialized()) return;

  // Activate screen and logfile output only for rank 0.
  if (comm->me == 0) {
    if (lmp->screen != NULL)
      interface.log.registerCFilePointer(&(lmp->screen));    
    if (lmp->logfile != NULL)
      interface.log.registerCFilePointer(&(lmp->logfile));    
  }

  // Initialize interface on all processors.
  interface.initialize(directory,
                       emap,
                       showew,
                       resetew,
                       showewsum,
                       maxew,
                       cflength,
                       cfenergy,
                       maxCutoffRadius,
                       atom->ntypes,
                       comm->me);

  // LAMMPS cutoff radius (given via pair_coeff) should not be smaller than
  // maximum symmetry function cutoff radius.
  if (maxCutoffRadius < interface.getMaxCutoffRadius())
    error->all(FLERR,"Inconsistent cutoff radius");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNNP::init_one(int i, int j)
{
  // TODO: Check how this actually works for different cutoffs.
  return maxCutoffRadius;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNNP::write_restart(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNNP::read_restart(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNNP::write_restart_settings(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNNP::read_restart_settings(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNNP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

void PairNNP::transferNeighborList()
{
  // Transfer neighbor list to NNP.
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
      if (d2 <= rc2) {
        interface.addNeighbor(i,j,atom->tag[j],atom->type[j],dx,dy,dz,d2);
      }
    }
  }
}

void PairNNP::handleExtrapolationWarnings()
{
  // Get number of extrapolation warnings for local atoms.
  // TODO: Is the conversion from std::size_t to long ok?
  long numCurrentEW = (long)interface.getNumExtrapolationWarnings();

  // Update (or set, resetew == true) total warnings counter.
  if (resetew) numExtrapolationWarningsTotal = numCurrentEW;
  else numExtrapolationWarningsTotal += numCurrentEW;

  // Update warnings summary counter.
  if(showewsum > 0) {
    numExtrapolationWarningsSummary += numCurrentEW;
  }

  // If requested write extrapolation warnings.
  // Requires communication of all symmetry functions statistics entries to
  // rank 0.
  if(showew > 0) {
    // First collect an overview of extrapolation warnings per process.
    long* numEWPerProc = NULL;
    if(comm->me == 0) numEWPerProc = new long[comm->nprocs];
    MPI_Gather(&numCurrentEW, 1, MPI_LONG, numEWPerProc, 1, MPI_LONG, 0, world);

    if(comm->me == 0) {
      for(int i=1;i<comm->nprocs;i++) {
        if(numEWPerProc[i] > 0) {
          long bs = 0;
          MPI_Status ms;
          // Get buffer size.
          MPI_Recv(&bs, 1, MPI_LONG, i, 0, world, &ms);
          char* buf = new char[bs];
          // Receive buffer.
          MPI_Recv(buf, bs, MPI_BYTE, i, 0, world, &ms);
          interface.extractEWBuffer(buf, bs);
          delete[] buf;
        }
      }
      interface.writeExtrapolationWarnings();
    }
    else if(numCurrentEW > 0) {
      // Get desired buffer length for all extrapolation warning entries.
      long bs = interface.getEWBufferSize();
      // Allocate and fill buffer.
      char* buf = new char[bs];
      interface.fillEWBuffer(buf, bs);
      // Send buffer size and buffer.
      MPI_Send(&bs, 1, MPI_LONG, 0, 0, world);
      MPI_Send(buf, bs, MPI_BYTE, 0, 0, world);
      delete[] buf;
    }

    if(comm->me == 0) delete[] numEWPerProc;
  }

  // If requested gather number of warnings to display summary.
  if(showewsum > 0 && update->ntimestep % showewsum == 0) {
    long globalEW = 0;
    // Communicate the sum over all processors to proc 0.
    MPI_Reduce(&numExtrapolationWarningsSummary,
               &globalEW, 1, MPI_LONG, MPI_SUM, 0, world);
    // Write to screen or logfile.
    if(comm->me == 0) {
      if(screen) {
        fprintf(screen,
                "### NNP EW SUMMARY ### TS: %10ld EW %10ld EWPERSTEP %10.3E\n",
                update->ntimestep,
                globalEW,
                double(globalEW) / showewsum);
      }
      if(logfile) {
        fprintf(logfile,
                "### NNP EW SUMMARY ### TS: %10ld EW %10ld EWPERSTEP %10.3E\n",
                update->ntimestep,
                globalEW,
                double(globalEW) / showewsum);
      }
    }
    // Reset summary counter.
    numExtrapolationWarningsSummary = 0;
  }

  // Stop if maximum number of extrapolation warnings is exceeded.
  if (numExtrapolationWarningsTotal > maxew) {
    error->one(FLERR,"Too many extrapolation warnings");
  }

  // Reset internal extrapolation warnings counters.
  interface.clearExtrapolationWarnings();
}
