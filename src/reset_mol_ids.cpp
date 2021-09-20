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
   Contributing author: Jacob Gissinger (jacob.gissinger@colorado.edu)
------------------------------------------------------------------------- */

#include "reset_mol_ids.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "compute_fragment_atom.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ResetMolIDs::ResetMolIDs(LAMMPS *lmp) : Command(lmp) {
  cfa = nullptr;
  cca = nullptr;

  // default settings

  compressflag = 1;
  singleflag = 0;
  offset = -1;

  idfrag.clear();
  idchunk.clear();
}

/* ---------------------------------------------------------------------- */

ResetMolIDs::~ResetMolIDs()
{
  if (!idfrag.empty()) modify->delete_compute(idfrag);
  if (compressflag && !idchunk.empty()) modify->delete_compute(idchunk);
}

/* ----------------------------------------------------------------------
   called as reset_mol_ids command in input script
------------------------------------------------------------------------- */

void ResetMolIDs::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_mol_ids command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_mol_ids unless atoms have IDs");
  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Can only use reset_mol_ids on molecular systems");

  // process args

  if (narg < 1) error->all(FLERR,"Illegal reset_mol_ids command");
  char *groupid = arg[0];

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal reset_mol_ids command");
      if (strcmp(arg[iarg+1],"yes") == 0) compressflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compressflag = 0;
      else error->all(FLERR,"Illegal reset_mol_ids command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"single") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal reset_mol_ids command");
      if (strcmp(arg[iarg+1],"yes") == 0) singleflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) singleflag = 0;
      else error->all(FLERR,"Illegal reset_mol_ids command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal reset_mol_ids command");
      offset = utils::tnumeric(FLERR,arg[iarg+1],true,lmp);
      if (offset < -1) error->all(FLERR,"Illegal reset_mol_ids command");
      iarg += 2;
    } else error->all(FLERR,"Illegal reset_mol_ids command");
  }

  if (comm->me == 0) utils::logmesg(lmp,"Resetting molecule IDs ...\n");

  // record wall time for resetting molecule IDs

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // initialize system since comm->borders() will be invoked

  lmp->init();

  // setup domain, communication
  // exchange will clear map, borders will reset
  // this is the map needed to lookup current global IDs for bond topology

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // create computes

  create_computes((char *) "COMMAND",groupid);

  // reset molecule IDs

  reset();

  // total time

  MPI_Barrier(world);

  if (comm->me == 0) {
    if (nchunk < 0)
      utils::logmesg(lmp,"  number of new molecule IDs = unknown\n");
    else
      utils::logmesg(lmp,"  number of new molecule IDs = {}\n",nchunk);
    utils::logmesg(lmp,"  reset_mol_ids CPU = {:.3f} seconds\n",
                   MPI_Wtime()-time1);
  }
}

/* ----------------------------------------------------------------------
   create computes used by reset_mol_ids
------------------------------------------------------------------------- */

void ResetMolIDs::create_computes(char *fixid, char *groupid)
{
  int igroup = group->find(groupid);
  if (igroup == -1) error->all(FLERR,"Could not find reset_mol_ids group ID");
  groupbit = group->bitmask[igroup];

  // create instances of compute fragment/atom, compute reduce (if needed),
  // and compute chunk/atom.  all use the group-ID for this command.
  // 'fixid' allows for creating independent instances of the computes

  idfrag = fmt::format("{}_reset_mol_ids_FRAGMENT_ATOM",fixid);
  if (singleflag)
    modify->add_compute(fmt::format("{} {} fragment/atom single yes",idfrag,groupid));
  else
    modify->add_compute(fmt::format("{} {} fragment/atom single no",idfrag,groupid));

  idchunk = fmt::format("{}_reset_mol_ids_CHUNK_ATOM",fixid);
  if (compressflag)
    modify->add_compute(fmt::format("{} {} chunk/atom molecule compress yes",
                                    idchunk,groupid));

  int icompute = modify->find_compute(idfrag);
  cfa = (ComputeFragmentAtom *) modify->compute[icompute];


  if (compressflag) {
    icompute = modify->find_compute(idchunk);
    cca = (ComputeChunkAtom *) modify->compute[icompute];
  }
}

/* ----------------------------------------------------------------------
   called from command() and directly from fixes that update molecules
------------------------------------------------------------------------- */

void ResetMolIDs::reset()
{
  // invoke peratom method of compute fragment/atom
  // walks bond connectivity and assigns each atom a fragment ID
  // if singleflag = 0, atoms w/out bonds will be assigned fragID = 0

  cfa->compute_peratom();
  double *fragIDs = cfa->vector_atom;

  // copy fragID to molecule ID for atoms in group

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molecule[i] = static_cast<tagint> (fragIDs[i]);

  // if compressflag = 0, done
  // set nchunk = -1 since cannot easily determine # of new molecule IDs

  nchunk = -1;

  // if compressflag = 1, invoke peratom method of compute chunk/atom
  // will compress new molecule IDs to be contiguous 1 to Nmol
  // NOTE: use of compute chunk/atom limits Nmol to a 32-bit int

  if (compressflag) {
    cca->compute_peratom();
    double *chunkIDs = cca->vector_atom;
    nchunk = cca->nchunk;

    // if singleflag = 0, check if any single (no-bond) atoms exist
    // if so, they have fragID = 0, and compression just set their chunkID = 1
    // singleexist = 0/1 if no/yes single atoms exist with chunkID = 1

    int singleexist = 0;
    if (!singleflag) {
      int mysingle = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          if (fragIDs[i] == 0.0) mysingle = 1;
      MPI_Allreduce(&mysingle,&singleexist,1,MPI_INT,MPI_MAX,world);
      if (singleexist) nchunk--;
    }

    // if offset < 0 (default), reset it
    // if group = all, offset = 0
    // else offset = largest molID of non-group atoms

    if (offset < 0) {
      if (groupbit != 1) {
        tagint mymol = 0;
        for (int i = 0; i < nlocal; i++)
          if (!(mask[i] & groupbit))
            mymol = MAX(mymol,molecule[i]);
        MPI_Allreduce(&mymol,&offset,1,MPI_LMP_TAGINT,MPI_MAX,world);
      } else offset = 0;
    }

    // reset molecule ID for all atoms in group
    // newID = chunkID + offset

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        tagint newid =  static_cast<tagint>(chunkIDs[i]);
        if (singleexist) {
          if (newid == 1) newid = 0;
          else newid += offset - 1;
        } else newid += offset;
        molecule[i] = newid;
      }
    }
  }
}
