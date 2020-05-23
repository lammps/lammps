/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing Author: Jacob Gissinger (jacob.gissinger@colorado.edu)
------------------------------------------------------------------------- */

#include "reset_mol_ids.h"
#include <mpi.h>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA_MSPEC 8 // max for molID first-neigh list

/* ---------------------------------------------------------------------- */

ResetMolIDs::ResetMolIDs(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ResetMolIDs::command(int narg, char ** /* arg */)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_mol_ids command before simulation box is defined");
  if (narg != 0) error->all(FLERR,"Illegal reset_mol_ids command");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_mol_ids unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Can only use reset_mol_ids on molecular systems");

  // NOTE: check if any fixes exist which store molecule IDs?
  // if so, this operation will mess up the fix

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Resetting molecule IDs ...\n");
    if (logfile) fprintf(logfile,"Resetting molecule IDs ...\n");
  }

  MPI_Comm_size(world,&nprocs);
  MPI_Comm_rank(world,&me);

  l2l_nlocal = 0;
  l2l_nglobal = 0;

  // create an atom map if one doesn't exist already

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // initialize system since comm->borders() will be invoked

  lmp->init();

  // setup domain, communication
  // acquire ghosts - is that really necessary?
  // exchange will clear map, borders will reset
  // this is the map needed to lookup current global IDs for bond topology

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int max_onetwo = 0;
  for (int i = 0; i < nlocal; i++)
    if (nspecial[i][0] > max_onetwo) max_onetwo = nspecial[i][0];
  MPI_Allreduce(MPI_IN_PLACE,&max_onetwo,1,MPI_INT,MPI_MAX,world);

  // xspecials are first neighs, also listed for ghosts (unlike specials)
  double **nxspecial; // N first-neighbors (includes ghosts)
  double **xspecial; // first-neighbor list (includes ghosts)
  memory->create(nxspecial,nall,1,"reset_mol_ids:nxspecial"); // only first neighs
  memory->create(xspecial,nall,max_onetwo,"reset_mol_ids:xspecial"); // only need space for first neighs
  for (int i = 0; i < nlocal; i++) {
    nxspecial[i][0] = nspecial[i][0];
    for (int j = 0; j < nspecial[i][0]; j++) {
      xspecial[i][j] = special[i][j];
    }
  }
  comm->forward_comm_array(1,nxspecial);
  comm->forward_comm_array(max_onetwo,xspecial); // only need first neighs

  // remapped == 1 if mol ID has already been reset
  int *remapped;
  memory->create(remapped,nall,"reset_mol_ids:remapped");
  for (int i = 0; i < nall; i++) {
    remapped[i] = 0;
  }

  // pioneers: atoms at edge of bondwalk
  // next_pioneers: holder for next set of pioneers
  tagint *pioneers;
  memory->create(pioneers,nall,"reset_mol_ids:pioneers");
  tagint *next_pioneers;
  memory->create(next_pioneers,nall,"reset_mol_ids:next_pioneers");

  // locally remapped molID
  tagint *localmolID;
  memory->create(localmolID,nall,"reset_mol_ids:localmolID");
  for (int i = 0; i < nall; i++) {
    localmolID[i] = 0;
  }

  // initialize algorithm with first atom
  int molIDincr = 1; // index local molIDs starting with 1
  int npioneers = 1;
  int next_npioneers = 0;
  int nmapped = 1; // number of atoms successfully remapped so far
  pioneers[0] = 0; // local atom indices starting with 0
  localmolID[0] = molIDincr;
  remapped[0] = 1;

  // walk bond list (via specials) to reassign mol IDs
  int localscan = 1; // min index start for unmapped atoms
  while (localscan < nlocal || npioneers > 0) { // only need nlocal

    // if reached end of molecule, npioneers will be zero
    // efficiency could be improved for single atoms...
    if (npioneers == 0) {
      for (int i = localscan; i < nlocal; i++) { // only need nlocal
        if (remapped[i] == 0) {
          npioneers = 1;
          pioneers[0] = i;
          localmolID[i] = ++molIDincr;
          remapped[i] = 1;
          nmapped++;
          localscan = i + 1;
          break;
        }
        if (i == nlocal-1) localscan = nlocal;
      }
    }

    // assumes first neighbors always listed in special (newton/special settings issues?)
    next_npioneers = 0;
    for (int i = 0; i < npioneers; i++) {
      for (int j = 0; j < nxspecial[(int) pioneers[i]][0]; j++) {
        int ispec = (int) atom->map(xspecial[(int) pioneers[i]][j]);
        if (ispec < 0) continue;
        if (remapped[ispec] == 0) {
          localmolID[ispec] = molIDincr;
          next_pioneers[next_npioneers++] = ispec;
          remapped[ispec] = 1;
          nmapped++;
        }
      }
    }

    npioneers = next_npioneers;
    for (int i = 0; i < next_npioneers; i++) {
      pioneers[i] = next_pioneers[i];
    }
  }

  // in serial, that's it
  if (nprocs == 1) {
    for (int i = 0; i < nlocal; i++)
      molecule[i] = localmolID[i]; // assumes forward_comm elsewhere...
  } else {

    // in parallel, there's more...
    // first step: each proc runs through its ghosts to determine equiv. molIDs
    // second step: proc 0 creates master list of equivalent molIDs
    // final step: broadcast master list, each proc updates its local atoms

    int maxID = molIDincr; // max local molID
    int *ghostly;
    memory->create(ghostly,maxID,"reset_mol_ids:ghostly");
    for (int i = 0; i < maxID; i++)
      ghostly[i] = 0;
    for (int i = nlocal; i < nall; i++) {
      // skip if we didn't assign atom a local molID
      if (localmolID[i] == 0) continue;
      ghostly[localmolID[i]-1] = 1;
    }

    int nghostly = 0; // # local molIDs that include ghost atoms
    for (int i = 0; i < maxID; i++)
      if (ghostly[i] == 1) nghostly++;

    // let's renumber all molecules such they are globally unique
    // we start with contiguous #s for all strictly-local molIDs
    int nstrict_local = maxID - nghostly;
    int *alln;
    memory->create(alln,nprocs,"reset_mol_ids:alln");
    for (int i = 0; i < nprocs; i++)
      alln[i] = 0;
    MPI_Allgather(&nstrict_local, 1, MPI_INT, alln, 1, MPI_INT, world);
    int total_nstrict = 0;
    for (int i = 0; i < nprocs; i++)
      total_nstrict = total_nstrict + alln[i];

    int goffset = total_nstrict+1; // first gostly molID

    int lstart = 1; // molID start for strictly-local IDs
    for (int i = 0; i < me; i++)
      lstart += alln[i];

    // then we start numbering ghostly molIDs, which need to be deduped eventually
    for (int i = 0; i < nprocs; i++)
      alln[i] = 0;
    MPI_Allgather(&nghostly, 1, MPI_INT, alln, 1, MPI_INT, world);
    int total_nghostly = 0;
    for (int i = 0; i < nprocs; i++)
      total_nghostly = total_nghostly + alln[i];

    int gstart = 0; // molID offset for ghostly IDs. add total_nstrict+1 later
    for (int i = 0; i < me; i++)
      gstart += alln[i];

    // go through and re-ID molecules (globally unique!) (but duplicates)
    int lincr = lstart; // incr new strictly-local molIDs
    int gincr = gstart + goffset; // incr new ghostly molIDs
    tagint *newID;
    memory->create(newID,maxID,"reset_mol_ids:newID");
    for (int i = 0; i < maxID; i++) {
      if (ghostly[i] == 1) newID[i] = gincr++;
      else newID[i] = lincr++;
    }
    for (int i = 0; i < nall; i++) {
      if (localmolID[i] == 0) continue;
      localmolID[i] = newID[localmolID[i]-1]; // can't use localmolID for indexing after this
    }

    // okay, now we have globally unique IDs
    // let's find equivs and send to root for walking

    double **ghostlymolID;
    memory->create(ghostlymolID,nall,1,"reset_mol_ids:ghostlymolID");
    for (int i = 0; i < nlocal; i++)
      ghostlymolID[i][0] = localmolID[i];
    comm->forward_comm_array(1,ghostlymolID);

    // run through ghosts
    int max_equiv = DELTA_MSPEC;
    int *nequiv_molIDs;
    memory->create(nequiv_molIDs,nghostly,"reset_mol_ids:nequiv_molIDs"); // diff. lengths on each proc
    tagint **equiv_molIDs;
    memory->create(equiv_molIDs,max_equiv,nghostly,"reset_mol_ids:equiv_molIDs"); // diff. lengths on each proc
    for (int i = 0; i < nghostly; i++)
      nequiv_molIDs[i] = 0;

    int addflag;
    int imolID;
    for (int i = nlocal; i < nall; i++) {
      // skip if we didn't assign atom a local molID
      if (localmolID[i] == 0) continue;
      // find other procs' equivalent to my ghostly atoms' molID
      imolID = localmolID[i]-gstart-goffset;
      addflag = 1;
      for (int j = 0; j < nequiv_molIDs[imolID]; j++) {
        if (static_cast<tagint> (ghostlymolID[i][0]) == equiv_molIDs[j][imolID]) { // double vs tagint?
          addflag = 0;
          break;
        }
      }

      if (addflag == 1) {
        if (nequiv_molIDs[imolID] == max_equiv) {
          max_equiv += DELTA_MSPEC;
          memory->grow(equiv_molIDs,max_equiv,nghostly,"reset_mol_ids:equiv_molIDs");
        }
        equiv_molIDs[nequiv_molIDs[imolID]++][imolID] = static_cast<tagint> (ghostlymolID[i][0]);
        l2l_nlocal++;
      }
    }

    // list of equiv molID pairs
    memory->create(local_l2l,2,l2l_nlocal,"reset_mol_ids:local_l2l");
    l2l_nlocal = 0;
    for (int i = 0; i < nghostly; i++) {
      for (int j = 0; j < nequiv_molIDs[i]; j++) {
        local_l2l[0][l2l_nlocal] = i+gstart+goffset; // local molID of locally-ghost
        local_l2l[1][l2l_nlocal++] = equiv_molIDs[j][i]; // local molID for owner
      }
    }

    // okay we need to send info to root (variable-gather)
    gather_molIDs();

    tagint *global_molIDs;
    if (me == 0) {
      // build special list (mspecials)
      int *nmspecial;
      tagint **mspecial;
      memory->create(nmspecial,total_nghostly,"reset_mol_ids:nmspecial");
      for (int i = 0; i < total_nghostly; i++)
        nmspecial[i] = 0;
      int max_mspec = DELTA_MSPEC;
      memory->create(mspecial,max_mspec,total_nghostly,"reset_mol_ids:mspecial");

      // to convert from molID to its index, subtract total_nstrict+1 (goffset)
      tagint mol1,mol2;
      int imol1,imol2;
      for (int i = 0; i < l2l_nglobal; i++) {
        mol1 = global_l2l[0][i];
        imol1 = mol1 - goffset;
        mol2 = global_l2l[1][i];
        imol2 = mol2 - goffset;
        // mol1
        addflag = 1;
        for (int j = 0; j < nmspecial[imol1]; j++) {
          if (mspecial[j][imol1] == mol2) {
            addflag = 0;
            break;
          }
        }
        if (addflag == 1) {
          if (nmspecial[imol1] == max_mspec) {
            max_mspec += DELTA_MSPEC;
            memory->grow(mspecial,max_mspec,total_nghostly,"reset_mol_ids:mspecial");
          }
          mspecial[nmspecial[imol1]++][imol1] = mol2;
        }
        // mol2
        addflag = 1;
        for (int j = 0; j < nmspecial[imol2]; j++) {
          if (mspecial[j][imol2] == mol1) {
            addflag = 0;
            break;
          }
        }
        if (addflag == 1) {
          if (nmspecial[imol2] == max_mspec) {
            max_mspec += DELTA_MSPEC;
            memory->grow(mspecial,max_mspec,total_nghostly,"reset_mol_ids:mspecial");
          }
          mspecial[nmspecial[imol2]++][imol2] = mol1;
        }
      }

      // walk mspecial
      memory->create(global_molIDs,total_nghostly,"reset_mol_ids:global_molIDs");
      for (int i = 0; i < total_nghostly; i++)
        global_molIDs[i] = 0;

      // remapped == 1 if mol ID has already been reset
      memory->destroy(remapped);
      memory->create(remapped,total_nghostly,"reset_mol_ids:remapped");
      for (int i = 0; i < total_nghostly; i++)
        remapped[i] = 0;

      // pioneers: atoms at edge of bondwalk
      // next_pioneers: holder for next set of pioneers
      int max_npioneers = pow(DELTA_MSPEC,2);
      memory->destroy(pioneers);
      memory->create(pioneers,max_npioneers,"reset_mol_ids:pioneers");
      memory->destroy(next_pioneers);
      memory->create(next_pioneers,max_npioneers,"reset_mol_ids:next_pioneers");

      // initialize algorithm with first molID
      molIDincr = goffset; // index global molIDs starting with nstrict_local+1
      npioneers = 1;
      next_npioneers = 0;
      nmapped = 1; // number of molIDs successfully remapped so far
      pioneers[0] = 0; // molID indices starting with 0
      global_molIDs[0] = molIDincr;
      remapped[0] = 1;

      // walk molID list (via mspecials) to reassign mol IDs
      localscan = 1;
      while (nmapped < total_nghostly) {

        // if reached end of molecule, npioneers will be zero
        // efficiency could be improved for single atoms...
        if (npioneers == 0) {
          for (int i = localscan; i < total_nghostly; i++) { // only need nlocal
            if (remapped[i] == 0) {
              npioneers = 1;
              pioneers[0] = i;
              global_molIDs[i] = ++molIDincr;
              remapped[i] = 1;
              nmapped++;
              localscan = i + 1;
              break;
            }
          }
        }


        next_npioneers = 0;
        for (int i = 0; i < npioneers; i++) {
          for (int j = 0; j < nmspecial[(int) pioneers[i]]; j++) {
            int ispec = (int) mspecial[j][(int) pioneers[i]] - goffset;
            if (remapped[ispec] == 0) {
              if (next_npioneers == max_npioneers) {
                max_npioneers += DELTA_MSPEC;
                memory->grow(pioneers,max_npioneers,"reset_mol_ids:pioneers");
                memory->grow(next_pioneers,max_npioneers,"reset_mol_ids:next_pioneers");
              }
              global_molIDs[ispec] = molIDincr;
              next_pioneers[next_npioneers++] = ispec;
              remapped[ispec] = 1;
              nmapped++;
            }
          }
        }

        npioneers = next_npioneers;
        for (int i = 0; i < next_npioneers; i++) {
          pioneers[i] = next_pioneers[i];
        }
      }

      // me cleanup
      memory->destroy(nmspecial);
      memory->destroy(mspecial);
    }

    // okay proc 0 now knows the new global molIDs
    // let's tell the others
    // reuse local_l2l locally for new global equivs to local_l2l
    memory->destroy(local_l2l);
    memory->create(local_l2l,nghostly,1,"reset_mol_ids:local_l2l");
    // recall gstart: tells us which undeduped molIDs belong to each proc
    int *allgstarts;
    memory->create(allgstarts,nprocs,"reset_mol_ids:allgstarts");
    MPI_Allgather(&gstart, 1, MPI_INT, allgstarts, 1, MPI_INT, world);
    MPI_Scatterv(&(global_molIDs[0]), alln, allgstarts, MPI_INT,
                 &(local_l2l[0][0]), nghostly, MPI_INT, 0, world);

    // finally, update localmolID
    for (int i = 0; i < nlocal; i++)
      if (localmolID[i] > total_nstrict)
        localmolID[i] = local_l2l[localmolID[i]-gstart-goffset][0];

    // done
    for (int i = 0; i < nlocal; i++)
      molecule[i] = localmolID[i]; // assumes forward_comm elsewhere...

    // parallel clean up
    memory->destroy(newID);
    memory->destroy(ghostlymolID);
    memory->destroy(nequiv_molIDs);
    memory->destroy(equiv_molIDs);
    memory->destroy(ghostly);
    memory->destroy(local_l2l);
    memory->destroy(global_l2l);
    memory->destroy(alln);
    memory->destroy(allgstarts);
    if (me == 0) memory->destroy(global_molIDs);
  }

  // delete temporary atom map

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  // clean up
  memory->destroy(nxspecial);
  memory->destroy(xspecial);
  memory->destroy(remapped);
  memory->destroy(pioneers);
  memory->destroy(next_pioneers);
  memory->destroy(localmolID);
}

/* ----------------------------------------------------------------------
gather locally-determined molID equivalences
------------------------------------------------------------------------- */

void ResetMolIDs::gather_molIDs()
{
#if !defined(MPI_STUBS)
  l2l_nglobal = 0;

  int *allncols;
  memory->create(allncols,nprocs,"reset_mol_ids:allncols");
  for (int i = 0; i < nprocs; i++)
    allncols[i] = 0;
  MPI_Allgather(&l2l_nlocal, 1, MPI_INT, allncols, 1, MPI_INT, world);
  for (int i = 0; i < nprocs; i++)
    l2l_nglobal = l2l_nglobal + allncols[i];

  if (l2l_nglobal == 0) {
    delete [] allncols;
    return;
  }

  int *allstarts; // indices for molID equivs sorted by processor
  memory->create(allstarts,nprocs,"reset_mol_ids:allstarts");

  int start = 0;
  for (int i = 0; i < me; i++)
    start += allncols[i];
  MPI_Allgather(&start, 1, MPI_INT, allstarts, 1, MPI_INT, world);
  MPI_Datatype columnunsized, column;
  int sizes[2]    = {2, l2l_nglobal};
  int subsizes[2] = {2, 1};
  int starts[2]   = {0,0};
  MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C,
                            MPI_LMP_TAGINT, &columnunsized);
  MPI_Type_create_resized (columnunsized, 0, sizeof(tagint), &column);
  MPI_Type_commit(&column);

  memory->create(global_l2l,2,l2l_nglobal,"reset_mol_ids:global_l2l");

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < l2l_nglobal; j++)
      global_l2l[i][j] = 0;

  if (l2l_nlocal > 0)
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < l2l_nlocal; j++)
        global_l2l[i][j+start] = local_l2l[i][j];

  // let's send to root
  if (me == 0) {
    MPI_Gatherv(MPI_IN_PLACE, l2l_nlocal, column, // Note: some values ignored for MPI_IN_PLACE
                &(global_l2l[0][0]), allncols, allstarts,
                column, 0, world);
  } else {
    MPI_Gatherv(&(global_l2l[0][start]), l2l_nlocal, column,
                &(global_l2l[0][0]), allncols, allstarts,
                column, 0, world);
  }

  memory->destroy(allncols);
  memory->destroy(allstarts);

  MPI_Type_free(&column);
  MPI_Type_free(&columnunsized);
#endif
}
