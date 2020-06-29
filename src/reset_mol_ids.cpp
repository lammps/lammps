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
#include "fmt/format.h"

using namespace LAMMPS_NS;

#define DELTA 100

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
    // second step: globally walk equivalent molIDs

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
    int max_equiv = DELTA;
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
          max_equiv += DELTA;
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

    // globally walk each molecule, via molIDs
    // all proc keep master list of 'growing' mol
    // if IDs added, send just added IDs out for next iteration

    // arrays used in gather_molIDs
    memory->create(allnpion,nprocs,"reset_mol_ids:allnpion");
    memory->create(allstarts,nprocs,"reset_mol_ids:allstarts");
    memory->create(global_pion,100,"reset_mol_ids:global_pion");
    
    maxmolID = maxnpion = 100;
    memory->create(pionIDs,maxnpion,"reset_mol_ids:pionIDs");
    memory->create(molIDlist,maxmolID,"reset_mol_ids:molIDlist");

    for (int thisproc = 0; thisproc < nprocs; thisproc++) {
      int quitflag = 0;
      while (1) {
        if (me == thisproc) {
          for (int i = 0; i < l2l_nlocal; i++) {
            if (local_l2l[0][i] == 0) {
              if (i < l2l_nlocal-1) continue;
              else {
                quitflag = 1;
                break;
              }
            }
            // else, initialize with first two equivs
            molIDlist[0] = local_l2l[0][i];
            molIDlist[1] = local_l2l[1][i];
            break;
          }
        }
        MPI_Bcast(&quitflag, 1, MPI_INT, thisproc, world);
        if (quitflag) break;
        nmolID = 2;
        MPI_Bcast(&molIDlist[0], nmolID, MPI_INT, thisproc, world);

        nadd = 2;
        while (nadd > 0) {
          npion = 0;
          for (int i = 0; i < l2l_nlocal; i++) {
            if (local_l2l[0][i] == 0) continue;
            for (int j = 0; j < 2; j++) {
              for (int k = 0; k < nadd; k++) {
                if (local_l2l[j][i] == molIDlist[nmolID-k-1]) {
                  if (npion == maxnpion) {
                    maxnpion += DELTA;
                    memory->grow(pionIDs,maxnpion,"reset_mol_ids:pionIDs");
                  }
                  pionIDs[npion++] = local_l2l[!j][i];
                }
              }
            }
          }

          // gather pioneers
          gather_molIDs();
        }

        total_nstrict++;
        for (int i = 0; i < nlocal; i++) {
          if (localmolID[i] < total_nstrict) continue;
          for (int j = 0; j < nmolID; j++) {
            if (localmolID[i] == molIDlist[j]) {
              localmolID[i] = total_nstrict; // wait, can we do this?
              break;
            }
          }
        }

        // zero out already remapped molIDs
        for (int i = 0; i < l2l_nlocal; i++) {
          if (local_l2l[0][i] == 0) continue;
          for (int j = 0; j < nmolID; j++) {
            if (local_l2l[0][i] == molIDlist[j]) {
              local_l2l[0][i] = 0;
              break;
            }
          }
        }
      }
    }

    // we're done
    for (int i = 0; i < nlocal; i++)
      molecule[i] = localmolID[i]; // assumes forward_comm elsewhere...

    // parallel clean up
    memory->destroy(newID);
    memory->destroy(ghostlymolID);
    memory->destroy(nequiv_molIDs);
    memory->destroy(equiv_molIDs);
    memory->destroy(ghostly);
    memory->destroy(local_l2l);
    memory->destroy(alln);
    memory->destroy(allnpion);
    memory->destroy(allstarts);
    memory->destroy(global_pion);
    memory->destroy(pionIDs);
    memory->destroy(molIDlist);
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
gather locally-determined molID equivalences, dedup, and add to molIDlist
------------------------------------------------------------------------- */

void ResetMolIDs::gather_molIDs()
{
  MPI_Allgather(&npion, 1, MPI_INT, allnpion, 1, MPI_INT, world);
  int npion_global = 0;
  for (int i = 0; i < nprocs; i++)
    npion_global += allnpion[i];

  if (npion_global == 0) return;

  allstarts[0] = 0;
  for (int i = 1; i < nprocs; i++)
    allstarts[i] = allstarts[i-1] + allnpion[i-1];

  memory->grow(global_pion,npion_global,"reset_mol_ids:global_pion");
  MPI_Allgatherv(&pionIDs[0], npion, MPI_INT, &global_pion[0],
                 allnpion, allstarts, MPI_INT, world);

  nadd = 0;
  int addflag;
  for (int i = 0; i < npion_global; i++) {
    addflag = 1;
    for (int j = 0; j < nmolID; j++) {
      if (global_pion[i] == molIDlist[j]) {
        addflag = 0;
        break;
      }
    }

    if (addflag == 1) {
      if (nmolID == maxmolID) {
        maxmolID += DELTA;
        memory->grow(molIDlist,maxmolID,"reset_mol_ids:molIDlist");
      }
      nadd++;
      molIDlist[nmolID++] = global_pion[i];
    }
  }
}
