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

#include <mpi.h>
#include <cstdio>
#include "special.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "accelerator_kokkos.h"
#include "hashlittle.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Special::Special(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  onetwo = onethree = onefour = NULL;
}

/* ---------------------------------------------------------------------- */

Special::~Special()
{
  memory->destroy(onetwo);
  memory->destroy(onethree);
  memory->destroy(onefour);
}

/* ----------------------------------------------------------------------
   create 1-2, 1-3, 1-4 lists of topology neighbors
   store in onetwo, onethree, onefour for each atom
   store 3 counters in nspecial[i]
------------------------------------------------------------------------- */

void Special::build()
{
  int i,j,k,m,n,size,proc;
  int max,maxall;
  char *buf;

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **nspecial = atom->nspecial;

  if (me == 0 && screen) {
    const double * const special_lj   = force->special_lj;
    const double * const special_coul = force->special_coul;
    fprintf(screen,"Finding 1-2 1-3 1-4 neighbors ...\n"
                   "  special bond factors lj:   %-10g %-10g %-10g\n"
                   "  special bond factors coul: %-10g %-10g %-10g\n",
                   special_lj[1],special_lj[2],special_lj[3],
                   special_coul[1],special_coul[2],special_coul[3]);
  }

  // initialize nspecial counters to 0

  for (i = 0; i < nlocal; i++) {
    nspecial[i][0] = 0;
    nspecial[i][1] = 0;
    nspecial[i][2] = 0;
  }

  // -----------------------------------------------------
  // compute nspecial[i][0] = # of 1-2 neighbors of atom i
  // create onetwo[i] = list of 1-2 neighbors for atom i
  // -----------------------------------------------------

  // ncount = # of my datums to send (newton or newton off)
  // include nlocal datums with owner of each atom

  int newton_bond = force->newton_bond;

  int ncount = 0;
  for (i = 0; i < nlocal; i++) ncount += num_bond[i];
  if (newton_bond) ncount *= 2;
  ncount += nlocal;

  int *proclist;
  memory->create(proclist,ncount,"special:proclist");
  InRvous *inbuf = (InRvous *) 
    memory->smalloc((bigint) ncount*sizeof(InRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // input datums = pairs of bonded atoms
  // owning proc for each datum = random hash of atomID
  // one datum for each owned atom: datum = owning proc, atomID
  // one datum for each bond partner: datum = atomID, bond partner ID
  //   add inverted datum when netwon_bond on

  m = 0;
  for (i = 0; i < nlocal; i++) {
    proc = hashlittle(&tag[i],sizeof(tagint),0) % nprocs;
    proclist[m] = proc;
    inbuf[m].me = me;
    inbuf[m].atomID = tag[i];
    inbuf[m].partnerID = 0;
    m++;

    for (j = 0; j < num_bond[i]; j++) {
      proclist[m] = proc;
      inbuf[m].me = -1;
      inbuf[m].atomID = tag[i];
      inbuf[m].partnerID = bond_atom[i][j];
      m++;
    }
    if (newton_bond) {
      for (j = 0; j < num_bond[i]; j++) {
	proclist[m] = hashlittle(&bond_atom[i][j],sizeof(tagint),0) % nprocs;
	inbuf[m].me = -1;
	inbuf[m].atomID = bond_atom[i][j];
	inbuf[m].partnerID = tag[i];
	m++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all info to form and return their onetwo lists

  int nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                                 rendezvous_1234,
                                 buf,sizeof(OutRvous),(void *) this);
  OutRvous *outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[0] and onetwo for all owned atoms based on output info
  // output datums = pairs of atoms that are 1-2 neighbors

  MPI_Allreduce(&max_rvous,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onetwo,nlocal,maxall,"special:onetwo");

  for (i = 0; i < nlocal; i++) nspecial[i][0] = 0;

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onetwo[i][nspecial[i][0]++] = outbuf[m].partnerID;
  }

  memory->sfree(outbuf);

  // compute and print max # of 1-2 neighbors

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-2 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-2 neighbors\n",maxall);
  }

  // -----------------------------------------------------
  // done if special_bond weights for 1-3, 1-4 are set to 1.0
  // -----------------------------------------------------

  if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0 &&
      force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    combine();
    fix_alteration();
    timer_output(time1);
    return;
  }

  // -----------------------------------------------------
  // compute nspecial[i][1] = # of 1-3 neighbors of atom i
  // create onethree[i] = list of 1-3 neighbors for atom i
  // -----------------------------------------------------

  // ncount = # of my datums to send
  // include nlocal datums with owner of each atom

  ncount = nlocal;
  for (i = 0; i < nlocal; i++) ncount += nspecial[i][0]*(nspecial[i][0]-1);

  memory->create(proclist,ncount,"special:proclist");
  inbuf = (InRvous *) 
    memory->smalloc((bigint) ncount*sizeof(InRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // input datums = all pairs of onetwo atoms (they are 1-3 neighbors)
  // owning proc for each datum = random hash of atomID
  // one datum for each owned atom: datum = owning proc, atomID
  // one datum for each onetwo pair: datum = atomID1, atomID2

  m = 0;
  for (i = 0; i < nlocal; i++) {
    proclist[m] = hashlittle(&tag[i],sizeof(tagint),0) % nprocs;
    inbuf[m].me = me;
    inbuf[m].atomID = tag[i];
    inbuf[m].partnerID = 0;
    m++;

    for (j = 0; j < nspecial[i][0]; j++) {
      proc = hashlittle(&onetwo[i][j],sizeof(tagint),0) % nprocs;
      for (k = 0; k < nspecial[i][0]; k++) {
        if (j == k) continue;
        proclist[m] = proc;
        inbuf[m].me = -1;
        inbuf[m].atomID = onetwo[i][j];
        inbuf[m].partnerID = onetwo[i][k];
        m++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all info to form and return their onethree lists

  nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                             rendezvous_1234,
                             buf,sizeof(OutRvous),(void *) this);
  outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[1] and onethree for all owned atoms based on output info
  // output datums = pairs of atoms that are 1-3 neighbors

  MPI_Allreduce(&max_rvous,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onethree,nlocal,maxall,"special:onethree");

  for (i = 0; i < nlocal; i++) nspecial[i][1] = 0;

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onethree[i][nspecial[i][1]++] = outbuf[m].partnerID;
  }

  memory->destroy(outbuf);

  // compute and print max # of 1-3 neighbors

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-3 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-3 neighbors\n",maxall);
  }

  // done if special_bond weights for 1-4 are set to 1.0

  if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    if (force->special_angle) angle_trim();
    combine();
    fix_alteration();
    timer_output(time1);
    return;
  }

  // -----------------------------------------------------
  // compute nspecial[i][2] = # of 1-4 neighbors of atom i
  // create onefour[i] = list of 1-4 neighbors for atom i
  // -----------------------------------------------------

  // ncount = # of my datums to send
  // include nlocal datums with owner of each atom

  ncount = nlocal;
  for (i = 0; i < nlocal; i++) ncount += nspecial[i][1]*nspecial[i][0];

  memory->create(proclist,ncount,"special:proclist");
  inbuf = (InRvous *) 
    memory->smalloc((bigint) ncount*sizeof(InRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // input datums = all pairs of onethree and onetwo atoms (they're 1-4 neighbors)
  // owning proc for each datum = random hash of atomID
  // one datum for each owned atom: datum = owning proc, atomID
  // one datum for each onethree/onetwo pair: datum = atomID1, atomID2

  m = 0;
  for (i = 0; i < nlocal; i++) {
    proclist[m] = hashlittle(&tag[i],sizeof(tagint),0) % nprocs;
    inbuf[m].me = me;
    inbuf[m].atomID = tag[i];
    inbuf[m].partnerID = 0;
    m++;

    for (j = 0; j < nspecial[i][1]; j++) {
      proc = hashlittle(&onethree[i][j],sizeof(tagint),0) % nprocs;
      for (k = 0; k < nspecial[i][0]; k++) {
        proclist[m] = proc;
        inbuf[m].me = -1;
        inbuf[m].atomID = onethree[i][j];
        inbuf[m].partnerID = onetwo[i][k];
        m++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all info to form and return their onefour lists

  nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                             rendezvous_1234,
                             buf,sizeof(OutRvous),(void *) this);
  outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[2] and onefour for all owned atoms based on output info
  // output datums = pairs of atoms that are 1-4 neighbors

  MPI_Allreduce(&max_rvous,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onefour,nlocal,maxall,"special:onefour");

  for (i = 0; i < nlocal; i++) nspecial[i][2] = 0;

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onefour[i][nspecial[i][2]++] = outbuf[m].partnerID;
  }

  memory->destroy(outbuf);

  // compute and print max # of 1-4 neighbors

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-4 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-4 neighbors\n",maxall);
  }

  // finish processing the onetwo, onethree, onefour lists

  dedup();
  if (force->special_angle) angle_trim();
  if (force->special_dihedral) dihedral_trim();
  combine();
  fix_alteration();
  timer_output(time1);
}

/* ----------------------------------------------------------------------
   remove duplicates within each of onetwo, onethree, onefour individually
------------------------------------------------------------------------- */

void Special::dedup()
{
  int i,j;
  tagint m;

  // clear map so it can be used as scratch space

  atom->map_clear();

  // use map to cull duplicates
  // exclude original atom explicitly
  // adjust onetwo, onethree, onefour values to reflect removed duplicates
  // must unset map for each atom

  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int unique;

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        onetwo[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][0] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onetwo[i][j],-1);
  }

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        onethree[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][1] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onethree[i][j],-1);
  }

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        onefour[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][2] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onefour[i][j],-1);
  }

  // re-create map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();
}

/* ----------------------------------------------------------------------
   concatenate onetwo, onethree, onefour into master atom->special list
   remove duplicates between 3 lists, leave dup in first list it appears in
   convert nspecial[0], nspecial[1], nspecial[2] into cumulative counters
------------------------------------------------------------------------- */

void Special::combine()
{
  int i,j;
  tagint m;

  int me;
  MPI_Comm_rank(world,&me);

  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // ----------------------------------------------------
  // compute culled maxspecial = max # of special neighs of any atom
  // ----------------------------------------------------

  // clear map so it can be used as scratch space

  atom->map_clear();

  // unique = # of unique nspecial neighbors of one atom
  // cull duplicates using map to check for them
  // exclude original atom explicitly
  // must unset map for each atom

  int unique;
  int maxspecial = 0;

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);

    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }
    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }
    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }

    maxspecial = MAX(maxspecial,unique);

    atom->map_one(tag[i],-1);
    for (j = 0; j < nspecial[i][0]; j++) atom->map_one(onetwo[i][j],-1);
    for (j = 0; j < nspecial[i][1]; j++) atom->map_one(onethree[i][j],-1);
    for (j = 0; j < nspecial[i][2]; j++) atom->map_one(onefour[i][j],-1);
  }

  // if atom->maxspecial has been updated before, make certain
  // we do not reset it to a smaller value. Since atom->maxspecial
  // is initialized to 1, this ensures that it is larger than zero.

  maxspecial = MAX(atom->maxspecial,maxspecial);

  // compute global maxspecial
  // add in extra factor from special_bonds command
  // allocate correct special array with same nmax, new maxspecial
  // previously allocated one must be destroyed
  // must make AtomVec class update its ptr to special

  MPI_Allreduce(&maxspecial,&atom->maxspecial,1,MPI_INT,MPI_MAX,world);
  atom->maxspecial += force->special_extra;

  // add force->special_extra only once

  force->special_extra = 0;

  if (me == 0) {
    if (screen)
      fprintf(screen,"  %d = max # of special neighbors\n",atom->maxspecial);
    if (logfile)
      fprintf(logfile,"  %d = max # of special neighbors\n",atom->maxspecial);
  }

  if (lmp->kokkos) {
    AtomKokkos* atomKK = (AtomKokkos*) atom;
    atomKK->modified(Host,SPECIAL_MASK);
    atomKK->sync(Device,SPECIAL_MASK);
    MemoryKokkos* memoryKK = (MemoryKokkos*) memory;
    memoryKK->grow_kokkos(atomKK->k_special,atom->special,
                        atom->nmax,atom->maxspecial,"atom:special");
    atomKK->modified(Device,SPECIAL_MASK);
    atomKK->sync(Host,SPECIAL_MASK);
  } else {
    memory->destroy(atom->special);
    memory->create(atom->special,atom->nmax,atom->maxspecial,"atom:special");
  }

  atom->avec->grow_reset();
  tagint **special = atom->special;

  // ----------------------------------------------------
  // fill special array with 1-2, 1-3, 1-4 neighs for each atom
  // ----------------------------------------------------

  // again use map to cull duplicates
  // exclude original atom explicitly
  // adjust nspecial[i] values to reflect removed duplicates
  // nspecial[i][1] and nspecial[i][2] now become cumulative counters

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);

    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][0] = unique;

    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][1] = unique;

    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][2] = unique;

    atom->map_one(tag[i],-1);
    for (j = 0; j < nspecial[i][2]; j++) atom->map_one(special[i][j],-1);
  }

  // re-create map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();
}

/* ----------------------------------------------------------------------
   trim list of 1-3 neighbors by checking defined angles
   delete a 1-3 neigh if they are not end atoms of a defined angle
     and if they are not 1,3 or 2,4 atoms of a defined dihedral
------------------------------------------------------------------------- */

void Special::angle_trim()
{
  int i,j,m,n,proc,index;

  int *num_angle = atom->num_angle;
  int *num_dihedral = atom->num_dihedral;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // stats on old 1-3 neighbor counts

  double onethreecount = 0.0;
  for (i = 0; i < nlocal; i++) onethreecount += nspecial[i][1];
  double allcount;
  MPI_Allreduce(&onethreecount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,
              "  %g = # of 1-3 neighbors before angle trim\n",allcount);
    if (logfile)
      fprintf(logfile,
              "  %g = # of 1-3 neighbors before angle trim\n",allcount);
  }

  // if angles or dihedrals are defined,
  // flag each 1-3 neigh if it appears in an angle or dihedral

  if ((num_angle && atom->nangles) || (num_dihedral && atom->ndihedrals)) {

    // ncount = # of my datums to send in 3 parts for each owned atom
    // proc owner, onethree list, angle end points
    // angle end points are from angle list and 1-3 and 2-4 pairs in dihedrals
    // latter is only for angles or dihedrlas where I own atom2

    int ncount = nlocal;
    for (i = 0; i < nlocal; i++) ncount += nspecial[i][1];
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < num_angle[i]; j++) {
        index = atom->map(angle_atom2[i][j]);
        if (index >= 0 && index < nlocal) ncount += 2;
      }
      for (j = 0; j < num_dihedral[i]; j++) {
        index = atom->map(dihedral_atom2[i][j]);
        if (index >= 0 && index < nlocal) ncount += 4;
      }
    }

    int *proclist;
    memory->create(proclist,ncount,"special:proclist");
    InRvous *inbuf = (InRvous *) 
      memory->smalloc((bigint) ncount*sizeof(InRvous),"special:inbuf");

    // setup input buf to rendezvous comm
    // one datum for each owned atom: datum = proc, atomID
    //   sent to owner of atomID
    // one datum for each 1-4 partner: datum = atomID, ID
    //   sent to owner of atomID
    // two datums for each dihedral 1-4 endatoms : datum = atomID, ID
    //   sent to owner of atomID

    m = 0;
    for (i = 0; i < nlocal; i++) {
      proc = hashlittle(&tag[i],sizeof(tagint),0) % nprocs;
      proclist[m] = proc;
      inbuf[m].me = me;
      inbuf[m].atomID = tag[i];
      inbuf[m].partnerID = 0;
      m++;

      for (j = 0; j < nspecial[i][1]; j++) {
        proclist[m] = proc;
        inbuf[m].me = -1;
        inbuf[m].atomID = tag[i];
        inbuf[m].partnerID = onethree[i][j];
        m++;
      }

      for (j = 0; j < num_angle[i]; j++) {
        index = atom->map(angle_atom2[i][j]);
        if (index < 0 || index >= nlocal) continue;

        proclist[m] = hashlittle(&angle_atom1[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = angle_atom1[i][j];
        inbuf[m].partnerID = angle_atom3[i][j];
        m++;

        proclist[m] = hashlittle(&angle_atom3[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = angle_atom3[i][j];
        inbuf[m].partnerID = angle_atom1[i][j];
        m++;
      }

      for (j = 0; j < num_dihedral[i]; j++) {
        index = atom->map(dihedral_atom2[i][j]);
        if (index < 0 || index >= nlocal) continue;

        proclist[m] = hashlittle(&dihedral_atom1[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom1[i][j];
        inbuf[m].partnerID = dihedral_atom3[i][j];
        m++;

        proclist[m] = hashlittle(&dihedral_atom2[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom2[i][j];
        inbuf[m].partnerID = dihedral_atom4[i][j];
        m++;

        proclist[m] = hashlittle(&dihedral_atom3[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom3[i][j];
        inbuf[m].partnerID = dihedral_atom1[i][j];
        m++;

        proclist[m] = hashlittle(&dihedral_atom4[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom4[i][j];
        inbuf[m].partnerID = dihedral_atom2[i][j];
        m++;
      }
    }

    // perform rendezvous operation
    // each proc owns random subset of atoms
    // func = compute bbox of each body, flag atom closest to geometric center
    // when done: each atom has atom ID of owning atom of its body
    
    char *buf;
    int nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                                   rendezvous_trim,
                                   buf,sizeof(OutRvous),(void *) this);
    OutRvous *outbuf = (OutRvous *) buf;

    memory->destroy(proclist);
    memory->sfree(inbuf);

    // reset nspecial[1] and onethree for all owned atoms based on output info

    for (i = 0; i < nlocal; i++) nspecial[i][1] = 0;

    for (m = 0; m < nreturn; m++) {
      i = atom->map(outbuf[m].atomID);
      onethree[i][nspecial[i][1]++] = outbuf[m].partnerID;
    }

    memory->destroy(outbuf);

  // if no angles or dihedrals are defined, delete all 1-3 neighs

  } else {
    for (i = 0; i < nlocal; i++) nspecial[i][1] = 0;
  }

  // stats on new 1-3 neighbor counts

  onethreecount = 0.0;
  for (i = 0; i < nlocal; i++) onethreecount += nspecial[i][1];
  MPI_Allreduce(&onethreecount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,
              "  %g = # of 1-3 neighbors after angle trim\n",allcount);
    if (logfile)
      fprintf(logfile,
              "  %g = # of 1-3 neighbors after angle trim\n",allcount);
  }
}

/* ----------------------------------------------------------------------
   trim list of 1-4 neighbors by checking defined dihedrals
   delete a 1-4 neigh if they are not end atoms of a defined dihedral
------------------------------------------------------------------------- */

void Special::dihedral_trim()
{
  int i,j,m,n,proc,index;

  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // stats on old 1-4 neighbor counts

  double onefourcount = 0.0;
  for (i = 0; i < nlocal; i++) onefourcount += nspecial[i][2];
  double allcount;
  MPI_Allreduce(&onefourcount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,
              "  %g = # of 1-4 neighbors before dihedral trim\n",allcount);
    if (logfile)
      fprintf(logfile,
              "  %g = # of 1-4 neighbors before dihedral trim\n",allcount);
  }

  // if dihedrals are defined, rendezvous onefour list with dihedral 1-4 pairs

  if (num_dihedral && atom->ndihedrals) {

    // ncount = # of my datums to send in 3 parts for each owned atom
    // onefour list, proc owner, dihedral end points
    // latter is only for dihedrals where I own atom2

    int ncount = nlocal;
    for (i = 0; i < nlocal; i++) ncount += nspecial[i][2];
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < num_dihedral[i]; j++) {
        index = atom->map(dihedral_atom2[i][j]);
        if (index >= 0 && index < nlocal) ncount += 2;
      }
    }

    int *proclist;
    memory->create(proclist,ncount,"special:proclist");
    InRvous *inbuf = (InRvous *) 
      memory->smalloc((bigint) ncount*sizeof(InRvous),"special:inbuf");

    // setup input buf to rendezvous comm
    // one datum for each owned atom: datum = proc, atomID
    //   sent to owner of atomID
    // one datum for each 1-4 partner: datum = atomID, ID
    //   sent to owner of atomID
    // two datums for each dihedral 1-4 endatoms : datum = atomID, ID
    //   sent to owner of atomID

    m = 0;
    for (i = 0; i < nlocal; i++) {
      proc = hashlittle(&tag[i],sizeof(tagint),0) % nprocs;
      proclist[m] = proc;
      inbuf[m].me = me;
      inbuf[m].atomID = tag[i];
      inbuf[m].partnerID = 0;
      m++;

      for (j = 0; j < nspecial[i][2]; j++) {
        proclist[m] = proc;
        inbuf[m].me = -1;
        inbuf[m].atomID = tag[i];
        inbuf[m].partnerID = onefour[i][j];
        m++;
      }

      for (j = 0; j < num_dihedral[i]; j++) {
        index = atom->map(dihedral_atom2[i][j]);
        if (index < 0 || index >= nlocal) continue;

        proclist[m] = hashlittle(&dihedral_atom1[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom1[i][j];
        inbuf[m].partnerID = dihedral_atom4[i][j];
        m++;

        proclist[m] = hashlittle(&dihedral_atom4[i][j],sizeof(tagint),0) % nprocs;
        inbuf[m].me = -2;
        inbuf[m].atomID = dihedral_atom4[i][j];
        inbuf[m].partnerID = dihedral_atom1[i][j];
        m++;
      }
    }

    // perform rendezvous operation
    // each proc owns random subset of atoms
    // func = compute bbox of each body, flag atom closest to geometric center
    // when done: each atom has atom ID of owning atom of its body
    
    char *buf;
    int nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                                   rendezvous_trim,
                                   buf,sizeof(OutRvous),(void *) this);
    OutRvous *outbuf = (OutRvous *) buf;

    memory->destroy(proclist);
    memory->sfree(inbuf);

    // reset nspecial[2] and onefour for all owned atoms based on output info

    for (i = 0; i < nlocal; i++) nspecial[i][2] = 0;

    for (m = 0; m < nreturn; m++) {
      i = atom->map(outbuf[m].atomID);
      onefour[i][nspecial[i][2]++] = outbuf[m].partnerID;
    }

    memory->destroy(outbuf);

  // if no dihedrals are defined, delete all 1-4 neighs

  } else {
    for (i = 0; i < nlocal; i++) nspecial[i][2] = 0;
  }

  // stats on new 1-4 neighbor counts

  onefourcount = 0.0;
  for (i = 0; i < nlocal; i++) onefourcount += nspecial[i][2];
  MPI_Allreduce(&onefourcount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,
              "  %g = # of 1-4 neighbors after dihedral trim\n",allcount);
    if (logfile)
      fprintf(logfile,
              "  %g = # of 1-4 neighbors after dihedral trim\n",allcount);
  }
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N InRvous datums
   create outbuf = list of Nout OutRvous datums
------------------------------------------------------------------------- */

int Special::rendezvous_1234(int n, char *inbuf,
                             int *&proclist, char *&outbuf,
                             void *ptr)
{
  int i,j,m;

  Special *sptr = (Special *) ptr;
  Atom *atom = sptr->atom;
  Memory *memory = sptr->memory;

  // clear atom map so it can be here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // initialize hash
  // ncount = number of atoms assigned to me
  // key = atom ID
  // value = index into Ncount-length data structure

  InRvous *in = (InRvous *) inbuf;
  //std::map<tagint,int> hash;
  tagint id;
  
  int ncount = 0;
  for (i = 0; i < n; i++)
    if (in[i].me >= 0)
      //hash[in[i].atomID] = ncount++;
      atom->map_one(in[i].atomID,ncount++);

  // procowner = caller proc that owns each atom
  // atomID = ID of each rendezvous atom I own

  int *procowner,*npartner;
  tagint *atomID;
  memory->create(procowner,ncount,"special:procowner");
  memory->create(atomID,ncount,"special:atomID");
  memory->create(npartner,ncount,"special:npartner");
  for (m = 0; m < ncount; m++) npartner[m] = 0;

  for (i = 0; i < n; i++) { 
    //m = hash.find(in[i].atomID)->second;
    m = atom->map(in[i].atomID);
    if (in[i].me >= 0) {
      procowner[m] = in[i].me;
      atomID[m] = in[i].atomID;
    } else npartner[m]++;
  }

  int max = 0;
  for (m = 0; m < ncount; m++)
    max = MAX(max,npartner[m]);
  sptr->max_rvous = max;

  int **partner;
  memory->create(partner,ncount,max,"special:partner");
  for (m = 0; m < ncount; m++) npartner[m] = 0;

  for (i = 0; i < n; i++) {
    if (in[i].me >= 0) continue;
    //m = hash.find(in[i].atomID)->second;
    m = atom->map(in[i].atomID);
    partner[m][npartner[m]++] = in[i].partnerID;
  }

  // pass list of OutRvous datums back to comm->rendezvous

  int nout = 0;
  for (m = 0; m < ncount; m++) nout += npartner[m];

  memory->create(proclist,nout,"special:proclist");
  OutRvous *out = (OutRvous *)
    memory->smalloc((bigint) nout*sizeof(OutRvous),"special:out");

  nout = 0;
  for (m = 0; m < ncount; m++)
    for (j = 0; j < npartner[m]; j++) {
      proclist[nout] = procowner[m];
      out[nout].atomID = atomID[m];
      out[nout].partnerID = partner[m][j];
      nout++;
    }

  outbuf = (char *) out;

  // clean up
  // Comm::rendezvous will delete proclist and out (outbuf)

  memory->destroy(procowner);
  memory->destroy(atomID);
  memory->destroy(npartner);
  memory->destroy(partner);

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  return nout;
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N InRvous datums
   create outbuf = list of Nout OutRvous datums
------------------------------------------------------------------------- */

int Special::rendezvous_trim(int n, char *inbuf,
                             int *&proclist, char *&outbuf,
                             void *ptr)
{
  int i,j,m;

  Special *sptr = (Special *) ptr;
  Atom *atom = sptr->atom;
  Memory *memory = sptr->memory;

  // clear atom map so it can be here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // initialize hash
  // ncount = number of atoms assigned to me
  // key = atom ID
  // value = index into Ncount-length data structure

  InRvous *in = (InRvous *) inbuf;
  //std::map<tagint,int> hash;
  tagint id;
  
  int ncount = 0;
  for (i = 0; i < n; i++)
    if (in[i].me >= 0)
      //hash[in[i].atomID] = ncount++;
      atom->map_one(in[i].atomID,ncount++);

  // procowner = caller proc that owns each atom
  // atomID = ID of each rendezvous atom I own
  // npartner = # of 1-3 partners for each atom I own

  int *procowner,*npartner;
  tagint *atomID;
  memory->create(procowner,ncount,"special:procowner");
  memory->create(atomID,ncount,"special:atomID");
  memory->create(npartner,ncount,"special:npartner");
  for (m = 0; m < ncount; m++) npartner[m] = 0;

  for (i = 0; i < n; i++) { 
    //m = hash.find(in[i].atomID)->second;
    m = atom->map(in[i].atomID);
    if (in[i].me >= 0) {
      procowner[m] = in[i].me;
      atomID[m] = in[i].atomID;
    } else if (in[i].me == -1) npartner[m]++;
  }

  int max = 0;
  for (m = 0; m < ncount; m++) max = MAX(max,npartner[m]);

  // partner = list of 1-3 or 1-4 partners for each atom I own

  int **partner;
  memory->create(partner,ncount,max,"special:partner");
  for (m = 0; m < ncount; m++) npartner[m] = 0;

  for (i = 0; i < n; i++) {
    if (in[i].me >= 0 || in[i].me == -2) continue;
    //m = hash.find(in[i].atomID)->second;
    m = atom->map(in[i].atomID);
    partner[m][npartner[m]++] = in[i].partnerID;
  }

  // flag = 1 if partner is in an actual angle or in a dihedral

  int **flag;
  memory->create(flag,ncount,max,"special:flag");
  
  for (i = 0; i < ncount; i++)
    for (j = 0; j < npartner[i]; j++)
      flag[i][j] = 0;

  tagint actual;
  for (i = 0; i < n; i++) {
    if (in[i].me != -2) continue;
    actual = in[i].partnerID;
    //m = hash.find(in[i].atomID)->second;
    m = atom->map(in[i].atomID);
    for (j = 0; j < npartner[m]; j++)
      if (partner[m][j] == actual) {
        flag[m][j] = 1;
        break;
      }
  }

  // pass list of OutRvous datums back to comm->rendezvous

  int nout = 0;
  for (m = 0; m < ncount; m++) nout += npartner[m];

  memory->create(proclist,nout,"special:proclist");
  OutRvous *out = (OutRvous *)
    memory->smalloc((bigint) nout*sizeof(OutRvous),"special:out");

  nout = 0;
  for (m = 0; m < ncount; m++)
    for (j = 0; j < npartner[m]; j++) {
      if (flag[m][j] == 0) continue;
      proclist[nout] = procowner[m];
      out[nout].atomID = atomID[m];
      out[nout].partnerID = partner[m][j];
      nout++;
    }

  outbuf = (char *) out;

  // clean up
  // Comm::rendezvous will delete proclist and out (outbuf)

  memory->destroy(procowner);
  memory->destroy(atomID);
  memory->destroy(npartner);
  memory->destroy(partner);
  memory->destroy(flag);

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  return nout;
}

/* ----------------------------------------------------------------------
   allow fixes to alter special list
   currently, only fix drude does this
     so that both the Drude core and electron are same level of neighbor
------------------------------------------------------------------------- */

void Special::fix_alteration()
{
  for (int ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix]->special_alter_flag)
      modify->fix[ifix]->rebuild_special();
}

/* ----------------------------------------------------------------------
   print timing output
------------------------------------------------------------------------- */

void Special::timer_output(double time1)
{
  double time2 = MPI_Wtime();
  if (comm->me == 0) {
    if (screen) fprintf(screen,"  special bonds CPU = %g secs\n",time2-time1);
    if (logfile) fprintf(logfile,"  special bonds CPU = %g secs\n",time2-time1);
  }
}
