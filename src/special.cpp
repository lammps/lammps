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

#include "mpi.h"
#include "stdio.h"
#include "special.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// allocate space for static class variable

Special *Special::sptr;

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
  int i,j,k,size;
  int max,maxall,nbuf;
  int *buf;

  MPI_Barrier(world);

  int nlocal = atom->nlocal;

  int *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **nspecial = atom->nspecial;

  if (me == 0 && screen) fprintf(screen,"Finding 1-2 1-3 1-4 neighbors ...\n");

  // initialize nspecial counters to 0

  for (i = 0; i < nlocal; i++) {
    nspecial[i][0] = 0;
    nspecial[i][1] = 0;
    nspecial[i][2] = 0;
  }

  // -----------------------------------------------------
  // compute nspecial[i][0] = # of 1-2 neighbors of atom i
  // -----------------------------------------------------

  // bond partners stored by atom itself

  for (i = 0; i < nlocal; i++) nspecial[i][0] = num_bond[i];

  // if newton_bond off, then done
  // else only counted 1/2 of all bonds, so count other half

  if (force->newton_bond) {

    // nbufmax = largest buffer needed to hold info from any proc
    // info for each atom = global tag of 2nd atom in each bond

    nbuf = 0;
    for (i = 0; i < nlocal; i++) nbuf += num_bond[i];
    memory->create(buf,nbuf,"special:buf");

    // fill buffer with global tags of bond partners of my atoms

    size = 0;
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++)
        buf[size++] = bond_atom[i][j];

    // cycle buffer around ring of procs back to self
    // when receive buffer, scan tags for atoms I own
    // when find one, increment nspecial count for that atom

    sptr = this;
    comm->ring(size,sizeof(int),buf,1,ring_one,NULL);

    memory->destroy(buf);
  }

  // ----------------------------------------------------
  // create onetwo[i] = list of 1-2 neighbors for atom i
  // ----------------------------------------------------

  max = 0;
  for (i = 0; i < nlocal; i++) max = MAX(max,nspecial[i][0]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-2 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-2 neighbors\n",maxall);
  }

  memory->create(onetwo,nlocal,maxall,"special:onetwo");

  // count = accumulating counter

  memory->create(count,nlocal,"special:count");
  for (i = 0; i < nlocal; i++) count[i] = 0;

  // add bond partners stored by atom to onetwo list

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++)
      onetwo[i][count[i]++] = bond_atom[i][j];

  // if newton_bond off, then done
  // else only stored 1/2 of all bonds, so store other half

  if (force->newton_bond) {

    // nbufmax = largest buffer needed to hold info from any proc
    // info for each atom = 2 global tags in each bond

    nbuf = 0;
    for (i = 0; i < nlocal; i++) nbuf += 2*num_bond[i];
    memory->create(buf,nbuf,"special:buf");

    // fill buffer with global tags of both atoms in bond

    size = 0;
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++) {
        buf[size++] = tag[i];
        buf[size++] = bond_atom[i][j];
      }

    // cycle buffer around ring of procs back to self
    // when receive buffer, scan 2nd-atom tags for atoms I own
    // when find one, add 1st-atom tag to onetwo list for 2nd atom

    sptr = this;
    comm->ring(size,sizeof(int),buf,2,ring_two,NULL);

    memory->destroy(buf);
  }

  memory->destroy(count);

  // -----------------------------------------------------
  // done if special_bonds for 1-3, 1-4 are set to 1.0
  // -----------------------------------------------------

  if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0 &&
      force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    combine();
    return;
  }

  // -----------------------------------------------------
  // compute nspecial[i][1] = # of 1-3 neighbors of atom i
  // -----------------------------------------------------

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 2 scalars + list of 1-2 neighbors

  nbuf = 0;
  for (i = 0; i < nlocal; i++) nbuf += 2 + nspecial[i][0];
  memory->create(buf,nbuf,"special:buf");

  // fill buffer with:
  // (1) = counter for 1-3 neighbors, initialized to 0
  // (2) = # of 1-2 neighbors
  // (3:N) = list of 1-2 neighbors

  size = 0;
  for (i = 0; i < nlocal; i++) {
    buf[size++] = 0;
    buf[size++] = nspecial[i][0];
    for (j = 0; j < nspecial[i][0]; j++) buf[size++] = onetwo[i][j];
  }

  // cycle buffer around ring of procs back to self
  // when receive buffer, scan list of 1-2 neighbors for atoms I own
  // when find one, increment 1-3 count by # of 1-2 neighbors of my atom,
  //   subtracting one since my list will contain original atom

  sptr = this;
  comm->ring(size,sizeof(int),buf,3,ring_three,buf);

  // extract count from buffer that has cycled back to me
  // nspecial[i][1] = # of 1-3 neighbors of atom i

  j = 0;
  for (i = 0; i < nlocal; i++) {
    nspecial[i][1] = buf[j];
    j += 2 + nspecial[i][0];
  }

  memory->destroy(buf);

  // ----------------------------------------------------
  // create onethree[i] = list of 1-3 neighbors for atom i
  // ----------------------------------------------------

  max = 0;
  for (i = 0; i < nlocal; i++) max = MAX(max,nspecial[i][1]);
  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-3 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-3 neighbors\n",maxall);
  }

  memory->create(onethree,nlocal,maxall,"special:onethree");

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 4 scalars + list of 1-2 neighs + list of 1-3 neighs

  nbuf = 0;
  for (i = 0; i < nlocal; i++) nbuf += 4 + nspecial[i][0] + nspecial[i][1];
  memory->create(buf,nbuf,"special:buf");

  // fill buffer with:
  // (1) = global tag of original atom
  // (2) = # of 1-2 neighbors
  // (3) = # of 1-3 neighbors
  // (4) = counter for 1-3 neighbors, initialized to 0
  // (5:N) = list of 1-2 neighbors
  // (N+1:2N) space for list of 1-3 neighbors

  size = 0;
  for (i = 0; i < nlocal; i++) {
    buf[size++] = tag[i];
    buf[size++] = nspecial[i][0];
    buf[size++] = nspecial[i][1];
    buf[size++] = 0;
    for (j = 0; j < nspecial[i][0]; j++) buf[size++] = onetwo[i][j];
    size += nspecial[i][1];
  }

  // cycle buffer around ring of procs back to self
  // when receive buffer, scan list of 1-2 neighbors for atoms I own
  // when find one, add its neighbors to 1-3 list
  //   increment the count in buf(i+4)
  //   exclude the atom whose tag = original
  //   this process may include duplicates but they will be culled later

  sptr = this;
  comm->ring(size,sizeof(int),buf,4,ring_four,buf);

  // fill onethree with buffer values that have been returned to me
  // sanity check: accumulated buf[i+3] count should equal
  //   nspecial[i][1] for each atom

  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (buf[j+3] != nspecial[i][1])
      error->one(FLERR,"1-3 bond count is inconsistent");
    j += 4 + nspecial[i][0];
    for (k = 0; k < nspecial[i][1]; k++)
      onethree[i][k] = buf[j++];
  }

  memory->destroy(buf);

  // done if special_bonds for 1-4 are set to 1.0

  if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    if (force->special_angle) angle_trim();
    combine();
    return;
  }

  // -----------------------------------------------------
  // compute nspecial[i][2] = # of 1-4 neighbors of atom i
  // -----------------------------------------------------

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 2 scalars + list of 1-3 neighbors

  nbuf = 0;
  for (i = 0; i < nlocal; i++) nbuf += 2 + nspecial[i][1];
  memory->create(buf,nbuf,"special:buf");

  // fill buffer with:
  // (1) = counter for 1-4 neighbors, initialized to 0
  // (2) = # of 1-3 neighbors
  // (3:N) = list of 1-3 neighbors

  size = 0;
  for (i = 0; i < nlocal; i++) {
    buf[size++] = 0;
    buf[size++] = nspecial[i][1];
    for (j = 0; j < nspecial[i][1]; j++) buf[size++] = onethree[i][j];
  }

  // cycle buffer around ring of procs back to self
  // when receive buffer, scan list of 1-3 neighbors for atoms I own
  // when find one, increment 1-4 count by # of 1-2 neighbors of my atom
  //   may include duplicates and original atom but they will be culled later

  sptr = this;
  comm->ring(size,sizeof(int),buf,5,ring_five,buf);

  // extract count from buffer that has cycled back to me
  // nspecial[i][2] = # of 1-4 neighbors of atom i

  j = 0;
  for (i = 0; i < nlocal; i++) {
    nspecial[i][2] = buf[j];
    j += 2 + nspecial[i][1];
  }

  memory->destroy(buf);

  // ----------------------------------------------------
  // create onefour[i] = list of 1-4 neighbors for atom i
  // ----------------------------------------------------

  max = 0;
  for (i = 0; i < nlocal; i++) max = MAX(max,nspecial[i][2]);
  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d = max # of 1-4 neighbors\n",maxall);
    if (logfile) fprintf(logfile,"  %d = max # of 1-4 neighbors\n",maxall);
  }

  memory->create(onefour,nlocal,maxall,"special:onefour");

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 3 scalars + list of 1-3 neighs + list of 1-4 neighs

  nbuf = 0;
  for (i = 0; i < nlocal; i++)
    nbuf += 3 + nspecial[i][1] + nspecial[i][2];
  memory->create(buf,nbuf,"special:buf");

  // fill buffer with:
  // (1) = # of 1-3 neighbors
  // (2) = # of 1-4 neighbors
  // (3) = counter for 1-4 neighbors, initialized to 0
  // (4:N) = list of 1-3 neighbors
  // (N+1:2N) space for list of 1-4 neighbors

  size = 0;
  for (i = 0; i < nlocal; i++) {
    buf[size++] = nspecial[i][1];
    buf[size++] = nspecial[i][2];
    buf[size++] = 0;
    for (j = 0; j < nspecial[i][1]; j++) buf[size++] = onethree[i][j];
    size += nspecial[i][2];
  }

  // cycle buffer around ring of procs back to self
  // when receive buffer, scan list of 1-3 neighbors for atoms I own
  // when find one, add its neighbors to 1-4 list
  //   incrementing the count in buf(i+4)
  //   this process may include duplicates but they will be culled later

  sptr = this;
  comm->ring(size,sizeof(int),buf,6,ring_six,buf);

  // fill onefour with buffer values that have been returned to me
  // sanity check: accumulated buf[i+2] count should equal
  //  nspecial[i][2] for each atom

  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (buf[j+2] != nspecial[i][2])
      error->one(FLERR,"1-4 bond count is inconsistent");
    j += 3 + nspecial[i][1];
    for (k = 0; k < nspecial[i][2]; k++)
      onefour[i][k] = buf[j++];
  }

  memory->destroy(buf);

  dedup();
  if (force->special_angle) angle_trim();
  if (force->special_dihedral) dihedral_trim();
  combine();
}

/* ----------------------------------------------------------------------
   remove duplicates within each of onetwo, onethree, onefour individually
------------------------------------------------------------------------- */

void Special::dedup()
{
  int i,j,m;

  // clear map so it can be used as scratch space

  atom->map_clear();

  // use map to cull duplicates
  // exclude original atom explicitly
  // adjust onetwo, onethree, onefour values to reflect removed duplicates
  // must unset map for each atom

  int **nspecial = atom->nspecial;
  int *tag = atom->tag;
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
  int i,j,m;

  int me;
  MPI_Comm_rank(world,&me);

  int **nspecial = atom->nspecial;
  int *tag = atom->tag;
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

  // compute global maxspecial, must be at least 1
  // add in extra factor from special_bonds command
  // allocate correct special array with same nmax, new maxspecial
  // previously allocated one must be destroyed
  // must make AtomVec class update its ptr to special

  MPI_Allreduce(&maxspecial,&atom->maxspecial,1,MPI_INT,MPI_MAX,world);
  atom->maxspecial += force->special_extra;
  atom->maxspecial = MAX(atom->maxspecial,1);

  if (me == 0) {
    if (screen)
      fprintf(screen,"  %d = max # of special neighbors\n",atom->maxspecial);
    if (logfile)
      fprintf(logfile,"  %d = max # of special neighbors\n",atom->maxspecial);
  }

  memory->destroy(atom->special);

  memory->create(atom->special,atom->nmax,atom->maxspecial,"atom:special");
  atom->avec->grow_reset();
  int **special = atom->special;

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
  int i,j,m,n;

  int *num_angle = atom->num_angle;
  int *num_dihedral = atom->num_dihedral;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom3 = atom->angle_atom3;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
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

    // dflag = flag for 1-3 neighs of all owned atoms

    int maxcount = 0;
    for (i = 0; i < nlocal; i++) maxcount = MAX(maxcount,nspecial[i][1]);
    memory->create(dflag,nlocal,maxcount,"special::dflag");

    for (i = 0; i < nlocal; i++) {
      n = nspecial[i][1];
      for (j = 0; j < n; j++) dflag[i][j] = 0;
    }

    // nbufmax = largest buffer needed to hold info from any proc
    // info for each atom = list of 1,3 atoms in each angle stored by atom
    //   and list of 1,3 and 2,4 atoms in each dihedral stored by atom

    int nbuf = 0;
    for (i = 0; i < nlocal; i++) nbuf += 2*num_angle[i] + 2*2*num_dihedral[i];
    int *buf;
    memory->create(buf,nbuf,"special:buf");

    // fill buffer with list of 1,3 atoms in each angle
    // and with list of 1,3 and 2,4 atoms in each dihedral

    int size = 0;
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_angle[i]; j++) {
        buf[size++] = angle_atom1[i][j];
        buf[size++] = angle_atom3[i][j];
      }
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++) {
        buf[size++] = dihedral_atom1[i][j];
        buf[size++] = dihedral_atom3[i][j];
        buf[size++] = dihedral_atom2[i][j];
        buf[size++] = dihedral_atom4[i][j];
      }

    // cycle buffer around ring of procs back to self
    // when receive buffer, scan list of 1,3 atoms looking for atoms I own
    // when find one, scan its 1-3 neigh list and mark I,J as in an angle

    sptr = this;
    comm->ring(size,sizeof(int),buf,7,ring_seven,NULL);

    // delete 1-3 neighbors if they are not flagged in dflag

    for (i = 0; i < nlocal; i++) {
      m = 0;
      for (j = 0; j < nspecial[i][1]; j++)
        if (dflag[i][j]) onethree[i][m++] = onethree[i][j];
      nspecial[i][1] = m;
    }

    // clean up

    memory->destroy(dflag);
    memory->destroy(buf);

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
  int i,j,m,n;

  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
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

  // if dihedrals are defined, flag each 1-4 neigh if it appears in a dihedral

  if (num_dihedral && atom->ndihedrals) {

    // dflag = flag for 1-4 neighs of all owned atoms

    int maxcount = 0;
    for (i = 0; i < nlocal; i++) maxcount = MAX(maxcount,nspecial[i][2]);
    memory->create(dflag,nlocal,maxcount,"special::dflag");

    for (i = 0; i < nlocal; i++) {
      n = nspecial[i][2];
      for (j = 0; j < n; j++) dflag[i][j] = 0;
    }

    // nbufmax = largest buffer needed to hold info from any proc
    // info for each atom = list of 1,4 atoms in each dihedral stored by atom

    int nbuf = 0;
    for (i = 0; i < nlocal; i++) nbuf += 2*num_dihedral[i];
    int *buf;
    memory->create(buf,nbuf,"special:buf");

    // fill buffer with list of 1,4 atoms in each dihedral

    int size = 0;
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++) {
        buf[size++] = dihedral_atom1[i][j];
        buf[size++] = dihedral_atom4[i][j];
      }

    // cycle buffer around ring of procs back to self
    // when receive buffer, scan list of 1,4 atoms looking for atoms I own
    // when find one, scan its 1-4 neigh list and mark I,J as in a dihedral

    sptr = this;
    comm->ring(size,sizeof(int),buf,8,ring_eight,NULL);

    // delete 1-4 neighbors if they are not flagged in dflag

    for (i = 0; i < nlocal; i++) {
      m = 0;
      for (j = 0; j < nspecial[i][2]; j++)
        if (dflag[i][j]) onefour[i][m++] = onefour[i][j];
      nspecial[i][2] = m;
    }

    // clean up

    memory->destroy(dflag);
    memory->destroy(buf);

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
   when receive buffer, scan tags for atoms I own
   when find one, increment nspecial count for that atom
------------------------------------------------------------------------- */

void Special::ring_one(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int *buf = (int *) cbuf;
  int m;

  for (int i = 0; i < ndatum; i++) {
    m = atom->map(buf[i]);
    if (m >= 0 && m < nlocal) nspecial[m][0]++;
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan 2nd-atom tags for atoms I own
   when find one, add 1st-atom tag to onetwo list for 2nd atom
------------------------------------------------------------------------- */

void Special::ring_two(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int nlocal = atom->nlocal;

  int **onetwo = sptr->onetwo;
  int *count = sptr->count;

  int *buf = (int *) cbuf;
  int m;

  for (int i = 1; i < ndatum; i += 2) {
    m = atom->map(buf[i]);
    if (m >= 0 && m < nlocal) onetwo[m][count[m]++] = buf[i-1];
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan list of 1-2 neighbors for atoms I own
   when find one, increment 1-3 count by # of 1-2 neighbors of my atom,
     subtracting one since my list will contain original atom
------------------------------------------------------------------------- */

void Special::ring_three(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int *buf = (int *) cbuf;
  int i,j,m,n,num12;

  i = 0;
  while (i < ndatum) {
    n = buf[i];
    num12 = buf[i+1];
    for (j = 0; j < num12; j++) {
      m = atom->map(buf[i+2+j]);
      if (m >= 0 && m < nlocal)
        n += nspecial[m][0] - 1;
    }
    buf[i] = n;
    i += 2 + num12;
  }
}

/* ----------------------------------------------------------------------
  when receive buffer, scan list of 1-2 neighbors for atoms I own
  when find one, add its neighbors to 1-3 list
    increment the count in buf(i+4)
    exclude the atom whose tag = original
    this process may include duplicates but they will be culled later
------------------------------------------------------------------------- */

void Special::ring_four(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int **onetwo = sptr->onetwo;

  int *buf = (int *) cbuf;
  int i,j,k,m,n,original,num12,num13;

  i = 0;
  while (i < ndatum) {
    original = buf[i];
    num12 = buf[i+1];
    num13 = buf[i+2];
    n = buf[i+3];
    for (j = 0; j < num12; j++) {
      m = atom->map(buf[i+4+j]);
      if (m >= 0 && m < nlocal)
        for (k = 0; k < nspecial[m][0]; k++)
          if (onetwo[m][k] != original)
            buf[i+4+num12+(n++)] = onetwo[m][k];
    }
    buf[i+3] = n;
    i += 4 + num12 + num13;
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan list of 1-3 neighbors for atoms I own
   when find one, increment 1-4 count by # of 1-2 neighbors of my atom
     may include duplicates and original atom but they will be culled later
------------------------------------------------------------------------- */

void Special::ring_five(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int *buf = (int *) cbuf;
  int i,j,m,n,num13;

  i = 0;
  while (i < ndatum) {
    n = buf[i];
    num13 = buf[i+1];
    for (j = 0; j < num13; j++) {
      m = atom->map(buf[i+2+j]);
      if (m >= 0 && m < nlocal) n += nspecial[m][0];
    }
      buf[i] = n;
      i += 2 + num13;
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan list of 1-3 neighbors for atoms I own
   when find one, add its neighbors to 1-4 list
     incrementing the count in buf(i+4)
     this process may include duplicates but they will be culled later
------------------------------------------------------------------------- */

void Special::ring_six(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int **onetwo = sptr->onetwo;

  int *buf = (int *) cbuf;
  int i,j,k,m,n,num13,num14;

  i = 0;
  while (i < ndatum) {
    num13 = buf[i];
    num14 = buf[i+1];
    n = buf[i+2];
    for (j = 0; j < num13; j++) {
      m = atom->map(buf[i+3+j]);
      if (m >= 0 && m < nlocal)
        for (k = 0; k < nspecial[m][0]; k++)
          buf[i+3+num13+(n++)] = onetwo[m][k];
    }
    buf[i+2] = n;
    i += 3 + num13 + num14;
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan list of 1,3 atoms looking for atoms I own
   when find one, scan its 1-3 neigh list and mark I,J as in an angle
------------------------------------------------------------------------- */

void Special::ring_seven(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int **onethree = sptr->onethree;
  int **dflag = sptr->dflag;

  int *buf = (int *) cbuf;
  int i,m,iglobal,jglobal,ilocal,jlocal;

  i = 0;
  while (i < ndatum) {
    iglobal = buf[i];
    jglobal = buf[i+1];
    ilocal = atom->map(iglobal);
    jlocal = atom->map(jglobal);
    if (ilocal >= 0 && ilocal < nlocal)
      for (m = 0; m < nspecial[ilocal][1]; m++)
        if (jglobal == onethree[ilocal][m]) {
          dflag[ilocal][m] = 1;
          break;
        }
    if (jlocal >= 0 && jlocal < nlocal)
      for (m = 0; m < nspecial[jlocal][1]; m++)
        if (iglobal == onethree[jlocal][m]) {
          dflag[jlocal][m] = 1;
          break;
        }
    i += 2;
  }
}

/* ----------------------------------------------------------------------
   when receive buffer, scan list of 1,4 atoms looking for atoms I own
   when find one, scan its 1-4 neigh list and mark I,J as in a dihedral
------------------------------------------------------------------------- */

void Special::ring_eight(int ndatum, char *cbuf)
{
  Atom *atom = sptr->atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int **onefour = sptr->onefour;
  int **dflag = sptr->dflag;

  int *buf = (int *) cbuf;
  int i,m,iglobal,jglobal,ilocal,jlocal;

  i = 0;
  while (i < ndatum) {
    iglobal = buf[i];
    jglobal = buf[i+1];
    ilocal = atom->map(iglobal);
    jlocal = atom->map(jglobal);
    if (ilocal >= 0 && ilocal < nlocal)
      for (m = 0; m < nspecial[ilocal][2]; m++)
        if (jglobal == onefour[ilocal][m]) {
          dflag[ilocal][m] = 1;
          break;
        }
    if (jlocal >= 0 && jlocal < nlocal)
      for (m = 0; m < nspecial[jlocal][2]; m++)
        if (iglobal == onefour[jlocal][m]) {
          dflag[jlocal][m] = 1;
          break;
        }
    i += 2;
  }
}
