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
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

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
  memory->destroy_2d_int_array(onetwo);
  memory->destroy_2d_int_array(onethree);
  memory->destroy_2d_int_array(onefour);
}

/* ----------------------------------------------------------------------
   create 1-2, 1-3, 1-4 lists of topology neighbors
   store in onetwo, onethree, onefour for each atom
   store 3 counters in nspecial[i] 
------------------------------------------------------------------------- */

void Special::build()
{
  int i,j,k,m,n,loop,size,original;
  int num12,num13,num14;
  int max,maxall,messtag,nbuf,nbufmax;
  int *buf,*bufcopy,*count;
  MPI_Request request;
  MPI_Status status;

  MPI_Barrier(world);

  int nlocal = atom->nlocal;

  int *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **nspecial = atom->nspecial;

  if (me == 0 && screen) fprintf(screen,"Finding 1-2 1-3 1-4 neighbors ...\n");

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

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
    MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

    buf = new int[nbufmax];
    bufcopy = new int[nbufmax];

    // fill buffer with global tags of bond partners of my atoms

    size = 0;
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++)
	buf[size++] = bond_atom[i][j];

    // cycle buffer around ring of procs back to self
    // when receive buffer, scan tags for atoms I own
    // when find one, increment nspecial count for that atom

    messtag = 1;
    for (loop = 0; loop < nprocs; loop++) {
      for (i = 0; i < size; i++) {
	m = atom->map(buf[i]);
	if (m >= 0 && m < nlocal) nspecial[m][0]++;
      }
      if (me != next) {
	MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
	MPI_Send(buf,size,MPI_INT,next,messtag,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&size);
	for (j = 0; j < size; j++) buf[j] = bufcopy[j];
      }
    }

    delete [] buf;
    delete [] bufcopy;
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

  onetwo = memory->create_2d_int_array(nlocal,maxall,"special:onetwo");

  // count = accumulating counter

  count = new int[nlocal];
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
    MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

    buf = new int[nbufmax];
    bufcopy = new int[nbufmax];

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

    messtag = 2;
    for (loop = 0; loop < nprocs; loop++) {
      for (i = 1; i < size; i += 2) {
	m = atom->map(buf[i]);
	if (m >= 0 && m < nlocal) onetwo[m][count[m]++] = buf[i-1];
      }
      if (me != next) {
	MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
	MPI_Send(buf,size,MPI_INT,next,messtag,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&size);
	for (j = 0; j < size; j++) buf[j] = bufcopy[j];
      }
    }

    delete [] buf;
    delete [] bufcopy;
  }

  delete [] count;

  // -----------------------------------------------------
  // done if special_bonds for 1-3, 1-4 are set to 1.0
  // -----------------------------------------------------

  if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0 &&
      force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
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
  MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

  buf = new int[nbufmax];
  bufcopy = new int[nbufmax];

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

  messtag = 3;
  for (loop = 0; loop < nprocs; loop++) {
    i = 0;
    while (i < size) {
      n = buf[i];
      num12 = buf[i+1];
      for (j = 0; j < num12; j++) {
	m = atom->map(buf[i+2+j]);
	if (m >= 0 && m < nlocal) n += nspecial[m][0] - 1;
      }
      buf[i] = n;
      i += 2 + num12;
    }
    if (me != next) {
      MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
      MPI_Send(buf,size,MPI_INT,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_INT,&size);
      for (j = 0; j < size; j++) buf[j] = bufcopy[j];
    }
  }

  // extract count from buffer that has cycled back to me
  // nspecial[i][1] = # of 1-3 neighbors of atom i

  j = 0;
  for (i = 0; i < nlocal; i++) {
    nspecial[i][1] = buf[j];
    j += 2 + nspecial[i][0];
  }

  delete [] buf;
  delete [] bufcopy;

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

  onethree = memory->create_2d_int_array(nlocal,maxall,"special:onethree");

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 4 scalars + list of 1-2 neighs + list of 1-3 neighs

  nbuf = 0;
  for (i = 0; i < nlocal; i++) nbuf += 4 + nspecial[i][0] + nspecial[i][1];
  MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

  buf = new int[nbufmax];
  bufcopy = new int[nbufmax];

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

  messtag = 4;
  for (loop = 0; loop < nprocs; loop++) {
    i = 0;
    while (i < size) {
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
    if (me != next) {
      MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
      MPI_Send(buf,size,MPI_INT,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_INT,&size);
      for (j = 0; j < size; j++) buf[j] = bufcopy[j];
    }
  }

  // fill onethree with buffer values that have been returned to me
  // sanity check: accumulated buf[i+3] count should equal
  //   nspecial[i][1] for each atom

  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (buf[j+3] != nspecial[i][1])
      error->one("1-3 bond count is inconsistent");
    j += 4 + nspecial[i][0];
    for (k = 0; k < nspecial[i][1]; k++) 
      onethree[i][k] = buf[j++];
  }

  delete [] buf;
  delete [] bufcopy;

  // done if special_bonds for 1-4 are set to 1.0

  if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
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
  MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

  buf = new int[nbufmax];
  bufcopy = new int[nbufmax];

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

  messtag = 5;
  for (loop = 0; loop < nprocs; loop++) {
    i = 0;
    while (i < size) {
      n = buf[i];
      num13 = buf[i+1];
      for (j = 0; j < num13; j++) {
	m = atom->map(buf[i+2+j]);
	if (m >= 0 && m < nlocal) n += nspecial[m][0];
      }
      buf[i] = n;
      i += 2 + num13;
    }
    if (me != next) {
      MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
      MPI_Send(buf,size,MPI_INT,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_INT,&size);
      for (j = 0; j < size; j++) buf[j] = bufcopy[j];
    }
  }

  // extract count from buffer that has cycled back to me
  // nspecial[i][2] = # of 1-4 neighbors of atom i

  j = 0;
  for (i = 0; i < nlocal; i++) {
    nspecial[i][2] = buf[j];
    j += 2 + nspecial[i][1];
  }

  delete [] buf;
  delete [] bufcopy;

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

  onefour = memory->create_2d_int_array(nlocal,maxall,"special:onefour");

  // nbufmax = largest buffer needed to hold info from any proc
  // info for each atom = 3 scalars + list of 1-3 neighs + list of 1-4 neighs

  nbuf = 0;
  for (i = 0; i < nlocal; i++) 
    nbuf += 3 + nspecial[i][1] + nspecial[i][2];
  MPI_Allreduce(&nbuf,&nbufmax,1,MPI_INT,MPI_MAX,world);

  buf = new int[nbufmax];
  bufcopy = new int[nbufmax];

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

  messtag = 6;
  for (loop = 0; loop < nprocs; loop++) {
    i = 0;
    while (i < size) {
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
    if (me != next) {
      MPI_Irecv(bufcopy,nbufmax,MPI_INT,prev,messtag,world,&request);
      MPI_Send(buf,size,MPI_INT,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_INT,&size);
      for (j = 0; j < size; j++) buf[j] = bufcopy[j];
    }
  }

  // fill onefour with buffer values that have been returned to me
  // sanity check: accumulated buf[i+2] count should equal
  //  nspecial[i][2] for each atom

  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (buf[j+2] != nspecial[i][2])
      error->one("1-4 bond count is inconsistent");
    j += 3 + nspecial[i][1];
    for (k = 0; k < nspecial[i][2]; k++) 
      onefour[i][k] = buf[j++];
  }

  delete [] buf;
  delete [] bufcopy;

  combine();
}

/* ----------------------------------------------------------------------
   concatenate onetwo, onethree, onefour into master atom->special list
   remove duplicates
   convert nspecial[0], nspecial[1], nspecial[2] into cumulative counters 
------------------------------------------------------------------------- */

void Special::combine()
{
  int i,j,m;

  int me;
  MPI_Comm_rank(world,&me);

  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  // ----------------------------------------------------
  // compute culled maxspecial = max # of special neighs of any atom
  // ----------------------------------------------------

  // clear map so it can be used as scratch space

  atom->map_clear();

  // unique = # of unique nspecial neighbors of one atom
  // cull duplicates using map to check for them
  // exclude original atom explicitly
  // must re-clear map for each atom

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
  // allocate correct special array with same nmax, new maxspecial
  // previously allocated one must be destroyed
  // must force AtomVec class to update its ptr to special

  MPI_Allreduce(&maxspecial,&atom->maxspecial,1,MPI_INT,MPI_MAX,world);
  atom->maxspecial = MAX(atom->maxspecial,1);

  if (me == 0) {
    if (screen)
      fprintf(screen,"  %d = max # of special neighbors\n",atom->maxspecial);
    if (logfile)
      fprintf(logfile,"  %d = max # of special neighbors\n",atom->maxspecial);
  }

  memory->destroy_2d_int_array(atom->special);

  atom->special = 
    memory->create_2d_int_array(atom->nmax,atom->maxspecial,"atom:special");
  atom->avec->reset_special();
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

  atom->map_set();
}
