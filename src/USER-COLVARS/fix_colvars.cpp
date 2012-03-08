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
   Contributing author:  Axel Kohlmeyer (TempleU)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>

#include "fix_colvars.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include "colvarproxy_lammps.h"

/* re-usable integer hash table code with static linkage. */

/** hash table top level data structure */
typedef struct inthash_t {
  struct inthash_node_t **bucket;        /* array of hash nodes */
  int size;                           /* size of the array */
  int entries;                        /* number of entries in table */
  int downshift;                      /* shift cound, used in hash function */
  int mask;                           /* used to select bits for hashing */
} inthash_t;

/** hash table node data structure */
typedef struct inthash_node_t {
  int data;                           /* data in hash node */
  int key;                            /* key for hash lookup */
  struct inthash_node_t *next;        /* next node in hash chain */
} inthash_node_t;

#define HASH_FAIL  -1
#define HASH_LIMIT  0.5

/* initialize new hash table  */
static void inthash_init(inthash_t *tptr, int buckets);
/* lookup entry in hash table */
static int inthash_lookup(void *tptr, int key);
/* insert an entry into hash table. */
static int inthash_insert(inthash_t *tptr, int key, int data);
/* delete the hash table */
static void inthash_destroy(inthash_t *tptr);
/* adapted sort for in-place sorting of map indices. */
static void id_sort(int *idmap, int left, int right);

/************************************************************************
 * integer hash code:
 ************************************************************************/

/* inthash() - Hash function returns a hash number for a given key.
 * tptr: Pointer to a hash table, key: The key to create a hash number for */
static int inthash(const inthash_t *tptr, int key) {
  int hashvalue;

  hashvalue = (((key*1103515249)>>tptr->downshift) & tptr->mask);
  if (hashvalue < 0) {
    hashvalue = 0;
  }    

  return hashvalue;
}

/*
 *  rebuild_table_int() - Create new hash table when old one fills up.
 *
 *  tptr: Pointer to a hash table
 */
static void rebuild_table_int(inthash_t *tptr) {
  inthash_node_t **old_bucket, *old_hash, *tmp;
  int old_size, h, i;

  old_bucket=tptr->bucket;
  old_size=tptr->size;

  /* create a new table and rehash old buckets */
  inthash_init(tptr, old_size<<1);
  for (i=0; i<old_size; i++) {
    old_hash=old_bucket[i];
    while(old_hash) {
      tmp=old_hash;
      old_hash=old_hash->next;
      h=inthash(tptr, tmp->key);
      tmp->next=tptr->bucket[h];
      tptr->bucket[h]=tmp;
      tptr->entries++;
    } /* while */
  } /* for */

  /* free memory used by old table */
  free(old_bucket);

  return;
}

/*
 *  inthash_init() - Initialize a new hash table.
 *
 *  tptr: Pointer to the hash table to initialize
 *  buckets: The number of initial buckets to create
 */
void inthash_init(inthash_t *tptr, int buckets) {

  /* make sure we allocate something */
  if (buckets==0)
    buckets=16;

  /* initialize the table */
  tptr->entries=0;
  tptr->size=2;
  tptr->mask=1;
  tptr->downshift=29;

  /* ensure buckets is a power of 2 */
  while (tptr->size<buckets) {
    tptr->size<<=1;
    tptr->mask=(tptr->mask<<1)+1;
    tptr->downshift--;
  } /* while */

  /* allocate memory for table */
  tptr->bucket=(inthash_node_t **) calloc(tptr->size, sizeof(inthash_node_t *));

  return;
}

/*
 *  inthash_lookup() - Lookup an entry in the hash table and return a pointer to
 *    it or HASH_FAIL if it wasn't found.
 *
 *  tptr: Pointer to the hash table
 *  key: The key to lookup
 */
int inthash_lookup(void *ptr, int key) {
  const inthash_t *tptr = (const inthash_t *) ptr;
  int h;
  inthash_node_t *node;


  /* find the entry in the hash table */
  h=inthash(tptr, key);
  for (node=tptr->bucket[h]; node!=NULL; node=node->next) {
    if (node->key == key)
      break;
  }

  /* return the entry if it exists, or HASH_FAIL */
  return(node ? node->data : HASH_FAIL);
}

/*
 *  inthash_insert() - Insert an entry into the hash table.  If the entry already
 *  exists return a pointer to it, otherwise return HASH_FAIL.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to insert into the hash table
 *  data: A pointer to the data to insert into the hash table
 */
int inthash_insert(inthash_t *tptr, int key, int data) {
  int tmp;
  inthash_node_t *node;
  int h;

  /* check to see if the entry exists */
  if ((tmp=inthash_lookup(tptr, key)) != HASH_FAIL)
    return(tmp);

  /* expand the table if needed */
  while (tptr->entries>=HASH_LIMIT*tptr->size)
    rebuild_table_int(tptr);

  /* insert the new entry */
  h=inthash(tptr, key);
  node=(struct inthash_node_t *) malloc(sizeof(inthash_node_t));
  node->data=data;
  node->key=key;
  node->next=tptr->bucket[h];
  tptr->bucket[h]=node;
  tptr->entries++;

  return HASH_FAIL;
}

/*
 * inthash_destroy() - Delete the entire table, and all remaining entries.
 * 
 */
void inthash_destroy(inthash_t *tptr) {
  inthash_node_t *node, *last;
  int i;

  for (i=0; i<tptr->size; i++) {
    node = tptr->bucket[i];
    while (node != NULL) { 
      last = node;   
      node = node->next;
      free(last);
    }
  }     

  /* free the entire array of buckets */
  if (tptr->bucket != NULL) {
    free(tptr->bucket);
    memset(tptr, 0, sizeof(inthash_t));
  }
}

/************************************************************************
 * integer list sort code:
 ************************************************************************/

/* sort for integer map. initial call  id_sort(idmap, 0, natoms - 1); */
static void id_sort(int *idmap, int left, int right)
{
  int pivot, l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = idmap[left];

  while (left < right) {
    while ((idmap[right] >= pivot) && (left < right))
      right--;
    if (left != right) {
      idmap[left] = idmap[right];
      left++;
    }
    while ((idmap[left] <= pivot) && (left < right))
      left++;
    if (left != right) {
      idmap[right] = idmap[left];
      right--;
    }
  }
  idmap[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;

  if (left < pivot)
    id_sort(idmap, left, pivot-1);
  if (right > pivot)
    id_sort(idmap, pivot+1, right);
}

/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;

// initialize static class members
int FixColvars::instances=0;

/***************************************************************
 create class and parse arguments in LAMMPS script. Syntax: 

 fix ID group-ID colvars <config_file> [optional flags...]

 optional keyword value pairs:

  input   <input prefix>    (for restarting/continuing, defaults to
                             NULL, but set to <output prefix> at end)
  output  <output prefix>   (defaults to 'out')
  seed    <integer>         (seed for RNG, defaults to '1966')
  tstat   <fix label>       (label of thermostatting fix)

 TODO: add (optional) arguments for RNG seed, temperature compute
 ***************************************************************/
FixColvars::FixColvars(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix colvars command: too few arguments");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix colvars without atom IDs defined");

  if (atom->rmass_flag)
    error->all(FLERR,"Cannot use fix colvars for atoms with rmass attribute");

  if (instances)
    error->all(FLERR,"Only one fix colvars can be active at a time");
  ++instances;

  scalar_flag = 1;
  global_freq = 1;
  nevery = 1;
  extscalar = 1;

  me = comm->me;
  
  conf_file = strdup(arg[3]);
  rng_seed = 1966;

  inp_name = NULL;
  out_name = NULL;
  tmp_name = NULL;

  /* parse optional arguments */
  int argsdone = 4;
  while (argsdone+1 < narg) {
    if (0 == strcmp(arg[argsdone], "input")) {
      inp_name = strdup(arg[argsdone+1]);
    } else if (0 == strcmp(arg[argsdone], "output")) {
      out_name = strdup(arg[argsdone+1]);
    } else if (0 == strcmp(arg[argsdone], "seed")) {
      rng_seed = atoi(arg[argsdone+1]);
    } else if (0 == strcmp(arg[argsdone], "tstat")) {
      tmp_name = strdup(arg[argsdone+1]);
    } else {
      error->all(FLERR,"Unknown fix imd parameter");
    }
    ++argsdone; ++argsdone;
  }

  if (!out_name) out_name = strdup("out");

  /* initialize various state variables. */
  tstat_id = -1;
  energy = 0.0;
  nlevels_respa = 0;
  num_coords = 0;
  coords = forces = oforce = comm_buf = NULL;
  proxy = NULL;
  idmap = NULL;

  /* storage required to communicate a single coordinate or force. */
  size_one = sizeof(struct commdata);
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixColvars::~FixColvars()
{
  memory->sfree(conf_file);
  memory->sfree(inp_name);
  memory->sfree(out_name);
  memory->sfree(tmp_name);
  deallocate();
  --instances;
}

/* ---------------------------------------------------------------------- */

void FixColvars::deallocate()
{
  memory->destroy(comm_buf);
  memory->destroy(forces);

  if (proxy) {
    delete proxy;
    inthash_t *hashtable = (inthash_t *)idmap;
    memory->sfree(coords);
    memory->sfree(oforce);
    inthash_destroy(hashtable);
    delete hashtable;
  }

  proxy = NULL;
  idmap = NULL;
  coords = NULL;
  forces = NULL;
  oforce = NULL;
  comm_buf = NULL;
}

/* ---------------------------------------------------------------------- */

void FixColvars::post_run()
{
  deallocate();
  memory->sfree(inp_name);
  inp_name = strdup(out_name);
}

/* ---------------------------------------------------------------------- */

int FixColvars::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= POST_RUN;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixColvars::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  return;
}

/* ---------------------------------------------------------------------- */

// initial setup of colvars run.

void FixColvars::setup(int)
{

  // nme:    number of atoms in group on this MPI task

  int i,nme;
  const int * const mask = atom->mask;
  const int * const tag  = atom->tag;
  const int * const type = atom->type;
  const double * const * const x = atom->x;
  const int nlocal = atom->nlocal;

  bigint n = group->count(igroup);
  if (n > MAXSMALLINT) error->all(FLERR,"Too many atoms for fix colvars");
  num_coords = static_cast<int> (n);

  nme=0;
  for (i=0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nme;

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  memory->create(comm_buf,nmax,"colvars:comm_buf");
  memory->create(forces,num_coords,"colvars:forces");
  memset(forces,0,num_coords*size_one);

  /* initialize and build hashtable. */
  inthash_t *hashtable=new inthash_t;
  inthash_init(hashtable, num_coords);
  idmap = (void *)hashtable;

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;

  if (me == 0) {
    memory->create(coords,num_coords,"colvars:forces");
    memory->create(oforce,num_coords,"colvars:oforce");

    int *taglist = new int[num_coords];

    // counter to map atom tags to a 0-based consecutive index list
    int numtag=0;
    
    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        taglist[numtag] = tag[i];
        ++numtag;
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv((void *)comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int j=0; j < ndata; ++j) {
        taglist[numtag] = comm_buf[j].tag;
        ++numtag;
      }
    }

    /* sort list of tags by value to have consistently the 
     * same list when running in parallel and build hash table. */
    id_sort(taglist, 0, num_coords-1);
    for (i=0; i < num_coords; ++i) {
      inthash_insert(hashtable, taglist[i], i);
    }
    delete[] taglist;

  } else {
    nme=0;
    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        comm_buf[nme].tag = tag[i];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  // after initializing the forward/backward tags we can do
  // a first coordinate communication.

  if (me == 0) {

    // store coordinate data in holding array, clear old forces

    for (i=0; i<nlocal; ++i) {
      if (mask[i] & groupbit) {
        const int j = inthash_lookup(idmap, tag[i]);
        if (j != HASH_FAIL) {
          oforce[j].tag  = coords[j].tag  = tag[i];
          oforce[j].type = coords[j].type = type[i];
          coords[j].x = x[i][0];
          coords[j].y = x[i][1];
          coords[j].z = x[i][2];
          oforce[j].x = oforce[j].y = oforce[j].z = 0.0;
        }
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {
        const int j = inthash_lookup(idmap, comm_buf[k].tag);
        if (j != HASH_FAIL) {
          oforce[j].tag  = coords[j].tag  = comm_buf[k].tag;
          oforce[j].type = coords[j].type = comm_buf[k].type;
          coords[j].x = comm_buf[k].x;
          coords[j].y = comm_buf[k].y;
          coords[j].z = comm_buf[k].z;
          oforce[j].x = oforce[j].y = oforce[j].z = 0.0;
        }
      }
    }

  } else { // me != 0
    /* copy coordinate data into communication buffer */
    nme = 0;
    for (i=0; i<nlocal; ++i) {
      if (mask[i] & groupbit) {
        comm_buf[nme].tag  = tag[i];
        comm_buf[nme].type = type[i];
        comm_buf[nme].x    = x[i][0];
        comm_buf[nme].y    = x[i][1];
        comm_buf[nme].z    = x[i][2];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  // now we are ready to initialize the colvar proxy
  if (me == 0) {
    double t_target = 0.0;
    if (tmp_name) {
      if (strcmp(tmp_name,"NULL") == 0)
        tstat_id = -1;
      else {
        tstat_id = modify->find_fix(tmp_name);
        if (tstat_id < 0) error->one(FLERR,"Could not find tstat fix ID");
        double *tt = (double*)modify->fix[tstat_id]->extract("t_target",tmp);
        if (tt) t_target = *tt;
      }
    }
    proxy = new colvarproxy_lammps(lmp,conf_file,inp_name,out_name,rng_seed,
                                   t_target,coords,forces,oforce,
                                   &inthash_lookup,idmap);
  }
  return;
}

/* ---------------------------------------------------------------------- */
/* Main colvars handler:
 * Send coodinates and add colvar forces to atoms. */
void FixColvars::post_force(int vflag)
{
  // some housekeeping: update status of the proxy as needed.
  if (me == 0) {
    if (tstat_id < 0) {
      proxy->set_temperature(0.0);
    } else {
      int tmp;
      // get thermostat target temperature from corresponding fix,
      // if the fix supports extraction.
      double *tt = (double *) modify->fix[tstat_id]->extract("t_target",tmp);
      if (tt) 
        proxy->set_temperature(*tt);
      else 
        proxy->set_temperature(0.0);
    }
  }

  const int * const mask  = atom->mask;
  const int * const tag = atom->tag;
  const int * const type = atom->type;
  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const int nlocal = atom->nlocal;

  /* check and potentially grow local communication buffers. */
  int i, k, nmax_new, nme=0;
  for (i=0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nme;

  MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
  if (nmax_new > nmax) {
    nmax = nmax_new;
    memory->grow(comm_buf,nmax,"colvars:comm_buf");
  }

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;
  
  if (me == 0) {

    // store coordinate data in holding array

    for (i=0; i<nlocal; ++i) {
      if (mask[i] & groupbit) {
        const int j = inthash_lookup(idmap, tag[i]);
        if (j != HASH_FAIL) {
          coords[j].tag  = tag[i];
          coords[j].type = type[i];
          coords[j].x = x[i][0];
          coords[j].y = x[i][1];
          coords[j].z = x[i][2];
        }
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {
        const int j = inthash_lookup(idmap, comm_buf[k].tag);
        if (j != HASH_FAIL) {
          coords[j].tag  = comm_buf[k].tag;
          coords[j].type = comm_buf[k].type;
          coords[j].x = comm_buf[k].x;
          coords[j].y = comm_buf[k].y;
          coords[j].z = comm_buf[k].z;
        }
      }
    }

  } else { // me != 0
    /* copy coordinate data into communication buffer */
    nme = 0;
    for (i=0; i<nlocal; ++i) {
      if (mask[i] & groupbit) {
        comm_buf[nme].tag  = tag[i];
        comm_buf[nme].type = type[i];
        comm_buf[nme].x    = x[i][0];
        comm_buf[nme].y    = x[i][1];
        comm_buf[nme].z    = x[i][2];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  ////////////////////////////////////////////////////////////////////////
  // call our workhorse and retrieve additional information.
  if (me == 0) {
    energy = proxy->compute();
    store_forces = proxy->need_system_forces();
  }
  ////////////////////////////////////////////////////////////////////////

  // broadcast force and energy data to all processors
  MPI_Bcast(forces, num_coords*size_one, MPI_BYTE, 0, world);
  MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&store_forces, 1, MPI_INT, 0, world);
  
  /* XXX. this is in principle O(N**2), i.e. not good. 
   * however we assume for now that the number of atoms 
   * that we manipulate via colvars will be small compared 
   * to the total system size, so we don't hurt too much. */
  for (int j=0; j < num_coords; ++j) {
    for (int i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        if (forces[j].tag == tag[i]) {
          f[i][0] += forces[j].x;
          f[i][1] += forces[j].y;
          f[i][2] += forces[j].z;
        }
      }
    }
  }
  return;
}

/* ---------------------------------------------------------------------- */
void FixColvars::post_force_respa(int vflag, int ilevel, int iloop)
{
  /* only process colvar forces on the outmost RESPA level. */
  if (ilevel == nlevels_respa-1) post_force(vflag);
  return;
}

/* ---------------------------------------------------------------------- */
void FixColvars::end_of_step()
{
  if (store_forces) {
    
    const int * const mask  = atom->mask;
    const int * const tag = atom->tag;
    const double * const * const f = atom->f;
    const int nlocal = atom->nlocal;

    /* check and potentially grow local communication buffers. */
    int i, k, nmax_new, nme=0;
    for (i=0; i < nlocal; ++i)
      if (mask[i] & groupbit) ++nme;

    MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
    if (nmax_new > nmax) {
      nmax = nmax_new;
      memory->grow(comm_buf,nmax,"colvars:comm_buf");
    }

    MPI_Status status;
    MPI_Request request;
    int tmp, ndata;
  
    if (me == 0) {

      // store old force data in holding array

      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          const int j = inthash_lookup(idmap, tag[i]);
          if (j != HASH_FAIL) {
            oforce[j].tag  = tag[i];
            oforce[j].type = 1;
            oforce[j].x = f[i][0];
            oforce[j].y = f[i][1];
            oforce[j].z = f[i][2];
          }
        }
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        int maxbuf = nmax*size_one;
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (int k=0; k<ndata; ++k) {
          const int j = inthash_lookup(idmap, comm_buf[k].tag);
          if (j != HASH_FAIL) {
            oforce[j].tag  = comm_buf[k].tag;
            oforce[j].type = comm_buf[k].type;
            oforce[j].x = comm_buf[k].x;
            oforce[j].y = comm_buf[k].y;
            oforce[j].z = comm_buf[k].z;
        }
      }
    }

  } else { // me != 0
    /* copy coordinate data into communication buffer */
    nme = 0;
    for (i=0; i<nlocal; ++i) {
      if (mask[i] & groupbit) {
        comm_buf[nme].tag  = tag[i];
        comm_buf[nme].type = 1;
        comm_buf[nme].x    = f[i][0];
        comm_buf[nme].y    = f[i][1];
        comm_buf[nme].z    = f[i][2];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixColvars::compute_scalar()
{
  return energy;
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixColvars::memory_usage(void)
{
  return static_cast<double>(2*num_coords+nmax)*size_one;
}

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
