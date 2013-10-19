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
   Contributing author:  Axel Kohlmeyer (ICTP)
------------------------------------------------------------------------- */

#include "fix_qmmm.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include "group.h"
#include "memory.h"

#include <stdlib.h>
#include <string.h>

#include "libqmmm.h"

// message tags for QM/MM inter communicator communication
// have to match with those from the QM code
enum {QMMM_TAG_OTHER=0, QMMM_TAG_SIZE=1, QMMM_TAG_COORD=2,QMMM_TAG_FORCE=3};

/****************************************************************************/

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
static int inthash_lookup(const inthash_t *tptr, int key);
/* generate list of keys for reverse lookups. */
static int *inthash_keys(inthash_t *tptr);
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
int inthash_lookup(const inthash_t *tptr, int key) {
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
 *  inthash_keys() - Return a list of keys.
 *  NOTE: the returned list must be freed with free(3).
 */
int *inthash_keys(inthash_t *tptr) {

  int *keys;
  inthash_node_t *node;

  keys = (int *)calloc(tptr->entries, sizeof(int));

  for (int i=0; i < tptr->size; ++i) {
    for (node=tptr->bucket[i]; node != NULL; node=node->next) {
      keys[node->data] = node->key;
    }
  }

  return keys;
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


using namespace LAMMPS_NS;
using namespace FixConst;

/* struct for packed data communication of coordinates and forces. */
struct commdata {
  int tag;
  float x,y,z;
};

/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID qmmm [couple <group-ID>]
 ***************************************************************/
FixQMMM::FixQMMM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  if (narg > 5)
    error->all(FLERR,"Illegal fix qmmm command");

  if (strcmp(update->unit_style,"metal") == 0) {
    qmmm_fscale = 1.0/23.0609;
  } else if (strcmp(update->unit_style,"real") == 0) {
    qmmm_fscale = 1.0;
  } else error->all(FLERR,"Fix qmmm requires real or metal units");

  /* define ec coupling group */
  mm_group = group->find("all");
  if ((narg == 5) && (0 == strcmp(arg[0],"couple"))) {
    mm_group = group->find(arg[4]);
  } else if (narg != 3) error->all(FLERR,"Illegal fix qmmm command");

  if (mm_group == -1)
    error->all(FLERR,"Could not find fix qmmm couple group ID");

  /* retrieve settings from global QM/MM configuration struct */
  comm_mode = qmmmcfg.comm_mode;
  qmmm_mode = qmmmcfg.qmmm_mode;
  qmmm_role = qmmmcfg.role;
  verbose   = qmmmcfg.verbose;

  if (comm_mode != QMMM_COMM_MPI)
    error->all(FLERR,"Only MPI communication mode is currently supported");

  if ((qmmm_role != QMMM_ROLE_MASTER) &&  (qmmm_role != QMMM_ROLE_SLAVE))
    error->all(FLERR,"LAMMPS can only function as MM master or MM slave");

  qm_comm = MPI_Comm_f2c(qmmmcfg.qm_comm);
  mm_comm = MPI_Comm_f2c(qmmmcfg.mm_comm);

  /* initialize storage */
  qm_idmap = mm_idmap = NULL;
  qm_remap = mm_remap = NULL;
  qm_coord = mm_coord = qm_force = mm_force = NULL;
  mm_type = NULL;

  do_init = 1;

  /* storage required to communicate a single coordinate or force. */
  size_one = sizeof(struct commdata);
  maxbuf = -1;
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixQMMM::~FixQMMM()
{

  if (qm_idmap) {
    inthash_t *hashtable = (inthash_t *)qm_idmap;
    inthash_destroy(hashtable);
    delete hashtable;
    free(qm_remap);
  }

  memory->destroy(comm_buf);
  memory->destroy(qm_coord);
  memory->destroy(qm_force);
}

/* ---------------------------------------------------------------------- */
int FixQMMM::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQMMM::exchange_positions()
{
  double **x = atom->x;
  const int * const mask  = atom->mask;
  const int * const tag  = atom->tag;
  const int nlocal = atom->nlocal;

  if ((comm->me == 0) && (verbose > 0)) {
    if (screen) fputs("QMMM: exchange positions\n",screen);
    if (logfile) fputs("QMMM: exchange positions\n",logfile);
  }

  if (qmmm_role == QMMM_ROLE_MASTER) {
    MPI_Status status;
    MPI_Request request;
    int i,tmp;

    /* check and potentially grow local communication buffers. */
    if (atom->nmax*size_one > maxbuf) {
      memory->destroy(comm_buf);
      maxbuf = atom->nmax*size_one;
      comm_buf = memory->smalloc(maxbuf,"qmmm:comm_buf");
    }
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);    

    if (comm->me == 0) {
      // insert local atoms into comm buffer
      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, tag[i]);
          if (j != 3*HASH_FAIL) {
            qm_coord[j]   = x[i][0];
            qm_coord[j+1] = x[i][1];
            qm_coord[j+2] = x[i][2];
          }
        }
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        int ndata;
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (int k=0; k<ndata; ++k) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, buf[k].tag);
          if (j != 3*HASH_FAIL) {
            qm_coord[j]   = buf[k].x;
            qm_coord[j+1] = buf[k].y;
            qm_coord[j+2] = buf[k].z;
          }
        }
      }

      /* done collecting coordinates, send it to dependent codes */
      if (comm->me == 0) {
        MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, qm_comm);
        MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, mm_comm);
      }

    } else {

      /* copy coordinate data into communication buffer */
      int ndata = 0;
      for (i=0; i<nlocal; ++i)
        if (mask[i] & groupbit) {
          buf[ndata].tag = tag[i];
          buf[ndata].x   = x[i][0];
          buf[ndata].y   = x[i][1];
          buf[ndata].z   = x[i][2];
          ++ndata;
        }

      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, ndata*size_one, MPI_BYTE, 0, 0, world);
    }
  } else if (qmmm_role == QMMM_ROLE_SLAVE) {

    MPI_Recv(qm_coord, 3*num_qm, MPI_DOUBLE, 0, QMMM_TAG_COORD,
             mm_comm, MPI_STATUS_IGNORE);
    // not needed for as long as we allow only one MPI task as slave
    MPI_Bcast(qm_coord, 3*num_qm, MPI_DOUBLE,0,world);

    /* update coordinates of (QM) atoms */
    for (int i=0; i < nlocal; ++i)
      if (mask[i] & groupbit)
        for (int j=0; j < num_qm; ++j)
          if (tag[i] == qm_remap[j]) {
            x[i][0] = qm_coord[3*j];
            x[i][1] = qm_coord[3*j+1];
            x[i][2] = qm_coord[3*j+2];
          }
  }
  return;
}

/* ---------------------------------------------------------------------- */

void FixQMMM::exchange_forces()
{
  double **f = atom->f;
  const int * const mask  = atom->mask;
  const int * const tag  = atom->tag;
  const int nlocal = atom->nlocal;

  if ((comm->me) == 0 && (verbose > 0)) {
    if (screen)  fputs("QMMM: exchange forces\n",screen);
    if (logfile) fputs("QMMM: exchange forces\n",logfile);
  }

  if (qmmm_role == QMMM_ROLE_MASTER) {
    MPI_Request req[2];
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);    

    if (comm->me == 0) {
      // receive QM forces from QE
      MPI_Irecv(qm_force,3*num_qm,MPI_DOUBLE,1,QMMM_TAG_FORCE,qm_comm,req);
      // receive MM forces from LAMMPS
      MPI_Irecv(qm_coord,3*num_qm,MPI_DOUBLE,1,QMMM_TAG_FORCE,mm_comm,req+1);
      MPI_Waitall(2,req,MPI_STATUS_IGNORE);

      // subtract MM forces from QM forces to get the delta
      // NOTE: QM forces are always sent in "real" units,
      // so we need to apply the scaling factor to get to the
      // supported internal units ("metal" or "real")
      for (int i=0; i < num_qm; ++i) {
        if  (verbose > 1) {
           const char fmt[] = "[%d]: QM(%g %g %g) MM(%g %g %g) /\\(%g %g %g)\n";
           if (screen) fprintf(screen, fmt, qm_remap[i],
                qmmm_fscale*qm_force[3*i+0], qmmm_fscale*qm_force[3*i+1],
                qmmm_fscale*qm_force[3*i+2], qm_coord[3*i+0], qm_coord[3*i+1],
                qm_coord[3*i+2],
                qmmm_fscale*qm_force[3*i+0] - qm_coord[3*i+0],
                qmmm_fscale*qm_force[3*i+1] - qm_coord[3*i+1],
                qmmm_fscale*qm_force[3*i+2] - qm_coord[3*i+2]);
           if (logfile) fprintf(logfile, fmt, qm_remap[i],
                qmmm_fscale*qm_force[3*i+0], qmmm_fscale*qm_force[3*i+1],
                qmmm_fscale*qm_force[3*i+2], qm_coord[3*i+0], qm_coord[3*i+1],
                qm_coord[3*i+2],
                qmmm_fscale*qm_force[3*i+0] - qm_coord[3*i+0],
                qmmm_fscale*qm_force[3*i+1] - qm_coord[3*i+1],
                qmmm_fscale*qm_force[3*i+2] - qm_coord[3*i+2]);
        }
        buf[i].tag = qm_remap[i];
        buf[i].x = qmmm_fscale*qm_force[3*i+0] - qm_coord[3*i+0];
        buf[i].y = qmmm_fscale*qm_force[3*i+1] - qm_coord[3*i+1];
        buf[i].z = qmmm_fscale*qm_force[3*i+2] - qm_coord[3*i+2];
    }
    MPI_Bcast(comm_buf,num_qm*size_one,MPI_BYTE,0,world);

    /* apply forces resulting from QM/MM coupling */
    for (int i=0; i < nlocal; ++i)
      if (mask[i] & groupbit)
        for (int j=0; j < num_qm; ++j)
          if (tag[i] == buf[j].tag) {
            f[i][0] += buf[j].x;
            f[i][0] += buf[j].y;
            f[i][0] += buf[j].z;
          }
    }

  } else if (qmmm_role == QMMM_ROLE_SLAVE) {

    memset(qm_force,0,3*num_qm*sizeof(double));
    for (int i=0; i < nlocal; ++i)
      if (mask[i] & groupbit) {
        const int j = 3*inthash_lookup((inthash_t *)qm_idmap, tag[i]);
        if (j != 3*HASH_FAIL) {
          qm_force[j]   = f[i][0];
          qm_force[j+1] = f[i][1];
          qm_force[j+2] = f[i][2];
        }
      }

    // collect and send MM slave forces to MM master
    // the reduction is not really needed with only one rank (for now)
    MPI_Reduce(qm_force,qm_coord,3*num_qm,MPI_DOUBLE,MPI_SUM,0,world);
    MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 0, QMMM_TAG_FORCE, mm_comm);
  }
  return;
}

/* ---------------------------------------------------------------------- */

void FixQMMM::init()
{
  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"Fix qmmm does not currently support r-RESPA");

  if (do_init) {
    int me = comm->me;
    do_init = 0;

    if (qmmm_role == QMMM_ROLE_MASTER) {
      MPI_Request req[2];
      int nat[2];

      if (me == 0) {
        // receive number of QM atoms from QE
        MPI_Irecv(nat, 1, MPI_INT, 1, QMMM_TAG_SIZE, qm_comm, req);
        // receive number of QM atoms from MM slave
        MPI_Irecv(nat+1, 1, MPI_INT, 1, QMMM_TAG_SIZE, mm_comm, req+1);
        MPI_Waitall(2,req,MPI_STATUS_IGNORE);
      }
      // broadcast across MM master processes
      MPI_Bcast(nat, 2, MPI_INT, 0, world);

      num_qm = group->count(igroup);
      num_mm = group->count(mm_group);

      // consistency check. the fix group and the QM and MM slave
      if ((num_qm != nat[0]) || (num_qm != nat[1]))
        error->all(FLERR,"Inconsistent number of QM/MM atoms");

      memory->create(qm_coord,3*num_qm,"qmmm:qm_coord");
      memory->create(qm_force,3*num_qm,"qmmm:qm_force");

      const char fmt[] = "Initializing QM/MM master with %d QM atoms\n";
      
      if (screen)  fprintf(screen,fmt,num_qm);
      if (logfile) fprintf(logfile,fmt,num_qm);

      if (qmmm_mode == QMMM_MODE_ELEC) {
        const char fm2[] = "Electrostatic coupling with %d atoms\n";
        if (screen)  fprintf(screen,fm2,num_mm);
        if (logfile) fprintf(logfile,fm2,num_mm);
      } 

    } else if (qmmm_role == QMMM_ROLE_SLAVE) {

      num_qm = group->count(igroup);

      if (me == 0) {
        /* send number of QM atoms to MM-master for confirmation */
        MPI_Send(&num_qm, 1, MPI_INT, 0, QMMM_TAG_SIZE, mm_comm);
      }
      memory->create(qm_coord,3*num_qm,"qmmm:qm_coord");
      memory->create(qm_force,3*num_qm,"qmmm:qm_force");

      const char fmt[] = "Initializing QM/MM slave with %d QM atoms\n";
      
      if (screen)  fprintf(screen,fmt,num_qm);
      if (logfile) fprintf(logfile,fmt,num_qm);

    }

    // communication buffer
    maxbuf = atom->nmax*size_one;
    comm_buf = (void *) memory->smalloc(maxbuf,"qmmm:comm_buf");

    /* initialize and build hashtable to map QM atoms */
    inthash_t *qm_hash=new inthash_t;
    inthash_init(qm_hash, num_qm);
    qm_idmap = (void *)qm_hash;

    MPI_Status status;
    MPI_Request request;
    const int nlocal = atom->nlocal;
    int i, j, tmp, ndata, qm_ntag;
    int *tag = atom->tag;
    int *mask  = atom->mask;
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);

    if (me == 0) {
      int *qm_taglist = new int[num_qm];
      qm_ntag = 0;
      for (i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit)
          qm_taglist[qm_ntag++] = tag[i];
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (j=0; j < ndata; ++j) {
          qm_taglist[qm_ntag++] = buf[j].tag;
        }
      }

      /* sort list of tags by value to have consistently the
       * same list when running in parallel and build hash table. */
      id_sort(qm_taglist, 0, num_qm-1);
      for (i=0; i < num_qm; ++i) {
        inthash_insert(qm_hash, qm_taglist[i], i);
      }
      delete[] qm_taglist;

      /* generate reverse index-to-tag map for communicating
       * qm/mm forces back to the proper atoms */
      qm_remap=inthash_keys(qm_hash);

      if (verbose > 1) {
        // print hashtable and reverse mapping
        for (i=0; i < num_qm; ++i) {
          if (screen) fprintf(screen, "qm_remap[%d]=%d  qm_hash[%d]=%d\n",
            i,qm_remap[i],qm_remap[i], inthash_lookup(qm_hash, qm_remap[i]));
          if (logfile) fprintf(logfile, "qm_remap[%d]=%d  qm_hash[%d]=%d\n",
            i,qm_remap[i],qm_remap[i], inthash_lookup(qm_hash, qm_remap[i]));
        }
      }

    } else {
      j = 0;
      for (i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          buf[j].x = -1.0;
          buf[j].tag = tag[i];
          ++j;
        }
      }
      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, j*size_one, MPI_BYTE, 0, 0, world);
    }

    // finally, after all is set up, do a first position synchronization
    exchange_positions();
  }
}

/* ---------------------------------------------------------------------- */

void FixQMMM::post_integrate()
{
  exchange_positions();
}

/* ---------------------------------------------------------------------- */

void FixQMMM::setup(int)
{
  exchange_forces();
}

/* ---------------------------------------------------------------------- */

void FixQMMM::post_force(int vflag)
{
  exchange_forces();
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixQMMM::memory_usage(void)
{
  double bytes;
  
  bytes = sizeof(FixQMMM);
  bytes += maxbuf;
  bytes += 6*num_qm*sizeof(double);

  return bytes;
}


// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
