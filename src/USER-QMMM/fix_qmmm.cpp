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
// #include "domain.h"
#include "force.h"
#include "error.h"
#include "group.h"
#include "memory.h"

// #include <math.h>
// #include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <errno.h>

#include "libqmmm.h"

// message tags for QM/MM inter communicator communication

enum {QMMM_TAG_OTHER=0, QMMM_TAG_SIZE=1, QMMM_TAG_COORD=2,QMMM_TAG_FORCE=3};


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
    qmmm_fscale = 23.0609;
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
  if (mm_idmap) {
    inthash_t *hashtable = (inthash_t *)mm_idmap;
    inthash_destroy(hashtable);
    delete hashtable;
    free(mm_remap);
  }

  memory->destroy(comm_buf);
  memory->destroy(qm_coord);
  memory->destroy(mm_coord);
  memory->destroy(qm_force);
  memory->destroy(mm_force);
  memory->destroy(mm_type);
}

/* ---------------------------------------------------------------------- */
int FixQMMM::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
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
        error->all(FLERR,"Inconsistant number of QM/MM atoms");

      memory->create(qm_coord,num_qm,"qmmm:qm_coord");
      memory->create(mm_coord,num_mm,"qmmm:mm_coord");

    } else if (qmmm_role == QMMM_ROLE_SLAVE) {

      num_qm = group->count(igroup);
      num_mm = group->count(mm_group);

      if (me == 0) {
        /* send number of QM atoms to MM-slave */
        MPI_Send(&num_qm, 1, MPI_INT, 0, QMMM_TAG_SIZE, mm_comm);
      }
      memory->create(qm_coord,num_qm,"qmmm:qm_coord");
      memory->create(mm_coord,num_mm,"qmmm:mm_coord");
    }

    // communication buffer
    maxbuf = atom->nmax*size_one;
    comm_buf = (void *) memory->smalloc(maxbuf,"qmmm:comm_buf");

    /* initialize and build hashtables. */
    inthash_t *qm_hash=new inthash_t;
    inthash_t *mm_hash=new inthash_t;
    inthash_init(qm_hash, num_qm);
    inthash_init(mm_hash, num_mm);
    qm_idmap = (void *)qm_hash;
    mm_idmap = (void *)mm_hash;
    mm_grbit = group->bitmask[mm_group];

    MPI_Status status;
    MPI_Request request;
    int i, j, tmp, ndata, qm_ntag, mm_ntag, nlocal = atom->nlocal;
    int *tag = atom->tag;
    int *mask  = atom->mask;
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);

    if (me == 0) {
      int *qm_taglist = new int[num_qm];
      int *mm_taglist = new int[num_mm];
      qm_ntag = mm_ntag = 0;
      for (i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit)
          qm_taglist[qm_ntag++] = tag[i];
        if (mask[i] & mm_grbit)
          mm_taglist[mm_ntag++] = tag[i];
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (j=0; j < ndata; ++j) {
          if (buf[j].x < 0.0)
            qm_taglist[qm_ntag++] = buf[j].tag;
          else
            mm_taglist[mm_ntag++] = buf[j].tag;
        }
      }

      /* sort list of tags by value to have consistently the
       * same list when running in parallel and build hash table. */
      id_sort(qm_taglist, 0, num_qm);
      id_sort(mm_taglist, 0, num_mm);
      for (i=0; i < num_qm; ++i) {
        inthash_insert(qm_hash, qm_taglist[i], i);
      }
      for (i=0; i < num_mm; ++i) {
        inthash_insert(mm_hash, mm_taglist[i], i);
      }
      delete[] qm_taglist;
      delete[] mm_taglist;

      /* generate reverse index-to-tag map for communicating
       * qm/mm forces back to the proper atoms */
      qm_remap=inthash_keys(qm_hash);
      mm_remap=inthash_keys(mm_hash);
    } else {
      j = 0;
      for (i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          buf[j].x = -1.0;
          buf[j].tag = tag[i];
          ++j;
        } else if (mask[i] & mm_grbit) {
          buf[j].x = 1.0;
          buf[j].tag = tag[i];
          ++j;
        }
      }
      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, j*size_one, MPI_BYTE, 0, 0, world);
    }
  }
}

void FixQMMM::exchange_positions()
{

  if (qmmm_role == QMMM_ROLE_MASTER) {
    int i,j,k;
    int *mask  = atom->mask;
    int *tag  = atom->tag;
    int nlocal = atom->nlocal;
    double **x = atom->x;

    /* check and potentially grow local communication buffers. */
    if (atom->nmax*size_one > maxbuf) {
      memory->destroy(comm_buf);
      maxbuf = atom->nmax*size_one;
      comm_buf = memory->smalloc(maxbuf,"qmmm:comm_buf");
    }
    
    MPI_Status status;
    MPI_Request request;
    int tmp, ndata;
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);

    if (comm->me == 0) {
      // insert local atoms into comm buffer
      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, tag[i]);
          if (j != HASH_FAIL) {
            qm_coord[j]   = x[i][0];
            qm_coord[j+1] = x[i][1];
            qm_coord[j+2] = x[i][2];
          }
        }
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (k=0; k<ndata; ++k) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, buf[k].tag);
          if (j != HASH_FAIL) {
            qm_coord[j]   = buf[k].x;
            qm_coord[j+1] = buf[k].y;
            qm_coord[j+2] = buf[k].z;
          }
        }
      }

      /* done collecting frame data now send it to QM and MM slave. */
      MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, qm_comm);
      printf("MM master: after send coords to QM\n");
      MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, mm_comm);
      printf("MM master: after send coords to MM\n");

    } else {

      /* copy coordinate data into communication buffer */
      ndata = 0;
      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          buf[ndata].tag = tag[i];
          buf[ndata].x   = x[i][0];
          buf[ndata].y   = x[i][1];
          buf[ndata].z   = x[i][2];
          ++ndata;
        }
      }

      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, ndata*size_one, MPI_BYTE, 0, 0, world);
    }
  } else if (qmmm_role == QMMM_ROLE_SLAVE) {
    MPI_Recv(qm_coord, num_qm, MPI_DOUBLE, 0, 0, mm_comm, MPI_STATUS_IGNORE);
    printf("MM slave: after received coords\n");
    // XXX: apply coordinates
  }
  return;
}



/* ---------------------------------------------------------------------- */
void FixQMMM::setup(int)
{

  // XXX: verify that the size of the groups has not changed
  // num_mm = group->count(mm_group);
  // num_qm = group->count(igroup);

  exchange_positions();

#if 0

  /* nme:    number of atoms in group on this MPI task
   * nmax:   max number of atoms in group across all MPI tasks
   * nlocal: all local atoms
   */
  int i,j;
  int nmax,nme,nlocal;
  int *mask  = atom->mask;
  int *tag  = atom->tag;
  nlocal = atom->nlocal;
  nme=0;
  for (i=0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nme;

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  memory->destroy(comm_buf);
  maxbuf = nmax*size_one;
  comm_buf = (void *) memory->smalloc(maxbuf,"imd:comm_buf");

  connect_msg = 1;
  reconnect();
  MPI_Bcast(&imd_inactive, 1, MPI_INT, 0, world);
  MPI_Bcast(&imd_terminate, 1, MPI_INT, 0, world);
  if (imd_terminate)
    error->all(FLERR,"LAMMPS terminated on error in setting up IMD connection.");

  /* initialize and build hashtable. */
  inthash_t *qm_hash=new inthash_t;
  inthash_init(qm_hash, num_coords);
  qm_idmap = (void *)qm_hash;

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;
  struct commdata *buf = static_cast<struct commdata *>(comm_buf);

  if (me == 0) {
    int *taglist = new int[num_coords];
    int numtag=0; /* counter to map atom tags to a 0-based consecutive index list */

    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        taglist[numtag] = tag[i];
        ++numtag;
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (j=0; j < ndata; ++j) {
        taglist[numtag] = buf[j].tag;
        ++numtag;
      }
    }

    /* sort list of tags by value to have consistently the
     * same list when running in parallel and build hash table. */
    id_sort(taglist, 0, num_coords-1);
    for (i=0; i < num_coords; ++i) {
      inthash_insert(qm_hash, taglist[i], i);
    }
    delete[] taglist;

    /* generate reverse index-to-tag map for communicating
     * IMD forces back to the proper atoms */
    qm_remap=inthash_keys(qm_hash);
  } else {
    nme=0;
    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        buf[nme].tag = tag[i];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  return;
#endif
}

#if 0

/* ---------------------------------------------------------------------- */
/* Main IMD protocol handler:
 * Send coodinates, energies, and add IMD forces to atoms. */
void FixQMMM::post_force(int vflag)
{
  /* check for reconnect */
  if (imd_inactive) {
    reconnect();
    MPI_Bcast(&imd_inactive, 1, MPI_INT, 0, world);
    MPI_Bcast(&imd_terminate, 1, MPI_INT, 0, world);
    if (imd_terminate)
      error->all(FLERR,"LAMMPS terminated on error in setting up IMD connection.");
    if (imd_inactive)
      return;     /* IMD client has detached and not yet come back. do nothing. */
  }

  int *tag = atom->tag;
  double **x = atom->x;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  int *mask  = atom->mask;
  struct commdata *buf;

  if (me == 0) {
    /* process all pending incoming data. */
    int imd_paused=0;
    while ((imdsock_selread(clientsock, 0) > 0) || imd_paused) {
      /* if something requested to turn off IMD while paused get out */
      if (imd_inactive) break;

      int32 length;
      int msg = imd_recv_header(clientsock, &length);

      switch(msg) {

        case IMD_GO:
          if (screen)
            fprintf(screen, "Ignoring unexpected IMD_GO message.\n");
          break;

        case IMD_IOERROR:
          if (screen)
            fprintf(screen, "IMD connection lost.\n");
          /* fallthrough */

        case IMD_DISCONNECT: {
          /* disconnect from client. wait for new connection. */
          imd_paused = 0;
          imd_forces = 0;
          memory->destroy(force_buf);
          force_buf = NULL;
          imdsock_destroy(clientsock);
          clientsock = NULL;
          if (screen)
            fprintf(screen, "IMD client detached. LAMMPS run continues.\n");

          connect_msg = 1;
          reconnect();
          if (imd_terminate) imd_inactive = 1;
          break;
        }

        case IMD_KILL:
          /* stop the simulation job and shutdown IMD */
          if (screen)
            fprintf(screen, "IMD client requested termination of run.\n");
          imd_inactive = 1;
          imd_terminate = 1;
          imd_paused = 0;
          imdsock_destroy(clientsock);
          clientsock = NULL;
          break;

        case IMD_PAUSE:
          /* pause the running simulation. wait for second IMD_PAUSE to continue. */
          if (imd_paused) {
            if (screen)
              fprintf(screen, "Continuing run on IMD client request.\n");
            imd_paused = 0;
          } else {
            if (screen)
              fprintf(screen, "Pausing run on IMD client request.\n");
            imd_paused = 1;
          }
          break;

        case IMD_TRATE:
          /* change the IMD transmission data rate */
          if (length > 0)
            imd_trate = length;
          if (screen)
            fprintf(screen, "IMD client requested change of transfer rate. Now it is %d.\n", imd_trate);
          break;

        case IMD_ENERGIES: {
          IMDEnergies dummy_energies;
          imd_recv_energies(clientsock, &dummy_energies);
          break;
        }

        case IMD_FCOORDS: {
          float *dummy_coords = new float[3*length];
          imd_recv_fcoords(clientsock, length, dummy_coords);
          delete[] dummy_coords;
          break;
        }

        case IMD_MDCOMM: {
          int32 *imd_tags = new int32[length];
          float *imd_fdat = new float[3*length];
          imd_recv_mdcomm(clientsock, length, imd_tags, imd_fdat);

          if (imd_forces < length) { /* grow holding space for forces, if needed. */
            memory->destroy(force_buf);
            force_buf = (void *) memory->smalloc(length*size_one,
                                                 "imd:force_buf");
          }
          imd_forces = length;
          buf = static_cast<struct commdata *>(force_buf);

          /* compare data to hash table */
          for (int ii=0; ii < length; ++ii) {
            buf[ii].tag = qm_remap[imd_tags[ii]];
            buf[ii].x   = imd_fdat[3*ii];
            buf[ii].y   = imd_fdat[3*ii+1];
            buf[ii].z   = imd_fdat[3*ii+2];
          }
          delete[] imd_tags;
          delete[] imd_fdat;
          break;
        }

        default:
          if (screen)
            fprintf(screen, "Unhandled incoming IMD message #%d. length=%d\n", msg, length);
          break;
      }
    }
  }

  /* update all tasks with current settings. */
  int old_imd_forces = imd_forces;
  MPI_Bcast(&imd_trate, 1, MPI_INT, 0, world);
  MPI_Bcast(&imd_inactive, 1, MPI_INT, 0, world);
  MPI_Bcast(&imd_forces, 1, MPI_INT, 0, world);
  MPI_Bcast(&imd_terminate, 1, MPI_INT, 0, world);
  if (imd_terminate)
    error->all(FLERR,"LAMMPS terminated on IMD request.");

  if (imd_forces > 0) {
    /* check if we need to readjust the forces comm buffer on the receiving nodes. */
    if (me != 0) {
      if (old_imd_forces < imd_forces) { /* grow holding space for forces, if needed. */
        if (force_buf != NULL)
          memory->sfree(force_buf);
        force_buf = memory->smalloc(imd_forces*size_one, "imd:force_buf");
      }
    }
    MPI_Bcast(force_buf, imd_forces*size_one, MPI_BYTE, 0, world);
  }

  /* Check if we need to communicate coordinates to the client.
   * Tuning imd_trate allows to keep the overhead for IMD low
   * at the expense of a more jumpy display. Rather than using
   * end_of_step() we do everything here in one go.
   *
   * If we don't communicate, only check if we have forces
   * stored away and apply them. */
  if (update->ntimestep % imd_trate) {
    if (imd_forces > 0) {
      double **f = atom->f;
      buf = static_cast<struct commdata *>(force_buf);

      /* XXX. this is in principle O(N**2) == not good.
       * however we assume for now that the number of atoms
       * that we manipulate via IMD will be small compared
       * to the total system size, so we don't hurt too much. */
      for (int j=0; j < imd_forces; ++j) {
        for (int i=0; i < nlocal; ++i) {
          if (mask[i] & groupbit) {
            if (buf[j].tag == tag[i]) {
              f[i][0] += imd_fscale*buf[j].x;
              f[i][1] += imd_fscale*buf[j].y;
              f[i][2] += imd_fscale*buf[j].z;
            }
          }
        }
      }
    }
    return;
  }

  /* check and potentially grow local communication buffers. */
  int i, k, nmax, nme=0;
  for (i=0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nme;

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  if (nmax*size_one > maxbuf) {
    memory->destroy(comm_buf);
    maxbuf = nmax*size_one;
    comm_buf = memory->smalloc(maxbuf,"imd:comm_buf");
  }

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;
  buf = static_cast<struct commdata *>(comm_buf);

  if (me == 0) {
    /* collect data into new array. we bypass the IMD API to save
     * us one extra copy of the data. */
    msglen = 3*sizeof(float)*num_coords+IMDHEADERSIZE;
    msgdata = new char[msglen];
    imd_fill_header((IMDheader *)msgdata, IMD_FCOORDS, num_coords);
    /* array pointer, to the offset where we receive the coordinates. */
    float *recvcoord = (float *) (msgdata+IMDHEADERSIZE);

    /* add local data */
    if (unwrap_flag) {
      double xprd = domain->xprd;
      double yprd = domain->yprd;
      double zprd = domain->zprd;
      double xy = domain->xy;
      double xz = domain->xz;
      double yz = domain->yz;

      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, tag[i]);
          if (j != HASH_FAIL) {
            int ix = (image[i] & IMGMASK) - IMGMAX;
            int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
            int iz = (image[i] >> IMG2BITS) - IMGMAX;

            if (domain->triclinic) {
              recvcoord[j]   = x[i][0] + ix * xprd + iy * xy + iz * xz;
              recvcoord[j+1] = x[i][1] + iy * yprd + iz * yz;
              recvcoord[j+2] = x[i][2] + iz * zprd;
            } else {
              recvcoord[j]   = x[i][0] + ix * xprd;
              recvcoord[j+1] = x[i][1] + iy * yprd;
              recvcoord[j+2] = x[i][2] + iz * zprd;
            }
          }
        }
      }
    } else {
      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          const int j = 3*inthash_lookup((inthash_t *)qm_idmap, tag[i]);
          if (j != HASH_FAIL) {
            recvcoord[j]   = x[i][0];
            recvcoord[j+1] = x[i][1];
            recvcoord[j+2] = x[i][2];
          }
        }
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (k=0; k<ndata; ++k) {
        const int j = 3*inthash_lookup((inthash_t *)qm_idmap, buf[k].tag);
        if (j != HASH_FAIL) {
          recvcoord[j]   = buf[k].x;
          recvcoord[j+1] = buf[k].y;
          recvcoord[j+2] = buf[k].z;
        }
      }
    }

    /* done collecting frame data now communicate with IMD client. */

    /* send coordinate data, if client is able to accept */
    if (clientsock && imdsock_selwrite(clientsock,0)) {
      imd_writen(clientsock, msgdata, msglen);
    }
    delete[] msgdata;

  } else {
    /* copy coordinate data into communication buffer */
    nme = 0;
    if (unwrap_flag) {
      double xprd = domain->xprd;
      double yprd = domain->yprd;
      double zprd = domain->zprd;
      double xy = domain->xy;
      double xz = domain->xz;
      double yz = domain->yz;

      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          int ix = (image[i] & IMGMASK) - IMGMAX;
          int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          int iz = (image[i] >> IMG2BITS) - IMGMAX;

          if (domain->triclinic) {
            buf[nme].tag = tag[i];
            buf[nme].x   = x[i][0] + ix * xprd + iy * xy + iz * xz;
            buf[nme].y   = x[i][1] + iy * yprd + iz * yz;
            buf[nme].z   = x[i][2] + iz * zprd;
          } else {
            buf[nme].tag = tag[i];
            buf[nme].x   = x[i][0] + ix * xprd;
            buf[nme].y   = x[i][1] + iy * yprd;
            buf[nme].z   = x[i][2] + iz * zprd;
          }
          ++nme;
        }
      }
    } else {
      for (i=0; i<nlocal; ++i) {
        if (mask[i] & groupbit) {
          buf[nme].tag = tag[i];
          buf[nme].x   = x[i][0];
          buf[nme].y   = x[i][1];
          buf[nme].z   = x[i][2];
          ++nme;
        }
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  return;
}
#endif

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixQMMM::memory_usage(void)
{
  double bytes;
  
  bytes = sizeof(FixQMMM);
  
  return bytes;
}


// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
