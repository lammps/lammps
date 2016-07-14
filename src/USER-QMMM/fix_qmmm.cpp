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
#include "domain.h"
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
enum {QMMM_TAG_OTHER=0, QMMM_TAG_SIZE=1, QMMM_TAG_COORD=2,
      QMMM_TAG_FORCE=3, QMMM_TAG_FORCE2=4, QMMM_TAG_CELL=5,
      QMMM_TAG_RADII=6, QMMM_TAG_CHARGE=7};

using namespace LAMMPS_NS;
using namespace FixConst;

// prototypes for local helper functions

static int match_element(double mass, int search_isotopes, double &delta);

/****************************************************************************/

/* re-usable integer hash table code with static linkage. */

/** hash table top level data structure */
typedef struct taginthash_t {
  struct taginthash_node_t **bucket;     /* array of hash nodes */
  tagint size;                           /* size of the array */
  tagint entries;                        /* number of entries in table */
  tagint downshift;                      /* shift cound, used in hash function */
  tagint mask;                           /* used to select bits for hashing */
} taginthash_t;

/** hash table node data structure */
typedef struct taginthash_node_t {
  tagint data;                           /* data in hash node */
  tagint key;                            /* key for hash lookup */
  struct taginthash_node_t *next;        /* next node in hash chain */
} taginthash_node_t;

#define HASH_FAIL  -1
#define HASH_LIMIT  0.5

/* initialize new hash table  */
static void taginthash_init(taginthash_t *tptr, tagint buckets);
/* lookup entry in hash table */
static tagint taginthash_lookup(const taginthash_t *tptr, tagint key);
/* generate list of keys for reverse lookups. */
static tagint *taginthash_keys(taginthash_t *tptr);
/* insert an entry into hash table. */
static tagint taginthash_insert(taginthash_t *tptr, tagint key, tagint data);
/* delete the hash table */
static void taginthash_destroy(taginthash_t *tptr);
/* adapted sort for in-place sorting of map indices. */
static void id_sort(tagint *idmap, tagint left, tagint right);

/************************************************************************
 * integer hash code:
 ************************************************************************/

/* taginthash() - Hash function returns a hash number for a given key.
 * tptr: Pointer to a hash table, key: The key to create a hash number for */
static tagint taginthash(const taginthash_t *tptr, tagint key) {
  tagint hashvalue;

  hashvalue = (((key*1103515249)>>tptr->downshift) & tptr->mask);
  if (hashvalue < 0) {
    hashvalue = 0;
  }

  return hashvalue;
}

/*
 *  rebuild_table_tagint() - Create new hash table when old one fills up.
 *
 *  tptr: Pointer to a hash table
 */
static void rebuild_table_tagint(taginthash_t *tptr) {
  taginthash_node_t **old_bucket, *old_hash, *tmp;
  tagint old_size, h, i;

  old_bucket=tptr->bucket;
  old_size=tptr->size;

  /* create a new table and rehash old buckets */
  taginthash_init(tptr, old_size<<1);
  for (i=0; i<old_size; i++) {
    old_hash=old_bucket[i];
    while(old_hash) {
      tmp=old_hash;
      old_hash=old_hash->next;
      h=taginthash(tptr, tmp->key);
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
 *  taginthash_init() - Initialize a new hash table.
 *
 *  tptr: Pointer to the hash table to initialize
 *  buckets: The number of initial buckets to create
 */
void taginthash_init(taginthash_t *tptr, tagint buckets) {

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
  tptr->bucket=(taginthash_node_t **) calloc(tptr->size, sizeof(taginthash_node_t *));

  return;
}

/*
 *  taginthash_lookup() - Lookup an entry in the hash table and return a pointer to
 *    it or HASH_FAIL if it wasn't found.
 *
 *  tptr: Pointer to the hash table
 *  key: The key to lookup
 */
tagint taginthash_lookup(const taginthash_t *tptr, tagint key) {
  tagint h;
  taginthash_node_t *node;


  /* find the entry in the hash table */
  h=taginthash(tptr, key);
  for (node=tptr->bucket[h]; node!=NULL; node=node->next) {
    if (node->key == key)
      break;
  }

  /* return the entry if it exists, or HASH_FAIL */
  return(node ? node->data : HASH_FAIL);
}


/*
 *  taginthash_keys() - Return a list of keys.
 *  NOTE: the returned list must be freed with free(3).
 */
tagint *taginthash_keys(taginthash_t *tptr) {

  tagint *keys;
  taginthash_node_t *node;

  keys = (tagint *)calloc(tptr->entries, sizeof(tagint));

  for (tagint i=0; i < tptr->size; ++i) {
    for (node=tptr->bucket[i]; node != NULL; node=node->next) {
      keys[node->data] = node->key;
    }
  }

  return keys;
}

/*
 *  taginthash_insert() - Insert an entry into the hash table.  If the entry already
 *  exists return a pointer to it, otherwise return HASH_FAIL.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to insert into the hash table
 *  data: A pointer to the data to insert into the hash table
 */
tagint taginthash_insert(taginthash_t *tptr, tagint key, tagint data) {
  tagint tmp;
  taginthash_node_t *node;
  tagint h;

  /* check to see if the entry exists */
  if ((tmp=taginthash_lookup(tptr, key)) != HASH_FAIL)
    return(tmp);

  /* expand the table if needed */
  while (tptr->entries>=HASH_LIMIT*tptr->size)
    rebuild_table_tagint(tptr);

  /* insert the new entry */
  h=taginthash(tptr, key);
  node=(struct taginthash_node_t *) malloc(sizeof(taginthash_node_t));
  node->data=data;
  node->key=key;
  node->next=tptr->bucket[h];
  tptr->bucket[h]=node;
  tptr->entries++;

  return HASH_FAIL;
}

/*
 * taginthash_destroy() - Delete the entire table, and all remaining entries.
 *
 */
void taginthash_destroy(taginthash_t *tptr) {
  taginthash_node_t *node, *last;
  tagint i;

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
    memset(tptr, 0, sizeof(taginthash_t));
  }
}

/************************************************************************
 * integer list sort code:
 ************************************************************************/

/* sort for integer map. initial call  id_sort(idmap, 0, natoms - 1); */
static void id_sort(tagint *idmap, tagint left, tagint right)
{
  tagint pivot, l_hold, r_hold;

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

/****************************************************************************/

/* struct for packed data communication of coordinates and forces. */
struct commdata {
  tagint tag;
  float x,y,z;
  float q;
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

  if (atom->tag_enable == 0)
    error->all(FLERR,"Fix qmmm requires atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Fix qmmm requires consecutive atom IDs");

  /* define ec coupling group */
  mm_group = group->find("all");
  if ((narg == 5) && (0 == strcmp(arg[0],"couple"))) {
#if 0  /* ignore ES coupling group for now */
    mm_group = group->find(arg[4]);
#endif
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
  qm_charge =NULL;
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
    taginthash_t *hashtable = (taginthash_t *)qm_idmap;
    taginthash_destroy(hashtable);
    delete hashtable;
    free(qm_remap);
  }

  memory->destroy(comm_buf);
  memory->destroy(qm_coord);
  memory->destroy(qm_force);
  memory->destroy(qm_charge);
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
  double *charge = atom->q;
  const int * const mask = atom->mask;
  const int * const type = atom->type;
  const tagint * const tag  = atom->tag;
  const int nlocal = atom->nlocal;
  const int natoms = (int) atom->natoms;

  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Too many MM atoms for fix qmmmm");

/*
 *  CHECK
 *  AK: this is not good. atom->natoms can be huge. instead:
 - collect list of atom tags (from local and remoted full MM procs) that
 match fix qmmm group
 - sort list -> master list of QM/MM interacting atoms
 - create hash table that connects master list index with tag
 - collect necessary data and get master list index via hash table
*/
  if (comm->me == 0) {  // AK: make aradii a per-atom property managed by
                        // fix qmmm so it migrates with the local atoms.
    double *aradii = (double *) calloc(natoms, sizeof(double));
    ec_fill_radii( aradii, natoms );

    // This array will pack all details about the cell
    double celldata[9];
    celldata[0] = domain->boxlo[0];
    celldata[1] = domain->boxlo[1];
    celldata[2] = domain->boxlo[2];
    celldata[3] = domain->boxhi[0];
    celldata[4] = domain->boxhi[1];
    celldata[5] = domain->boxhi[2];
    celldata[6] = domain->xy;
    celldata[7] = domain->xz;
    celldata[8] = domain->yz;

    if (qmmm_role == QMMM_ROLE_MASTER) {
      int isend_buf[3];
      isend_buf[0] = natoms;
      isend_buf[1] = num_qm;
      isend_buf[2] =num_mm;
      MPI_Send( isend_buf, 3,   MPI_INTEGER,1, QMMM_TAG_SIZE,  qm_comm );
      MPI_Send( celldata,  9,   MPI_DOUBLE, 1, QMMM_TAG_CELL,  qm_comm );
      MPI_Send( aradii, natoms, MPI_DOUBLE, 1, QMMM_TAG_RADII, qm_comm );
      free(aradii);
    }
    if (verbose > 0) {
      if (screen) fputs("QMMM: exchange positions\n",screen);
      if (logfile) fputs("QMMM: exchange positions\n",logfile);
    }
  }

  if (qmmm_role == QMMM_ROLE_MASTER) {
    int i;
    double *mm_coord_all = (double *) calloc(3*natoms, sizeof(double));
    double *mm_charge_all = (double *) calloc(natoms, sizeof(double));
    int *mm_mask_all = (int *) calloc(natoms, sizeof(int));

    /* check and potentially grow local communication buffers. */
    if (atom->nmax*size_one > maxbuf) {
      memory->destroy(comm_buf);
      maxbuf = atom->nmax*size_one;
      comm_buf = memory->smalloc(maxbuf,"qmmm:comm_buf");
    }

    if (comm->me == 0) {

      // insert local atoms into comm buffer and global arrays
      for (i=0; i<nlocal; ++i) {
        const int k = (int) tag[i]-1;
        mm_mask_all[ k ] = -1;
        if (mask[i] & groupbit) {
          const int j = 3*taginthash_lookup((taginthash_t *)qm_idmap, tag[i]);
          if (j != 3*HASH_FAIL) {
            qm_coord[j]   = x[i][0];
            qm_coord[j+1] = x[i][1];
            qm_coord[j+2] = x[i][2];
            qm_charge[j/3] = charge[i];
            mm_mask_all[k] = type[i];
          }
        }
        mm_coord_all[3*k + 0] = x[i][0];
        mm_coord_all[3*k + 1] = x[i][1];
        mm_coord_all[3*k + 2] = x[i][2];
        mm_charge_all[k] = charge[i];
      }

      /* done collecting coordinates, send it to dependent codes */
      /* to QM code */
      MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, qm_comm);
      MPI_Send(qm_charge, num_qm, MPI_DOUBLE, 1, QMMM_TAG_CHARGE, qm_comm);
      MPI_Send(mm_charge_all, natoms, MPI_DOUBLE, 1, QMMM_TAG_COORD, qm_comm);
      MPI_Send(mm_coord_all, 3*natoms, MPI_DOUBLE, 1, QMMM_TAG_COORD, qm_comm);
      MPI_Send(mm_mask_all, natoms, MPI_INT, 1, QMMM_TAG_COORD, qm_comm);

      /* to MM slave code */
      MPI_Send(qm_coord, 3*num_qm, MPI_DOUBLE, 1, QMMM_TAG_COORD, mm_comm);

      free(mm_coord_all);
      free(mm_charge_all);
      free(mm_mask_all);

    } else {
      error->one(FLERR,"Cannot handle parallel MM (yet)");
    }
  } else if (qmmm_role == QMMM_ROLE_SLAVE) {

    MPI_Recv(qm_coord, 3*num_qm, MPI_DOUBLE, 0, QMMM_TAG_COORD, mm_comm, MPI_STATUS_IGNORE);
    // not needed for as long as we allow only one MPI task as slave
    MPI_Bcast(qm_coord, 3*num_qm, MPI_DOUBLE,0,world);

    /* update coordinates of (QM) atoms */
    for (int i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        for (int j=0; j < num_qm; ++j) {
          if (tag[i] == qm_remap[j]) {
            x[i][0] = qm_coord[3*j];
            x[i][1] = qm_coord[3*j+1];
            x[i][2] = qm_coord[3*j+2];
          }
        }
      }
    }
  }

  return;
}

/* ---------------------------------------------------------------------- */

void FixQMMM::exchange_forces()
{
  double **f = atom->f;
  const int * const mask  = atom->mask;
  const tagint * const tag  = atom->tag;
  const int nlocal = atom->nlocal;
  const int natoms = (int) atom->natoms;

  if ((comm->me) == 0 && (verbose > 0)) {
    if (screen)  fputs("QMMM: exchange forces\n",screen);
    if (logfile) fputs("QMMM: exchange forces\n",logfile);
  }

  if (qmmm_role == QMMM_ROLE_MASTER) {
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);

    double *mm_force_all = (double *) calloc(natoms*3, sizeof(double));
    double *mm_force_on_qm_atoms = qm_coord; // use qm_coord as a buffer

    if (comm->me == 0) {
      // receive QM forces from QE
      MPI_Recv(qm_force,3*num_qm,MPI_DOUBLE,1,QMMM_TAG_FORCE,qm_comm,MPI_STATUS_IGNORE);
      // receive ec contribution to MM forces from QE
      MPI_Recv(mm_force_all,3*natoms,MPI_DOUBLE,1,QMMM_TAG_FORCE2,qm_comm,MPI_STATUS_IGNORE);
      // receive MM forces from LAMMPS
      MPI_Recv( mm_force_on_qm_atoms, 3*num_qm,MPI_DOUBLE,1,QMMM_TAG_FORCE,mm_comm,MPI_STATUS_IGNORE);

      // subtract MM forces from QM forces to get the delta
      // NOTE: QM forces are always sent in "real" units,
      // so we need to apply the scaling factor to get to the
      // supported internal units ("metal" or "real")
      for (int i=0; i < num_qm; ++i) {
        if  (verbose > 1) {
           const char fmt[] = "[" TAGINT_FORMAT "]: QM(%g %g %g) MM(%g %g %g) /\\(%g %g %g)\n";
           if (screen) fprintf(screen, fmt, qm_remap[i],
                qmmm_fscale*qm_force[3*i+0], qmmm_fscale*qm_force[3*i+1], qmmm_fscale*qm_force[3*i+2],
                mm_force_on_qm_atoms[3*i+0], mm_force_on_qm_atoms[3*i+1], mm_force_on_qm_atoms[3*i+2],
                qmmm_fscale*qm_force[3*i+0] - mm_force_on_qm_atoms[3*i+0],
                qmmm_fscale*qm_force[3*i+1] - mm_force_on_qm_atoms[3*i+1],
                qmmm_fscale*qm_force[3*i+2] - mm_force_on_qm_atoms[3*i+2]);
           if (logfile) fprintf(logfile, fmt, qm_remap[i],
                qmmm_fscale*qm_force[3*i+0], qmmm_fscale*qm_force[3*i+1], qmmm_fscale*qm_force[3*i+2],
                mm_force_on_qm_atoms[3*i+0], mm_force_on_qm_atoms[3*i+1], mm_force_on_qm_atoms[3*i+2],
                qmmm_fscale*qm_force[3*i+0] - mm_force_on_qm_atoms[3*i+0],
                qmmm_fscale*qm_force[3*i+1] - mm_force_on_qm_atoms[3*i+1],
                qmmm_fscale*qm_force[3*i+2] - mm_force_on_qm_atoms[3*i+2]);
        }
        buf[i].tag = qm_remap[i];
        buf[i].x = qmmm_fscale*qm_force[3*i+0] - mm_force_on_qm_atoms[3*i+0];
        buf[i].y = qmmm_fscale*qm_force[3*i+1] - mm_force_on_qm_atoms[3*i+1];
        buf[i].z = qmmm_fscale*qm_force[3*i+2] - mm_force_on_qm_atoms[3*i+2];
      }
    }
    MPI_Bcast(comm_buf,num_qm*size_one,MPI_BYTE,0,world);

    // Inefficient... use buffers.
    MPI_Bcast(mm_force_all,natoms*3,MPI_DOUBLE,0,world);

    /* apply forces resulting from QM/MM coupling */
    if (qmmm_mode == QMMM_MODE_MECH) {
      for (int i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          for (int j=0; j < num_qm; ++j) {
            if (tag[i] == buf[j].tag) {
              f[i][0] += buf[j].x;
              f[i][1] += buf[j].y;
              f[i][2] += buf[j].z;
            }
          }
        }
      }
    } else if (qmmm_mode == QMMM_MODE_ELEC) {
      for (int i=0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          for (int j=0; j < num_qm; ++j) {
            if (tag[i] == buf[j].tag) {
              f[i][0] += buf[j].x;
              f[i][1] += buf[j].y;
              f[i][2] += buf[j].z;
            }
          }
        } else {
          const int k = (int) tag[i]-1;
          f[i][0] += qmmm_fscale * mm_force_all[ 3*k + 0 ];
          f[i][1] += qmmm_fscale * mm_force_all[ 3*k + 1 ];
          f[i][2] += qmmm_fscale * mm_force_all[ 3*k + 2 ];
        }
      }
    }

    free(mm_force_all);

  } else if (qmmm_role == QMMM_ROLE_SLAVE) {

    // use qm_force and qm_coord as communication buffer
    double * mm_force_on_qm_atoms = qm_force;
    double * reduced_mm_force_on_qm_atoms = qm_coord;

    memset( mm_force_on_qm_atoms, 0, 3*num_qm*sizeof(double) );
    for (int i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        const int j = 3*taginthash_lookup((taginthash_t *)qm_idmap, tag[i]);
        if (j != 3*HASH_FAIL) {
          mm_force_on_qm_atoms[j]   = f[i][0];
          mm_force_on_qm_atoms[j+1] = f[i][1];
          mm_force_on_qm_atoms[j+2] = f[i][2];
        }
      }
    }

    // collect and send MM slave forces to MM master
    // the reduction is not really needed with only one rank (for now)
    MPI_Reduce(mm_force_on_qm_atoms, reduced_mm_force_on_qm_atoms, 3*num_qm, MPI_DOUBLE, MPI_SUM, 0, world);
    // use qm_coord array as a communication buffer
    MPI_Send(reduced_mm_force_on_qm_atoms, 3*num_qm, MPI_DOUBLE, 0, QMMM_TAG_FORCE, mm_comm);
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
      memory->create(qm_charge,num_qm,"qmmm:qm_charge");
      memory->create(qm_force,3*num_qm,"qmmm:qm_force");

      const char fmt1[] = "Initializing QM/MM master with %d QM atoms\n";
      const char fmt2[] = "Initializing QM/MM master with %d MM atoms\n";
      const char fmt3[] = "Electrostatic coupling with %d atoms\n";

      if (screen)  {
        fprintf(screen,fmt1,num_qm);
        fprintf(screen,fmt2,num_mm);
        if (qmmm_mode == QMMM_MODE_ELEC) fprintf(screen,fmt3,num_mm-num_qm);
      }

      if (logfile) {
        fprintf(logfile,fmt1,num_qm);
        fprintf(logfile,fmt2,num_mm);
        if (qmmm_mode == QMMM_MODE_ELEC) fprintf(logfile,fmt3,num_mm-num_qm);
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
    taginthash_t *qm_hash=new taginthash_t;
    taginthash_init(qm_hash, num_qm);
    qm_idmap = (void *)qm_hash;

    const int nlocal = atom->nlocal;
    int i, j, tmp, ndata, qm_ntag;
    tagint *tag = atom->tag;
    int *mask  = atom->mask;
    struct commdata *buf = static_cast<struct commdata *>(comm_buf);

    if (me == 0) {
      MPI_Status status;
      MPI_Request request;
      tagint *qm_taglist = new tagint[num_qm];
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
        taginthash_insert(qm_hash, qm_taglist[i], i);
      }
      delete[] qm_taglist;

      /* generate reverse index-to-tag map for communicating
       * qm/mm forces back to the proper atoms */
      qm_remap=taginthash_keys(qm_hash);

      if (verbose > 1) {
        const char fmt[] = "qm_remap[%d]=" TAGINT_FORMAT
                           "  qm_hash[" TAGINT_FORMAT "]=" TAGINT_FORMAT "\n";
        // print hashtable and reverse mapping
        for (i=0; i < num_qm; ++i) {
          if (screen) fprintf(screen,fmt,i,qm_remap[i],qm_remap[i],
                              taginthash_lookup(qm_hash, qm_remap[i]));
          if (logfile) fprintf(logfile,fmt,i,qm_remap[i],qm_remap[i],
                               taginthash_lookup(qm_hash, qm_remap[i]));
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


/* ---------------------------------------------------------------------- */
/* Manage the atomic number */

#define QMMM_ISOTOPES 351

static const int FixQMMM_Z[QMMM_ISOTOPES] = {
  1, 1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9,
  10, 10, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16, 16, 16,
  16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 20,
  21, 22, 22, 22, 22, 22, 23, 23, 24, 24, 24, 24, 25, 26, 26,
  26, 26, 27, 28, 28, 28, 28, 28, 29, 29, 30, 30, 30, 30, 30,
  31, 31, 32, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 34, 35,
  35, 36, 36, 36, 36, 36, 36, 37, 37, 38, 38, 38, 38, 39, 40,
  40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43,
  44, 44, 44, 44, 44, 44, 44, 45, 46, 46, 46, 46, 46, 46, 47,
  47, 48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 50, 50, 50, 50,
  50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52, 52, 52, 52, 52,
  52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 56, 56, 56,
  56, 56, 56, 56, 57, 57, 58, 58, 58, 58, 59, 60, 60, 60, 60,
  60, 60, 60, 61, 61, 62, 62, 62, 62, 62, 62, 62, 63, 63, 64,
  64, 64, 64, 64, 64, 64, 65, 66, 66, 66, 66, 66, 66, 66, 67,
  68, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71,
  71, 72, 72, 72, 72, 72, 72, 73, 73, 74, 74, 74, 74, 74, 75,
  75, 76, 76, 76, 76, 76, 76, 76, 77, 77, 78, 78, 78, 78, 78,
  78, 79, 80, 80, 80, 80, 80, 80, 80, 81, 81, 82, 82, 82, 82,
  83, 84, 84, 85, 85, 86, 86, 86, 87, 88, 88, 88, 88, 89, 90,
  90, 91, 92, 92, 92, 92, 92, 93, 93, 94, 94, 94, 94, 94, 94,
  95, 95, 96, 96, 96, 96, 96, 96, 97, 97, 98, 98, 98, 98, 99,
  100, 101, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
  111, 112, 113, 114, 115};

static const char FixQMMM_EL[QMMM_ISOTOPES][4] = {
  "H", "H", "H", "He", "He", "Li", "Li", "Be", "B", "B",
  "C", "C", "C", "N", "N", "O", "O", "O", "F", "Ne", "Ne",
  "Ne", "Na", "Mg", "Mg", "Mg", "Al", "Si", "Si", "Si",
  "P", "S", "S", "S", "S", "Cl", "Cl", "Ar", "Ar", "Ar",
  "K", "K", "K", "Ca", "Ca", "Ca", "Ca", "Ca", "Ca", "Sc",
  "Ti", "Ti", "Ti", "Ti", "Ti", "V", "V", "Cr", "Cr", "Cr",
  "Cr", "Mn", "Fe", "Fe", "Fe", "Fe", "Co", "Ni", "Ni",
  "Ni", "Ni", "Ni", "Cu", "Cu", "Zn", "Zn", "Zn", "Zn", "Zn",
  "Ga", "Ga", "Ge", "Ge", "Ge", "Ge", "Ge", "As", "Se",
  "Se", "Se", "Se", "Se", "Se", "Br", "Br", "Kr", "Kr",
  "Kr", "Kr", "Kr", "Kr", "Rb", "Rb", "Sr", "Sr", "Sr",
  "Sr", "Y", "Zr", "Zr", "Zr", "Zr", "Zr", "Nb", "Mo", "Mo",
  "Mo", "Mo", "Mo", "Mo", "Mo", "Tc", "Tc", "Tc", "Ru", "Ru",
  "Ru", "Ru", "Ru", "Ru", "Ru", "Rh", "Pd", "Pd", "Pd", "Pd",
  "Pd", "Pd", "Ag", "Ag", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd",
  "Cd", "Cd", "In", "In", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn",
  "Sn", "Sn", "Sn", "Sn", "Sb", "Sb", "Te", "Te", "Te", "Te",
  "Te", "Te", "Te", "Te", "I", "Xe", "Xe", "Xe", "Xe", "Xe",
  "Xe", "Xe", "Xe", "Xe", "Cs", "Ba", "Ba", "Ba", "Ba", "Ba",
  "Ba", "Ba", "La", "La", "Ce", "Ce", "Ce", "Ce", "Pr", "Nd",
  "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Pm", "Pm", "Sm", "Sm",
  "Sm", "Sm", "Sm", "Sm", "Sm", "Eu", "Eu", "Gd", "Gd", "Gd",
  "Gd", "Gd", "Gd", "Gd", "Tb", "Dy", "Dy", "Dy", "Dy", "Dy",
  "Dy", "Dy", "Ho", "Er", "Er", "Er", "Er", "Er", "Er", "Tm",
  "Yb", "Yb", "Yb", "Yb", "Yb", "Yb", "Yb", "Lu", "Lu", "Hf",
  "Hf", "Hf", "Hf", "Hf", "Hf", "Ta", "Ta", "W", "W", "W", "W",
  "W", "Re", "Re", "Os", "Os", "Os", "Os", "Os", "Os", "Os",
  "Ir", "Ir", "Pt", "Pt", "Pt", "Pt", "Pt", "Pt", "Au", "Hg",
  "Hg", "Hg", "Hg", "Hg", "Hg", "Hg", "Tl", "Tl", "Pb", "Pb",
  "Pb", "Pb", "Bi", "Po", "Po", "At", "At", "Rn", "Rn", "Rn",
  "Fr", "Ra", "Ra", "Ra", "Ra", "Ac", "Th", "Th", "Pa", "U",
  "U", "U", "U", "U", "Np", "Np", "Pu", "Pu", "Pu", "Pu", "Pu",
  "Pu", "Am", "Am", "Cm", "Cm", "Cm", "Cm", "Cm", "Cm", "Bk",
  "Bk", "Cf", "Cf", "Cf", "Cf", "Es", "Fm", "Md", "Md", "No",
  "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
  "Uut", "Uuq", "Uup"};

static const int FixQMMM_A[QMMM_ISOTOPES] = {
  1, 2, 3, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14, 14, 15, 16, 17, 18,
  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
  36, 35, 37, 36, 38, 40, 39, 40, 41, 40, 42, 43, 44, 46, 48, 45,
  46, 47, 48, 49, 50, 50, 51, 50, 52, 53, 54, 55, 54, 56, 57, 58,
  59, 58, 60, 61, 62, 64, 63, 65, 64, 66, 67, 68, 70, 69, 71, 70,
  72, 73, 74, 76, 75, 74, 76, 77, 78, 80, 82, 79, 81, 78, 80, 82,
  83, 84, 86, 85, 87, 84, 86, 87, 88, 89, 90, 91, 92, 94, 96, 93,
  92, 94, 95, 96, 97, 98, 100, 97, 98, 99, 96, 98, 99, 100, 101,
  102, 104, 103, 102, 104, 105, 106, 108, 110, 107, 109, 106, 108,
  110, 111, 112, 113, 114, 116, 113, 115, 112, 114, 115, 116, 117,
  118, 119, 120, 122, 124, 121, 123, 120, 122, 123, 124, 125, 126,
  128, 130, 127, 124, 126, 128, 129, 130, 131, 132, 134, 136, 133,
  130, 132, 134, 135, 136, 137, 138, 138, 139, 136, 138, 140, 142,
  141, 142, 143, 144, 145, 146, 148, 150, 145, 147, 144, 147, 148,
  149, 150, 152, 154, 151, 153, 152, 154, 155, 156, 157, 158, 160,
  159, 156, 158, 160, 161, 162, 163, 164, 165, 162, 164, 166, 167,
  168, 170, 169, 168, 170, 171, 172, 173, 174, 176, 175, 176, 174,
  176, 177, 178, 179, 180, 180, 181, 180, 182, 183, 184, 186, 185,
  187, 184, 186, 187, 188, 189, 190, 192, 191, 193, 190, 192, 194,
  195, 196, 198, 197, 196, 198, 199, 200, 201, 202, 204, 203, 205,
  204, 206, 207, 208, 209, 209, 210, 210, 211, 211, 220, 222, 223,
  223, 224, 226, 228, 227, 230, 232, 231, 233, 234, 235, 236, 238,
  236, 237, 238, 239, 240, 241, 242, 244, 241, 243, 243, 244, 245,
  246, 247, 248, 247, 249, 249, 250, 251, 252, 252, 257, 258, 260,
  259, 262, 265, 268, 271, 272, 270, 276, 281, 280, 285, 284, 289,
  288};

static const double FixQMMM_MASS[QMMM_ISOTOPES] = {
  1.00782503207, 2.0141017778, 3.0160492777, 3.0160293191, 4.00260325415, 6.015122795, 7.01600455, 9.0121822,
  10.0129370, 11.0093054, 12.0000000, 13.0033548378, 14.003241989, 14.0030740048, 15.0001088982, 15.99491461956,
  16.99913170, 17.9991610, 18.99840322, 19.9924401754, 20.99384668, 21.991385114, 22.9897692809, 23.985041700,
  24.98583692, 25.982592929, 26.98153863, 27.9769265325, 28.976494700, 29.97377017, 30.97376163, 31.97207100,
  32.97145876, 33.96786690, 35.96708076, 34.96885268, 36.96590259, 35.967545106, 37.9627324, 39.9623831225,
  38.96370668, 39.96399848, 40.96182576, 39.96259098, 41.95861801, 42.9587666, 43.9554818, 45.9536926,
  47.952534, 44.9559119, 45.9526316, 46.9517631, 47.9479463, 48.9478700, 49.9447912, 49.9471585,
  50.9439595, 49.9460442, 51.9405075, 52.9406494, 53.9388804, 54.9380451, 53.9396105, 55.9349375,
  56.9353940, 57.9332756, 58.9331950, 57.9353429, 59.9307864, 60.9310560, 61.9283451, 63.9279660,
  62.9295975, 64.9277895, 63.9291422, 65.9260334, 66.9271273, 67.9248442, 69.9253193, 68.9255736,
  70.9247013, 69.9242474, 71.9220758, 72.9234589, 73.9211778, 75.9214026, 74.9215965, 73.9224764,
  75.9192136, 76.9199140, 77.9173091, 79.9165213, 81.9166994, 78.9183371, 80.9162906, 77.9203648,
  79.9163790, 81.9134836, 82.914136, 83.911507, 85.91061073, 84.911789738, 86.909180527, 83.913425,
  85.9092602, 86.9088771, 87.9056121, 88.9058483, 89.9047044, 90.9056458, 91.9050408, 93.9063152,
  95.9082734, 92.9063781, 91.906811, 93.9050883, 94.9058421, 95.9046795, 96.9060215, 97.9054082,
  99.907477, 96.906365, 97.907216, 98.9062547, 95.907598, 97.905287, 98.9059393, 99.9042195,
  100.9055821, 101.9043493, 103.905433, 102.905504, 101.905609, 103.904036, 104.905085, 105.903486,
  107.903892, 109.905153, 106.905097, 108.904752, 105.906459, 107.904184, 109.9030021, 110.9041781,
  111.9027578, 112.9044017, 113.9033585, 115.904756, 112.904058, 114.903878, 111.904818, 113.902779,
  114.903342, 115.901741, 116.902952, 117.901603, 118.903308, 119.9021947, 121.9034390, 123.9052739,
  120.9038157, 122.9042140, 119.904020, 121.9030439, 122.9042700, 123.9028179, 124.9044307, 125.9033117,
  127.9044631, 129.9062244, 126.904473, 123.9058930, 125.904274, 127.9035313, 128.9047794, 129.9035080,
  130.9050824, 131.9041535, 133.9053945, 135.907219, 132.905451933, 129.9063208, 131.9050613, 133.9045084,
  134.9056886, 135.9045759, 136.9058274, 137.9052472, 137.907112, 138.9063533, 135.907172, 137.905991,
  139.9054387, 141.909244, 140.9076528, 141.9077233, 142.9098143, 143.9100873, 144.9125736, 145.9131169,
  147.916893, 149.920891, 144.912749, 146.9151385, 143.911999, 146.9148979, 147.9148227, 148.9171847,
  149.9172755, 151.9197324, 153.9222093, 150.9198502, 152.9212303, 151.9197910, 153.9208656, 154.9226220,
  155.9221227, 156.9239601, 157.9241039, 159.9270541, 158.9253468, 155.924283, 157.924409, 159.9251975,
  160.9269334, 161.9267984, 162.9287312, 163.9291748, 164.9303221, 161.928778, 163.929200, 165.9302931,
  166.9320482, 167.9323702, 169.9354643, 168.9342133, 167.933897, 169.9347618, 170.9363258, 171.9363815,
  172.9382108, 173.9388621, 175.9425717, 174.9407718, 175.9426863, 173.940046, 175.9414086, 176.9432207,
  177.9436988, 178.9458161, 179.9465500, 179.9474648, 180.9479958, 179.946704, 181.9482042, 182.9502230,
  183.9509312, 185.9543641, 184.9529550, 186.9557531, 183.9524891, 185.9538382, 186.9557505, 187.9558382,
  188.9581475, 189.9584470, 191.9614807, 190.9605940, 192.9629264, 189.959932, 191.9610380, 193.9626803,
  194.9647911, 195.9649515, 197.967893, 196.9665687, 195.965833, 197.9667690, 198.9682799, 199.9683260,
  200.9703023, 201.9706430, 203.9734939, 202.9723442, 204.9744275, 203.9730436, 205.9744653, 206.9758969,
  207.9766521, 208.9803987, 208.9824304, 209.9828737, 209.987148, 210.9874963, 210.990601, 220.0113940,
  222.0175777, 223.0197359, 223.0185022, 224.0202118, 226.0254098, 228.0310703, 227.0277521, 230.0331338,
  232.0380553, 231.0358840, 233.0396352, 234.0409521, 235.0439299, 236.0455680, 238.0507882, 236.046570,
  237.0481734, 238.0495599, 239.0521634, 240.0538135, 241.0568515, 242.0587426, 244.064204, 241.0568291,
  243.0613811, 243.0613891, 244.0627526, 245.0654912, 246.0672237, 247.070354, 248.072349, 247.070307,
  249.0749867, 249.0748535, 250.0764061, 251.079587, 252.081626, 252.082980, 257.095105, 258.098431,
  260.10365, 259.10103, 262.10963, 265.11670, 268.12545, 271.13347, 272.13803, 270.13465,
  276.15116, 281.16206, 280.16447, 285.17411, 284.17808, 289.18728, 288.19249};

/*
 * This function matches the element which has the least absolute
 * difference between the masses. delta is set to the difference
 * between the mass provided and the best match.
 *
 * THIS FUNCTION RELIES ON THE CORRECT ORDER OF THE DATA TABLES
 * IN ORDER TO SKIP ISOTOPES FROM THE MATCH
 *
 *                     */
int match_element(double mass, int search_isotopes, double &delta)
{
  int i;
  int bestcandidate;
  double diff, mindiff;
  int lastz;

  mindiff =  1e6;
  bestcandidate = -1;
  lastz = -1;
  for(i=0; i<QMMM_ISOTOPES; i++) {
    if(!search_isotopes && lastz == FixQMMM_Z[i])
      continue;
    diff = fabs(FixQMMM_MASS[i] - mass);
    if(diff < mindiff) {
      mindiff = diff;
      bestcandidate = i;
    }
    lastz = FixQMMM_Z[i];
  }
  delta = mindiff;
  return bestcandidate;
}

/* ---------------------------------------------------------------------- */
//! Atomic radii of all atoms. Master/QM only
//
//

/*
 *     Source: wikipedia
 *         [ http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) ]
 *
 *             Atomic radii for the various elements in picometers.
 *
 *                 A value of -1 has the meaning of "data not available".
 *
 *                 */

/// Number of elements in the arrays
#define EC_ELEMENTS 116

static const double ec_r_covalent[EC_ELEMENTS+1] = {
  -1 /* Z=0 */,
  0.38, 0.32, 1.34, 0.90, 0.82, 0.77, 0.75, 0.73, 0.71, 0.69, 1.54, 1.30, 1.18, 1.11, 1.06, 1.02, 0.99, 0.97, 1.96, 1.74,
  1.44, 1.36, 1.25, 1.27, 1.39, 1.25, 1.26, 1.21, 1.38, 1.31, 1.26, 1.22, 1.19, 1.16, 1.14, 1.10, 2.11, 1.92, 1.62, 1.48,
  1.37, 1.45, 1.56, 1.26, 1.35, 1.31, 1.53, 1.48, 1.44, 1.41, 1.38, 1.35, 1.33, 1.30, 2.25, 1.98, 1.69,   -1,   -1,   -1,
  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, 1.60, 1.50, 1.38, 1.46, 1.59, 1.28, 1.37, 1.28, 1.44, 1.49,
  1.48, 1.47, 1.46,   -1,   -1, 1.45,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1
};

int FixQMMM::ec_fill_radii( double *aradii, int ndata )
{
  int i, el;
  char errmsg[150];
  int *type2element;

  memory->create(type2element,atom->ntypes+1,"qmmm:type2element");
  type2element[0] = 0;

  for (i=1; i <= atom->ntypes; ++i) {
    double delta;
    el = match_element(atom->mass[i], 0, delta);
    sprintf(errmsg,"FixQMMM: type %2d (mass: %8g) matches %2s with:"
            " Z = %-3d A = %-3d r_cov = %5.2f (error = %-8.2g -> %-.2g%%)\n",
            i, atom->mass[i], FixQMMM_EL[el], FixQMMM_Z[el], FixQMMM_A[el],
            ec_r_covalent[FixQMMM_Z[el]], delta, delta/atom->mass[i] * 100.0);
    type2element[i] = FixQMMM_Z[el];
    if(screen) fprintf(screen, "%s", errmsg);
    if(logfile) fprintf(logfile, "%s", errmsg);
  }

  for (i=0; i<ndata; i++) {
    el = type2element[atom->type[i]];
    if (el < 0 || el > EC_ELEMENTS) {
      sprintf(errmsg,"Unable to find element Z=%d in table of covalent radii", el);
      error->one(FLERR,errmsg);
    }
    aradii[i] = ec_r_covalent[el];
    if (ec_r_covalent[el] < 0.0) {
      sprintf(errmsg, "Covalent radius for atom of element Z=%d not availabe", el);
      error->one(FLERR,errmsg);
    }
  }
  memory->destroy(type2element);
  return 0;
}

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
