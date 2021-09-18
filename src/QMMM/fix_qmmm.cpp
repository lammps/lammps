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
   Contributing author:  Axel Kohlmeyer (ICTP)
------------------------------------------------------------------------- */

#include "fix_qmmm.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "update.h"

#include <cstring>

#include "libqmmm.h"

// message tags for QM/MM inter communicator communication
// have to match with those from the QM code
enum {QMMM_TAG_OTHER=0, QMMM_TAG_SIZE=1, QMMM_TAG_COORD=2,
      QMMM_TAG_FORCE=3, QMMM_TAG_FORCE2=4, QMMM_TAG_CELL=5,
      QMMM_TAG_RADII=6, QMMM_TAG_CHARGE=7, QMMM_TAG_TYPE=8,
      QMMM_TAG_MASS=9};

using namespace LAMMPS_NS;
using namespace FixConst;

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
    while (old_hash) {
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
  for (node=tptr->bucket[h]; node!=nullptr; node=node->next) {
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
    for (node=tptr->bucket[i]; node != nullptr; node=node->next) {
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
    while (node != nullptr) {
      last = node;
      node = node->next;
      free(last);
    }
  }

  /* free the entire array of buckets */
  if (tptr->bucket != nullptr) {
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
  qm_idmap = mm_idmap = nullptr;
  qm_remap = mm_remap = nullptr;
  qm_coord = mm_coord = qm_force = mm_force = nullptr;
  qm_charge =nullptr;
  mm_type = nullptr;

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
  int * const mask = atom->mask;
  int * const type = atom->type;
  double * const mass = atom->mass;
  tagint * const tag  = atom->tag;
  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
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
  if (comm->me == 0) {
    // Pack various cell dimension properties into one chunk.
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
      int isend_buf[4];
      isend_buf[0] = natoms;
      isend_buf[1] = num_qm;
      isend_buf[2] = num_mm;
      isend_buf[3] = ntypes;
      MPI_Send(isend_buf, 4, MPI_INT,    1, QMMM_TAG_SIZE, qm_comm);
      MPI_Send(celldata,  9, MPI_DOUBLE, 1, QMMM_TAG_CELL, qm_comm);
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
      MPI_Send(type, natoms, MPI_INT, 1, QMMM_TAG_TYPE, qm_comm);
      MPI_Send(mass, ntypes+1, MPI_DOUBLE, 1, QMMM_TAG_MASS, qm_comm);

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
  if (utils::strmatch(update->integrate_style,"^respa"))
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
  bytes += (double)maxbuf;
  bytes += (double)6*num_qm*sizeof(double);

  return bytes;
}


// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
