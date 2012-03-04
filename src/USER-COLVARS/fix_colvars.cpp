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
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"

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

/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;

/* struct for packed data communication of coordinates and forces. */
struct commdata { 
  int tag; 
  float x,y,z; 
};

/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax: 
 * fix ID group-ID colvars <config_file>
 ***************************************************************/
FixColvars::FixColvars(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 2) 
    error->all(FLERR,"Illegal fix colvars command");

#if 0
  imd_port = atoi(arg[3]); 
  if (imd_port < 1024)
    error->all(FLERR,"Illegal fix imd parameter: port < 1024");

  /* default values for optional flags */
  unwrap_flag = 0;
  nowait_flag = 0;
  connect_msg = 1;
  imd_fscale = 1.0;
  imd_trate = 1;
  
  /* parse optional arguments */
  int argsdone = 4;
  while (argsdone+1 < narg) {
    if (0 == strcmp(arg[argsdone], "unwrap")) {
      if (0 == strcmp(arg[argsdone+1], "on")) {  
        unwrap_flag = 1;
      } else {
        unwrap_flag = 0;
      }
    } else if (0 == strcmp(arg[argsdone], "nowait")) {
      if (0 == strcmp(arg[argsdone+1], "on")) {  
        nowait_flag = 1;
      } else {
        nowait_flag = 0;
      }
    } else if (0 == strcmp(arg[argsdone], "fscale")) {
      imd_fscale = atof(arg[argsdone+1]);
    } else if (0 == strcmp(arg[argsdone], "trate")) {
      imd_trate = atoi(arg[argsdone+1]);
    } else {
      error->all(FLERR,"Unknown fix imd parameter");
    }
    ++argsdone; ++argsdone;
  }

  /* sanity check on parameters */
  if (imd_trate < 1)
    error->all(FLERR,"Illegal fix imd parameter. trate < 1.");

  bigint n = group->count(igroup);
  if (n > MAXSMALLINT) error->all(FLERR,"Too many atoms for fix imd");
  num_coords = static_cast<int> (n);

  MPI_Comm_rank(world,&me);

  /* initialize various imd state variables. */
  clientsock = NULL;
  localsock  = NULL;
  nlevels_respa = 0;
  imd_inactive = 0;
  imd_terminate = 0;
  imd_forces = 0;
  force_buf = NULL;
  maxbuf = 0;
  msgdata = NULL;
  msglen = 0;
  comm_buf = NULL;
  idmap = NULL;
  rev_idmap = NULL;
  
  if (me == 0) {
    /* set up incoming socket on MPI rank 0. */
    imdsock_init();
    localsock = imdsock_create();
    clientsock = NULL;
    if (imdsock_bind(localsock,imd_port)) {
      perror("bind to socket failed");
      imdsock_destroy(localsock);
      imd_terminate = 1;
    } else {
      imdsock_listen(localsock);
    }
  }
  MPI_Bcast(&imd_terminate, 1, MPI_INT, 0, world);
  if (imd_terminate)
    error->all(FLERR,"LAMMPS Terminated on error in IMD.");
    
  /* storage required to communicate a single coordinate or force. */
  size_one = sizeof(struct commdata);

#if defined(LAMMPS_ASYNC_IMD)
  /* set up for i/o worker thread on MPI rank 0.*/
  if (me == 0) {
    if (screen)
      fputs("Using fix imd with asynchronous I/O.\n",screen);
    if (logfile)
      fputs("Using fix imd with asynchronous I/O.\n",logfile);

    /* set up mutex and condition variable for i/o thread */
    /* hold mutex before creating i/o thread to keep it waiting. */
    pthread_mutex_init(&read_mutex, NULL);
    pthread_mutex_init(&write_mutex, NULL);
    pthread_cond_init(&write_cond, NULL);

    pthread_mutex_lock(&write_mutex);
    buf_has_data=0;
    pthread_mutex_unlock(&write_mutex);

    /* set up and launch i/o thread */
    pthread_attr_init(&iot_attr);
    pthread_attr_setdetachstate(&iot_attr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&iothread, &iot_attr, &fix_imd_ioworker, this);
  }
#endif

#endif
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixColvars::~FixColvars()
{

#if 0
  inthash_t *hashtable = (inthash_t *)idmap;
  memory->sfree(comm_buf);
  memory->sfree(force_buf);
  inthash_destroy(hashtable);
  delete hashtable;
  free(rev_idmap);
  // close sockets
  imdsock_shutdown(clientsock);
  imdsock_destroy(clientsock);
  imdsock_shutdown(localsock);
  imdsock_destroy(localsock);
  clientsock=NULL;
  localsock=NULL;
  return;
#endif
}

/* ---------------------------------------------------------------------- */
int FixColvars::setmask()
{
  int mask = 0;
#if 0
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
#endif
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
/* wait for IMD client (e.g. VMD) to respond, initialize communication 
 * buffers and collect tag/id maps. */
void FixColvars::setup(int)
{

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
  maxbuf = nmax*size_one;
  comm_buf = (void *) memory->smalloc(maxbuf,"imd:comm_buf");

  connect_msg = 1;
  reconnect();
  MPI_Bcast(&imd_inactive, 1, MPI_INT, 0, world);
  MPI_Bcast(&imd_terminate, 1, MPI_INT, 0, world);
  if (imd_terminate)
    error->all(FLERR,"LAMMPS terminated on error in setting up IMD connection.");

  /* initialize and build hashtable. */
  inthash_t *hashtable=new inthash_t;
  inthash_init(hashtable, num_coords);
  idmap = (void *)hashtable;
  
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
      inthash_insert(hashtable, taglist[i], i);
    }
    delete[] taglist;

    /* generate reverse index-to-tag map for communicating 
     * IMD forces back to the proper atoms */
    rev_idmap=inthash_keys(hashtable);
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
void FixColvars::post_force(int vflag)
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
  int *image = atom->image;
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
            buf[ii].tag = rev_idmap[imd_tags[ii]];
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
    memory->sfree(comm_buf);
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
          const int j = 3*inthash_lookup((inthash_t *)idmap, tag[i]);
          if (j != HASH_FAIL) {
            int ix = (image[i] & 1023) - 512;
            int iy = (image[i] >> 10 & 1023) - 512;
            int iz = (image[i] >> 20) - 512;

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
          const int j = 3*inthash_lookup((inthash_t *)idmap, tag[i]);
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
        const int j = 3*inthash_lookup((inthash_t *)idmap, buf[k].tag);
        if (j != HASH_FAIL) {
          recvcoord[j]   = buf[k].x;
          recvcoord[j+1] = buf[k].y;
          recvcoord[j+2] = buf[k].z;
        }
      }
    }

    /* done collecting frame data now communicate with IMD client. */

#if defined(LAMMPS_ASYNC_IMD)
    /* wake up i/o worker thread and release lock on i/o buffer
     * we can go back to our MD and let the i/o thread do the rest */
    buf_has_data=1;
    pthread_cond_signal(&write_cond);
    pthread_mutex_unlock(&write_mutex);
#else
    /* send coordinate data, if client is able to accept */
    if (clientsock && imdsock_selwrite(clientsock,0)) {
      imd_writen(clientsock, msgdata, msglen);
    }
    delete[] msgdata;
#endif

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
          int ix = (image[i] & 1023) - 512;
          int iy = (image[i] >> 10 & 1023) - 512;
          int iz = (image[i] >> 20) - 512;

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

/* ---------------------------------------------------------------------- */
void FixColvars::post_force_respa(int vflag, int ilevel, int iloop)
{
  /* only process IMD on the outmost RESPA level. */
  if (ilevel == nlevels_respa-1) post_force(vflag);
  return;
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixColvars::memory_usage(void)
{
  return static_cast<double>(num_coords+maxbuf+imd_forces)*size_one;
}

#endif
/* End of FixColvars class implementation. */

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:

