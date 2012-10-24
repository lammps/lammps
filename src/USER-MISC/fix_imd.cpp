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
   The FixIMD class contains code from VMD and NAMD which is copyrighted
   by the Board of Trustees of the University of Illinois and is free to
   use with LAMMPS according to point 2 of the UIUC license (10% clause):

" Licensee may, at its own expense, create and freely distribute
complimentary works that interoperate with the Software, directing others to
the TCBG server to license and obtain the Software itself. Licensee may, at
its own expense, modify the Software to make derivative works.  Except as
explicitly provided below, this License shall apply to any derivative work
as it does to the original Software distributed by Illinois.  Any derivative
work should be clearly marked and renamed to notify users that it is a
modified version and not the original Software distributed by Illinois.
Licensee agrees to reproduce the copyright notice and other proprietary
markings on any derivative work and to include in the documentation of such
work the acknowledgement:

 "This software includes code developed by the Theoretical and Computational
  Biophysics Group in the Beckman Institute for Advanced Science and
  Technology at the University of Illinois at Urbana-Champaign."

Licensee may redistribute without restriction works with up to 1/2 of their
non-comment source code derived from at most 1/10 of the non-comment source
code developed by Illinois and contained in the Software, provided that the
above directions for notice and acknowledgement are observed.  Any other
distribution of the Software or any derivative work requires a separate
license with Illinois.  Licensee may contact Illinois (vmd@ks.uiuc.edu) to
negotiate an appropriate license for such distribution."
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:  Axel Kohlmeyer (TempleU)
   IMD API, hash, and socket code written by: John E. Stone,
   Justin Gullingsrud, and James Phillips, (TCBG, Beckman Institute, UIUC)
------------------------------------------------------------------------- */

#include "fix_imd.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
#include <winsock2.h>
#else
#include <arpa/inet.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/file.h>
#endif

#include <errno.h>

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

/********** API definitions of the VMD/NAMD code ************************
 * This code was taken and adapted from VMD-1.8.7/NAMD-2.7 in Sep 2009. *
 * If there are any bugs or problems, please contact akohlmey@gmail.com *
 ************************************************************************/

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/* part 1: Interactive MD (IMD) API */

#include <limits.h>

#if ( INT_MAX == 2147483647 )
typedef int     int32;
#else
typedef short   int32;
#endif

typedef struct {
  int32 type;
  int32 length;
} IMDheader;

#define IMDHEADERSIZE 8
#define IMDVERSION 2

typedef enum IMDType_t {
  IMD_DISCONNECT,   /**< close IMD connection, leaving sim running */
  IMD_ENERGIES,     /**< energy data block                         */
  IMD_FCOORDS,      /**< atom coordinates                          */
  IMD_GO,           /**< start the simulation                      */
  IMD_HANDSHAKE,    /**< endianism and version check message       */
  IMD_KILL,         /**< kill the simulation job, shutdown IMD     */
  IMD_MDCOMM,       /**< MDComm style force data                   */
  IMD_PAUSE,        /**< pause the running simulation              */
  IMD_TRATE,        /**< set IMD update transmission rate          */
  IMD_IOERROR       /**< indicate an I/O error                     */
} IMDType;          /**< IMD command message type enumerations */

typedef struct {
  int32 tstep;      /**< integer timestep index                    */
  float T;          /**< Temperature in degrees Kelvin             */
  float Etot;       /**< Total energy, in Kcal/mol                 */
  float Epot;       /**< Potential energy, in Kcal/mol             */
  float Evdw;       /**< Van der Waals energy, in Kcal/mol         */
  float Eelec;      /**< Electrostatic energy, in Kcal/mol         */
  float Ebond;      /**< Bond energy, Kcal/mol                     */
  float Eangle;     /**< Angle energy, Kcal/mol                    */
  float Edihe;      /**< Dihedral energy, Kcal/mol                 */
  float Eimpr;      /**< Improper energy, Kcal/mol                 */
} IMDEnergies;      /**< IMD simulation energy report structure    */

/** Send control messages - these consist of a header with no subsequent data */
static int imd_handshake(void *);    /**< check endianness, version compat */
/** Receive header and data */
static IMDType imd_recv_header(void *, int32 *);
/** Receive MDComm-style forces, units are Kcal/mol/angstrom */
static int imd_recv_mdcomm(void *, int32, int32 *, float *);
/** Receive energies */
static int imd_recv_energies(void *, IMDEnergies *);
/** Receive atom coordinates. */
static int imd_recv_fcoords(void *, int32, float *);
/** Prepare IMD data packet header */
static void imd_fill_header(IMDheader *header, IMDType type, int32 length);
/** Write data to socket */
static int32 imd_writen(void *s, const char *ptr, int32 n);

/* part 2: abstracts platform-dependent routines/APIs for using sockets */

typedef struct {
  struct sockaddr_in addr; /* address of socket provided by bind() */
  int addrlen;             /* size of the addr struct */
  int sd;                  /* socket file descriptor */
} imdsocket;

static int   imdsock_init(void);
static void *imdsock_create(void);
static int   imdsock_bind(void *, int);
static int   imdsock_listen(void *);
static void *imdsock_accept(void *);  /* return new socket */
static int   imdsock_write(void *, const void *, int);
static int   imdsock_read(void *, void *, int);
static int   imdsock_selread(void *, int);
static int   imdsock_selwrite(void *, int);
static void  imdsock_shutdown(void *);
static void  imdsock_destroy(void *);

/***************************************************************
 * End of API definitions of the VMD/NAMD code.                *
 * The implementation follows at the end of the file.          *
 ***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;

/* struct for packed data communication of coordinates and forces. */
struct commdata {
  int tag;
  float x,y,z;
};

/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID imd <imd_trate> <imd_port> [unwrap (on|off)] [fscale <imd_fscale>]
 ***************************************************************/
FixIMD::FixIMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix imd command");

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
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixIMD::~FixIMD()
{

#if defined(LAMMPS_ASYNC_IMD)
  if (me == 0) {
    pthread_mutex_lock(&write_mutex);
    buf_has_data=-1;
    pthread_cond_signal(&write_cond);
    pthread_mutex_unlock(&write_mutex);
    pthread_join(iothread, NULL);

    /* cleanup */
    pthread_attr_destroy(&iot_attr);
    pthread_mutex_destroy(&write_mutex);
    pthread_cond_destroy(&write_cond);
  }
#endif

  inthash_t *hashtable = (inthash_t *)idmap;
  memory->destroy(comm_buf);
  memory->destroy(force_buf);
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
}

/* ---------------------------------------------------------------------- */
int FixIMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixIMD::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  return;
}

/* ---------------------------------------------------------------------- */

/* (re-)connect to an IMD client (e.g. VMD). return 1 if
   new connection was made, 0 if not. */
int FixIMD::reconnect()
{
  /* set up IMD communication, but only if needed. */
  imd_inactive = 0;
  imd_terminate = 0;

  if (me == 0) {
    if (clientsock) return 1;
    if (screen && connect_msg)
      if (nowait_flag)
        fprintf(screen,"Listening for IMD connection on port %d. Transfer rate %d.\n",imd_port, imd_trate);
      else
        fprintf(screen,"Waiting for IMD connection on port %d. Transfer rate %d.\n",imd_port, imd_trate);

    connect_msg = 0;
    clientsock = NULL;
    if (nowait_flag) {
      int retval = imdsock_selread(localsock,0);
      if (retval > 0) {
        clientsock = imdsock_accept(localsock);
      } else {
        imd_inactive = 1;
        return 0;
      }
    } else {
      int retval=0;
      do {
        retval = imdsock_selread(localsock, 60);
      } while (retval <= 0);
      clientsock = imdsock_accept(localsock);
    }

    if (!imd_inactive && !clientsock) {
      if (screen)
        fprintf(screen, "IMD socket accept error. Dropping connection.\n");
      imd_terminate = 1;
      return 0;
    } else {
      /* check endianness and IMD protocol version. */
      if (imd_handshake(clientsock)) {
        if (screen)
          fprintf(screen, "IMD handshake error. Dropping connection.\n");
        imdsock_destroy(clientsock);
        imd_terminate = 1;
        return 0;
      } else {
        int32 length;
        if (imdsock_selread(clientsock, 1) != 1 ||
            imd_recv_header(clientsock, &length) != IMD_GO) {
          if (screen)
            fprintf(screen, "Incompatible IMD client version? Dropping connection.\n");
          imdsock_destroy(clientsock);
          imd_terminate = 1;
          return 0;
        } else {
          return 1;
        }
      }
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */
/* wait for IMD client (e.g. VMD) to respond, initialize communication
 * buffers and collect tag/id maps. */
void FixIMD::setup(int)
{
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
}

/* worker threads for asynchronous i/o */
#if defined(LAMMPS_ASYNC_IMD)
/* c bindings wrapper */
void *fix_imd_ioworker(void *t)
{
  FixIMD *imd=(FixIMD *)t;
  imd->ioworker();
  return NULL;
}

/* the real i/o worker thread */
void FixIMD::ioworker()
{
  while (1) {
    pthread_mutex_lock(&write_mutex);
    if (buf_has_data < 0) {
      /* master told us to go away */
      fprintf(screen,"Asynchronous I/O thread is exiting.\n");
      buf_has_data=0;
      pthread_mutex_unlock(&write_mutex);
      pthread_exit(NULL);
    } else if (buf_has_data > 0) {
      /* send coordinate data, if client is able to accept */
      if (clientsock && imdsock_selwrite(clientsock,0)) {
        imd_writen(clientsock, msgdata, msglen);
      }
      delete[] msgdata;
      buf_has_data=0;
      pthread_mutex_unlock(&write_mutex);
    } else {
      /* nothing to write out yet. wait on condition. */
      pthread_cond_wait(&write_cond, &write_mutex);
      pthread_mutex_unlock(&write_mutex);
    }
  }
}
#endif

/* ---------------------------------------------------------------------- */
/* Main IMD protocol handler:
 * Send coodinates, energies, and add IMD forces to atoms. */
void FixIMD::post_force(int vflag)
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
          const int j = 3*inthash_lookup((inthash_t *)idmap, tag[i]);
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

/* ---------------------------------------------------------------------- */
void FixIMD::post_force_respa(int vflag, int ilevel, int iloop)
{
  /* only process IMD on the outmost RESPA level. */
  if (ilevel == nlevels_respa-1) post_force(vflag);
  return;
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixIMD::memory_usage(void)
{
  return static_cast<double>(num_coords+maxbuf+imd_forces)*size_one;
}

/* End of FixIMD class implementation. */

/***************************************************************************/

/* NOTE: the following code is the based on the example implementation
 * of the IMD protocol API from VMD and NAMD. The UIUC license allows
 * to re-use up to 10% of a project's code to be used in other software */

/***************************************************************************
 * DESCRIPTION:
 *   Socket interface, abstracts machine dependent APIs/routines.
 ***************************************************************************/

int imdsock_init(void) {
#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
  int rc = 0;
  static int initialized=0;

  if (!initialized) {
    WSADATA wsdata;
    rc = WSAStartup(MAKEWORD(1,1), &wsdata);
    if (rc == 0)
      initialized = 1;
  }

  return rc;
#else
  return 0;
#endif
}


void * imdsock_create(void) {
  imdsocket * s;

  s = (imdsocket *) malloc(sizeof(imdsocket));
  if (s != NULL)
    memset(s, 0, sizeof(imdsocket));

  if ((s->sd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
    printf("Failed to open socket.");
    free(s);
    return NULL;
  }

  return (void *) s;
}

int imdsock_bind(void * v, int port) {
  imdsocket *s = (imdsocket *) v;
  memset(&(s->addr), 0, sizeof(s->addr));
  s->addr.sin_family = PF_INET;
  s->addr.sin_port = htons(port);

  return bind(s->sd, (struct sockaddr *) &s->addr, sizeof(s->addr));
}

int imdsock_listen(void * v) {
  imdsocket *s = (imdsocket *) v;
  return listen(s->sd, 5);
}

void *imdsock_accept(void * v) {
  int rc;
  imdsocket *new_s = NULL, *s = (imdsocket *) v;
#if defined(ARCH_AIX5) || defined(ARCH_AIX5_64) || defined(ARCH_AIX6_64)
  unsigned int len;
#define _SOCKLEN_TYPE unsigned int
#elif defined(SOCKLEN_T)
  SOCKLEN_T len;
#define _SOCKLEN_TYPE SOCKLEN_T
#elif defined(_POSIX_SOURCE) || (defined(__APPLE__) && defined(__MACH__)) || defined(__linux)
  socklen_t len;
#define _SOCKLEN_TYPE socklen_t
#else
#define _SOCKLEN_TYPE int
  int len;
#endif

  len = sizeof(s->addr);
  rc = accept(s->sd, (struct sockaddr *) &s->addr, ( _SOCKLEN_TYPE * ) &len);
  if (rc >= 0) {
    new_s = (imdsocket *) malloc(sizeof(imdsocket));
    if (new_s != NULL) {
      *new_s = *s;
      new_s->sd = rc;
    }
  }
  return (void *)new_s;
}

int  imdsock_write(void * v, const void *buf, int len) {
  imdsocket *s = (imdsocket *) v;
#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
  return send(s->sd, (const char*) buf, len, 0);  /* windows lacks the write() call */
#else
  return write(s->sd, buf, len);
#endif
}

int  imdsock_read(void * v, void *buf, int len) {
  imdsocket *s = (imdsocket *) v;
#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
  return recv(s->sd, (char*) buf, len, 0); /* windows lacks the read() call */
#else
  return read(s->sd, buf, len);
#endif

}

void imdsock_shutdown(void *v) {
  imdsocket * s = (imdsocket *) v;
  if (s == NULL)
    return;

#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
  shutdown(s->sd, SD_SEND);
#else
  shutdown(s->sd, 1);  /* complete sends and send FIN */
#endif
}

void imdsock_destroy(void * v) {
  imdsocket * s = (imdsocket *) v;
  if (s == NULL)
    return;

#if defined(_MSC_VER) || defined(__MINGW32_VERSION)
  closesocket(s->sd);
#else
  close(s->sd);
#endif
  free(s);
}

int imdsock_selread(void *v, int sec) {
  imdsocket *s = (imdsocket *)v;
  fd_set rfd;
  struct timeval tv;
  int rc;

  if (v == NULL) return 0;

  FD_ZERO(&rfd);
  FD_SET(s->sd, &rfd);
  memset((void *)&tv, 0, sizeof(struct timeval));
  tv.tv_sec = sec;
  do {
    rc = select(s->sd+1, &rfd, NULL, NULL, &tv);
  } while (rc < 0 && errno == EINTR);
  return rc;

}

int imdsock_selwrite(void *v, int sec) {
  imdsocket *s = (imdsocket *)v;
  fd_set wfd;
  struct timeval tv;
  int rc;

  if (v == NULL) return 0;

  FD_ZERO(&wfd);
  FD_SET(s->sd, &wfd);
  memset((void *)&tv, 0, sizeof(struct timeval));
  tv.tv_sec = sec;
  do {
    rc = select(s->sd + 1, NULL, &wfd, NULL, &tv);
  } while (rc < 0 && errno == EINTR);
  return rc;
}

/* end of socket code. */
/*************************************************************************/

/*************************************************************************/
/* start of imd API code. */
/* Only works with aligned 4-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap4_aligned(void *v, long ndata) {
  int *data = (int *) v;
  long i;
  int *N;
  for (i=0; i<ndata; i++) {
    N = data + i;
    *N=(((*N>>24)&0xff) | ((*N&0xff)<<24) |
        ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
  }
}


/** structure used to perform byte swapping operations */
typedef union {
  int32 i;
  struct {
    unsigned int highest : 8;
    unsigned int high    : 8;
    unsigned int low     : 8;
    unsigned int lowest  : 8;
  } b;
} netint;

static int32 imd_htonl(int32 h) {
  netint n;
  n.b.highest = h >> 24;
  n.b.high    = h >> 16;
  n.b.low     = h >> 8;
  n.b.lowest  = h;
  return n.i;
}

static int32 imd_ntohl(int32 n) {
  netint u;
  u.i = n;
  return (u.b.highest << 24 | u.b.high << 16 | u.b.low << 8 | u.b.lowest);
}

static void imd_fill_header(IMDheader *header, IMDType type, int32 length) {
  header->type = imd_htonl((int32)type);
  header->length = imd_htonl(length);
}

static void swap_header(IMDheader *header) {
  header->type = imd_ntohl(header->type);
  header->length= imd_ntohl(header->length);
}

static int32 imd_readn(void *s, char *ptr, int32 n) {
  int32 nleft;
  int32 nread;

  nleft = n;
  while (nleft > 0) {
    if ((nread = imdsock_read(s, ptr, nleft)) < 0) {
      if (errno == EINTR)
        nread = 0;         /* and call read() again */
      else
        return -1;
    } else if (nread == 0)
      break;               /* EOF */
    nleft -= nread;
    ptr += nread;
  }
  return n-nleft;
}

static int32 imd_writen(void *s, const char *ptr, int32 n) {
  int32 nleft;
  int32 nwritten;

  nleft = n;
  while (nleft > 0) {
    if ((nwritten = imdsock_write(s, ptr, nleft)) <= 0) {
      if (errno == EINTR)
        nwritten = 0;
      else
        return -1;
    }
    nleft -= nwritten;
    ptr += nwritten;
  }
  return n;
}

int imd_disconnect(void *s) {
  IMDheader header;
  imd_fill_header(&header, IMD_DISCONNECT, 0);
  return (imd_writen(s, (char *)&header, IMDHEADERSIZE) != IMDHEADERSIZE);
}

int imd_handshake(void *s) {
  IMDheader header;
  imd_fill_header(&header, IMD_HANDSHAKE, 1);
  header.length = IMDVERSION;   /* Not byteswapped! */
  return (imd_writen(s, (char *)&header, IMDHEADERSIZE) != IMDHEADERSIZE);
}

/* The IMD receive functions */

IMDType imd_recv_header(void *s, int32 *length) {
  IMDheader header;
  if (imd_readn(s, (char *)&header, IMDHEADERSIZE) != IMDHEADERSIZE)
    return IMD_IOERROR;
  swap_header(&header);
  *length = header.length;
  return IMDType(header.type);
}

int imd_recv_mdcomm(void *s, int32 n, int32 *indices, float *forces) {
  if (imd_readn(s, (char *)indices, 4*n) != 4*n) return 1;
  if (imd_readn(s, (char *)forces, 12*n) != 12*n) return 1;
  return 0;
}

int imd_recv_energies(void *s, IMDEnergies *energies) {
  return (imd_readn(s, (char *)energies, sizeof(IMDEnergies))
          != sizeof(IMDEnergies));
}

int imd_recv_fcoords(void *s, int32 n, float *coords) {
  return (imd_readn(s, (char *)coords, 12*n) != 12*n);
}

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
