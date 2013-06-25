/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#if !defined( _CLOG_COMMSET )
#define _CLOG_COMMSET

#include "clog_inttypes.h"
#include "clog_uuid.h"
#include "clog_const.h"

#define CLOG_CommGID_t            CLOG_Uuid_t
#define CLOG_CommLID_t            CLOG_int32_t
#define CLOG_ThreadLID_t          int
/*
typedef CLOG_Uuid_t  CLOG_CommGID_t
typedef CLOG_int32_t CLOG_CommLID_t
typedef int          CLOG_ThreadLID_t
*/

/* Define CLOG communicator event types used in log_mpi_xxxx.c */
#define CLOG_COMM_WORLD_CREATE    0    /* MPI_COMM_WORLD creation */
#define CLOG_COMM_SELF_CREATE     1    /* MPI_COMM_SELF  creation */
#define CLOG_COMM_FREE            10   /* MPI_Comm_free() */
#define CLOG_COMM_INTRA_CREATE    100  /* intracomm creation */
#define CLOG_COMM_INTRA_LOCAL     101  /* local intracomm of intercomm */
#define CLOG_COMM_INTRA_REMOTE    102  /* remote intracomm of intercomm */
#define CLOG_COMM_INTER_CREATE    1000 /* intercomm creation */

/* Define CLOG communicator kind */
#define CLOG_COMM_KIND_UNKNOWN   -1    /* is UNKNOWN */
#define CLOG_COMM_KIND_INTER      0    /* is intercommunicator */
#define CLOG_COMM_KIND_INTRA      1    /* is intracommunicator */
#define CLOG_COMM_KIND_LOCAL      2    /* is local  intracommunicator */
#define CLOG_COMM_KIND_REMOTE     3    /* is remote intracommunicator */

/* Some CLOG communicator internal constant */
#define CLOG_COMM_TABLE_INCRE     10
#define CLOG_COMM_LID_NULL       -999999999
#define CLOG_COMM_RANK_NULL      -1
#define CLOG_COMM_WRANK_NULL     -1
#define CLOG_COMM_TAG_START       100000

#if !defined( CLOG_NOMPI )
#include "mpi.h"
#else 
/*
    To avoid mpi_null.h from being exposed to the user's include_dir,
    the definition of MPI_Comm has to match that in mpi_null.h.
    The "#if !defined( _MPI_NULL_MPI_COMM )" is to avoid duplicated
    definition of MPI_Comm when both mpi_null.h and this .h are used
    in the same .c file.
*/
#if !defined( _MPI_NULL_MPI_COMM )
#define _MPI_NULL_MPI_COMM
typedef int  MPI_Comm;
#endif
#endif /* Endof if !defined( CLOG_NOMPI ) */

typedef struct CLOG_CommIDs_t_ {
          CLOG_CommGID_t   global_ID;  /* global comm ID */
          CLOG_CommLID_t   local_ID;   /* local comm ID */
          CLOG_int32_t     kind;       /* value = CLOG_COMM_KIND_xxxx */
          int              world_rank; /* rank in MPI_COMM_WORLD */
          int              comm_rank;  /* rank of comm labelled by global_ID */
          MPI_Comm         comm;
   struct CLOG_CommIDs_t_ *next;       /* related CLOG_CommIDs_t* */
} CLOG_CommIDs_t;

/* Definition for CLOG Communicator control data structure */
typedef struct {
    int                 LID_key;
    int                 world_size;    /* size returned by MPI_Comm_size */
    int                 world_rank;    /* rank returned by MPI_Comm_rank */
    unsigned int        max;
    CLOG_CommLID_t      count;
    CLOG_CommIDs_t     *table;
    CLOG_CommIDs_t     *IDs4world;
    CLOG_CommIDs_t     *IDs4self;
} CLOG_CommSet_t;

#if defined(__cplusplus)
extern "C" {
#endif

CLOG_CommSet_t* CLOG_CommSet_create( void );

void CLOG_CommSet_free( CLOG_CommSet_t **comm_handle );

#if defined( CLOG_IMPL )
void CLOG_CommSet_add_GID(       CLOG_CommSet_t *commset,
                           const CLOG_CommGID_t  commgid );

void CLOG_CommSet_append_GIDs(       CLOG_CommSet_t *parent_commset,
                                     int             child_table_count,
                               const CLOG_CommIDs_t *child_table );

CLOG_BOOL_T CLOG_CommSet_sync_IDs(       CLOG_CommSet_t *parent_commset,
                                         int             child_table_count,
                                   const CLOG_CommIDs_t *child_table );
#endif

void CLOG_CommSet_init( CLOG_CommSet_t *commset );

const CLOG_CommIDs_t* CLOG_CommSet_add_intracomm( CLOG_CommSet_t *commset,
                                                  MPI_Comm comm );

const CLOG_CommIDs_t*
CLOG_CommSet_add_intercomm(       CLOG_CommSet_t *commset,
                                  MPI_Comm        intercomm,
                            const CLOG_CommIDs_t *intracommIDs );

CLOG_CommLID_t CLOG_CommSet_get_LID( CLOG_CommSet_t *commset, MPI_Comm comm );

const CLOG_CommIDs_t* CLOG_CommSet_get_IDs( CLOG_CommSet_t *commset,
                                            MPI_Comm comm );

void CLOG_CommSet_merge( CLOG_CommSet_t *commset );

#if defined( CLOG_IMPL )
int CLOG_CommSet_write( const CLOG_CommSet_t *commset, int fd,
                              CLOG_BOOL_T     do_byte_swap );

int CLOG_CommSet_read( CLOG_CommSet_t *commset, int fd,
                       CLOG_BOOL_T     do_byte_swap );

void CLOG_CommSet_print( CLOG_CommSet_t *commset, FILE *stream );
#endif

#if defined(__cplusplus)
}
#endif

#endif /* of _CLOG_COMMSET */
