/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#if !defined( _CLOG_UUID )
#define _CLOG_UUID

#include "clog_inttypes.h"
/*
   Only CLOG_UUID_NAME_SIZE-1 is useful.
   The last char in CLOG_UUID_NAME_SIZE is NULL char,
   the terminating character.
#define CLOG_UUID_NAME_SIZE 17
#define CLOG_UUID_SIZE sizeof(CLOG_int32_t)+sizeof(double)+CLOG_UUID_NAME_SIZE-1

typedef struct {
  CLOG_int32_t  rand;
  double        time;
  char          name[CLOG_UUID_NAME_SIZE];
} CLOG_Uuid_t;
*/
#define CLOG_UUID_NAME_SIZE   20
#define CLOG_UUID_SIZE        sizeof(CLOG_int32_t) + sizeof(double) \
                            + CLOG_UUID_NAME_SIZE

/* size of string representation of CLOG_Uuit_t */
#define CLOG_UUID_STR_SIZE   80

#if defined(__cplusplus)
extern "C" {
#endif

typedef char  CLOG_Uuid_t[ CLOG_UUID_SIZE ];

void CLOG_Uuid_init( void );

void CLOG_Uuid_finalize( void );

void CLOG_Uuid_generate( CLOG_Uuid_t uuid );

void CLOG_Uuid_sprint( CLOG_Uuid_t uuid, char *str );

int  CLOG_Uuid_is_equal( const CLOG_Uuid_t uuid1, const CLOG_Uuid_t uuid2 );

int  CLOG_Uuid_compare( const void *obj1, const void *obj2 );

void CLOG_Uuid_copy( const CLOG_Uuid_t src_uuid, CLOG_Uuid_t dest_uuid );

void CLOG_Uuid_swap_bytes( CLOG_Uuid_t uuid );

#if defined(__cplusplus)
}
#endif

#endif /* of _CLOG_UUID */
