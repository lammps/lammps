/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#if !defined( _CLOG_CONST )
#define _CLOG_CONST

/*
   CLOG_FILE_TYPE determines default CLOG2 file extension, i.e. ".clog2"
   as well as the prefix for the local temporary clog2 file.
*/
#define  CLOG_FILE_TYPE      "clog2"

#define  CLOG_BOOL_T           int
#define  CLOG_BOOL_NULL       -1
#define  CLOG_BOOL_FALSE       0
#define  CLOG_BOOL_TRUE        1

#define  CLOG_PATH_STRLEN    256
#define  CLOG_ERR_STRLEN     256

#define  CLOG_PROCID_NULL     -1
#define  CLOG_NULL_FILE       -5

#define  CLOG_MAXTIME          100000000.0       /* later than all times */

#endif /* of _CLOG_CONST */
