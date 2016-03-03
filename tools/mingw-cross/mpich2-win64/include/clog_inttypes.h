/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
/*
   MPE Logging's generated header file for fixed size integer types
   and their printf format specifiers used in CLOG2.
*/
#if !defined( _CLOG_INTTYPES )
#define  _CLOG_INTTYPES

/* #include <stdint.h> */
/* #include <inttypes.h> */

typedef __int8         CLOG_int8_t;
typedef __int16        CLOG_int16_t;
typedef __int32        CLOG_int32_t;
typedef __int64        CLOG_int64_t;

/* 
   Define address-sized integer for the use of MPI_Comm_xxx_attr
   in clog2_commset.c.
*/
typedef __int32          CLOG_Pint;

#define i8fmt        "%hhd"
#define i16fmt       "%hd"
#define i32fmt       "%ld"
#define i64fmt       "%lld"

#endif
