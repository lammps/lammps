/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
/*
   MPE Miscellaneous routine headers
*/

#ifndef _MPE_MISC
#define _MPE_MISC


#ifndef MPE_NOMPI 
#include "mpi.h"

#if defined(__cplusplus)
extern "C" {
#endif

void MPE_Seq_begin ( MPI_Comm, int );
void MPE_Seq_end   ( MPI_Comm, int );

int  MPE_DelTag     ( MPI_Comm, int, void *, void * );
int  MPE_GetTags    ( MPI_Comm, int, MPI_Comm *, int * );
int  MPE_ReturnTags ( MPI_Comm, int, int );
int  MPE_TagsEnd    ( void );

void MPE_IO_Stdout_to_file ( char *, int );

void MPE_GetHostName       ( char *, int );

void MPE_Start_debugger ( void );
void MPE_Errors_call_debugger ( char *, char *, char ** );
void MPE_Errors_call_xdbx     ( char *, char * );
void MPE_Errors_call_dbx_in_xterm ( char *, char * );
void MPE_Signals_call_debugger ( void );

int  MPE_Decomp1d ( int, int, int, int *, int * );

void MPE_Comm_global_rank ( MPI_Comm, int, int * );

#if defined(__cplusplus)
}
#endif

#if (defined(__STDC__) || defined(__cplusplus))
void MPE_Errors_to_dbx ( MPI_Comm *, int *, ... );
#else
void MPE_Errors_to_dbx ( MPI_Comm *, int *, char *, char *, int * );
#endif

#else

#if defined(__cplusplus)
extern "C" {
#endif
void MPE_GetHostName       ( char *, int );
#if defined(__cplusplus)
}
#endif

#endif

#endif
