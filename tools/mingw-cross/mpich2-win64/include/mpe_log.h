/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#ifndef _MPE_LOG_H_
#define _MPE_LOG_H_

#ifdef HAVE_WINDOWS_H
#ifdef USE_MPE_STATIC_LIBRARY
# define MPEU_DLL_SPEC
#else
# ifdef MPE_EXPORTS
#  define MPEU_DLL_SPEC __declspec(dllexport)
# else
#  define MPEU_DLL_SPEC __declspec(dllimport)
# endif
#endif
#else
# define MPEU_DLL_SPEC
#endif

/*
  Constants, MPE_Log_XXX, are for backward compatibility reasons.
  MPE is currently returning MPE_LOG_XXX status.
*/
/* function return values */
#define MPE_LOG_OK                0
#define MPE_Log_OK                MPE_LOG_OK
  /* no problems */
#define MPE_LOG_LOCKED_OUT        1
#define MPE_Log_LOCKED_OUT        MPE_LOG_LOCKED_OUT
  /* logs are being worked on, cannot insert any new entries */
#define MPE_LOG_NO_MEMORY         2
#define MPE_Log_NO_MEMORY         MPE_LOG_NO_MEMORY
  /* could not allocate memory for logging data */
#define MPE_LOG_FILE_PROB         3
#define MPE_Log_FILE_PROB         MPE_LOG_FILE_PROB
  /* cound not open file for writing out the logged info */
#define MPE_LOG_NOT_INITIALIZED   4
#define MPE_Log_NOT_INITIALIZED   MPE_LOG_NOT_INITIALIZED
  /* logging not initialized */
#define MPE_LOG_PACK_FAIL         5
#define MPE_Log_PACK_FAIL         MPE_LOG_PACK_FAIL

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

#include "clog_commset.h"

#if defined(__cplusplus)
extern "C" {
#endif

/* call before calling any other logging functions */
int MPE_Init_log( void );

int MPE_Initialized_logging( void );

/* create state with byte info data description lines in MPI_Comm */
int MPE_Describe_comm_state( MPI_Comm comm,
                             int state_startID, int state_finalID,
                             const char *name, const char *color,
                             const char *format );

/* create state with byte info data description lines in MPI_COMM_WORLD */
int MPE_Describe_info_state( int state_startID, int state_finalID,
                             const char *name, const char *color,
                             const char *format );

/* Internal MPE routine for MPI logging in MPI_Comm */
int MPE_Describe_known_state( const CLOG_CommIDs_t *commIDs, int local_thread,
                              int stateID, int state_startID, int state_finalID,
                              const char *name, const char *color,
                              const char *format );

/* create state description lines in MPI_COMM_WORLD */
int MPE_Describe_state( int state_startID, int state_finalID,
                        const char *name, const char *color );

/* create event with byte info data description lines in MPI_comm */
int MPE_Describe_comm_event( MPI_Comm comm, int eventID,
                             const char *name, const char *color,
                             const char *format );

/* create event with byte info data description lines in MPI_COMM_WORLD */
int MPE_Describe_info_event( int eventID,
                             const char *name, const char *color,
                             const char *format );

/* Internal MPE routine for MPI logging in MPI_Eventm */
int MPE_Describe_known_event( const CLOG_CommIDs_t *commIDs, int local_thread,
                              int eventID,
                              const char *name, const char *color,
                              const char *format );

/* create event description lines in MPI_COMM_WORLD */
int MPE_Describe_event( int eventID,
                        const char *name, const char *color );

/* fetch new user-space event number(s) */
int MPE_Log_get_event_number( void );
int MPE_Log_get_state_eventIDs( int *statedef_startID,
                                int *statedef_finalID );
int MPE_Log_get_solo_eventID( int *eventdef_eventID );

/* Internal MPE routine for MPI logging */

/* fetch a new known event ID of state for MPE Logging system */
int MPE_Log_get_known_eventID( void );
/* fetch a new known state ID for MPE Logging system */
int MPE_Log_get_known_stateID( void );
/* fetch a new known solo event ID for MPE Logging system */
int MPE_Log_get_known_solo_eventID( void );

/* log a MPI_Comm creation/destruction event */
int MPE_Log_commIDs_intracomm( const CLOG_CommIDs_t *commIDs, int local_thread,
                               int comm_etype,
                               const CLOG_CommIDs_t *intracommIDs );
int MPE_Log_commIDs_nullcomm( const CLOG_CommIDs_t *commIDs, int local_thread,
                              int comm_etype );
int MPE_Log_commIDs_intercomm( const CLOG_CommIDs_t *commIDs, int local_thread,
                               int comm_etype,
                               const CLOG_CommIDs_t *intercommIDs );

/* log the sending of a message in MPI_Comm */
int MPE_Log_commIDs_send( const CLOG_CommIDs_t *commIDs, int local_thread,
                          int receiver, int tag, int size );
int MPE_Log_comm_send( MPI_Comm comm, int receiver, int tag, int size );

/* log the sending of a message in MPI_COMM_WORLD */
int MPE_Log_send( int receiver, int tag, int size );

/* log the receiving of a message in MPI_Comm */
int MPE_Log_commIDs_receive( const CLOG_CommIDs_t *commIDs, int local_thread,
                             int sender, int tag, int size );
int MPE_Log_comm_receive( MPI_Comm comm, int sender, int tag, int size );

/* log the receiving of a message in MPI_COMM_WORLD */
int MPE_Log_receive( int sender, int tag, int size );

/* MPE_LOG_BYTESIZE is equal to sizeof(CLOG_BYTES) defined in clog.h */
#define MPE_LOG_BYTESIZE  (4 * sizeof(double))
typedef char MPE_LOG_BYTES[ MPE_LOG_BYTESIZE ];

int MPE_Log_pack( MPE_LOG_BYTES bytebuf, int *position,
                  char tokentype, int count, const void *data );

int MPE_Log_commIDs_event( const CLOG_CommIDs_t *commIDs, int local_thread,
                           int event, const char *bytebuf );

/* log an event in MPI_Comm */
int MPE_Log_comm_event( MPI_Comm comm, int event, const char *bytebuf );


/* log an event in MPI_COMM_WORLD */
int MPE_Log_event( int event, int data, const char *bytebuf );

/* log a bare event in MPI_COMM_WORLD */
int MPE_Log_bare_event( int event );

/* log an infomational event in MPI_COMM_WORLD */
int MPE_Log_info_event( int event, const char *bytebuf );

int MPE_Log_sync_clocks( void );

/* start logging events */
int MPE_Start_log( void );

/* stop logging events */
int MPE_Stop_log( void );

/* Synchronize thread related data to all processes */
void MPE_Log_thread_sync( int local_thread_count );

/* write out data to a file */
int MPE_Finish_log( const char *filename );

/* get the immutable merged logfile name */
const char* MPE_Log_merged_logfilename( void );

#if defined(__cplusplus)
}
#endif

#endif



/*
   The following documentation has little resemblance to current CLOG-2
   format.  They are kept here for historical reference purpose.
   2/17/2005         Anthony Chan 
*/
/*
The format:

Each line:
  type process task data cycle timestamp [comment]

    type - nonnegative integer representing a user-defined event type
    process - an integer representing the process in which the event occurred
    task - an integer representing a different notion of task.  Usually 
           ignored.
    data - an integer representing user data for the event
    cycle - an integer representing a time cycle, used to distinguish
            between time returned by a timer that "rolls over" during
            the run
    timestamp - an integer representing (when considered in conjuction
                with the cycle number) a time for the event.  Upshot treats
                the units as microseconds
    comment - an optional character string representing user data.  Currently
              12 character maximum, will soon hopefully be any length (really!)

All events from -100 to -1 are reserved header information events.  When
a log is produced, all [-100,-1] events will be moved to the top of the
logfile and have their timestamps set to 0.

All event from -101 and below are reserved system events.  This is to
provide some standardization for the logfiles, so various interpreting
programs can glean similar data from the same logfile.  All [-101,...)
events will have valid timestamps and will be left in time-sorted
order in the logfile.

Formats for reserved types:

  -1 Creation data                *not used*
     Comment: Creator and date

  -2 Number of events in the logfile   *not used*
     Data: number of events

  -3 Number of processors in the run
     Data: number of processes

  -4 Number of tasks used in the run  *not used*
     Task: number of tasks

  -5 Number of event types used        *not used*
     Data: number event types

  -6 Start time of the run
     Timestamp: start time

  -7 End time of the run
     Timestamp: end time

  -8 Number of times the timer cycled
     For example, if the timer's units are in microseconds, and it has a
     range of 0 - 2^32, and a run lasts 3 hours (range=4294 seconds, 3 hours=
     10800 seconds), the timer would have cycled at least twice.
     Data: number of timer cycles

  -9 Decription of event types     *not used*
     Data: event type
     Comment: Description

  -10 printf string for event types   *not used*
      Data: event type
      Comment: printf string

  -11 Rollover point
      The point at which the timer values 'rollover'
      Timestamp: rollover point

  -13 State definition
      Define a state based on the events that signal the beginning and end
      of the state.  Also, define what to call the state and what color/
      stipple pattern to give it in a graphical visualization tool.
      Task: start event
      Data: end event
      Comment: color:bitmap state name

      example:  -13 0 3 4 0 0 Green:boxes Rhode Island
      An event with type 3 will signify the entrance into a 'Rhode Island'
      state.  An event wil type 4 will signify the exit of the 'Rhode Island'
      state.

      States may be overlapped (enter a 'Rhode Island' state while in a
      'Wisconsin' state while in a 'Nevada' state), and the state name may
      have whitspace in it.
      
   -100 Synchronization event
        Sync events are used internally to sychronize timers on the various
        processes.  They do not appear in the logfiles.

   -101 Send message
        Represents the sending of a message
	Data: process ID of the receiving process
	Comment: <message-type tag of message> <size of the message, in bytes>

   -102 Receive message
        Represents the receiving of a message
	Data: process ID of the sending process
	Comment: <message-type tag of message> <size of the message, in bytes>

*/
