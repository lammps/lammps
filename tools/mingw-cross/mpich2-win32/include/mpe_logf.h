!      /* -*- Mode: Fortran; -*- */
!
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!
!
!
!  MPE Logging Return Codes
!
      integer MPE_LOG_OK, MPE_LOG_LOCKED_OUT
      parameter ( MPE_LOG_OK = 0, MPE_LOG_LOCKED_OUT = 1 )
      integer MPE_LOG_NO_MEMORY, MPE_LOG_FILE_PROB
      parameter ( MPE_LOG_NO_MEMORY = 2,  MPE_LOG_FILE_PROB = 3 )
      integer MPE_LOG_NOT_INITIALIZED, MPE_LOG_PACK_FAIL
      parameter ( MPE_LOG_NOT_INITIALIZED = 4, MPE_LOG_PACK_FAIL = 5 )
!
!  MPE Logging Function Prototypes
!
      integer  MPE_Init_log
      external MPE_Init_log
      integer  MPE_Initialized_logging
      external MPE_Initialized_logging
      integer  MPE_Describe_info_state
      external MPE_Describe_info_state
      integer  MPE_Describe_state
      external MPE_Describe_state
      integer  MPE_Describe_info_event
      external MPE_Describe_info_event
      integer  MPE_Describe_event
      external MPE_Describe_event
      integer  MPE_Log_get_event_number
      external MPE_Log_get_event_number
      integer  MPE_Log_get_state_eventIDs
      external MPE_Log_get_state_eventIDs
      integer  MPE_Log_get_solo_eventID
      external MPE_Log_get_solo_eventID
      integer  MPE_Start_log
      external MPE_Start_log
      integer  MPE_Log_send
      external MPE_Log_send
      integer  MPE_Log_receive
      external MPE_Log_receive
      integer  MPE_Log_pack
      external MPE_Log_pack
      integer  MPE_Log_event
      external MPE_Log_event
      integer  MPE_Log_bare_event
      external MPE_Log_bare_event
      integer  MPE_Log_info_event
      external MPE_Log_info_event
      external MPE_Log_sync_clocks
      integer  MPE_Log_sync_clocks
      integer  MPE_Stop_log
      external MPE_Stop_log
      integer  MPE_Finish_log
      external MPE_Finish_log
