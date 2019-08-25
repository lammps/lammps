//# ifdef USE_STDAFX
//# include "stdafx.h"
//# endif

# include "logexc.h"

message_logger std_log;

message_logger &message_logger::global(){
  if(!glogp){
    std_log.set_global(true);
  }
  return *glogp;
}

message_logger *message_logger::glogp=NULL;
stdfile_logger default_log("",0,stdout,stderr,vblALLBAD|vblMESS1,vblFATAL,1);

const char *fmt(const char *format,...){
  va_list args;
  va_start(args,format);
  static char buff[1024];
  vsnprintf(buff,1024,format,args);
  return buff;
}
