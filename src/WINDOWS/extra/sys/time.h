#ifndef GETTIMEOFDAY_H
#define GETTIMEOFDAY_H

//#include <config.h>
#include <sys/timeb.h>
//#include "../include/time.h"


struct timeval 
{
    time_t tv_sec;
    time_t tv_usec;
};


__inline int gettimeofday(struct timeval *tp, void *tzp)
{
    
    struct _timeb timebuffer;
    
    _ftime(&timebuffer);
    tp->tv_sec = timebuffer.time;
    tp->tv_usec = timebuffer.millitm * 1000;
    
    return 0;
}

#endif /* GETTIMEOFDAY_H */