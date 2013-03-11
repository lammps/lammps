#include "timer.h"

Timer::Timer()
{
  flag = 0;
  start();
return;
}

void Timer::start()
{
 t1 = clock();
 flag |= 1;

return;
}

void Timer::stop()
{
  if ( flag&1 ) {
    t2 = clock();
    flag |= 2;
  }
return;
}

void Timer::print()
{
  if ( (flag&3) != 3) return;

  cpu_time_used = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
  printf("Total CPU time used: %g seconds.\n", cpu_time_used);

return;
}

double Timer::elapse()
{
  if ( (flag&3) != 3) return 0.;
  else return ((double) (t2 - t1)) / CLOCKS_PER_SEC;
}
