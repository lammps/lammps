#include "timer.h"

#include <cstdio>

/* -----------------------------------------------------------------------------
 * Initialization of time
 * -------------------------------------------------------------------------- */
Timer::Timer()
{
   flag = 0;
   start();

   return;
}

/* -----------------------------------------------------------------------------
 * public function, start the timer
 * -------------------------------------------------------------------------- */
void Timer::start()
{
   t1 = clock();
   flag |= 1;
   return;
}

/* -----------------------------------------------------------------------------
 * public function, stop the timer
 * -------------------------------------------------------------------------- */
void Timer::stop()
{
   if ( flag&1 ) {
      t2 = clock();
      flag |= 2;
   }
   return;
}

/* -----------------------------------------------------------------------------
 * public function, print the total time used after timer stops
 * -------------------------------------------------------------------------- */
void Timer::print()
{
   if ( (flag&3) != 3) return;

   cpu_time_used = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
   printf("Total CPU time used: %g seconds.\n", cpu_time_used);

   return;
}

/* -----------------------------------------------------------------------------
 * public function, return the total time used up to now, in seconds
 * -------------------------------------------------------------------------- */
double Timer::elapse()
{
   if ( (flag&3) != 3) return 0.;
   else return ((double) (t2 - t1)) / CLOCKS_PER_SEC;
}
