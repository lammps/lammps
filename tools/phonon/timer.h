#ifndef TIMER_H
#define TIMER_H

#include "stdio.h"
#include "stdlib.h"
#include "time.h"

class Timer {
public:
  Timer();

  void start();
  void stop();
  void print();
  double elapse();

private:
  clock_t t1, t2;
  double cpu_time_used;

  int flag;
};

#endif
