#include <iso646.h>
#if !defined(__MINGW32__)
#include "erf.h"
#endif
#include "direct.h"
#include "math.h"
// LAMMPS uses usleep with 100 ms arguments, no microsecond precision needed
#if !defined(__MINGW32__)
#include "sleep.h"
#endif

// some symbols have different names in Windows

#undef ATOBIGINT
#define ATOBIGINT _atoi64

#define pclose _pclose
#define __restrict__ __restrict

// the following functions ared defined to get rid of
// 'ambiguous call to overloaded function' error in VSS for mismathched type arguments

#if defined(__MINGW32__)
static inline double pow(const double &x, const int n) {
  double yy,ww;

  if (x == 0.0) return 0.0;
  int nn = (n > 0) ? n : -n;
  ww = x;

  for (yy = 1.0; nn != 0; nn >>= 1, ww *=ww)
    if (nn & 1) yy *= ww;

  return (n > 0) ? yy : 1.0/yy;
}

static inline int pow(const int x, const int n) {
int yy,ww,nn;

  if ((x == 0) || (x == 1)) return x;
  if (n < 0) return 0;
  ww = x;
  nn = n;

  for (yy = 1; nn != 0; nn >>= 1, ww *=ww)
    if (nn & 1) yy *= ww;

  return yy;
}

#else
inline double pow(int i, int j){
  return pow((double)i,j);
}
#endif

inline double sqrt(int i){
  return sqrt((double) i);
}

inline double fabs(int i){
  return fabs((double) i);
}

inline double trunc(double x) {
  return x > 0 ? floor(x) : ceil(x);
}

// Windows version of mkdir function does not have permission flags
# define S_IRWXU 0
# define S_IRGRP 0
# define S_IXGRP 0
inline int mkdir(const char *path, int){
  return _mkdir(path);
}
