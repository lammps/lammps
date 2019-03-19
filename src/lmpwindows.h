#include <ciso646>
#if !defined(__MINGW32__)
#include "erf.h"
#endif
#include <direct.h>
#include <cmath>
#include <cstring>
// LAMMPS uses usleep with 100 ms arguments, no microsecond precision needed
#if !defined(__MINGW32__)
#include "sleep.h"
#endif

// some symbols have different names in Windows

#undef ATOBIGINT
#define ATOBIGINT _atoi64

#define pclose _pclose
#define strdup _strdup

// the following functions ared defined to get rid of
// 'ambiguous call to overloaded function' error in VSS for mismathched type arguments
#if !defined(__MINGW32__)
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
#ifndef S_IRWXU
# define S_IRWXU 0
#endif
#ifndef S_IRGRP
# define S_IRGRP 0
#endif
#ifndef S_IXGRP
# define S_IXGRP 0
#endif
inline int mkdir(const char *path, int){
  return _mkdir(path);
}

