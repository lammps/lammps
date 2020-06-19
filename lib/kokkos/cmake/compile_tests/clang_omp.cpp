#include <omp.h>

int main(int, char**) {
  int thr = omp_get_num_threads();
  if (thr > 0)
    return thr;
  else
    return 0;
}
