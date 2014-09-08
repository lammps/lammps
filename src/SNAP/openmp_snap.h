
#ifndef LMP_OPENMP_SNAP_H
#define LMP_OPENMP_SNAP_H

#if defined(_OPENMP)
#include <omp.h>
#else
enum omp_sched_t {omp_sched_static, omp_sched_dynamic, omp_sched_guided, omp_sched_auto};
inline int omp_get_thread_num() { return 0;}
inline int omp_set_num_threads(int num_threads) {return 1;}
/* inline int __sync_fetch_and_add(int* ptr, int value) {int tmp = *ptr; ptr[0]+=value; return tmp;} */
inline void omp_set_schedule(omp_sched_t schedule,int modifier=1) {}
inline int omp_in_parallel() {return 0;}
#endif

#endif
