#include <pthread.h>

void* kokkos_test(void* args) { return args; }

int main(void) {
  pthread_t thread;
  /* Use NULL to avoid C++11. Some compilers
     do not have C++11 by default.  Forcing C++11
     in the compile tests can be done, but is unnecessary
  */
  pthread_create(&thread, NULL, kokkos_test, NULL);
  pthread_join(thread, NULL);
  return 0;
}
