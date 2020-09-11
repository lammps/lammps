#include <pthread.h>

void* kokkos_test(void* args) { return args; }

int main(void) {
  pthread_t thread;
  /* Use nullptr to avoid C++11. Some compilers
     do not have C++11 by default.  Forcing C++11
     in the compile tests can be done, but is unnecessary
  */
  pthread_create(&thread, nullptr, kokkos_test, nullptr);
  pthread_join(thread, nullptr);
  return 0;
}
