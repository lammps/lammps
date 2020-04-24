#include <pthread.h>

void* kokkos_test(void* args) { return args; }

int main(void) {
  pthread_t thread;
  pthread_create(&thread, NULL, kokkos_test, NULL);
  pthread_join(thread, NULL);
  return 0;
}
