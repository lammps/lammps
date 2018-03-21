#include <cstdio>
#include <cuda_runtime_api.h>
int main()
{
	cudaDeviceProp prop;
  const cudaError_t err_code = cudaGetDeviceProperties(&prop, 0);
  if (cudaSuccess != err_code) {
		fprintf(stderr,"cudaGetDeviceProperties failed: %s\n", cudaGetErrorString(err_code));
		return -1;
	}
  switch (prop.major) {
    case 3:
      printf("Kepler"); break;
    case 5:
      printf("Maxwell"); break;
    case 6:
      printf("Pascal"); break;
    default:
      fprintf(stderr, "Unspported Device %d%d\n", (int)prop.major, (int)prop.minor);
      return -1;
  }
  printf("%d%d\n", (int)prop.major, (int)prop.minor);
  return 0;
}
