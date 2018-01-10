#include<Kokkos_Macros.hpp>
#if ( CUDA_VERSION < 9000 )
#define KOKKOS_IMPL_CUDA_BALLOT(x) __ballot(x)
#define KOKKOS_IMPL_CUDA_SHFL(x,y,z) __shfl(x,y,z)
#define KOKKOS_IMPL_CUDA_SHFL_UP(x,y,z) __shfl_up(x,y,z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN(x,y,z) __shfl_down(x,y,z)
#else
#define KOKKOS_IMPL_CUDA_BALLOT(x) __ballot_sync(0xffffffff,x)
#define KOKKOS_IMPL_CUDA_SHFL(x,y,z) __shfl_sync(0xffffffff,x,y,z)
#define KOKKOS_IMPL_CUDA_SHFL_UP(x,y,z) __shfl_up_sync(0xffffffff,x,y,z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN(x,y,z) __shfl_down_sync(0xffffffff,x,y,z)
#endif 
