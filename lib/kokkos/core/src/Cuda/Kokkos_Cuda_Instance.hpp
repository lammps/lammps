#ifndef KOKKOS_CUDA_INSTANCE_HPP_
#define KOKKOS_CUDA_INSTANCE_HPP_

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct CudaTraits {
  enum { WarpSize       = 32      /* 0x0020 */ };
  enum { WarpIndexMask  = 0x001f  /* Mask for warpindex */ };
  enum { WarpIndexShift = 5       /* WarpSize == 1 << WarpShift */ };

  enum { ConstantMemoryUsage    = 0x008000 /* 32k bytes */ };
  enum { ConstantMemoryCache    = 0x002000 /*  8k bytes */ };
  enum { KernelArgumentLimit    = 0x001000 /*  4k bytes */ };

  typedef unsigned long
    ConstantGlobalBufferType[ ConstantMemoryUsage / sizeof(unsigned long) ];

#if defined(KOKKOS_ARCH_VOLTA) || \
    defined(KOKKOS_ARCH_PASCAL)
  enum { ConstantMemoryUseThreshold = 0x000200 /* 0 bytes -> always use constant (or global)*/ };
#else
  enum { ConstantMemoryUseThreshold = 0x000200 /* 512 bytes */ };
#endif

  KOKKOS_INLINE_FUNCTION static
  CudaSpace::size_type warp_count( CudaSpace::size_type i )
    { return ( i + WarpIndexMask ) >> WarpIndexShift ; }

  KOKKOS_INLINE_FUNCTION static
  CudaSpace::size_type warp_align( CudaSpace::size_type i )
    {
      enum { Mask = ~CudaSpace::size_type( WarpIndexMask ) };
      return ( i + WarpIndexMask ) & Mask ;
    }
};

//----------------------------------------------------------------------------

CudaSpace::size_type cuda_internal_multiprocessor_count();
CudaSpace::size_type cuda_internal_maximum_warp_count();
CudaSpace::size_type cuda_internal_maximum_grid_count();
CudaSpace::size_type cuda_internal_maximum_shared_words();

CudaSpace::size_type cuda_internal_maximum_concurrent_block_count();

CudaSpace::size_type * cuda_internal_scratch_flags( const Cuda&, const CudaSpace::size_type size );
CudaSpace::size_type * cuda_internal_scratch_space( const Cuda&, const CudaSpace::size_type size );
CudaSpace::size_type * cuda_internal_scratch_unified( const Cuda&, const CudaSpace::size_type size );

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

class CudaInternal {
private:

  CudaInternal( const CudaInternal & );
  CudaInternal & operator = ( const CudaInternal & );


public:

  typedef Cuda::size_type size_type ;

  int         m_cudaDev ;

  // Device Properties
  int         m_cudaArch ;
  unsigned    m_multiProcCount ;
  unsigned    m_maxWarpCount ;
  unsigned    m_maxBlock ;
  unsigned    m_maxSharedWords ;
  uint32_t    m_maxConcurrency ;
  int         m_shmemPerSM ;
  int         m_maxShmemPerBlock ;
  int         m_regsPerSM ;
  int         m_maxBlocksPerSM ;
  int         m_maxThreadsPerSM ;
  int         m_maxThreadsPerBlock ;

  mutable size_type   m_scratchSpaceCount ;
  mutable size_type   m_scratchFlagsCount ;
  mutable size_type   m_scratchUnifiedCount ;
  mutable size_type   m_scratchFunctorSize ;
  size_type   m_scratchUnifiedSupported ;
  size_type   m_streamCount ;
  mutable size_type * m_scratchSpace ;
  mutable size_type * m_scratchFlags ;
  mutable size_type * m_scratchUnified ;
  mutable size_type * m_scratchFunctor ;
  uint32_t  * m_scratchConcurrentBitset ;
  cudaStream_t m_stream ;

  static int was_initialized;
  static int was_finalized;

  static CudaInternal & singleton();

  int verify_is_initialized( const char * const label ) const ;

  int is_initialized() const
    { return 0 != m_scratchSpace && 0 != m_scratchFlags ; }

  void initialize( int cuda_device_id , cudaStream_t stream = 0 );
  void finalize();

  void print_configuration( std::ostream & ) const ;

  void fence() const ;

  ~CudaInternal();

  CudaInternal()
    : m_cudaDev( -1 )
    , m_cudaArch( -1 )
    , m_multiProcCount( 0 )
    , m_maxWarpCount( 0 )
    , m_maxBlock( 0 )
    , m_maxSharedWords( 0 )
    , m_maxConcurrency( 0 )
    , m_shmemPerSM( 0 )
    , m_maxShmemPerBlock( 0 )
    , m_regsPerSM( 0 )
    , m_maxBlocksPerSM( 0 )
    , m_maxThreadsPerSM( 0 )
    , m_maxThreadsPerBlock( 0 )
    , m_scratchSpaceCount( 0 )
    , m_scratchFlagsCount( 0 )
    , m_scratchUnifiedCount( 0 )
    , m_scratchFunctorSize( 0 )
    , m_scratchUnifiedSupported( 0 )
    , m_streamCount( 0 )
    , m_scratchSpace( 0 )
    , m_scratchFlags( 0 )
    , m_scratchUnified( 0 )
    , m_scratchFunctor( 0 )
    , m_scratchConcurrentBitset( 0 )
    , m_stream( 0 )
    {}

  size_type * scratch_space( const size_type size ) const ;
  size_type * scratch_flags( const size_type size ) const ;
  size_type * scratch_unified( const size_type size ) const ;
  size_type * scratch_functor( const size_type size ) const ;
};

} // Namespace Impl
} // Namespace Kokkos
#endif
