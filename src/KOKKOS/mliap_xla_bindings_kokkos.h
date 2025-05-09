#ifdef WITH_JAX

#include "pair_mliap_kokkos.h"
#include "xla/ffi/api/c_api.h"    // https://github.com/openxla/xla/blob/main/xla/ffi/api/c_api.h
#include "xla/ffi/api/ffi.h"      // https://github.com/openxla/xla/blob/main/xla/ffi/api/ffi.h
namespace ffi = xla::ffi;

//! TODO: add FP64 support

template <ffi::DataType T>
ffi::Error forward_host(cudaStream_t stream, long pair_handle, ffi::Buffer<T> copy_from,
                        ffi::Result<ffi::Buffer<T>> copy_to)
{
  auto dims = copy_from.dimensions();
  if (dims.size() == 0) { return ffi::Error::InvalidArgument("copy_from has no dimensions"); }
  int vec_len = copy_from.element_count() / dims[0];

  Kokkos::fence();

  LAMMPS_NS::PairMLIAPKokkos<LMPDeviceType> *data =
      (LAMMPS_NS::PairMLIAPKokkos<LMPDeviceType> *) pair_handle;
  data->forward_comm(copy_from.typed_data(), copy_to->typed_data(), vec_len);

  Kokkos::fence();

  return ffi::Error::Success();
}
XLA_FFI_DEFINE_HANDLER_SYMBOL(xla_forward_exchange_fp32, forward_host<ffi::DataType::F32>,
                              ffi::Ffi::Bind()
                                  .Ctx<ffi::PlatformStream<cudaStream_t>>()
                                  .Attr<long>("pair_handle")
                                  .Arg<ffi::Buffer<ffi::DataType::F32>>()
                                  .Ret<ffi::Buffer<ffi::DataType::F32>>());
XLA_FFI_DEFINE_HANDLER_SYMBOL(xla_forward_exchange_fp64, forward_host<ffi::DataType::F64>,
                              ffi::Ffi::Bind()
                                  .Ctx<ffi::PlatformStream<cudaStream_t>>()
                                  .Attr<long>("pair_handle")
                                  .Arg<ffi::Buffer<ffi::DataType::F64>>()
                                  .Ret<ffi::Buffer<ffi::DataType::F64>>());

template <ffi::DataType T>
ffi::Error reverse_host(cudaStream_t stream, long pair_handle, ffi::Buffer<T> copy_from,
                        ffi::Result<ffi::Buffer<T>> copy_to)
{
  auto dims = copy_from.dimensions();
  if (dims.size() == 0) { return ffi::Error::InvalidArgument("copy_from has no dimensions"); }
  int vec_len = copy_from.element_count() / dims[0];

  Kokkos::fence();

  LAMMPS_NS::PairMLIAPKokkos<LMPDeviceType> *data =
      (LAMMPS_NS::PairMLIAPKokkos<LMPDeviceType> *) pair_handle;
  data->reverse_comm(copy_from.typed_data(), copy_to->typed_data(), vec_len);

  Kokkos::fence();

  return ffi::Error::Success();
}
XLA_FFI_DEFINE_HANDLER_SYMBOL(xla_reverse_exchange_fp32, reverse_host<ffi::DataType::F32>,
                              ffi::Ffi::Bind()
                                  .Ctx<ffi::PlatformStream<cudaStream_t>>()
                                  .Attr<long>("pair_handle")
                                  .Arg<ffi::Buffer<ffi::DataType::F32>>()
                                  .Ret<ffi::Buffer<ffi::DataType::F32>>());
XLA_FFI_DEFINE_HANDLER_SYMBOL(xla_reverse_exchange_fp64, reverse_host<ffi::DataType::F64>,
                              ffi::Ffi::Bind()
                                  .Ctx<ffi::PlatformStream<cudaStream_t>>()
                                  .Attr<long>("pair_handle")
                                  .Arg<ffi::Buffer<ffi::DataType::F64>>()
                                  .Ret<ffi::Buffer<ffi::DataType::F64>>());

#define EXPORT_XLA_BINDINGS(T) \
  extern "C" void *T##_ptr()   \
  {                            \
    return (void *) T;         \
  }

#else

#define EXPORT_XLA_BINDINGS(T) \
  extern "C" void *T##_ptr()   \
  {                            \
    return nullptr;            \
  }

#endif    // WITH_JAX

EXPORT_XLA_BINDINGS(xla_forward_exchange_fp32)
EXPORT_XLA_BINDINGS(xla_reverse_exchange_fp32)
EXPORT_XLA_BINDINGS(xla_forward_exchange_fp64)
EXPORT_XLA_BINDINGS(xla_reverse_exchange_fp64)
