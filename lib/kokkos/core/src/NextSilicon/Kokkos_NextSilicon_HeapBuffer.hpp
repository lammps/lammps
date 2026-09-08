// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_HEAPBUFFER_HPP
#define KOKKOS_NEXTSILICON_HEAPBUFFER_HPP

#include <NextSilicon/Kokkos_NextSiliconSpace.hpp>

namespace Kokkos::Experimental::Impl {

class NextSiliconHeapBuffer {
  size_t data_size_ = 0;
  void *data_       = nullptr;

  void allocate(const char *label, size_t requested) {
    Kokkos::Experimental::NextSiliconSharedSpace space;
    data_      = space.allocate(label, requested);
    data_size_ = requested;
  }

  void deallocate() {
    if (data_) {
      Kokkos::Experimental::NextSiliconSharedSpace space;
      space.deallocate(data_, data_size_);
      data_      = nullptr;
      data_size_ = 0;
    }
  }

 public:
  NextSiliconHeapBuffer()                                         = default;
  NextSiliconHeapBuffer(const NextSiliconHeapBuffer &)            = delete;
  NextSiliconHeapBuffer &operator=(const NextSiliconHeapBuffer &) = delete;
  NextSiliconHeapBuffer(NextSiliconHeapBuffer &&other) noexcept {
    deallocate();
    data_      = std::exchange(other.data_, nullptr);
    data_size_ = std::exchange(other.data_size_, 0);
  }
  NextSiliconHeapBuffer &operator=(NextSiliconHeapBuffer &&other) noexcept {
    if (this != &other) {
      deallocate();
      data_      = std::exchange(other.data_, nullptr);
      data_size_ = std::exchange(other.data_size_, 0);
    }
    return *this;
  }
  ~NextSiliconHeapBuffer() { deallocate(); }

  std::byte *ensure(const char *label, size_t requested) {
    if (requested > data_size_) {
      deallocate();
      allocate(label, requested);
    }
    return static_cast<std::byte *>(data_);
  }
};
}  // namespace Kokkos::Experimental::Impl

#endif

// HeapBuffer -> NextSiliconSpace
