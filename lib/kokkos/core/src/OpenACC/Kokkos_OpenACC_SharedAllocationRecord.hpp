//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_OPENACC_SHARED_ALLOCATION_RECORD_HPP
#define KOKKOS_OPENACC_SHARED_ALLOCATION_RECORD_HPP

#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <openacc.h>

template <>
class Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace,
                                           void>
    : public HostInaccessibleSharedAllocationRecordCommon<
          Kokkos::Experimental::OpenACCSpace> {
 private:
  friend class HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenACCSpace>;
  friend class SharedAllocationRecordCommon<Kokkos::Experimental::OpenACCSpace>;
  friend Kokkos::Experimental::OpenACCSpace;

  using base_t = HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenACCSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  /**\brief  Root record for tracked allocations from this OpenACCSpace
   * instance */
  static RecordBase s_root_record;

  const Kokkos::Experimental::OpenACCSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  template <typename ExecutionSpace>
  SharedAllocationRecord(
      const ExecutionSpace& /*exec_space*/,
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate)
      : SharedAllocationRecord(arg_space, arg_label, arg_alloc_size,
                               arg_dealloc) {}

  SharedAllocationRecord(
      const Kokkos::Experimental::OpenACC& exec_space,
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

  SharedAllocationRecord(
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);
};

#endif
