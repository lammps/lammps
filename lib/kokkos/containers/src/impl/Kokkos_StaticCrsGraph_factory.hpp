/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_STATICCRSGRAPH_FACTORY_HPP
#define KOKKOS_IMPL_STATICCRSGRAPH_FACTORY_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
template <class DataType, class Arg1Type, class Arg2Type, typename SizeType,
          class Arg3Type>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType,
                               Arg3Type>::HostMirror
create_mirror_view(
    const StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType, Arg3Type>&
        view,
    typename Impl::enable_if<ViewTraits<DataType, Arg1Type, Arg2Type,
                                        Arg3Type>::is_hostspace>::type* = 0) {
  return view;
}
#else
template <class DataType, class Arg1Type, class Arg2Type, class Arg3Type,
          typename SizeType>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type,
                               SizeType>::HostMirror
create_mirror_view(
    const StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type, SizeType>&
        view,
    typename Impl::enable_if<ViewTraits<DataType, Arg1Type, Arg2Type,
                                        Arg3Type>::is_hostspace>::type* = 0) {
  return view;
}
#endif

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
template <class DataType, class Arg1Type, class Arg2Type, typename SizeType,
          class Arg3Type>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType,
                               Arg3Type>::HostMirror
create_mirror(const StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType,
                                   Arg3Type>& view) {
  // Force copy:
  // typedef Impl::ViewAssignment< Impl::ViewDefault > alloc ; // unused
  typedef StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType, Arg3Type>
      staticcrsgraph_type;
#else
template <class DataType, class Arg1Type, class Arg2Type, class Arg3Type,
          typename SizeType>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type,
                               SizeType>::HostMirror
create_mirror(const StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type,
                                   SizeType>& view) {
  // Force copy:
  // typedef Impl::ViewAssignment< Impl::ViewDefault > alloc ; // unused
  typedef StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type, SizeType>
      staticcrsgraph_type;
#endif

  typename staticcrsgraph_type::HostMirror tmp;
  typename staticcrsgraph_type::row_map_type::HostMirror tmp_row_map =
      create_mirror(view.row_map);
  typename staticcrsgraph_type::row_block_type::HostMirror
      tmp_row_block_offsets = create_mirror(view.row_block_offsets);

  // Allocation to match:
  tmp.row_map = tmp_row_map;  // Assignment of 'const' from 'non-const'
  tmp.entries = create_mirror(view.entries);
  tmp.row_block_offsets =
      tmp_row_block_offsets;  // Assignment of 'const' from 'non-const'

  // Deep copy:
  deep_copy(tmp_row_map, view.row_map);
  deep_copy(tmp.entries, view.entries);
  deep_copy(tmp_row_block_offsets, view.row_block_offsets);

  return tmp;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
template <class DataType, class Arg1Type, class Arg2Type, typename SizeType,
          class Arg3Type>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType,
                               Arg3Type>::HostMirror
create_mirror_view(
    const StaticCrsGraph<DataType, Arg1Type, Arg2Type, SizeType, Arg3Type>&
        view,
    typename Impl::enable_if<!ViewTraits<DataType, Arg1Type, Arg2Type,
                                         Arg3Type>::is_hostspace>::type* = 0)
#else
template <class DataType, class Arg1Type, class Arg2Type, class Arg3Type,
          typename SizeType>
inline typename StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type,
                               SizeType>::HostMirror
create_mirror_view(
    const StaticCrsGraph<DataType, Arg1Type, Arg2Type, Arg3Type, SizeType>&
        view,
    typename Impl::enable_if<!ViewTraits<DataType, Arg1Type, Arg2Type,
                                         Arg3Type>::is_hostspace>::type* = 0)
#endif
{
  return create_mirror(view);
}
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template <class StaticCrsGraphType, class InputSizeType>
inline typename StaticCrsGraphType::staticcrsgraph_type create_staticcrsgraph(
    const std::string& label, const std::vector<InputSizeType>& input) {
  typedef StaticCrsGraphType output_type;
  // typedef std::vector< InputSizeType >  input_type ; // unused

  typedef typename output_type::entries_type entries_type;

  typedef View<typename output_type::size_type[],
               typename output_type::array_layout,
               typename output_type::execution_space,
               typename output_type::memory_traits>
      work_type;

  output_type output;

  // Create the row map:

  const size_t length = input.size();

  {
    work_type row_work("tmp", length + 1);

    typename work_type::HostMirror row_work_host = create_mirror_view(row_work);

    size_t sum       = 0;
    row_work_host[0] = 0;
    for (size_t i = 0; i < length; ++i) {
      row_work_host[i + 1] = sum += input[i];
    }

    deep_copy(row_work, row_work_host);

    output.entries = entries_type(label, sum);
    output.row_map = row_work;
  }

  return output;
}

//----------------------------------------------------------------------------

template <class StaticCrsGraphType, class InputSizeType>
inline typename StaticCrsGraphType::staticcrsgraph_type create_staticcrsgraph(
    const std::string& label,
    const std::vector<std::vector<InputSizeType> >& input) {
  typedef StaticCrsGraphType output_type;
  typedef typename output_type::entries_type entries_type;

  static_assert(entries_type::rank == 1, "Graph entries view must be rank one");

  typedef View<typename output_type::size_type[],
               typename output_type::array_layout,
               typename output_type::execution_space,
               typename output_type::memory_traits>
      work_type;

  output_type output;

  // Create the row map:

  const size_t length = input.size();

  {
    work_type row_work("tmp", length + 1);

    typename work_type::HostMirror row_work_host = create_mirror_view(row_work);

    size_t sum       = 0;
    row_work_host[0] = 0;
    for (size_t i = 0; i < length; ++i) {
      row_work_host[i + 1] = sum += input[i].size();
    }

    deep_copy(row_work, row_work_host);

    output.entries = entries_type(label, sum);
    output.row_map = row_work;
  }

  // Fill in the entries:
  {
    typename entries_type::HostMirror host_entries =
        create_mirror_view(output.entries);

    size_t sum = 0;
    for (size_t i = 0; i < length; ++i) {
      for (size_t j = 0; j < input[i].size(); ++j, ++sum) {
        host_entries(sum) = input[i][j];
      }
    }

    deep_copy(output.entries, host_entries);
  }

  return output;
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSARRAY_FACTORY_HPP */
