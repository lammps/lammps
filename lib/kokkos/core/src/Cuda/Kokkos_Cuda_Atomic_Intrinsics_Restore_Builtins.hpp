/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifdef KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED

#undef __atomic_load_n
#undef __atomic_load
#undef __atomic_store_n
#undef __atomic_store
#undef __atomic_exchange_n
#undef __atomic_exchange
#undef __atomic_compare_exchange_n
#undef __atomic_compare_exchange
#undef __atomic_fetch_add
#undef __atomic_fetch_sub
#undef __atomic_fetch_and
#undef __atomic_fetch_xor
#undef __atomic_fetch_or
#undef __atomic_test_and_set
#undef __atomic_clear
#undef __atomic_always_lock_free
#undef __atomic_is_lock_free
#undef __atomic_thread_fence
#undef __atomic_signal_fence

#undef KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED

#endif // KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED
