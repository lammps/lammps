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

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_THREADS)

#include <Kokkos_Core_fwd.hpp>

/* Standard C++ libraries */

#include <cstdlib>
#include <string>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <mutex>

#include <Kokkos_Threads.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

std::mutex host_internal_cppthread_mutex;

// std::thread compatible driver.
// Recovery from an exception would require constant intra-thread health
// verification; which would negatively impact runtime.  As such simply
// abort the process.

void internal_cppthread_driver() {
  try {
    ThreadsExec::driver();
  } catch (const std::exception& x) {
    std::cerr << "Exception thrown from worker thread: " << x.what()
              << std::endl;
    std::cerr.flush();
    std::abort();
  } catch (...) {
    std::cerr << "Exception thrown from worker thread" << std::endl;
    std::cerr.flush();
    std::abort();
  }
}

}  // namespace

//----------------------------------------------------------------------------
// Spawn a thread

void ThreadsExec::spawn() {
  std::thread t(internal_cppthread_driver);
  t.detach();
}

//----------------------------------------------------------------------------

bool ThreadsExec::is_process() {
  static const std::thread::id master_pid = std::this_thread::get_id();

  return master_pid == std::this_thread::get_id();
}

void ThreadsExec::global_lock() { host_internal_cppthread_mutex.lock(); }

void ThreadsExec::global_unlock() { host_internal_cppthread_mutex.unlock(); }

//----------------------------------------------------------------------------

void ThreadsExec::wait_yield(volatile int& flag, const int value) {
  while (value == flag) {
    std::this_thread::yield();
  }
}

}  // namespace Impl
}  // namespace Kokkos

#else
void KOKKOS_CORE_SRC_THREADS_EXEC_BASE_PREVENT_LINK_ERROR() {}
#endif /* end #if defined( KOKKOS_ENABLE_THREADS ) */
