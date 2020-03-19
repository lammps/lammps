/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <ostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <impl/Kokkos_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void host_abort( const char * const message )
{
  fwrite(message,1,strlen(message),stderr);
  fflush(stderr);
  ::abort();
}

void throw_runtime_exception( const std::string & msg )
{
  std::ostringstream o ;
  o << msg ;
  traceback_callstack( o );
  throw std::runtime_error( o.str() );
}


std::string human_memory_size(size_t arg_bytes)
{
  double bytes = arg_bytes;
  const double K = 1024;
  const double M = K*1024;
  const double G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << std::setprecision(4) << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << std::setprecision(4) << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << std::setprecision(4) << bytes << " M";
  } else {
    bytes /= G;
    out << std::setprecision(4) << bytes << " G";
  }
  return out.str();
}

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __GNUC__ ) && defined( ENABLE_TRACEBACK )

/*  This is only known to work with GNU C++
 *  Must be compiled with '-rdynamic'
 *  Must be linked with   '-ldl'
 */

/* Print call stack into an error stream,
 * so one knows in which function the error occurred.
 *
 * Code copied from:
 *   http://stupefydeveloper.blogspot.com/2008/10/cc-call-stack.html
 *
 * License on this site:
 *   This blog is licensed under a
 *   Creative Commons Attribution-Share Alike 3.0 Unported License.
 *
 *   http://creativecommons.org/licenses/by-sa/3.0/
 *
 * Modified to output to std::ostream.
 */
#include <signal.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <dlfcn.h>

#include <cstdlib>

namespace Kokkos {
namespace Impl {

void traceback_callstack( std::ostream & msg )
{
  using namespace abi;

  enum { MAX_DEPTH = 32 };

  void *trace[MAX_DEPTH];
  Dl_info dlinfo;

  int status;

  int trace_size = backtrace(trace, MAX_DEPTH);

  msg << std::endl << "Call stack {" << std::endl ;

  for (int i=1; i<trace_size; ++i)
  {
    if(!dladdr(trace[i], &dlinfo))
        continue;

    const char * symname = dlinfo.dli_sname;

    char * demangled = __cxa_demangle(symname, NULL, 0, &status);

    if ( status == 0 && demangled ) {
      symname = demangled;
    }

    if ( symname && *symname != 0 ) {
      msg << "  object: " << dlinfo.dli_fname
          << " function: " << symname
          << std::endl ;
    }

    if ( demangled ) {
        free(demangled);
    }
  }
  msg << "}" ;
}

}
}

#else

namespace Kokkos {
namespace Impl {

void traceback_callstack( std::ostream & msg )
{
  msg << std::endl << "Traceback functionality not available" << std::endl ;
}

}
}

#endif

