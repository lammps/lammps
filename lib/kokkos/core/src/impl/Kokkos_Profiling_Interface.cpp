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
 // Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
 //
 // ************************************************************************
 //@HEADER
 */

#include <impl/Kokkos_Profiling_Interface.hpp>

#ifdef KOKKOSP_ENABLE_PROFILING
#include <string.h>

namespace Kokkos {
  namespace Experimental {
    bool profileLibraryLoaded() {
       	return (NULL != initProfileLibrary);
    }

    void beginParallelFor(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
        if(NULL != beginForCallee) {
            Kokkos::fence();
            (*beginForCallee)(kernelPrefix.c_str(), devID, kernelID);
        }
    };

    void endParallelFor(const uint64_t kernelID) {
        if(NULL != endForCallee) {
            Kokkos::fence();
            (*endForCallee)(kernelID);
        }
    };

    void beginParallelScan(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
        if(NULL != beginScanCallee) {
            Kokkos::fence();
            (*beginScanCallee)(kernelPrefix.c_str(), devID, kernelID);
        }
    };

    void endParallelScan(const uint64_t kernelID) {
        if(NULL != endScanCallee) {
            Kokkos::fence();
            (*endScanCallee)(kernelID);
        }
    };
    
    void beginParallelReduce(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
        if(NULL != beginReduceCallee) {
            Kokkos::fence();
            (*beginReduceCallee)(kernelPrefix.c_str(), devID, kernelID);
        }
    };
    
    void endParallelReduce(const uint64_t kernelID) {
        if(NULL != endReduceCallee) {
            Kokkos::fence();
            (*endReduceCallee)(kernelID);
        }
    };
    
    void initialize() {
        void* firstProfileLibrary;

        char* envProfileLibrary  = getenv("KOKKOS_PROFILE_LIBRARY");

	// If we do not find a profiling library in the environment then exit
	// early.
	if( NULL == envProfileLibrary ) {
		return ;
	}

	char* profileLibraryName = strtok(envProfileLibrary, ";");

        if( (NULL != profileLibraryName) && (strcmp(profileLibraryName, "") != 0) ) {
            firstProfileLibrary = dlopen(profileLibraryName, RTLD_NOW | RTLD_GLOBAL);

            if(NULL == firstProfileLibrary) {
                std::cerr << "Error: Unable to load KokkosP library: " <<
                profileLibraryName << std::endl;
            } else {
                std::cout << "KokkosP: Library Loaded: " << profileLibraryName << std::endl;

                beginForCallee = (beginFunction) dlsym(firstProfileLibrary, "kokkosp_begin_parallel_for");
                beginScanCallee = (beginFunction) dlsym(firstProfileLibrary, "kokkosp_begin_parallel_scan");
                beginReduceCallee = (beginFunction) dlsym(firstProfileLibrary, "kokkosp_begin_parallel_reduce");

                endScanCallee = (endFunction) dlsym(firstProfileLibrary, "kokkosp_end_parallel_scan");
                endForCallee = (endFunction) dlsym(firstProfileLibrary, "kokkosp_end_parallel_for");
                endReduceCallee = (endFunction) dlsym(firstProfileLibrary, "kokkosp_end_parallel_reduce");

                initProfileLibrary = (initFunction) dlsym(firstProfileLibrary, "kokkosp_init_library");
                finalizeProfileLibrary = (finalizeFunction) dlsym(firstProfileLibrary, "kokkosp_finalize_library");
            }
        }

        if(NULL != initProfileLibrary) {
            (*initProfileLibrary)(0,
		(uint64_t) KOKKOSP_INTERFACE_VERSION,
		(uint32_t) 0,
		NULL);
        }
    };

    void finalize() {
      if(NULL != finalizeProfileLibrary) {
        (*finalizeProfileLibrary)();

        // Set all profile hooks to NULL to prevent
        // any additional calls. Once we are told to
        // finalize, we mean it
        beginForCallee = NULL;
        beginScanCallee = NULL;
        beginReduceCallee = NULL;
        endScanCallee = NULL;
        endForCallee = NULL;
        endReduceCallee = NULL;
        initProfileLibrary = NULL;
        finalizeProfileLibrary = NULL;
      }
    };
  }
}

#endif
