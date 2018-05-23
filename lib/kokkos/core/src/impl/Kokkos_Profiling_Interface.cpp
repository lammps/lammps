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

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_PROFILING)

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <cstring>

namespace Kokkos {
namespace Profiling {

static initFunction initProfileLibrary = nullptr;
static finalizeFunction finalizeProfileLibrary = nullptr;

static beginFunction beginForCallee = nullptr;
static beginFunction beginScanCallee = nullptr;
static beginFunction beginReduceCallee = nullptr;
static endFunction endForCallee = nullptr;
static endFunction endScanCallee = nullptr;
static endFunction endReduceCallee = nullptr;

static pushFunction pushRegionCallee = nullptr;
static popFunction popRegionCallee = nullptr;

static allocateDataFunction allocateDataCallee = nullptr;
static deallocateDataFunction deallocateDataCallee = nullptr;

static beginDeepCopyFunction beginDeepCopyCallee = nullptr;
static endDeepCopyFunction endDeepCopyCallee = nullptr;

static createProfileSectionFunction createSectionCallee = nullptr;
static startProfileSectionFunction startSectionCallee = nullptr;
static stopProfileSectionFunction stopSectionCallee = nullptr;
static destroyProfileSectionFunction destroySectionCallee = nullptr;

static profileEventFunction profileEventCallee = nullptr;

SpaceHandle::SpaceHandle(const char* space_name) {
  strncpy(name,space_name,64);
}

bool profileLibraryLoaded() {
  return (nullptr != initProfileLibrary);
}

void beginParallelFor(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
  if(nullptr != beginForCallee) {
    Kokkos::fence();
    (*beginForCallee)(kernelPrefix.c_str(), devID, kernelID);
  }
}

void endParallelFor(const uint64_t kernelID) {
  if(nullptr != endForCallee) {
    Kokkos::fence();
    (*endForCallee)(kernelID);
  }
}

void beginParallelScan(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
  if(nullptr != beginScanCallee) {
    Kokkos::fence();
    (*beginScanCallee)(kernelPrefix.c_str(), devID, kernelID);
  }
}

void endParallelScan(const uint64_t kernelID) {
  if(nullptr != endScanCallee) {
    Kokkos::fence();
    (*endScanCallee)(kernelID);
  }
}

void beginParallelReduce(const std::string& kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
  if(nullptr != beginReduceCallee) {
    Kokkos::fence();
    (*beginReduceCallee)(kernelPrefix.c_str(), devID, kernelID);
  }
}

void endParallelReduce(const uint64_t kernelID) {
  if(nullptr != endReduceCallee) {
    Kokkos::fence();
    (*endReduceCallee)(kernelID);
  }
}


void pushRegion(const std::string& kName) {
  if( nullptr != pushRegionCallee ) {
    Kokkos::fence();
    (*pushRegionCallee)(kName.c_str());
  }
}

void popRegion() {
  if( nullptr != popRegionCallee ) {
    Kokkos::fence();
    (*popRegionCallee)();
  }
}

void allocateData(const SpaceHandle space, const std::string label, const void* ptr, const uint64_t size) {
  if(nullptr != allocateDataCallee) {
    (*allocateDataCallee)(space,label.c_str(),ptr,size);
  }
}

void deallocateData(const SpaceHandle space, const std::string label, const void* ptr, const uint64_t size) {
  if(nullptr != deallocateDataCallee) {
    (*deallocateDataCallee)(space,label.c_str(),ptr,size);
  }
}

void beginDeepCopy(const SpaceHandle dst_space, const std::string dst_label, const void* dst_ptr,
    const SpaceHandle src_space, const std::string src_label, const void* src_ptr,
    const uint64_t size) {
  if(nullptr != beginDeepCopyCallee) {
    (*beginDeepCopyCallee)(dst_space, dst_label.c_str(), dst_ptr,
                      src_space, src_label.c_str(), src_ptr,
                      size);
  }
}

void endDeepCopy() {
  if(nullptr != endDeepCopyCallee) {
    (*endDeepCopyCallee)();
  }
}

void createProfileSection(const std::string& sectionName, uint32_t* secID) {

	if(nullptr != createSectionCallee) {
		(*createSectionCallee)(sectionName.c_str(), secID);
	}
}

void startSection(const uint32_t secID) {
	if(nullptr != startSectionCallee) {
		(*startSectionCallee)(secID);
	}
}

void stopSection(const uint32_t secID) {
	if(nullptr != stopSectionCallee) {
		(*stopSectionCallee)(secID);
	}
}

void destroyProfileSection(const uint32_t secID) {
	if(nullptr != destroySectionCallee) {
		(*destroySectionCallee)(secID);
	}
}

void markEvent(const std::string& eventName) {
	if(nullptr != profileEventCallee) {
		(*profileEventCallee)(eventName.c_str());
	}
}

void initialize() {

  // Make sure initialize calls happens only once
  static int is_initialized = 0;
  if(is_initialized) return;
  is_initialized = 1;

  void* firstProfileLibrary;

  char* envProfileLibrary  = getenv("KOKKOS_PROFILE_LIBRARY");

  // If we do not find a profiling library in the environment then exit
  // early.
  if( nullptr == envProfileLibrary ) {
    return ;
  }

  char* envProfileCopy = (char*) malloc(sizeof(char) * (strlen(envProfileLibrary) + 1));
  sprintf(envProfileCopy, "%s", envProfileLibrary);

  char* profileLibraryName = strtok(envProfileCopy, ";");

  if( (nullptr != profileLibraryName) && (strcmp(profileLibraryName, "") != 0) ) {
    firstProfileLibrary = dlopen(profileLibraryName, RTLD_NOW | RTLD_GLOBAL);

    if(nullptr == firstProfileLibrary) {
      std::cerr << "Error: Unable to load KokkosP library: " <<
        profileLibraryName << std::endl;
      std::cerr << "dlopen(" << profileLibraryName << ", RTLD_NOW | RTLD_GLOBAL) failed with "
        << dlerror() << '\n';
    } else {
#ifdef KOKKOS_ENABLE_PROFILING_LOAD_PRINT
      std::cout << "KokkosP: Library Loaded: " << profileLibraryName << std::endl;
#endif

      // dlsym returns a pointer to an object, while we want to assign to pointer to function
      // A direct cast will give warnings hence, we have to workaround the issue by casting pointer to pointers.
      auto p1 = dlsym(firstProfileLibrary, "kokkosp_begin_parallel_for");
      beginForCallee = *((beginFunction*) &p1);
      auto p2 = dlsym(firstProfileLibrary, "kokkosp_begin_parallel_scan");
      beginScanCallee = *((beginFunction*) &p2);
      auto p3 = dlsym(firstProfileLibrary, "kokkosp_begin_parallel_reduce");
      beginReduceCallee = *((beginFunction*) &p3);

      auto p4 = dlsym(firstProfileLibrary, "kokkosp_end_parallel_scan");
      endScanCallee = *((endFunction*) &p4);
      auto p5 = dlsym(firstProfileLibrary, "kokkosp_end_parallel_for");
      endForCallee = *((endFunction*) &p5);
      auto p6 = dlsym(firstProfileLibrary, "kokkosp_end_parallel_reduce");
      endReduceCallee = *((endFunction*) &p6);

      auto p7 = dlsym(firstProfileLibrary, "kokkosp_init_library");
      initProfileLibrary = *((initFunction*) &p7);
      auto p8 = dlsym(firstProfileLibrary, "kokkosp_finalize_library");
      finalizeProfileLibrary = *((finalizeFunction*) &p8);

      auto p9 = dlsym(firstProfileLibrary, "kokkosp_push_profile_region");
      pushRegionCallee = *((pushFunction*) &p9);
      auto p10 = dlsym(firstProfileLibrary, "kokkosp_pop_profile_region");
      popRegionCallee = *((popFunction*) &p10);

      auto p11 = dlsym(firstProfileLibrary, "kokkosp_allocate_data");
      allocateDataCallee = *((allocateDataFunction*) &p11);
      auto p12 = dlsym(firstProfileLibrary, "kokkosp_deallocate_data");
      deallocateDataCallee = *((deallocateDataFunction*) &p12);

      auto p13 = dlsym(firstProfileLibrary, "kokkosp_begin_deep_copy");
      beginDeepCopyCallee = *((beginDeepCopyFunction*) &p13);
      auto p14 = dlsym(firstProfileLibrary, "kokkosp_end_deep_copy");
      endDeepCopyCallee = *((endDeepCopyFunction*) &p14);
      
      auto p15 = dlsym(firstProfileLibrary, "kokkosp_create_profile_section");
      createSectionCallee = *((createProfileSectionFunction*) &p15);
      auto p16 = dlsym(firstProfileLibrary, "kokkosp_start_profile_section");
      startSectionCallee = *((startProfileSectionFunction*) &p16);
      auto p17 = dlsym(firstProfileLibrary, "kokkosp_stop_profile_section");
      stopSectionCallee = *((stopProfileSectionFunction*) &p17);      
      auto p18 = dlsym(firstProfileLibrary, "kokkosp_destroy_profile_section");
      destroySectionCallee = *((destroyProfileSectionFunction*) &p18);
      
      auto p19 = dlsym(firstProfileLibrary, "kokkosp_profile_event");
      profileEventCallee = *((profileEventFunction*) &p19);
    }
  }

  if(nullptr != initProfileLibrary) {
    (*initProfileLibrary)(0,
        (uint64_t) KOKKOSP_INTERFACE_VERSION,
        (uint32_t) 0,
        nullptr);
  }

  free(envProfileCopy);
}

void finalize() {
  // Make sure finalize calls happens only once
  static int is_finalized = 0;
  if(is_finalized) return;
  is_finalized = 1;

  if(nullptr != finalizeProfileLibrary) {
    (*finalizeProfileLibrary)();

    // Set all profile hooks to nullptr to prevent
    // any additional calls. Once we are told to
    // finalize, we mean it
    initProfileLibrary = nullptr;
    finalizeProfileLibrary = nullptr;

    beginForCallee = nullptr;
    beginScanCallee = nullptr;
    beginReduceCallee = nullptr;
    endScanCallee = nullptr;
    endForCallee = nullptr;
    endReduceCallee = nullptr;

    pushRegionCallee = nullptr;
    popRegionCallee = nullptr;

    allocateDataCallee = nullptr;
    deallocateDataCallee = nullptr;

    beginDeepCopyCallee = nullptr;
    endDeepCopyCallee = nullptr;
    
    createSectionCallee = nullptr;
	startSectionCallee = nullptr;
	stopSectionCallee = nullptr;
	destroySectionCallee = nullptr;

	profileEventCallee = nullptr;
  }
}
}
}

#else

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <cstring>

namespace Kokkos {
namespace Profiling {

bool profileLibraryLoaded() { return false; }


void beginParallelFor(const std::string& , const uint32_t , uint64_t* ) {}
void endParallelFor(const uint64_t ) {}
void beginParallelScan(const std::string& , const uint32_t , uint64_t* ) {}
void endParallelScan(const uint64_t ) {}
void beginParallelReduce(const std::string& , const uint32_t , uint64_t* ) {}
void endParallelReduce(const uint64_t ) {}

void pushRegion(const std::string& ) {}
void popRegion() {}
void createProfileSection(const std::string& , uint32_t* ) {}
void startSection(const uint32_t ) {}
void stopSection(const uint32_t ) {}
void destroyProfileSection(const uint32_t ) {}

void markEvent(const std::string& ) {}

void allocateData(const SpaceHandle , const std::string , const void* , const uint64_t ) {}
void deallocateData(const SpaceHandle , const std::string , const void* , const uint64_t ) {}

void beginDeepCopy(const SpaceHandle , const std::string , const void* , 
    const SpaceHandle , const std::string , const void* ,
    const uint64_t ) {}
void endDeepCopy() {}

void initialize() {}
void finalize() {}

}} // end namespace Kokkos::Profiling

#endif
