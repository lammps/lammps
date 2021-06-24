
#include <inttypes.h>
#include <iostream>

struct Kokkos_Profiling_KokkosPDeviceInfo;

// just get the basename for print_help/parse_args
std::string get_basename(char* cmd, int idx = 0) {
  if (idx > 0) return cmd;
  std::string _cmd = cmd;
  auto _pos        = _cmd.find_last_of('/');
  if (_pos != std::string::npos) return _cmd.substr(_pos + 1);
  return _cmd;
}

struct SpaceHandle {
  char name[64];
};

const int parallel_for_id    = 0;
const int parallel_reduce_id = 1;
const int parallel_scan_id   = 2;

extern "C" void kokkosp_init_library(
    const int /*loadSeq*/, const uint64_t /*interfaceVer*/,
    const uint32_t /*devInfoCount*/,
    Kokkos_Profiling_KokkosPDeviceInfo* /* deviceInfo */) {
  std::cout << "kokkosp_init_library::";
}

extern "C" void kokkosp_finalize_library() {
  std::cout << "kokkosp_finalize_library::";
}

extern "C" void kokkosp_print_help(char* exe) {
  std::cout << "kokkosp_print_help:" << get_basename(exe) << "::";
}

extern "C" void kokkosp_parse_args(int argc, char** argv) {
  std::cout << "kokkosp_parse_args:" << argc;
  for (int i = 0; i < argc; ++i) std::cout << ":" << get_basename(argv[i], i);
  std::cout << "::";
}

extern "C" void kokkosp_begin_parallel_for(const char* name,
                                           const uint32_t devID,
                                           uint64_t* kID) {
  *kID = parallel_for_id;
  std::cout << "kokkosp_begin_parallel_for:" << name << ":" << devID << ":"
            << *kID << "::";
}

extern "C" void kokkosp_end_parallel_for(const uint64_t kID) {
  std::cout << "kokkosp_end_parallel_for:" << kID << "::";
}

extern "C" void kokkosp_begin_parallel_scan(const char* name,
                                            const uint32_t devID,
                                            uint64_t* kID) {
  *kID = parallel_scan_id;
  std::cout << "kokkosp_begin_parallel_scan:" << name << ":" << devID << ":"
            << *kID << "::";
}

extern "C" void kokkosp_end_parallel_scan(const uint64_t kID) {
  std::cout << "kokkosp_end_parallel_scan:" << kID << "::";
}

extern "C" void kokkosp_begin_parallel_reduce(const char* name,
                                              const uint32_t devID,
                                              uint64_t* kID) {
  *kID = parallel_reduce_id;
  std::cout << "kokkosp_begin_parallel_reduce:" << name << ":" << devID << ":"
            << *kID << "::";
}

extern "C" void kokkosp_end_parallel_reduce(const uint64_t kID) {
  std::cout << "kokkosp_end_parallel_reduce:" << kID << "::";
}

extern "C" void kokkosp_push_profile_region(char* regionName) {
  std::cout << "kokkosp_push_profile_region:" << regionName << "::";
}

extern "C" void kokkosp_pop_profile_region() {
  std::cout << "kokkosp_pop_profile_region::";
}

extern "C" void kokkosp_allocate_data(SpaceHandle handle, const char* name,
                                      void* ptr, uint64_t size) {
  std::cout << "kokkosp_allocate_data:" << handle.name << ":" << name << ":"
            << ptr << ":" << size << "::";
}

extern "C" void kokkosp_deallocate_data(SpaceHandle handle, const char* name,
                                        void* ptr, uint64_t size) {
  std::cout << "kokkosp_deallocate_data:" << handle.name << ":" << name << ":"
            << ptr << ":" << size << "::";
}

extern "C" void kokkosp_begin_deep_copy(SpaceHandle dst_handle,
                                        const char* dst_name,
                                        const void* dst_ptr,
                                        SpaceHandle src_handle,
                                        const char* src_name,
                                        const void* src_ptr, uint64_t size) {
  std::cout << "kokkosp_begin_deep_copy:" << dst_handle.name << ":" << dst_name
            << ":" << dst_ptr << ":" << src_handle.name << ":" << src_name
            << ":" << src_ptr << ":" << size << "::";
}

extern "C" void kokkosp_end_deep_copy() {
  std::cout << "kokkosp_end_deep_copy::";
}

uint32_t section_id = 3;
extern "C" void kokkosp_create_profile_section(const char* name,
                                               uint32_t* sec_id) {
  *sec_id = section_id;
  std::cout << "kokkosp_create_profile_section:" << name << ":" << *sec_id
            << "::";
}

extern "C" void kokkosp_start_profile_section(uint32_t sec_id) {
  std::cout << "kokkosp_start_profile_section:" << sec_id << "::";
}

extern "C" void kokkosp_stop_profile_section(uint32_t sec_id) {
  std::cout << "kokkosp_stop_profile_section:" << sec_id << "::";
}
extern "C" void kokkosp_destroy_profile_section(uint32_t sec_id) {
  std::cout << "kokkosp_destroy_profile_section:" << sec_id << "::";
}

extern "C" void kokkosp_profile_event(const char* name) {
  std::cout << "kokkosp_profile_event:" << name << "::";
}
extern "C" void kokkosp_declare_metadata(const char* key, const char* value) {
  std::cout << "kokkosp_declare_metadata:" << key << ":" << value << "::";
}
