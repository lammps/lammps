set(temp "#ifndef LMP_GIT_VERSION_H\n#define LMP_GIT_VERSION_H\n")
set(temp_git_commit "(unknown)")
set(temp_git_branch "(unknown)")
set(temp_git_describe "(unknown)")
set(temp_git_info "false")

message(STATUS "Git Directory: ${LAMMPS_DIR}/.git")
if(GIT_FOUND AND EXISTS ${LAMMPS_DIR}/.git)
  set(temp_git_info "true")
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    OUTPUT_VARIABLE temp_git_commit
    ERROR_QUIET
    WORKING_DIRECTORY ${LAMMPS_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE temp_git_branch
    ERROR_QUIET
    WORKING_DIRECTORY ${LAMMPS_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --dirty=-modified --always
    OUTPUT_VARIABLE temp_git_describe
    ERROR_QUIET
    WORKING_DIRECTORY ${LAMMPS_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(temp "${temp}bool LAMMPS_NS::LAMMPS::has_git_info() { return ${temp_git_info}; }\n")
set(temp "${temp}const char *LAMMPS_NS::LAMMPS::git_commit() { return \"${temp_git_commit}\"; }\n")
set(temp "${temp}const char *LAMMPS_NS::LAMMPS::git_branch() { return \"${temp_git_branch}\"; }\n")
set(temp "${temp}const char *LAMMPS_NS::LAMMPS::git_descriptor() { return \"${temp_git_describe}\"; }\n")
set(temp "${temp}#endif\n\n")

message(STATUS "Generating lmpgitversion.h...")

string(REPLACE "\\ " " " LAMMPS_GIT_HEADER "${LAMMPS_STYLE_HEADERS_DIR}/lmpgitversion.h")
file(WRITE "${LAMMPS_GIT_HEADER}.tmp" "${temp}" )
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${LAMMPS_GIT_HEADER}.tmp" "${LAMMPS_GIT_HEADER}")
