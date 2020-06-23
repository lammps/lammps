# Find clang-format
find_program(ClangFormat_EXECUTABLE NAMES clang-format
                                          clang-format-10.0
                                          clang-format-9.0
                                          clang-format-8.0
                                          clang-format-7.0
                                          clang-format-6.0
                                     DOC "clang-format executable")
mark_as_advanced(ClangFormat_EXECUTABLE)

if(ClangFormat_EXECUTABLE)
  # find version
  execute_process(COMMAND ${ClangFormat_EXECUTABLE} --version
                  OUTPUT_VARIABLE clang_format_version
                  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)


  if(clang_format_version MATCHES "^clang-format version .*")
    # Arch Linux
    # clang-format version 10.0.0

    # Ubuntu 18.04 LTS Output
    # clang-format version 6.0.0-1ubuntu2 (tags/RELEASE_600/final)
    string(REGEX REPLACE "clang-format version ([0-9.]+).*"
                         "\\1"
                         ClangFormat_VERSION
                         "${clang_format_version}")
  elseif(clang_format_version MATCHES ".*LLVM version .*")
    # CentOS 7 Output
    # LLVM (http://llvm.org/):
    #   LLVM version 3.4.2
    #   Optimized build.
    #   Built Nov  1 2018 (15:06:24).
    #   Default target: x86_64-redhat-linux-gnu
    #   Host CPU: x86-64
    string(REGEX REPLACE ".*LLVM version ([0-9.]+).*"
                         "\\1"
                         ClangFormat_VERSION
                         "${clang_format_version}")
  else()
    set(ClangFormat_VERSION "0.0")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ClangFormat REQUIRED_VARS ClangFormat_EXECUTABLE VERSION_VAR ClangFormat_VERSION)
