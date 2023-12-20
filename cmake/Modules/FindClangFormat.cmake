# Find clang-format
find_program(ClangFormat_EXECUTABLE NAMES clang-format
                                          clang-format-17.0
                                          clang-format-16.0
                                          clang-format-15.0
                                          clang-format-14.0
                                          clang-format-13.0
                                          clang-format-12.0
                                          clang-format-11.0
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

  if(clang_format_version MATCHES "^(Ubuntu |Debian |)clang-format version .*")
    # Arch Linux output:
    # clang-format version 10.0.0
    #
    # Ubuntu 18.04 LTS output:
    # clang-format version 6.0.0-1ubuntu2 (tags/RELEASE_600/final)
    #
    # Ubuntu 20.04 LTS output:
    # clang-format version 10.0.0-4ubuntu1
    #
    # Ubuntu 22.04 LTS output:
    # Ubuntu clang-format version 14.0.0-1ubuntu1
    #
    # Debian 11 output:
    # Debian clang-format version 11.0.1-2
    #
    # Debian 12 output:
    # Debian clang-format version 14.0.6
    #
    # Fedora 36 output:
    # clang-format version 14.0.5 (Fedora 14.0.5-1.fc36)
    string(REGEX REPLACE "^(Ubuntu |Debian |)clang-format version ([0-9.]+).*"
                         "\\2"
                         ClangFormat_VERSION
                         "${clang_format_version}")
  elseif(clang_format_version MATCHES ".*LLVM version .*")
    # CentOS 7 output:
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
