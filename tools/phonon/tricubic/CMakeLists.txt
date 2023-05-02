
# Support Linux from Ubuntu 20.04LTS onward, CentOS 7.x (with EPEL),
# macOS, MSVC 2019 (=Version 16)
cmake_minimum_required(VERSION 3.10)

# set up project
project(tricubic VERSION 1.1 DESCRIPTION "Tricubic library" LANGUAGES CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# hacks for MSVC to prevent lots of pointless warnings about "unsafe" functions
if(MSVC)
  add_compile_options(/wd4244)
  add_compile_options(/wd4267)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

add_library(tricubic STATIC tricubic.cpp)
target_include_directories(tricubic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
