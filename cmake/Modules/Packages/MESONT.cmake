

set(MESONT_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/MESONT)
include(StyleHeaderUtils)
include(CheckLanguage)

# always include C++-only sources
file(GLOB MESONT_SOURCES ${CONFIGURE_DEPENDS} ${MESONT_SOURCES_DIR}/[^.]*mesocnt*.cpp)
file(GLOB MESONT_HEADERS ${CONFIGURE_DEPENDS} ${MESONT_SOURCES_DIR}/[^.]*mesocnt*.h)
# remove derived class when base class is not available
if(NOT PKG_MOLECULE)
  list(REMOVE_ITEM MESONT_SOURCES ${MESONT_SOURCES_DIR}/bond_mesocnt.cpp)
  list(REMOVE_ITEM MESONT_HEADERS ${MESONT_SOURCES_DIR}/bond_mesocnt.h)
endif()

# include styles dependent on Fortran library only when Fortran is available.
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
  file(GLOB MESONT_LIB_SOURCES ${CONFIGURE_DEPENDS} ${LAMMPS_LIB_SOURCE_DIR}/mesont/[^.]*.f90)
  add_library(mesont STATIC ${MESONT_LIB_SOURCES})
  set_target_properties(mesont PROPERTIES OUTPUT_NAME lammps_mesont${LAMMPS_MACHINE})
  target_link_libraries(lammps PRIVATE mesont)

  list(APPEND MESONT_SOURCES ${MESONT_SOURCES_DIR}/pair_mesont_tpm.cpp)
  list(APPEND MESONT_SOURCES ${MESONT_SOURCES_DIR}/atom_vec_mesont.cpp)
  list(APPEND MESONT_SOURCES ${MESONT_SOURCES_DIR}/compute_mesont.cpp)
  list(APPEND MESONT_HEADERS ${MESONT_SOURCES_DIR}/pair_mesont_tpm.h)
  list(APPEND MESONT_HEADERS ${MESONT_SOURCES_DIR}/atom_vec_mesont.h)
  list(APPEND MESONT_HEADERS ${MESONT_SOURCES_DIR}/compute_mesont.h)
endif()

# check for package files in src directory due to old make system
DetectBuildSystemConflict(${LAMMPS_SOURCE_DIR} ${MESONT_SOURCES} ${MESONT_HEADERS})

# manually register style headers
get_property(alist GLOBAL PROPERTY ANGLE)
get_property(blist GLOBAL PROPERTY BOND)
get_property(clist GLOBAL PROPERTY COMPUTE)
get_property(plist GLOBAL PROPERTY PAIR)
get_property(vlist GLOBAL PROPERTY ATOM_VEC)
foreach(fname ${MESONT_HEADERS})
  file(STRINGS ${fname} is_style LIMIT_COUNT 1 REGEX ANGLE_CLASS)
  if(is_style)
    list(APPEND alist ${fname})
  endif()
  file(STRINGS ${fname} is_style LIMIT_COUNT 1 REGEX BOND_CLASS)
  if(is_style)
    list(APPEND blist ${fname})
  endif()
  file(STRINGS ${fname} is_style LIMIT_COUNT 1 REGEX COMPUTE_CLASS)
  if(is_style)
    list(APPEND clist ${fname})
  endif()
  file(STRINGS ${fname} is_style LIMIT_COUNT 1 REGEX PAIR_CLASS)
  if(is_style)
    list(APPEND plist ${fname})
  endif()
  file(STRINGS ${fname} is_style LIMIT_COUNT 1 REGEX ATOM_CLASS)
  if(is_style)
    list(APPEND vlist ${fname})
  endif()
endforeach()
set_property(GLOBAL PROPERTY ANGLE "${alist}")
set_property(GLOBAL PROPERTY BOND "${blist}")
set_property(GLOBAL PROPERTY COMPUTE "${clist}")
set_property(GLOBAL PROPERTY PAIR "${plist}")
set_property(GLOBAL PROPERTY ATOM_VEC "${vlist}")

target_sources(lammps PRIVATE ${MESONT_SOURCES})
target_include_directories(lammps PRIVATE ${MESONT_SOURCES_DIR})

RegisterPackages(${MESONT_SOURCES_DIR})
