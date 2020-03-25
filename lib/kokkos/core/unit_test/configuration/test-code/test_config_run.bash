
SRC_DIR=${KOKKOS_PATH}/core/unit_test/configuration/test-code

# List of parallel device types 
MakeDevices=$1
CMakeDevices=$2
MakeArch=$3
CMakeArch=$4
MakeOptions=$5
CMakeOptions=$6

cd gnu-make
rm -rf *
make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX KokkosCore_config.h &>out
make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX print-cxx-flags &> cxxflags

cd ../cmake
rm -rf *
cmake -DKokkos_SKIP_VALIDATION=ON \
      -DCMAKE_CXX_COMPILER=$CXX \
      $CMakeDevices \
      $CMakeArch \
      $CMakeOptions \
      $SRC_DIR &> config_out
cd ..
grep define gnu-make/KokkosCore_config.h | sort -u &> make_config_defines
grep define cmake/kokkos/KokkosCore_config.h | sort -u &> cmake_config_defines

diff make_config_defines cmake_config_defines &> config_defines_diff
diff_exists=`cat config_defines_diff | wc -l`
if [ $diff_exists -gt 0 ]
then
  echo ""
  echo ""
  echo "Failed #define test"
  echo Make: "make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX KokkosCore_config.h"
  echo CMake: "cmake -DCMAKE_CXX_COMPILER=$CXX $CMakeDevices $CMakeArch $CMakeOptions $SRC_DIR"
  cat config_defines_diff
  echo "Sleeping for 3 seconds if you want to stop and explore..."
  echo ""
  sleep 3
else
  echo ""
  echo ""
  echo "Passed #define test"
  echo Make: "make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX KokkosCore_config.h"
  echo CMake: "cmake -DCMAKE_CXX_COMPILER=$CXX $CMakeDevices $CMakeArch $CMakeOptions $SRC_DIR"
fi

#find because it goes in different locations
#grep out compiler warnings
#head multiple matches
#sed a bunch of stuff to clean up cmake garbage
#awk trim whitespace
#awk print each on new line
#grep remove empty lines
#grep don't consider -std flags in the comparison
#sort and print unique flags
find cmake/kokkos -name KokkosTargets.cmake -exec grep -h INTERFACE_COMPILE_OPTIONS {} \; \
  | grep -v skew \
  | head -n 1 \
  | sed 's/INTERFACE_COMPILE_OPTIONS//g' \
  | sed 's/;/ /g' \
  | sed 's/"//g' \
  | sed 's/\\$<\\$<//g' \
  | sed 's/COMPILE_LANGUAGE:CXX>://g' \
  | sed 's/> / /g' \
  | sed 's/>$//g' \
  | awk '{$1=$1;print}' \
  | awk -v RS=" " '{print}' \
  | grep -v -e '^$' \
  | grep -v '\-std' \
  | sort | uniq > cmake_cxx_flags

#-I flags and -std= flags are not part of CMake's compile options
#that's fine, let's ignore thse below
#redundant lines - tail the last one
#awk print each on new line
#grep out blank lines
#grep out include flags
#grep out -std flags
#sort and print unique flags
tail -n 1 gnu-make/cxxflags \
  | awk -v RS=" " '{print}' \
  | grep -v -e '^$' \
  | grep -v '\-I' \
  | grep -v '\-std=' \
  | grep -v 'gcc-toolchain' \
  | sort | uniq > gnu_make_cxx_flags
diff gnu_make_cxx_flags cmake_cxx_flags &> config_cxxflags_diff
diff_exists=`cat config_cxxflags_diff | wc -l`

if [ $diff_exists -gt 0 ]
then
  echo ""
  echo ""
  echo "Failed CXXFLAGS test"
  echo Make: "make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX KokkosCore_config.h"
  echo CMake: "cmake -DCMAKE_CXX_COMPILER=$CXX $CMakeDevices $CMakeArch $CMakeOptions $SRC_DIR"
  cat config_cxxflags_diff
  echo "Sleeping for 3 seconds if you want to stop and explore..."
  echo ""
  sleep 3
else
  echo ""
  echo ""
  echo "Passed CXXFLAGS test"
  echo Make: "make -f ${SRC_DIR}/Makefile KOKKOS_DEVICES=$MakeDevices KOKKOS_ARCH=$MakeArch $MakeOptions CXX=$CXX KokkosCore_config.h"
  echo CMake: "cmake -DCMAKE_CXX_COMPILER=$CXX $CMakeDevices $CMakeArch $CMakeOptions $SRC_DIR"
fi

