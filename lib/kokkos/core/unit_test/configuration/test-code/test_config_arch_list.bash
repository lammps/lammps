
# List of parallel device types 
HostArch=(SNB HSW SKX KNL)
DeviceArch=(Kepler35 Kepler37 Pascal60 Pascal61 Volta70)
if [ ! -z "$KOKKOS_HOST_ARCH_TEST" ]; then
  export KOKKOS_ARCH_TEST=1
  HostArch=(WSM SNB HSW SKX WSM AMDAVX ARMv80 ARMv81 BDW KNC KNL BGQ Power7 Power8 Power9 Ryzen EPYC ARMv8_ThunderX ARMv8_ThunderX2)
  DeviceArch=()
fi

if [ ! -z "$KOKKOS_DEVICE_ARCH_TEST" ]; then
  export KOKKOS_ARCH_TEST=1
  HostArch=(SNB)
  DeviceArch=(Kepler30 Kepler32 Kepler35 Kepler37 Maxwell50 Maxwell52 Maxwell53 Pascal60 Pascal61 Volta70 Volta72)
fi

MakeDevices=$1
CMakeDevices=$2

SRC_DIR=${KOKKOS_PATH}/core/unit_test/configuration/test-code

for harch in "${HostArch[@]}"
do
  harch_up=`echo $harch | tr a-z A-Z`
  CMAKE_HARCH="-DKokkos_ARCH_${harch_up}=ON"

  if [ "$harch" == "ARMv8_ThunderX2" ]; then
    harch="ARMv8-TX2"
  elif [ "$harch" == "ARMv8_ThunderX" ]; then
    harch="ARMv8-ThunderX"
  fi

  if [ ! -z "$DeviceArch" ]
  then
    for darch in "${DeviceArch[@]}"
    do
      darch_up=`echo $darch | tr a-z A-Z`
      CMAKE_DARCH="-DKokkos_ARCH_${darch_up}=ON"
      ${SRC_DIR}/test_config_options_list.bash "$MakeDevices" "$CMakeDevices" "$harch,$darch" "${CMAKE_HARCH} ${CMAKE_DARCH}"
    done
  else
    ${SRC_DIR}/test_config_options_list.bash "$MakeDevices" "$CMakeDevices" "$harch" "${CMAKE_HARCH}"
  fi 
done

