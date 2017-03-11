#!/bin/bash 

. /etc/profile.d/modules.sh

echo "build-dir $1"
echo "backend $2"
echo "module $3"
echo "compiler $4"
echo "cxxflags $5"
echo "architecrure $6"
echo "debug $7"
echo "kokkos-options $8"
echo "kokkos-cuda-options $9"
echo "hwloc $9"

NOW=`date "+%Y%m%d%H%M%S"`
BASEDIR="$1-$NOW"

mkdir $BASEDIR
cd $BASEDIR

module load $2

if [ $9 == "yes" ]; then
if [ $7 == "debug" ]; then
  ../generate_makefile.sh --with-devices=$2 \
  	--compiler=$4 \
	--cxxflags=$5 \
        --arch=$6 \
        --debug \
	--with-options=$8 \
        --with-cuda-options=$9
        --with-hwloc=${HWLOC_ROOT}
else
    ../generate_makefile.sh --with-devices=$2 \
        --compiler=$4 \
        --cxxflags=$5 \
        --arch=$6 \
        --debug \
        --with-options=$8 \
        --with-cuda-options=$9 
        --with-hwloc=${HWLOC_ROOT}
fi
else
if [ $7 == "debug" ]; then
  ../generate_makefile.sh --with-devices=$2 \
        --compiler=$4 \
        --cxxflags=$5 \
        --arch=$6 \
        --debug \
        --with-options=$8 \
        --with-cuda-options=$9
else
    ../generate_makefile.sh --with-devices=$2 \
        --compiler=$4 \
        --cxxflags=$5 \
        --arch=$6 \
        --debug \
        --with-options=$8 \
        --with-cuda-options=$9 
fi
fi


make test
return $?
