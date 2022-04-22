DEEPMD_ROOT=/home/yifanl/.conda/envs/dpdev4
TENSORFLOW_INCLUDE_DIRS="/home/yifanl/.conda/envs/dpdev4/include;/home/yifanl/.conda/envs/dpdev4/include"
TENSORFLOW_LIBRARY_PATH="/home/yifanl/.conda/envs/dpdev4/lib;/home/yifanl/.conda/envs/dpdev4/lib"

TF_INCLUDE_DIRS=`echo $TENSORFLOW_INCLUDE_DIRS | sed "s/;/ -I/g"`
TF_LIBRARY_PATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -L/g"`
TF_RPATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -Wl,-rpath=/g"`

NNP_INC=" -std=c++14 -DLOW_PREC  -DLAMMPS_VERSION_NUMBER=20210929 -I$TF_INCLUDE_DIRS -I$DEEPMD_ROOT/include/ "
NNP_PATH=" -L$TF_LIBRARY_PATH -L$DEEPMD_ROOT/lib"
NNP_LIB=" -Wl,--no-as-needed -ldeepmd_cc_low -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$TF_RPATH -Wl,-rpath=$DEEPMD_ROOT/lib"
