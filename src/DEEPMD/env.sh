DEEPMD_ROOT=/scratch/gpfs/yifanl/usr/licensed/anaconda3/2020.7/dp-dev-4
TENSORFLOW_INCLUDE_DIRS="/scratch/gpfs/yifanl/usr/licensed/anaconda3/2020.7/dp-dev-4/include;/scratch/gpfs/yifanl/usr/licensed/anaconda3/2020.7/dp-dev-4/include"
TENSORFLOW_LIBRARY_PATH="/scratch/gpfs/yifanl/usr/licensed/anaconda3/2020.7/dp-dev-4/lib;/scratch/gpfs/yifanl/usr/licensed/anaconda3/2020.7/dp-dev-4/lib"

TF_INCLUDE_DIRS=`echo $TENSORFLOW_INCLUDE_DIRS | sed "s/;/ -I/g"`
TF_LIBRARY_PATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -L/g"`
TF_RPATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -Wl,-rpath=/g"`

NNP_INC=" -std=c++11 -DHIGH_PREC   -I$TF_INCLUDE_DIRS -I$DEEPMD_ROOT/include/ "
NNP_PATH=" -L$TF_LIBRARY_PATH -L$DEEPMD_ROOT/lib"
NNP_LIB=" -Wl,--no-as-needed -ldeepmd_op_cuda -ldeepmd_op -ldeepmd_cc -ldeepmd -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$TF_RPATH -Wl,-rpath=$DEEPMD_ROOT/lib"
