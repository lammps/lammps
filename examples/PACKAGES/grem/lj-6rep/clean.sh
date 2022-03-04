#/bin/bash

for i in $(ls -d [0-9]*)
do
    rm -f $i/final* 
    rm -f $i/log*
    rm -f $i/ent*
    rm -f $i/output
    cp $i/restart.init $i/restart_file
done

echo 1 > lastexchange
cp walker.bkp lastwalker

exit 0
