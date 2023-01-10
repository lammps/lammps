# run in a loop to test for rare segfaults
for (( ; ; ))
do 
    ~/lammps/build-msmeam/lmp -in test3.in
    # terminate loop if seg fault
    if [[ $? -eq 139 ]]; then 
        exit 1
    fi
done