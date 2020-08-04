now="$(date)" 
printf "Start at %s\n" "$now"

valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file="valgrind.cgheptane.cubic200.log" /mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.cubic200.lmp

now="$(date)" 
printf "Finish at %s\n" "$now"

 