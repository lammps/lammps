#!/bin/sh

# first argument is required and the LAMMPS executable
if [ "$#" -lt 1 ]
then \
    echo "Usage: ${0} <lammps executable> [input files]"
    exit 0;
fi
export OMP_NUM_THREADS=1
lmp="${1}"
shift
# the list of inputs is optional. fall back to find all files with in.* in the current folder
# pattern for files
pattern='in.*'
inputs=''
if [ "$#" -lt 1 ]
then \
    n=0
    # remove backup files
    rm -vf ${pattern}~
    for arg in ${pattern}
    do \
        [ "${arg}" = "${pattern}" ] && break
        [ ! -f "${arg}" ] && continue
        inputs="${inputs} ${arg}"
        n=$((n + 1))
    done
else
    n=0
    inputs=""
    for arg in "$@"
    do \
        [ ! -f "${arg}" ] && continue
        inputs="${inputs} ${arg}"
        n=$((n + 1))
    done
fi

if [ -z "${inputs}" ]
then \
    echo "No inputs. Nothing to do"
    exit 0
fi
echo "Processing inputs: ${inputs}"
for in in ${inputs}
do \
    echo "Input: ${in}"
    dir=$(dirname "${in}")
    stem=$(basename "${in}" | sed -e 's/^in\.//' -e 's/\.lmp$//' -e 's/~$//')
    date=$(date +%d%b%Y)
    for n in 1 4
    do \
        mpirun -np "${n}" "${lmp}" -in "${in}" -log "${dir}/log.${date}.${stem}.g++.${n}"
    done
done
