#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: $0 <tag or branch>"
    exit 1
fi

# check out desired LAMMPS version
git checkout $1

# make sure we are in the LAMMPS root directory
if [ ! -e README ] || [ ! -e LICENSE ] || [ ! -e SECURITY.md ] || [ ! -e CITATION.cff ]
then
    echo "Must be in the LAMMPS root folder to run this script"
    exit 2
fi

# get version for filename
lammps_release_tag=$(git describe --tags --abbrev=0 | cut -d_ -f2)

# update docenv
make -C doc clean-all
make -C doc html || exit 3
make -C doc pdf || exit 4

# extract initial release tarball from git sources
tarball=lammps-src-${lammps_release_tag}.tar
rm -vf ${tarball} ${tarball}.gz
git archive --output=${tarball} --prefix=lammps-${lammps_release_tag}/ HEAD
mkdir -p lammps-${lammps_release_tag}/doc
mv -v doc/html doc/Manual.pdf lammps-${lammps_release_tag}/doc
tar -rf ${tarball} lammps-${lammps_release_tag}/doc
gzip -9v ${tarball}
rm -r lammps-${lammps_release_tag}
