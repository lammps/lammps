#!/bin/bash
set -euo pipefail

APP_NAME=lammps
DOCDIR="$1"
DESTDIR="${PWD}"

echo "Delete tar files, if they exist"
rm -rvf "${PWD}"/lammps-src-*.tar.gz
pushd "${DOCDIR}/.."
VERSION=$(grep LAMMPS_VERSION src/version.h | sed 's/^.*LAMMPS_VERSION //' | tr -d \" | tr -d \ )
TARNAME=lammps-src-${VERSION}.tar
TARPATH="${DESTDIR}/${TARNAME}"
git archive --output="${TARPATH}" --prefix=lammps-${VERSION}/ HEAD

cd "${DOCDIR}"
make clean-all
make html pdf
rm html/index.html
cp html/Manual.html html/index.html
tar -rf "${TARPATH}" --owner=root --group=root --transform "s,^,lammps-${VERSION}/doc/," html Manual.pdf
popd
gzip -9v "${TARNAME}"
exit 0
