#!/bin/bash -vx

APP_NAME=lammps-gui
DESTDIR=${PWD}/../LAMMPS_GUI

echo "Delete old files, if they exist"
rm -rvf ${DESTDIR} LAMMPS-Linux-amd64.zip

echo "Create staging area for deployment and populate"
DESTDIR=${DESTDIR} cmake --install .  --prefix "/"

echo "Remove libc related files"
rm -v ${DESTDIR}/lib ld*.so ${DESTDIR}/lib/ld*.so.[0-9]
rm -v ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}.so.?
rm -v ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}-[0-9].[0-9]*.so

# get qt dir
QTDIR=$(ldd ${DESTDIR}/bin/lammps-gui | grep libQt5Core | sed -e 's/^.*=> *//' -e 's/libQt5Core.so.*$/qt5/')

echo "Add required dependencies for Qt"
for dir in styles platforms imageformats
do \
    mkdir -p ${DESTDIR}/lib/${dir}
    cp -r  ${QTDIR}/plugins/${dir} ${DESTDIR}/lib/
done

pushd ..
zip -9rv LAMMPS-Linux-amd64.zip LAMMPS_GUI
popd
exit 0
