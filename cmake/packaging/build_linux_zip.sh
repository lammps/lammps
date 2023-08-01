#!/bin/bash -vx

APP_NAME=lammps-gui
DESTDIR=${PWD}/../LAMMPS_GUI

echo "Delete old files, if they exist"
rm -rvf ${DESTDIR} LAMMPS-Linux-amd64.zip

echo "Create staging area for deployment and populate"
DESTDIR=${DESTDIR} cmake --install .  --prefix "/"

echo "Remove debug info"
for s in ${DESTDIR}/bin/* ${DESTDIR}/lib/liblammps*
do \
        strip --strip-debug $s
done

echo "Remove libc, gcc, and X11 related files"
rm -v ${DESTDIR}/lib/ld*.so ${DESTDIR}/lib/ld*.so.[0-9]
rm -v ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}.so.?
rm -v ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}-[0-9].[0-9]*.so
rm -v ${DESTDIR}/lib/libX* ${DESTDIR}/lib/libxcb*

# get qt dir
QTDIR=$(ldd ${DESTDIR}/bin/lammps-gui | grep libQt5Core | sed -e 's/^.*=> *//' -e 's/libQt5Core.so.*$/qt5/')
cat > ${DESTDIR}/bin/qt.conf <<EOF
[Paths]
Plugins = ../qt5plugins
EOF

echo "Add required dependencies for Qt"
for dir in styles platforms imageformats
do \
    mkdir -p ${DESTDIR}/qt5plugins
    cp -r  ${QTDIR}/plugins/${dir} ${DESTDIR}/qt5plugins/
done

pushd ..
zip -9rv LAMMPS-Linux-amd64.zip LAMMPS_GUI
popd
exit 0
