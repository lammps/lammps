#!/bin/bash -vx

APP_NAME=lammps-gui
DESTDIR=${PWD}/../LAMMPS_GUI

echo "Delete old files, if they exist"
rm -rvf ${DESTDIR} LAMMPS-Win10-amd64.zip

echo "Create staging area for deployment and populate"
DESTDIR=${DESTDIR} cmake --install .  --prefix "/"

echo "Add required dependencies for Qt"
for dll in Qt5Core.dll Qt5Gui.dll Qt5Widgets.dll
do \
    cp /usr/x86_64-w64-mingw32/sys-root/mingw/bin/${dll} ${DESTDIR}/bin/
done
for dir in styles platforms imageformats
do \
    mkdir -p ${DESTDIR}/${dir}
    cp -r /usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt5/plugins/${dir}/*.dll ${DESTDIR}/${dir}
done

pushd ..
zip -9rv LAMMPS-Win10-amd64.zip LAMMPS_GUI
popd
exit 0
