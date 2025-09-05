#!/bin/bash

APP_NAME=lammps-gui
DESTDIR=${PWD}/../LAMMPS_GUI
VERSION="$1"

echo "Delete old files, if they exist"
rm -rf ${DESTDIR} LAMMPS-Linux-x86_64-GUI-*.tar.gz

echo "Create staging area for deployment and populate"
DESTDIR=${DESTDIR} cmake --install .  --prefix "/"
cp lammps-gui_build-prefix/bin/lammps-gui ${DESTDIR}/bin/

echo "Remove debug info"
for s in ${DESTDIR}/bin/* ${DESTDIR}/lib/liblammps*
do \
        test -f $s && strip --strip-debug $s
done

echo "Remove libc, gcc, and X11 related shared libs"
rm -f ${DESTDIR}/lib/ld*.so ${DESTDIR}/lib/ld*.so.[0-9]
rm -f ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}.so.?
rm -f ${DESTDIR}/lib/lib{c,dl,rt,m,pthread}-[0-9].[0-9]*.so
rm -f ${DESTDIR}/lib/libX* ${DESTDIR}/lib/libxcb*
rm -f ${DESTDIR}/lib/libgcc_s*
rm -f ${DESTDIR}/lib/libstdc++*

# get qt dir
QTDIR=$(ldd ${DESTDIR}/bin/lammps-gui | grep libQt.Core | sed -e 's/^.*=> *//' -e 's/libQt\(.\)Core.so.*$/qt\1/')
cat > ${DESTDIR}/bin/qt.conf <<EOF
[Paths]
Plugins = ../qtplugins
EOF

# platform plugin
mkdir -p ${DESTDIR}/qtplugins/platforms
cp ${QTDIR}/plugins/platforms/libqxcb.so ${DESTDIR}/qtplugins/platforms

# get platform plugin dependencies
QTDEPS=$(LD_LIBRARY_PATH=${DESTDIR}/lib ldd ${QTDIR}/plugins/platforms/libqxcb.so | grep -v ${DESTDIR} | grep libQt[56] | sed -e 's/^.*=> *//' -e 's/\(libQt[56].*.so.*\) .*$/\1/')
for dep in ${QTDEPS}
do \
    cp ${dep} ${DESTDIR}/lib
done

echo "Add additional plugins for Qt"
for dir in styles imageformats
do \
    cp -r  ${QTDIR}/plugins/${dir} ${DESTDIR}/qtplugins/
done

# get imageplugin dependencies
for s in ${DESTDIR}/qtplugins/imageformats/*.so
do \
    QTDEPS=$(LD_LIBRARY_PATH=${DESTDIR}/lib ldd $s | grep -v ${DESTDIR} | grep -E '(libQt.|jpeg)' | sed -e 's/^.*=> *//' -e 's/\(lib.*.so.*\) .*$/\1/')
    for dep in ${QTDEPS}
    do \
        cp ${dep} ${DESTDIR}/lib
    done
done

echo "Set up wrapper script"
MYDIR=$(dirname "$0")
cp ${MYDIR}/xdg-open ${DESTDIR}/bin
cp ${MYDIR}/linux_wrapper.sh ${DESTDIR}/bin
for s in ${DESTDIR}/bin/*
do \
        EXE=$(basename $s)
        test ${EXE} = linux_wrapper.sh && continue
        test ${EXE} = qt.conf && continue
        test ${EXE} = xdg-open && continue
        ln -s bin/linux_wrapper.sh ${DESTDIR}/${EXE}
done

pushd ..
tar -czvvf LAMMPS-Linux-x86_64-GUI-${VERSION}.tar.gz LAMMPS_GUI
popd
mv -v ../LAMMPS-Linux-x86_64-GUI-${VERSION}.tar.gz .

echo "Cleanup dir"
rm -r ${DESTDIR}
exit 0
