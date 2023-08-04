#!/bin/bash -vx

APP_NAME=lammps-gui
DESTDIR=${PWD}/LAMMPS_GUI
SYSROOT="$1"

echo "Delete old files, if they exist"
rm -rvf ${DESTDIR}/LAMMPS_GUI ${DESTDIR}/LAMMPS-Win10-amd64.zip

echo "Create staging area for deployment and populate"
DESTDIR=${DESTDIR} cmake --install .  --prefix "/"

echo "Copying required DLL files"
for dll in $(objdump -p *.exe *.dll | sed -n -e '/DLL Name:/s/^.*DLL Name: *//p' | sort | uniq)
do \
    doskip=0
    for skip in ADVAPI32 CFGMGR32 KERNEL32 msvcrt ole32 PSAPI WSOCK32 liblammps
    do \
        test ${dll} = ${skip}.dll && doskip=1
    done
    test ${doskip} -eq 1 && continue
    test -f ${DESTDIR}/bin/${dll} || cp -v ${SYSROOT}/bin/${dll} ${DESTDIR}/bin
done
    
echo "Copy required Qt plugins"
mkdir -p ${DESTDIR}/qt5plugins
for plugin in imageformats platforms styles
do \
    cp -r ${SYSROOT}/lib/qt5/plugins/${plugin} ${DESTDIR}/qt5plugins/
done

cat > ${DESTDIR}/bin/qt.conf <<EOF
[Paths]
Plugins = ../qt5plugins
EOF
zip -9rv LAMMPS-Win10-amd64.zip LAMMPS_GUI
exit 0
