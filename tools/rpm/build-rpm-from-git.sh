#!/bin/sh
# automated build script
if ! which rpmbuild > /dev/null 2>&1
then
    echo "Cannot build rpms without 'rpmbuild' installed"
    exit 1
fi

MY_RPM_BUILD_DIR=$(rpmbuild --showrc 2>&1 | sed -n -e 's/^.*: _topdir[	 ]\+\(.*\)$/\1/p')
echo Building RPMs in ${MY_RPM_BUILD_DIR}

for d in "${PWD}" "${PWD%/src}" "${PWD%/tools/rpm}" "$1"
do \
    if test -d "${d}/.git"
    then
        LAMMPS_PATH="${d}"
        break
    fi
done

if test -z "${LAMMPS_PATH}"
then
    echo "'${PWD}' is not a suitable working directory"
    exit 1
fi
pushd "${LAMMPS_PATH}"

for d in SPECS SOURCES RPMS SRPMS BUILD BUILDROOT
do \
  dir="${MY_RPM_BUILD_DIR}/${d}"
  test -d "${dir}" || mkdir -p "${dir}" || exit 2
done

datestr=$(date +%Y%m%d)
sed -e "/^Version/s/\(Version:[ 	]\+\)[0-9].*$/\1${datestr}/" tools/rpm/lammps.spec > ${MY_RPM_BUILD_DIR}/SPECS/lammps.spec

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc/Manual.pdf src lib python bench potentials \
    tools/*.cpp tools/*.f \
    | gzip -9c - > "${MY_RPM_BUILD_DIR}/SOURCES/lammps-current.tar.gz"

rpmbuild --clean --rmsource --rmspec -ba "${MY_RPM_BUILD_DIR}/SPECS/lammps.spec"

popd
