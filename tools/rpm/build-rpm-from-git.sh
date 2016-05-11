#!/bin/sh
# automated build script
if ! which rpmbuild > /dev/null 2>&1
then
    echo "Cannot build rpms without 'rpmbuild' installed"
    exit 1
fi

MYRPM_BUILD_DIR=$(rpmbuild --showrc 2>&1 | sed -n -e 's/^.*: _topdir[	 ]\+\(.*\)$/\1/p')
echo Building RPMs in ${MYRPM_BUILD_DIR}

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

for d in SPECS SOURCES RPMS SRPMS BUILD BUILDROOT
do \
  dir="${MYRPM_BUILD_DIR}/${d}"
  test -d "${dir}" || mkdir -p "${dir}" || exit 2
done

for d in cache debug drpms repodata
do \
  dir="${MYRPM_BUILD_DIR}/RPMS/${d}"
  test -d "${dir}" || mkdir -p "${dir}" || exit 2
done

pushd "${LAMMPS_PATH}"
datestr=$(date +%Y%m%d)
sed -e "/^Version/s/\(Version:[ 	]\+\)[0-9].*$/\1${datestr}/" tools/rpm/lammps.spec > ${MYRPM_BUILD_DIR}/SPECS/lammps.spec
cp -pv tools/rpm/lammps.sh ${MYRPM_BUILD_DIR}/SOURCES/
cp -pv tools/rpm/lammps.csh ${MYRPM_BUILD_DIR}/SOURCES/

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc/Manual.pdf doc/src/PDF src lib python  \
    examples/{README,ASPHERE,KAPPA,MC,VISCOSITY,dipole,peri,hugoniostat,colloid,crack,friction,msst,obstacle,body,sputter,pour,ELASTIC,neb,ellipse,flow,meam,min,indent,deposit,micelle,shear,srd,dreiding,eim,prd,rigid,COUPLE,peptide,melt,comb,tad,reax,balance,snap,USER/{awpmd,misc,phonon,cg-cmm,fep}} \
    bench potentials tools/*.cpp tools/*.f tools/msi2lmp tools/xmgrace tools/createatoms tools/colvars \
    | gzip -9c - > "${MYRPM_BUILD_DIR}/SOURCES/lammps-current.tar.gz"
popd

# we use the source rpm on suse, since they do not have a "dist" tag.
if grep -q -i suse\ 12 /etc/issue && test $(uname -i) = i386
then
    rpmbuild --clean --rmsource --rmspec -ba "${MYRPM_BUILD_DIR}/SPECS/lammps.spec" || exit 1
else
    rpmbuild --clean --rmsource --rmspec -bb "${MYRPM_BUILD_DIR}/SPECS/lammps.spec" || exit 1
fi

# clean up any unwanted subdirs in the RPMS folder
for ext in i386 i586 i686 x86_64
do \
    if [ -d ${MYRPM_BUILD_DIR}/RPMS/${ext} ]
    then
        mv -v ${MYRPM_BUILD_DIR}/RPMS/${ext}/* ${MYRPM_BUILD_DIR}/RPMS
        rmdir ${MYRPM_BUILD_DIR}/RPMS/${ext}
    fi
done
