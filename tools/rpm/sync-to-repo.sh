#!/bin/sh
# sync data from build directory to repository

if ! which rpmbuild > /dev/null 2>&1
then
    echo "Cannot build rpms without 'rpmbuild' installed"
    exit 1
fi

MYRPM_BUILD_DIR=$(rpmbuild --showrc 2>&1 | sed -n -e 's/^.*: _topdir[	 ]\+\(.*\)$/\1/p')

if [ -n "${MYRPM_REPO_USER}" ] \
&& [ -n "${MYRPM_REPO_HOST}" ] \
&& [ -n "${MYRPM_REPO_DIR}" ]
then
    pushd ${MYRPM_BUILD_DIR}/RPMS
    mv -v *-debuginfo-*.rpm debug/
    createrepo -v --deltas --num-deltas 5 .
    rsync -arpv debug repodata drpms *.rpm \
        ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST}:${MYRPM_REPO_DIR}/
    popd
fi
