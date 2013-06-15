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
    # move debuginfo rpms, if present
    for f in *-debuginfo-*.rpm
    do \
        test -f $f && mv -v $f debug/
    done
    createrepo -v --deltas --num-deltas 2 --max-delta-rpm-size 30000000 .
    rsync -arpv debug repodata drpms *.rpm \
        ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST}:${MYRPM_REPO_DIR}/
    # we use the source rpm on suse, since they do not have a "dist" tag.
    if grep -q -i suse /etc/issue && test $(uname -i) = i386
    then
        cd ../SRPMS
        rsync -arpv *.rpm \
            ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST}:${MYRPM_REPO_DIR}/../../../source/
    fi
    ssh ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST} "cd ${MYRPM_REPO_DIR}/../../../; ./mkindexhtml.sh"
    popd
else
    cat <<EOF

Required environment variables to determine the target repository
account, server and location are not fully configured.
MYRPM_REPO_USER=${MYRPM_REPO_USER}
MYRPM_REPO_HOST=${MYRPM_REPO_HOST}
MYRPM_REPO_DIR=${MYRPM_REPO_DIR}

EOF
    exit 1
fi
