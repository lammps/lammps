#!/bin/sh
# sync windows installers to repository


if [ -n "${MINGW_REPO_USER}" ] \
&& [ -n "${MINGW_REPO_HOST}" ] \
&& [ -n "${MINGW_REPO_DIR}" ]
then
    pushd ${HOME}/mingw-cross

    rsync -arpv --delete lammps-32bit-20??????.exe \
        ${MINGW_REPO_USER}@${MINGW_REPO_HOST}:${MINGW_REPO_DIR}/32bit
    rsync -arpv --delete lammps-32bit-mpich2-20??????.exe \
        ${MINGW_REPO_USER}@${MINGW_REPO_HOST}:${MINGW_REPO_DIR}/32bit-mpi
    rsync -arpv --delete lammps-64bit-20??????.exe \
        ${MINGW_REPO_USER}@${MINGW_REPO_HOST}:${MINGW_REPO_DIR}/64bit
    rsync -arpv --delete lammps-64bit-mpich2-20??????.exe \
        ${MINGW_REPO_USER}@${MINGW_REPO_HOST}:${MINGW_REPO_DIR}/64bit-mpi

    ssh ${MINGW_REPO_USER}@${MINGW_REPO_HOST} "cd ${MINGW_REPO_DIR}/../; ./mkhtmlindex.sh *"
    popd
else
    cat <<EOF

Required environment variables to determine the target repository
account, server and location are not fully configured.
MINGW_REPO_USER=${MINGW_REPO_USER}
MINGW_REPO_HOST=${MINGW_REPO_HOST}
MINGW_REPO_DIR=${MINGW_REPO_DIR}

EOF
    exit 1
fi
