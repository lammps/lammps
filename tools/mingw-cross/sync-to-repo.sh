#!/bin/sh
# sync windows installers to repository

# function to selectively remove older packages
# and make a symlink to the latest version
prune_exe () {
  dir=$1
  sub=$2
  ref=$(date +%s)
  old=999999999

  for exe in ${dir}/lammps${sub}-20[0-9][0-9]*.exe
  do \
      [ -f ${exe} ] || continue
      # re-set symbolic link to latest entry
      p=$(echo ${exe} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(\.exe\)$@\1@')
      r=$(echo ${exe} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(\.exe\)$@\2@')
      t=$(echo ${exe} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(\.exe\)$@\3@')
      e=$(echo ${exe} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(\.exe\)$@\4@')

      # compute age difference in days
      y=$(echo ${t} | cut -c 1-4)
      m=$(echo ${t} | cut -c 5-6)
      d=$(echo ${t} | cut -c 7-8)
      s=$(date +%s -d "${m}/${d}/${y}")
      age=$(expr \( $ref - $s \) / 86400)

      if [ $age -lt $old ] 
      then
          old=$age
          sym="${r}${t}${e}"
          sto="${p}${r}latest${e}"
      fi

      # after one weeks we keep only one per week.
      unset tmp
      if [ $age -gt 7 ]
      then
          w=$(expr $age / 7)
          eval tmp=\$week${w}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          fi
          eval week${w}=1
      fi

      # after about three months we keep only one per month
      unset tmp
      if [ $age -gt 84 ]
      then
          m=$(expr $age / 28)
          eval tmp=\$month${m}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          fi
          eval month${m}=1
      fi
      unset tmp

      # after about one year we keep only one per year
      if [ $age -gt 336 ]
      then
          y=$(expr $age / 336)
          eval tmp=\$year${y}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          fi
          eval year${y}=1
      fi
  done
  rm -f ${sto}
  ln -s ${sym} ${sto}
}

if [ -n "${MINGW_REPO_USER}" ] \
&& [ -n "${MINGW_REPO_HOST}" ] \
&& [ -n "${MINGW_REPO_DIR}" ]
then
    pushd ${HOME}/mingw-cross

    rsync -arpv --delete 32bit 64bit \
        ${MINGW_REPO_USER}@${MINGW_REPO_HOST}:${MINGW_REPO_DIR}/

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
