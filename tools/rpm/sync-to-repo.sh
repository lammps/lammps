#!/bin/sh
# sync data from build directory to repository

# function to selectively remove older rpms and make a symlink to the latest version
prune_rpm () {
  dir=$1
  sub=$2
  ref=$(date +%s)
  old=999999999

  for rpm in ${dir}/lammps${sub}-20[0-9][0-9]*.rpm
  do \
      [ -f ${rpm} ] || continue
      # re-set symbolic link to latest entry
      p=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[0-9].*\)\(\..*\.rpm\)$@\1@')
      r=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[0-9].*\)\(\..*\.rpm\)$@\2@')
      t=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[0-9].*\)\(\..*\.rpm\)$@\3@')
      v=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[0-9].*\)\(\..*\.rpm\)$@\4@')
      e=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[0-9].*\)\(\..*\.rpm\)$@\5@')

      # compute age difference in days
      y=$(echo ${t} | cut -c 1-4)
      m=$(echo ${t} | cut -c 5-6)
      d=$(echo ${t} | cut -c 7-8)
      s=$(date +%s -d "${m}/${d}/${y}")
      age=$(expr \( $ref - $s \) / 86400)

      if [ $age -lt $old ] 
      then
          old=$age
          sym="${r}${t}${v}${e}"
          sto="${p}${r}latest${e}"
      fi

      # NOTE: to simplify the math, for the following we
      # define one month to have 28 days and a year of 
      # 12 months to have correspondingly only 336 days.

      # after about one year we keep only one per year
      unset tmp
      if [ $age -gt 336 ]
      then
          y=$(expr $age / 336)
          eval tmp=\$year${sub#-}${y}
          if [ -n "$tmp" ]
          then
              rm -vf ${rpm}
          else
              echo "first in year $y $rpm"
          fi
          export year${sub#-}${y}=1
      fi

      # after about three months we keep only one per month
      unset tmp
      if [ $age -gt 84 ] && [ $age -lt 336 ]
      then
          m=$(expr $age / 28)
          eval tmp=\$month${sub#-}${m}
          if [ -n "$tmp" ]
          then
              rm -vf ${rpm}
          else
              echo "first in month $m $rpm"
          fi
          export month${sub#-}${m}=1
      fi

      # after one week we keep only one per week.
      unset tmp
      if [ $age -gt 7 ] && [ $age -lt 84 ]
      then
          w=$(expr $age / 7)
          eval tmp=\$week${sub#-}${w}
          if [ -n "$tmp" ]
          then
              rm -vf ${rpm}
          else
              echo "first in week $w $rpm"
          fi
          export week${sub#-}${w}=1
      fi

  done

  if [ -n "${sto}" ]
  then
      rm -f ${sto}
      ln -s ${sym} ${sto}
  fi
}


# function to build delta rpms
make_drpms () {
  dir=$1
  sub=$2
  prv=

  test -d "${p}drpms" || mkdir -p "${p}drpms"

  for rpm in ${dir}/lammps${sub}-20[0-9][0-9]*.rpm
  do \
      [ -f ${rpm} ] || continue
      # re-set symbolic link to latest entry
      p=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[-a-z0-9]*\)\(\..*\.rpm\)$@\1@')
      r=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[-a-z0-9]*\)\(\..*\.rpm\)$@\2@')
      t=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[-a-z0-9]*\)\(\..*\.rpm\)$@\3@')
      v=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[-a-z0-9]*\)\(\..*\.rpm\)$@\4@')
      e=$(echo ${rpm} | sed -e 's@^\(.*/\)\(lammps-.*\)\(20[0-9][0-9]\+\)\(-[-a-z0-9]*\)\(\..*\)\.rpm$@\5.drpm@')

      [ -n "$prv" ] && makedeltarpm -v ${prv} ${rpm} "${p}drpms/${r}${tp}${vp}_${t}${v}${e}"

      prv=${rpm}
      tp=${t}
      vp=${v}
  done
}


if ! which rpmbuild > /dev/null 2>&1
then
    echo "Cannot build rpms without 'rpmbuild' installed"
    exit 1
fi

MYRPM_BUILD_DIR=$(rpmbuild --showrc 2>&1 | sed -n -e 's/^.*: _topdir[	 ]\+\(.*\)$/\1/p')
datestr=$(date +%Y%m%d)

if [ -n "${MYRPM_REPO_USER}" ] \
&& [ -n "${MYRPM_REPO_HOST}" ] \
&& [ -n "${MYRPM_REPO_DIR}" ]
then
    pushd ${MYRPM_BUILD_DIR}
    # move debuginfo rpms, if present
    for f in RPMS/*-debuginfo-*.rpm
    do \
        test -f $f && mv -v $f RPMS/debug/
    done

    prune_rpm RPMS/debug -debuginfo
    prune_rpm RPMS
    prune_rpm RPMS -common
    prune_rpm RPMS -doc
    prune_rpm RPMS -openmpi
    prune_rpm RPMS -mpich2
    prune_rpm RPMS -mpich
    prune_rpm RPMS -python

    # remove all previous delta rpms of this group. we will recreate all
    # of them and do not want to have to prune here, too.
    rm -vf RPMS/drpms/*.drpm

    #make_drpms RPMS/debug -debuginfo
    make_drpms RPMS -common
    make_drpms RPMS -doc
    make_drpms RPMS -openmpi
    make_drpms RPMS -mpich2
    make_drpms RPMS -mpich
    make_drpms RPMS -python

    createrepo --deltas -x \*latest\* -v RPMS
    rsync -arpv --exclude cache --delete RPMS/ \
        ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST}:${MYRPM_REPO_DIR}
    # we use the source rpm on suse, since they do not have a "dist" tag.
    if grep -q -i suse\ 12 /etc/issue && test $(uname -i) = i386
    then
        prune_rpm SRPMS
        rsync -arpv SRPMS/ \
            ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST}:${MYRPM_REPO_DIR}/../../../source
    fi
    ssh ${MYRPM_REPO_USER}@${MYRPM_REPO_HOST} "cd ${MYRPM_REPO_DIR}/../../../; ./mkhtmlindex.sh *"
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
