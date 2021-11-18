#!/bin/bash
#helper script to update URLs from CMake folder

CMAKE_FOLDER=../../../cmake

function extract_setting()
{
  export $1=$(grep -Rh "set($1" ../../../cmake/ | cut -d\" -f2)
}

function update_setting()
{
  echo Setting $1 to $2
  sed -i "/^$1=/c$1=\"$2\"" init_http_cache.sh
}

DETECTED_URLS=$(grep -PoRh "\w+_URL" ../../../cmake/ | sort | uniq | grep -v ^LAMMPS)
KNOWN_URLS=$(grep -Rh "_URL=" init_http_cache.sh | grep -v ^LAMMPS | grep -v SCAFACOS_FIX | cut -d= -f1)

# check if init_http_cache.sh contains all URLs
for URL in $DETECTED_URLS
do
  grep -q ^$URL= init_http_cache.sh
  if [ $? -ne 0 ]
  then
    FILENAME_VAR="${URL/_URL/_FILENAME}"
    echo $URL is not known. Please update 'init_http_cache.sh' as follows:
    echo
    echo 1. add the following line:
    echo
    echo $URL=""
    echo
    echo 2. Define a new $FILENAME_VAR if necessary
    echo
    echo $FILENAME_VAR="pkgname-0.0.0.tar.gz"
    echo
    echo 3. extend TARBALLS with $URL
    echo
    echo TARBALLS=\(
    echo "    ..."
    echo "    $URL"
    echo \)
    echo
    echo 4. Rerun this script
    echo
    exit 1
  fi
done

# update URLs by grabbing the latest ones from cmake files
for URL in $KNOWN_URLS
do
  extract_setting "$URL"
  update_setting "$URL" ${!URL}
done
