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


URLS=$(grep -Rh "_URL=" init_http_cache.sh | grep -v ^LAMMPS | grep -v SCAFACOS_FIX | cut -d= -f1)

# update URLs by grabbing the latest ones from cmake files
for URL in $URLS
do
  extract_setting "$URL"
  update_setting "$URL" ${!URL}
done
