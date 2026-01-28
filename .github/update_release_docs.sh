#!/bin/bash

# check out LAMMPS release version
git checkout release || exit 1

# make sure we are in the LAMMPS root directory
if [ ! -e README ] || [ ! -e LICENSE ] || [ ! -e SECURITY.md ] || [ ! -e CITATION.cff ]
then
    echo "Must be in the LAMMPS root folder to run this script"
    exit 2
fi

export LC_ALL=C
export webserver=publish@www.lammps.org
export webroot=/var/www/lammps/download
export docroot=/var/www/lammps/docs

# build documentation using its own virtual environment
export LAMMPS_WEBSITE_BUILD=1
export LAMMPS_WEBSITE_BUILD_VERSION=release
export LAMMPS_WEBSITE_BASEURL="https://docs.lammps.org/release/"
make -C doc clean
make -C doc upgrade
make -C doc html WEB_SEARCH=YES
make -C doc pdf
rsync -arp doc/html/* ${webserver}:${docroot}/release/
rsync -p doc/Manual.pdf  ${webserver}:${docroot}/release/
