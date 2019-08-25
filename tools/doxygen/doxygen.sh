#!/bin/bash

LAMMPSBASE=${PWD} ;
LAMMPSSRC="${LAMMPSBASE}/src" ;
LAMMPSDOC="${LAMMPSBASE}/doc" ;
LAMMPSDOCHTML="${LAMMPSDOC}/html" ;
LAMMPSDEVELOPERDOC="${LAMMPSBASE}/doc/src/Developer" ;
LAMMPSDOCTMP="/tmp/lammps-docs-*" ;
LAMMPSDOCMAKE="${LAMMPSDOC}/Makefile" ;
LAMMPSDOXYGEN="${LAMMPSBASE}/tools/doxygen" ;
LAMMPSDOXYGENDOC="${LAMMPSDOXYGEN}/doc" ;
LAMMPSDOXYGENDOCMANUAL="${LAMMPSDOXYGEN}/doc/html/Manual" ;

LAMMPSDOXYFILE=${LAMMPSDOXYGEN}/Doxyfile.lammps ;
DOXYFILE=${LAMMPSDOXYGEN}/Doxyfile ;

LAMMPSDEVELOPERDOXFILE=${LAMMPSDOXYGEN}/Developer.dox.lammps ;
DEVELOPERDOXFILE=${LAMMPSDOXYGEN}/Developer.dox ;

FIG2DEV=`which fig2dev` ;
DOXYGEN=`which doxygen` ;


if [[ -d ${LAMMPSSRC} && -d ${LAMMPSDEVELOPERDOC} && -d ${LAMMPSDOC} && -f ${LAMMPSDOCMAKE} && -d ${LAMMPSDOXYGEN} && -f ${LAMMPSDOXYFILE} ]] ;
 then
   cd ${LAMMPSSRC} ;
   make no-all ;
   cd ${LAMMPSBASE} ;
 else
   echo "Cannot configure LAMMPS sources - Please run doxygen.sh from the LAMMPS root directory." >&2 ;
   exit 1 ;
 fi ;


if [ ! -x ${FIG2DEV} ] ;
  then
    echo "${FIG2DEV} not found - Please install ${FIG2DEV} for Your operating system." >&2 ;
    exit 1 ;
  fi ;


if [ -d ${LAMMPSDEVELOPERDOC} ] ;
  then
    ${FIG2DEV} -L png ${LAMMPSDEVELOPERDOC}/classes.fig > tools/doxygen/classes.png ;
    ${FIG2DEV} -L eps ${LAMMPSDEVELOPERDOC}/classes.fig > tools/doxygen/classes.eps ;
  else
    echo "LAMMPS developer documentation not found - Please control Your LAMMPS installation." >&2 ;
    exit 1 ;
  fi ;


if [ ! -x ${DOXYGEN} ] ;
  then
    echo "doxygen not found - Please install doxygen for Your operating system." >&2 ;
    exit 1 ;
  fi ;


if [[ -d ${LAMMPSSRC} && -f ${LAMMPSDOXYFILE} && -f ${LAMMPSDEVELOPERDOXFILE} ]] ;
  then
    ICS=' ';
    read -r -a array < ${LAMMPSSRC}/version.h ;
    version=`echo ${array[2]} ${array[3]} ${array[4]} | sed s/\"//g`;
    cp -av ${LAMMPSDOXYFILE} ${DOXYFILE} ;
    sed -i "s/LAMMPS_VERSION/${version}/g"  ${DOXYFILE} ;
    cp -av ${LAMMPSDEVELOPERDOXFILE} ${DEVELOPERDOXFILE} ;
    sed -i "s/LAMMPS_VERSION/${version}/g"  ${DEVELOPERDOXFILE} ;
    ${DOXYGEN} ${DOXYFILE} ;
  else
    echo "Doxyfile prototype not found - Please control Your LAMMPS installation." >&2 ;
    exit 1 ;
  fi ;


if [[ -d ${LAMMPSDOC} && -f ${LAMMPSDOCMAKE} ]] ;
  then
    cd ${LAMMPSDOC} ;
    for d in ${LAMMPSDOCTMP} ;
      do
        rm -vrf $d;
      done ;
    make clean ;
    make html ;
    make pdf ;
    cd ${LAMMPSBASE} ;
    echo "${LAMMPSDOC} exists." ;
  else
    echo "Cannot build LAMMPS native documentation - Please run doxygen.sh from the LAMMPS root directory."i >&2 ;
    exit 1 ;
  fi ;


if [ -d ${LAMMPSDOCHTML} ] ;
  then
    if [ -d ${LAMMPSDOXYGENDOCMANUAL} ] ;
      then
        rm -vrf ${LAMMPSDOXYGENDOCMANUAL} ;
        echo "${LAMMPSDOXYGENDOCMANUAL} removed." ;
      fi ;
    cp -av ${LAMMPSDOCHTML} ${LAMMPSDOXYGENDOCMANUAL} ;
    cp -av ${LAMMPSDOC}/*.pdf ${LAMMPSDOXYGENDOCMANUAL} ;
    echo "${LAMMPSDOXYGENDOCMANUAL} copied." ;
  else
    echo "Cannot include LAMMPS native documentation into the doxygen documentation - Please run doxygen.sh from the LAMMPS root directory." >&2 ;
    exit 1 ;
  fi ;

