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


error=0;


if [[ -d ${LAMMPSSRC} && -d ${LAMMPSDEVELOPERDOC} && -d ${LAMMPSDOC} && -f ${LAMMPSDOCMAKE} && -d ${LAMMPSDOXYGEN} && -f ${LAMMPSDOXYFILE} ]] ;
 then
   cd ${LAMMPSSRC} ;
   make no-all ;
   cd ${LAMMPSBASE} ;
 else
  error=1 ;
  echo "Cannot configure LAMMPS sources - Please run doxygen.sh from the LAMMPS root directory."
 fi ;


if [ ${error} == 0 ] ;
  then
    if [[ ! -x ${FIG2DEV} ]] ;
      then
        error=1 ;
        echo "${FIG2DEV} not found - Please install ${FIG2DEV} for Your operating system." ;
      fi ;
  fi ;


if [ ${error} == 0 ] ;
  then
    if [[ -d ${LAMMPSDEVELOPERDOC} ]] ;
      then
        ${FIG2DEV} -L png ${LAMMPSDEVELOPERDOC}/classes.fig > tools/doxygen/classes.png ;
        ${FIG2DEV} -L eps ${LAMMPSDEVELOPERDOC}/classes.fig > tools/doxygen/classes.eps ;
      else
        error=1 ;
        echo "LAMMPS developer documentation not found - Please control Your LAMMPS installation." ;
      fi;
  fi;


if [ ${error} == 0 ] ;
  then
    if [[ ! -x ${DOXYGEN} ]] ;
      then
        error=1;
        echo "doxygen not found - Please install doxygen for Your operating system." ;
      fi ;
  fi ;


if [ ${error} == 0 ] ;
  then
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
        echo "Doxyfile prototype not found - Please control Your LAMMPS installation." ;
      fi ;
  fi ;


if [ ${error} == 0 ] ;
  then
    if [[ -d ${LAMMPSDOC} && -f ${LAMMPSDOCMAKE} ]] ;
      then
        cd ${LAMMPSDOC} ;
        if [ -d ${LAMMPSDOCTMP} ] ;
         then
          rm -vrf ${LAMMPSDOCTMP} ;
         fi ;
        make clean ;
        make html ;
        make pdf ;
        cd ${LAMMPSBASE} ;
      echo "${LAMMPSDOC} exists." ;
    else
      error=1 ;
      echo "Cannot build LAMMPS native documentation - Please run doxygen.sh from the LAMMPS root directory." ;
    fi ;
  fi ;


if [ ${error} == 0 ] ;
  then
    if [[ ${error} == 0 && -d ${LAMMPSDOCHTML} ]] ;
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
          error=1 ;
          echo "Cannot include LAMMPS native documentation into the doxygen documentation - Please run doxygen.sh from the LAMMPS root directory." ;
    fi ;
  fi ;

