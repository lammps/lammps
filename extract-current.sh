#!/bin/sh

git archive -v --format=tar --prefix=lammps-current/  HEAD  README LICENSE doc/Manual.pdf src lib tools/*.cpp tools/*.f | gzip -9c - > ../lammps-current.tar.gz

