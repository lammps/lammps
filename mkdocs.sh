#!/bin/sh

# update equations

for s in doc/Eqs/*.tex 
do \
	n=`basename $s .tex`
	if [ ! -e doc/Eqs/${n}.jpg ] || [ $s -nt doc/Eqs/${n}.jpg ]
	then
		pdflatex $s
		pdftoppm ${n}.pdf | pnmcrop | pnmtojpeg > doc/Eqs/${n}.jpg
        rm -f ${n}.aux ${n}.log ${n}.pdf
	fi
done

make -C txt2html

./txt2html/txt2html doc/*.txt
htmldoc --batch lammps.book

