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

# build txt2html tool
make -C txt2html

# convert all .txt files to html
./txt2html/txt2html -b doc/*.txt

# check if we have any new html files,
# that are not yet in the book file.
for s in `echo doc/*.txt | sed -e 's,\.txt,\.html,g'`
do \
  grep -q $s lammps.book || echo doc file $s missing in lammps.book
done
htmldoc --batch lammps.book

