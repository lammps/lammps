#!/bin/sh

# update equations

for s in src/Eqs/*.tex 
do \
	n=`basename $s .tex`
	if [ ! -e src/Eqs/${n}.jpg ] || [ $s -nt src/Eqs/${n}.jpg ]
	then
		pdflatex $s
		pdftoppm ${n}.pdf | pnmcrop | pnmtojpeg > src/Eqs/${n}.jpg
        rm -f ${n}.aux ${n}.log ${n}.pdf
	fi
done

# convert all .txt files to html
make html

# check if we have any new html files,
# that are not yet in the book file.
for s in `echo src/*.txt | sed -e 's,src/,,g' -e 's,\.txt,\.html,g'`
do \
  grep -q $s lammps.book || echo doc file $s missing in lammps.book
done
( cd html; htmldoc --batch ../lammps.book )

git co -- html src

