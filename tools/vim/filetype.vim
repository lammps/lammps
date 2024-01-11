augroup filetypedetect
 au! BufRead,BufNewFile in.*           setfiletype lammps
 au! BufRead,BufNewFile *.lmp          setfiletype lammps
 au! BufRead,BufNewFile *.in           setfiletype lammps
augroup END
