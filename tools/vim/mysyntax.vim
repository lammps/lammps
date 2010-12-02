augroup syntax
au  BufNewFile,BufReadPost *.lmp so ~/.vim/lammps.vim
au  BufNewFile,BufReadPost in.* so ~/.vim/lammps.vim
augroup END
