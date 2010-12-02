=== Vim Syntax Highlighting ===
===============================

The files provided in this directory will enable syntax highlighting
for the lammps script syntax in vim. The simulation scripts have to
end on *.lmp or start with in.* (see mysyntax.vim).  
By far not all commands are included
in the syntax file (lammps.vim). You can easily add new ones.

=To enable the highlighting:
============================
(0)   Create a ~/.vimrc
      You can have a look in /usr/local/share/vim/current/vimrc.example
(1)   Insert in ~/.vimrc
         let mysyntaxfile = "~/.vim/mysyntax.vim"
      just before
         syntax on
(2)   Create directory ~/.vim and place mysyntax.vim and lammps.vim there

=Here is an alternate method for VIM Version 7.2 and later:
===========================================================

(0) Create/edit ~/.vimrc to contain:
         syntax on
(1) Create directories ~/.vim and ~/.vim/syntax
(2) Copy lammps.vim to ~/.vim/syntax/lammps.vim
(3) Create/edit ~/.vim/filetype.vim to contain

" vim syntax highlight customizations
if exists("did_load_filetypes")
 finish
endif

augroup filetypedetect
 au! BufRead,BufNewFile in.*           setfiletype lammps
 au! BufRead,BufNewFile *.lmp          setfiletype lammps
augroup END
(4) the end


Gerolf Ziegenhain <gerolf@ziegenhain.com> 2007

---------------

updated by Sam Bateman, 11/2010
updated by Aidan Thompson, 12/2010

Sam Bateman
Naval Research Laboratory
Code 7434
1005 Balch Blvd.
Stennis Space Center, MS 39529
Phone: (228) 688-4328
Email: sam.bateman@nrlssc.navy.mil
