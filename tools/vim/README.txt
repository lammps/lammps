=== Vim Syntax Highlighting ===
===============================

The files provided in this directory will enable syntax highlighting
for the lammps script syntax in vim. The simulation scripts have to
end on *.lmp or start with in.* (see mysyntax.vim).  
By far not all commands are included
in the syntax file (lammps.vim). You can easily add new ones.

=To enable the highlighting (compatible with old versions of VIM):
============================
(0)   Create a ~/.vimrc
      You can have a look in /usr/share/vim/vim*/vimrc_example.vim
(1)   Insert in ~/.vimrc
         let mysyntaxfile = "~/.vim/mysyntax.vim"
      just before
         syntax on
(2)   Create directory ~/.vim and place mysyntax.vim and lammps.vim there

=Here is an alternate method for VIM Version 7.2 and later:
===========================================================

(0) Create/edit ~/.vimrc to contain:
         syntax on
(1) Create directories ~/.vim/syntax and ~/.vim/ftdetect
(2) Copy lammps.vim to ~/.vim/syntax/lammps.vim
(3) Copy filetype.vim to ~/.vim/ftdetect/lammps.vim


Gerolf Ziegenhain <gerolf@ziegenhain.com> 2007

Distribution Packaging guidelines:
==================================

(0) Copy lammps.vim to ${VIMFILES}/syntax/lammps.vim
(1) Copy filetype.vim as ${VIMFILES}/ftdetect/lammps.vim

${VIMFILES} is typically /usr/share/vim/vimfiles
Consult your packaging guidlines for exact location.

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
