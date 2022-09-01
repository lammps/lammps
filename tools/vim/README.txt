=== Vim Syntax Highlighting ===
===============================

The files provided in this directory will enable syntax highlighting
for the lammps script syntax in vim. The simulation scripts have to
end on .lmp or .in or start with in. (see mysyntax.vim or filetype.vim).
By far not all commands are included in the syntax file (lammps.vim).
You can easily add new ones.

=To enable the highlighting (compatible with old versions of VIM):
============================
(1)   Create a ~/.vimrc
      You can have a look in /usr/share/vim/vim*/vimrc_example.vim
(2)   Insert in ~/.vimrc
         let mysyntaxfile = "~/.vim/mysyntax.vim"
      just before
         syntax on
(3)   Create directory ~/.vim and place mysyntax.vim and lammps.vim there

=Here is an alternate method for VIM Version 7.2 and later:
===========================================================

(1) Create/edit ~/.vimrc to contain:
         syntax on
(2) Create directories ~/.vim/syntax and ~/.vim/ftdetect
(3) Copy lammps.vim to ~/.vim/syntax/lammps.vim
(4) Copy filetype.vim to ~/.vim/ftdetect/lammps.vim

Distribution Packaging guidelines:
==================================

(1) Copy lammps.vim to ${VIMFILES}/syntax/lammps.vim
(2) Copy filetype.vim as ${VIMFILES}/ftdetect/lammps.vim

${VIMFILES} is typically /usr/share/vim/vimfiles
Consult your packaging guidlines for exact location.

Gerolf Ziegenhain <gerolf@ziegenhain.com> 2007

=Disable insertion of <TAB> characters
========================================

(1) The LAMMPS developers strongly discourage writing
   files (input and C++ source code) with <TAB> characters.
   To change the default behavior of vim use in your vimrc:
       set expandtab

---------------

updated by Sam Bateman, 11/2010
updated by Aidan Thompson, 12/2010
updated by Axel Kohlmeyer, 08/2022
