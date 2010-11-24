=== Vim Syntax Highlighting ===
===============================

The files provided in this directory will enable syntax highlighting
for the lammps script syntax in vim. The simulation scripts have to
end on .lmp (see mysyntax.vim).  By far not all commands are included
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

Gerolf Ziegenhain <gerolf@ziegenhain.com> 2007

---------------

updated by Sam Bateman, 11/2010

Sam Bateman
Naval Research Laboratory
Code 7434
1005 Balch Blvd.
Stennis Space Center, MS 39529
Phone: (228) 688-4328
Email: sam.bateman@nrlssc.navy.mil
