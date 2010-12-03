=== Emacs Syntax Highlighting ===
Created by Aidan Thompson 12/2010
===============================

The lammps.el file provided in this directory will enable syntax 
highlighting for the lammps script syntax in emacs. The groupings
of commands were copied from tools/vim. The simulation scripts have to
end on *.lmp or start with in.* (see lammps.el). By far not all 
commands are included in the syntax file (lammps.el). 
You can easily add new ones to the existing classes.
'lammps-mode' is derived from 'shell-script-mode' which provides
some basic syntax highlighting of strings, comments, etc.

=To enable the highlighting:
============================
(0) Create/edit the emacs init file ~/.emacs to contain:
   
(load "~/.emacs.d/lammps")

This file may also be called ~/.emacs.el, or ~/.emacs.d/init.el

(1) Copy lammps.el to the directory ~/.emacs.d

