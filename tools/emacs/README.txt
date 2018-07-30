=== Emacs Syntax Highlighting ===
Created by Aidan Thompson 12/2010
===============================
Updated by Roit Goswami Mon Jul 30 2018 

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
   
(load "~/.emacs.d/lammps-mode.el")

This file may also be called ~/.emacs.el, or ~/.emacs.d/init.el

(1) Copy lammps-mode.el to the directory ~/.emacs.d

=Update:
========

The package may now also be installed by a MELPA style recipe, namely:

```lisp
(lammps-mode :fetcher github :repo "HaoZeke/lammps-mode")
```

For a simpler installation with `use-package` simply add:

```lisp
(use-package lammps-mode)
```

The latest version of the package will be kept in sync as a squashed update on
the lammps repository as well.

It is advisable to use the MELPA installation methods listed here:
https://melpa.org/#/getting-started

For autoloading and auto-recognizing "in.*" and "*.lmp" files add the following
to `.emacs`:

```lisp
(autoload 'lammps-mode "lammps-mode.el" "LAMMPS mode." t)
(setq auto-mode-alist (append auto-mode-alist
                              '(("in\\." . lammps-mode))
                              '(("\\.lmp\\'" . lammps-mode))
                              ))
```
