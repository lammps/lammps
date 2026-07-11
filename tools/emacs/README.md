# GNU Emacs Syntax Highlighting

> Copyright (C) 2010-2018  Aidan Thompson <athomps at sandia.gov>
> Copyright (C) 2018  Rohit Goswami <r95g10 at gmail.com>

The `lammps-mode.el` file provided in this directory will enable syntax
highlighting for the lammps script syntax in GNU Emacs. The groupings of
commands were originally copied from `tools/vim`.

## Installation
**Requirements: GNU Emacs 24.\***

### Obtaining the Package

#### MELPA

The easiest installation method is via MELPA and it is advisable to use one of
the many [MELPA installation methods](https://melpa.org/#/getting-started).

For example, with [use-package](https://github.com/jwiegley/use-package) one can
simply use the following:

``` emacs-lisp
(use-package lammps-mode)
```

#### Manually

Assuming for some reason you have downloaded the file to `~/.emacs.d/lisp` you
would do the following (kanged [from here](http://ergoemacs.org/emacs/emacs_installing_packages.html)):

``` emacs-lisp
;; Tell emacs where is your personal elisp lib dir
(add-to-list 'load-path "~/.emacs.d/lisp/")

;; load the package.
(load "lammps-mode")
```

### Autoloading \& Auto-recognition

To automatically turn on the LAMMPS mode for editing your input scripts,
use the following line as the **first** line of your script:
```
# -*- lammps -*-
```

For automatically switching on the LAMMPS mode based on filename patterns,
e.g. for `in.*` and `*.lmp` files, add the following code to your `.emacs`:

``` emacs-lisp
(autoload 'lammps-mode "lammps-mode.el" "LAMMPS mode." t)
(setq auto-mode-alist (append auto-mode-alist
                              '(("in\\." . lammps-mode))
                              '(("\\.lmp\\'" . lammps-mode))
                              ))
```

## Status

By far not all commands are included in the syntax file (lammps-mode.el). You
can easily add new ones to the existing classes.

## Implementation Details

`lammps-mode` is derived from `shell-script-mode` which provides some basic
syntax highlighting of strings, comments, etc.

The MELPA recipe used for this package is simply:

``` emacs-lisp
(lammps-mode :fetcher github :repo "HaoZeke/lammps-mode")
```

## Caveats

* Does not work with Xemacs [See [this comment](https://github.com/lammps/lammps/pull/1022#issuecomment-408871233)]

## License

[GNU GPL v2](https://github.com/HaoZeke/lammps-mode/blob/master/LICENSE).
Check the file for more details.
