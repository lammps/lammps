# Outline of LAMMPS documentation file conventions

This purpose of this document is to provide a point of reference
for LAMMPS developers and contributors as to what conventions
should be used to structure and format files in the LAMMPS manual.

Last change: 2020-04-23

## File format and tools

In fall 2019, the LAMMPS documentation file format has changed from
a home grown minimal markup designed to generate HTML format files
from a mostly plain text format to using the reStructuredText file
format.  For a transition period all files in the old .txt format
were transparently converted to .rst and then processed.  The txt2rst
tool is still included in the distribution to obtain an initial .rst
file for integration into the manual.  Since the transition to
reStructured text as source format, many of the artifacts or the
translation have been removed though and parts of the documentation
refactored and expanded to take advantage of the capabilities
reStructuredText and associated tools.  The conversion from the
source to the final formats (HTML, PDF, and optionally e-book
reader formats ePUB and MOBI) is mostly automated and controlled
by a Makefile in the `doc` folder. This makefile assumes that the
processing is done on a Unix-like machine and Python 3.5 or later
and a matching virtualenv module are available.  Additional Python
packages (like the Sphinx tool and several extensions) are
transparently installed into a virtual environment over the
internet using the `pip` package manager.  Further requirements
and details are discussed in the manual.

## Work in progress

The refactoring and improving of the documentation is an ongoing
process, so statements in this document may not always be fully
up-to-date.  If in doubt, contact the LAMMPS developers.

## General structure

The layout and formatting of added files should follow the example
of the existing files.  Since those are directly derived from their
former .txt format versions and the manual has been maintained in
that format for many years, there is a large degree of consistency
already, so comparision with similar files should give you a good
idea what kind of information and sections are needed.

## Formatting conventions

Filenames, folders, paths, (shell) commands, definitions, makefile
settings and similar should be formatted as "literals" with
double backward quotes bracketing the item: \`\`path/to/some/file\`\`

Keywords and options are formatted in italics:  \*option\*

Mathematical expressions, equations, symbols are typeset using
either a `.. math:`` block or the `:math:` role.

Groups of shell commands or LAMMPS input script or C/C++ source
code should be typeset into a `.. code-block::` section. A syntax
highlighting extension for LAMMPS input scripts is provided, so
`LAMMPS` can be used to indicate the language in the code block
in addition to `bash`, `c`, or `python`.  When no syntax style
is indicated, no syntax highlighting is performed.

As an alternative, e.g. to typeset the syntax of file formats
a `.. parsed-literal::` block can be used, which allows some
formatting directives, which means that related characters need
to be escaped with a preceding backslash: `\*`.

Special remarks can be highlighted with a `.. note::` block and
strong warnings can be put into a `.. warning::` block.

## Required steps when adding a custom style to LAMMPS

When adding a new style (e.g. pair style or a compute or a fix)
or a new command, it is **required** to include the corresponding
documentation.  Those are often new files that need to be added.
In order to be included in the documentation, those new files
need to be reference in a `.. toctree::` block.  Most of those
use patterns with wildcards, so the addition will be automatic.
However, those additions also need to be added to some lists of
styles or commands.  The `make style\_check` command will perform
a test and report any missing entries and list the affected files.
Any references defined with `.. \_refname:` have to be unique
across all documentation files and this can be checked for with
`make anchor\_check`.  Finally, a spell-check should be done,
which is triggered via `make spelling`.  Any offenses need to
be corrected and false positives should be added to the file
`utils/sphinx-config/false\_positives.txt`.

## Required additional steps when adding a new package to LAMMPS

TODO
