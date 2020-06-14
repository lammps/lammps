.. index:: pair_coeff

pair_coeff command
==================

Syntax
""""""

.. code-block:: LAMMPS

   pair_coeff I J args

* I,J = atom types (see asterisk form below)
* args = coefficients for one or more pairs of atom types

Examples
""""""""

.. code-block:: LAMMPS

   pair_coeff 1 2 1.0 1.0 2.5
   pair_coeff 2 * 1.0 1.0
   pair_coeff 3* 1*2 1.0 1.0 2.5
   pair_coeff * * 1.0 1.0
   pair_coeff * * nialhjea 1 1 2
   pair_coeff * 3 morse.table ENTRY1
   pair_coeff 1 2 lj/cut 1.0 1.0 2.5 (for pair_style hybrid)

Description
"""""""""""

Specify the pairwise force field coefficients for one or more pairs of
atom types.  The number and meaning of the coefficients depends on the
pair style.  Pair coefficients can also be set in the data file read
by the :doc:`read_data <read_data>` command or in a restart file.

I and J can be specified in one of two ways.  Explicit numeric values
can be used for each, as in the first example above.  I <= J is
required.  LAMMPS sets the coefficients for the symmetric J,I
interaction to the same values.

A wildcard asterisk can be used in place of or in conjunction with the
I,J arguments to set the coefficients for multiple pairs of atom
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of atom types, then an asterisk with no numeric values means all
types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).  Note that only type pairs with I <= J are considered; if
asterisks imply type pairs where J < I, they are ignored.

Note that a pair_coeff command can override a previous setting for the
same I,J pair.  For example, these commands set the coeffs for all I,J
pairs, then overwrite the coeffs for just the I,J = 2,3 pair:

.. code-block:: LAMMPS

   pair_coeff * * 1.0 1.0 2.5
   pair_coeff 2 3 2.0 1.0 1.12

A line in a data file that specifies pair coefficients uses the exact
same format as the arguments of the pair_coeff command in an input
script, with the exception of the I,J type arguments.  In each line of
the "Pair Coeffs" section of a data file, only a single type I is
specified, which sets the coefficients for type I interacting with
type I.  This is because the section has exactly N lines, where N =
the number of atom types.  For this reason, the wild-card asterisk
should also not be used as part of the I argument.  Thus in a data
file, the line corresponding to the first example above would be listed
as

.. parsed-literal::

   2 1.0 1.0 2.5

For many potentials, if coefficients for type pairs with I != J are
not set explicitly by a pair_coeff command, the values are inferred
from the I,I and J,J settings by mixing rules; see the
:doc:`pair_modify <pair_modify>` command for a discussion.  Details on
this option as it pertains to individual potentials are described on
the doc page for the potential.

Many pair styles, typically for many-body potentials, use tabulated
potential files as input, when specifying the pair_coeff command.
Potential files provided with LAMMPS are in the potentials directory
of the distribution.  For some potentials, such as EAM, other archives
of suitable files can be found on the Web.  They can be used with
LAMMPS so long as they are in the format LAMMPS expects, as discussed
on the individual doc pages.

When a pair_coeff command using a potential file is specified, LAMMPS
looks for the potential file in 2 places.  First it looks in the
location specified.  E.g. if the file is specified as "niu3.eam", it
is looked for in the current working directory.  If it is specified as
"../potentials/niu3.eam", then it is looked for in the potentials
directory, assuming it is a sister directory of the current working
directory.  If the file is not found, it is then looked for in the
directory specified by the LAMMPS_POTENTIALS environment variable.
Thus if this is set to the potentials directory in the LAMMPS distribution,
then you can use those files from anywhere on your system, without
copying them into your working directory.  Environment variables are
set in different ways for different shells.  Here are example settings
for

csh, tcsh:

.. parsed-literal::

   % setenv LAMMPS_POTENTIALS /path/to/lammps/potentials

bash:

.. parsed-literal::

   % export LAMMPS_POTENTIALS=/path/to/lammps/potentials

Windows:

.. parsed-literal::

   % set LAMMPS_POTENTIALS="C:\\Path to LAMMPS\\Potentials"

----------

The alphabetic list of pair styles defined in LAMMPS is given on the
:doc:`pair_style <pair_style>` doc page.  They are also listed in more
compact form on the :doc:`Commands pair <Commands_pair>` doc page.

Click on the style to display the formula it computes and its
coefficients as specified by the associated pair_coeff command.

----------

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

Related commands
""""""""""""""""

:doc:`pair_style <pair_style>`, :doc:`pair_modify <pair_modify>`,
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`,
:doc:`pair_write <pair_write>`

**Default:** none
