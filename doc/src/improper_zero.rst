.. index:: improper_style zero

improper_style zero command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style zero [nocoeff]

Examples
""""""""

.. code-block:: LAMMPS

   improper_style zero
   improper_style zero nocoeff
   improper_coeff *

Description
"""""""""""

Using an improper style of zero means improper forces and energies are
not computed, but the geometry of improper quadruplets is still
accessible to other commands.

As an example, the :doc:`compute improper/local <compute_improper_local>` command can be used to
compute the chi values for the list of quadruplets of improper atoms
listed in the data file read by the :doc:`read_data <read_data>`
command.  If no improper style is defined, this command cannot be
used.

The optional *nocoeff* flag allows to read data files with a ImproperCoeff
section for any improper style. Similarly, any improper_coeff commands
will only be checked for the improper type number and the rest ignored.

Note that the :doc:`improper_coeff <improper_coeff>` command must be
used for all improper types, though no additional values are
specified.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

none


:doc:`improper_style none <improper_none>`

Default
"""""""

none
