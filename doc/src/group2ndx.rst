.. index:: group2ndx
.. index:: ndx2group

group2ndx command
=================

ndx2group command
=================

Syntax
""""""

.. code-block:: LAMMPS

   group2ndx file args
   ndx2group file args

* file = name of index file to write out or read in
* args = zero or more group IDs may be appended

Examples
""""""""

.. code-block:: LAMMPS

   group2ndx allindex.ndx
   group2ndx someindex.ndx upper lower mobile
   ndx2group someindex.ndx
   ndx2group someindex.ndx mobile

Description
"""""""""""

Write or read a Gromacs style index file in text format that associates
atom IDs with the corresponding group definitions. This index file can be
used with in combination with Gromacs analysis tools or to import group
definitions into the :doc:`fix colvars <fix_colvars>` input file.

It can also be used to save and restore group definitions for static groups
using the individual atom IDs. This may be important if the original
group definition depends on a region or otherwise on the geometry and thus
cannot be easily recreated.

Another application would be to import atom groups defined for Gromacs
simulation into LAMMPS.  When translating Gromacs topology and geometry
data to LAMMPS.

The *group2ndx* command will write group definitions to an index file.
Without specifying any group IDs, all groups will be written to the
index file.  When specifying group IDs, only those groups will be
written to the index file.  In order to follow the Gromacs conventions,
the group *all* will be renamed to *System* in the index file.

The *ndx2group* command will create of update group definitions from
those stored in an index file.  Without specifying any group IDs, all
groups except *System* will be read from the index file and the
corresponding groups recreated.  If a group of the same name already
exists, it will be completely reset.  When specifying group IDs, those
groups, if present, will be read from the index file and restored.

File Format
"""""""""""

The file format is equivalent and compatible with what is produced by
the `Gromacs make_ndx command <https://manual.gromacs.org/current/onlinehelp/gmx-make_ndx.html>`_.
and follows the `Gromacs definition of an ndx file <https://manual.gromacs.org/current/reference-manual/file-formats.html#ndx>`_

Each group definition begins with the group name in square brackets with
blanks, e.g. \[ water \] and is then followed by the list of atom
indices, which may be spread over multiple lines.  Here is a small
example file:

.. code-block:: ini

   [ Oxygen ]
   1  4  7
   [ Hydrogen ]
   2  3  5  6
   8  9
   [ Water ]
   1 2 3 4 5 6 7 8 9

The index file defines 3 groups: Oxygen, Hydrogen, and Water and the
latter happens to be the union of the first two.

----------

Restrictions
""""""""""""

These commands require that atoms have atom IDs, since this is the
information that is written to the index file.

These commands are part of the EXTRA-COMMAND package.  They are only
enabled if LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`group <group>`, :doc:`dump <dump>`, :doc:`fix colvars <fix_colvars>`

Default
"""""""

none
