.. index:: group2ndx

group2ndx command
=================

ndx2group command
=================

Syntax
""""""

.. parsed-literal::

   group2ndx file group-ID ...
   ndx2group file group-ID ...

* file = name of index file to write out or read in
* zero or more group IDs may be appended

Examples
""""""""

.. parsed-literal::

   group2ndx allindex.ndx
   group2ndx someindex.ndx upper lower mobile
   ndx2group someindex.ndx
   ndx2group someindex.ndx mobile

Description
"""""""""""

Write or read a Gromacs style index file in text format that associates
atom IDs with the corresponding group definitions. This index file can be
used with in combination with Gromacs analysis tools or to import group
definitions into the :doc:`fix colvars <fix_colvars>` input file. It can
also be used to save and restore group definitions for static groups.

The *group2ndx* command will write group definitions to an index file.
Without specifying any group IDs, all groups will be written to the index
file. When specifying group IDs, only those groups will be written to the
index file. In order to follow the Gromacs conventions, the group *all*
will be renamed to *System* in the index file.

The *ndx2group* command will create of update group definitions from those
stored in an index file. Without specifying any group IDs, all groups except
*System* will be read from the index file and the corresponding groups
recreated. If a group of the same name already exists, it will be completely
reset. When specifying group IDs, those groups, if present, will be read
from the index file and restored.

----------

Restrictions
""""""""""""

This command requires that atoms have atom IDs, since this is the
information that is written to the index file.

These commands are part of the USER-COLVARS package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`group <group>`, :doc:`dump <dump>`, :doc:`fix colvars <fix_colvars>`

**Default:** none
