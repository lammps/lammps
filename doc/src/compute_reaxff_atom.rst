.. index:: compute reaxff/atom
.. index:: compute reaxff/atom/kk

compute reaxff/atom command
===========================

Accelerator Variants: *reaxff/atom/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID reaxff/atom attribute args ... keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* reaxff/atom = name of this compute command
* attribute = *pair*

  .. parsed-literal::

       *pair* args = nsub
         nsub = *n*-instance of a sub-style, if a pair style is used multiple times in a hybrid style

* keyword = *bonds*

  .. parsed-literal::

       *bonds* value = *no* or *yes*
         *no* = ignore list of local bonds
         *yes* = include list of local bonds

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all reaxff/atom bonds yes

Description
"""""""""""

.. versionadded:: 7Feb2024

Define a computation that extracts bond information computed by the ReaxFF
potential specified by :doc:`pair_style reaxff <pair_reaxff>`.

By default, it produces per-atom data that includes the following columns:

* abo = atom bond order (sum of all bonds)
* nlp = number of lone pairs
* nb = number of bonds

Bonds will only be included if its atoms are in the group.

In addition, if ``bonds`` is set to ``yes``, the compute will also produce a
local array of all bonds on the current processor whose atoms are in the group.
The columns of each entry of this local array are:

* id_i = atom i id of bond
* id_j = atom j id of bond
* bo = bond order of bond

Output info
"""""""""""

This compute calculates a per-atom array and local array depending on the
number of keywords. The number of rows in the local array is the number of
bonds as described above. Both per-atom and local array have 3 columns.

The arrays can be accessed by any command that uses local and per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` page for
an overview of LAMMPS output options.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The compute reaxff/atom command requires that the :doc:`pair_style reaxff
<pair_reaxff>` is invoked.  This fix is part of the REAXFF package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`

Default
"""""""

The option defaults are *bonds* = *no*.
