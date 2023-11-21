.. index:: fix reaxff/atom
.. index:: fix reaxff/atom/kk

fix reaxff/atom command
=======================

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

TODO

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
