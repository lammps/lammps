
Internal Styles
---------------

LAMMPS has a number of styles that are not meant to be used in an input
file and thus are not documented in the :doc:`LAMMPS command
documentation <Commands_all>`.  The differentiation between user
commands and internal commands is through the case of the command name:
user commands and styles are all lower case, internal styles are all
upper case.  Internal styles are not called from the input file, but
their classes are instantiated by other styles.  Often they are
created by other styles to store internal data or to perform actions
regularly at specific steps of the simulation.

The paragraphs below document some of those styles that have general
utility and may be used to avoid redundant implementation.

DEPRECATED Styles
^^^^^^^^^^^^^^^^^

The styles called DEPRECATED (e.g. pair, bond, fix, compute, region, etc.)
have the purpose to inform users that a specific style has been removed
or renamed.  This is achieved by creating an alias for the deprecated
style to the corresponding class.  For example, the fix style DEPRECATED
is aliased to fix style ave/spatial and fix style ave/spatial/sphere with
the following code:

.. code-block:: c++

   FixStyle(DEPRECATED,FixDeprecated);
   FixStyle(ave/spatial,FixDeprecated);
   FixStyle(ave/spatial/sphere,FixDeprecated);

The individual class will then determine based on the style name
what action to perform:

- inform that the style has been removed and what style replaces it, if any, and then error out
- inform that the style has been renamed and then either execute the replacement or error out
- inform that the style is no longer required, and it is thus ignored and continue

There is also a section in the user's guide for :doc:`removed commands
and packages <Commands_removed>` with additional explanations.

Internal fix styles
^^^^^^^^^^^^^^^^^^^

fix DUMMY
"""""""""

Most fix classes cannot be instantiated before the simulation box has
been created since they access data that is only available then.
However, in some cases it is required that a fix must be at or close to
the top of the list of all fixes.  In those cases an instance of the
DUMMY fix style may be created by calling ``Modify::add_fix()`` and then
later replaced by calling ``Modify::replace_fix()``.

fix STORE/ATOM
""""""""""""""

Fix STORE/ATOM can be used as persistent storage of per-atom data.

**Syntax**

.. code-block:: LAMMPS

   fix ID group-ID STORE/ATOM N1 N2 gflag rflag

* ID, group-ID are documented in :doc:`fix <fix>` command
* STORE/ATOM = style name of this fix command
* N1 = 1, N2 = 0 : data is per-atom vector = single value per atom
* N1 > 1, N2 = 0 : data is per-atom array = N1 values per atom
* N1 > 0, N2 > 0 : data is per-atom tensor = N1xN2 values per atom
* gflag = 1 communicate per-atom values with ghost atoms, 0 do not update ghost atom data
* rflag = 1 store per-atom value in restart file, 0 do not store data in restart

Similar functionality is also available through using custom per-atom
properties with :doc:`fix property/atom <fix_property_atom>`.  The
choice between the two fixes should be based on whether the user should
be able to access this per-atom data: if yes, then fix property/atom is
preferred, otherwise fix STORE/ATOM.

fix STORE/GLOBAL
""""""""""""""""

Fix STORE/GLOBAL can be used as persistent storage of global data with support for restarts

**Syntax**

.. code-block:: LAMMPS

   fix ID group-ID STORE/GLOBAL N1 N2

* ID, group-ID are documented in :doc:`fix <fix>` command
* STORE/GLOBAL = style name of this fix command
* N1 >=1 : number of global items to store
* N2 = 1 : data is global vector of length N1
* N2 > 1 : data is global N1xN2 array

fix STORE/LOCAL
"""""""""""""""

Fix STORE/LOCAL can be used as persistent storage for local data

**Syntax**

.. code-block:: LAMMPS

   fix ID group-ID STORE/LOCAL Nreset Nvalues

* ID, group-ID are documented in :doc:`fix <fix>` command
* STORE/LOCAL = style name of this fix command
* Nreset = frequency at which local data is available
* Nvalues = number of values per local item, that is the number of columns
