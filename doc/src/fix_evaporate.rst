.. index:: fix evaporate

fix evaporate command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID evaporate N M region-ID seed

* ID, group-ID are documented in :doc:`fix <fix>` command
* evaporate = style name of this fix command
* N = delete atoms every this many timesteps
* M = number of atoms to delete each time
* region-ID = ID of region within which to perform deletions
* seed = random number seed to use for choosing atoms to delete
* zero or more keyword/value pairs may be appended
  
  .. parsed-literal::
  
     keyword = *molecule*
       *molecule* value = *no* or *yes*



Examples
""""""""


.. parsed-literal::

   fix 1 solvent evaporate 1000 10 surface 49892
   fix 1 solvent evaporate 1000 10 surface 38277 molecule yes

Description
"""""""""""

Remove M atoms from the simulation every N steps.  This can be used,
for example, to model evaporation of solvent particles or molecules
(i.e. drying) of a system.  Every N steps, the number of atoms in the
fix group and within the specified region are counted.  M of these are
chosen at random and deleted.  If there are less than M eligible
particles, then all of them are deleted.

If the setting for the *molecule* keyword is *no*\ , then only single
atoms are deleted.  In this case, you should insure you do not delete
only a portion of a molecule (only some of its atoms), or LAMMPS will
soon generate an error when it tries to find those atoms.  LAMMPS will
warn you if any of the atoms eligible for deletion have a non-zero
molecule ID, but does not check for this at the time of deletion.

If the setting for the *molecule* keyword is *yes*\ , then when an atom
is chosen for deletion, the entire molecule it is part of is deleted.
The count of deleted atoms is incremented by the number of atoms in
the molecule, which may make it exceed *M*\ .  If the molecule ID of the
chosen atom is 0, then it is assumed to not be part of a molecule, and
just the single atom is deleted.

As an example, if you wish to delete 10 water molecules every *N*
steps, you should set *M* to 30.  If only the water's oxygen atoms
were in the fix group, then two hydrogen atoms would be deleted when
an oxygen atom is selected for deletion, whether the hydrogen atoms
are inside the evaporation region or not.

Note that neighbor lists are re-built on timesteps that atoms are
removed.  Thus you should not remove atoms too frequently or you will
incur overhead due to the cost of building neighbor lists.

.. note::

   If you are monitoring the temperature of a system where the atom
   count is changing due to evaporation, you typically should use the
   :doc:`compute_modify dynamic yes <compute_modify>` command for the
   temperature compute you are using.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global scalar, which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
number of deleted atoms.  The scalar value calculated by this fix is
"intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix deposit <fix_deposit>`

Default
"""""""

The option defaults are molecule = no.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
