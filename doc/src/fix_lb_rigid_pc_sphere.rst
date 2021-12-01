.. index:: fix lb/rigid/pc/sphere

fix lb/rigid/pc/sphere command
==============================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID lb/rigid/pc/sphere bodystyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/rigid/pc/sphere = style name of this fix command
* bodystyle = *single* or *molecule* or *group*

  .. parsed-literal::

       *single* args = none
       *molecule* args = none
       *group* args = N groupID1 groupID2 ...
         N = # of groups

* zero or more keyword/value pairs may be appended
* keyword = *force* or *torque* or *innerNodes*

  .. parsed-literal::

       *force* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass force is active
       *torque* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass torque is active
       *innerNodes* values = innergroup-ID
         innergroup-ID = ID of the atom group which does not experience a hydrodynamic force from the lattice-Boltzmann fluid

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 spheres lb/rigid/pc/sphere
   fix 1 all lb/rigid/pc/sphere force 1 0 0 innerNodes ForceAtoms

Description
"""""""""""

This fix was part of the old LATBOLTZ package and is now defunct.  LAMMPS standard :doc:`fix rigid <fix_rigid>`  can be used in its place.


