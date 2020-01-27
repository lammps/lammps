.. index:: fix lb/momentum

fix lb/momentum command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID lb/momentum nevery keyword values ...

* ID, group-ID are documented in the :doc:`fix <fix>` command
* lb/momentum = style name of this fix command
* nevery = adjust the momentum every this many timesteps
* zero or more keyword/value pairs may be appended
* keyword = *linear*
  
  .. parsed-literal::
  
       *linear* values = xflag yflag zflag
         xflag,yflag,zflag = 0/1 to exclude/include each dimension.



Examples
""""""""


.. parsed-literal::

   fix 1 sphere lb/momentum
   fix 1 all lb/momentum linear 1 1 0

Description
"""""""""""

This fix is based on the :doc:`fix momentum <fix_momentum>` command, and
was created to be used in place of that command, when a
lattice-Boltzmann fluid is present.

Zero the total linear momentum of the system, including both the atoms
specified by group-ID and the lattice-Boltzmann fluid every nevery
timesteps.  This is accomplished by adjusting the particle velocities
and the fluid velocities at each lattice site.

.. note::

   This fix only considers the linear momentum of the system.

By default, the subtraction is performed for each dimension.  This can
be changed by specifying the keyword *linear*\ , along with a set of
three flags set to 0/1 in order to exclude/ include the corresponding
dimension.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


Can only be used if a lattice-Boltzmann fluid has been created via the
:doc:`fix lb/fluid <fix_lb_fluid>` command, and must come after this
command.

This fix is part of the USER-LB package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix momentum <fix_momentum>`, :doc:`fix lb/fluid <fix_lb_fluid>`

Default
"""""""

Zeros the total system linear momentum in each dimension.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
