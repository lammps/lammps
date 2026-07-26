.. index:: fix uvt

fix uvt command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID uvt temp Tstart Tstop Tdamp mu Mustart Mustop Mudamp ne Ne0 dedn source [keyword value ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ``temp`` values = Tstart Tstop Tdamp
* ``mu`` values = Mustart Mustop Mudamp
* ``ne`` value = initial electronic coordinate
* ``dedn`` value = global scalar or vector reference used as dE/dN

Optional keywords are the Nose-Hoover keywords supported by
:doc:`fix nvt <fix_nh>`, plus ``ne_velocity`` and ``dedn_defer``.

Examples
""""""""

.. code-block:: LAMMPS

   variable dEdN equal 5.0*(f_cp[13]-1.0)
   fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 2.0 0.5 ne 1.8 dedn v_dEdN

Description
"""""""""""

``fix uvt`` performs classical molecular dynamics at constant chemical
potential by adding one global electronic coordinate, ``Ne``, to the
standard Nose-Hoover NVT equations of motion.  The electronic coordinate
is propagated with an extended-mass variable and is driven by
``-dE/dN + mu``.

The ``dedn`` source may be an equal-style variable (``v_name``), a
global compute (``c_ID``), or a global fix (``f_ID``).  If the source
provides a global vector, an entry can be selected with the usual
``[index]`` syntax.

The global fix vector appends six entries after the regular
Nose-Hoover vector entries: ``Ne``, ``Ne_dot``, ``dEdN``, ``mu``,
the electronic kinetic energy, and the electronic potential contribution
``-mu*Ne``.

The fix writes restart data for the Nose-Hoover state and electronic
degree of freedom, so simulations may be continued with
:doc:`read_restart <read_restart>` followed by the same ``fix uvt``
command.

Restrictions
""""""""""""

This fix supports temperature control only.  Pressure control keywords
are not allowed.

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix pimd/uvt <fix_pimd_uvt>`,
:doc:`read_restart <read_restart>`
