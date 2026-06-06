fix uvt command
================

.. versionadded:: 6Jun2026

Syntax
""""""

.. parsed-literal::

   fix ID group-ID uvt temp Tstart Tstop Tdamp [keyword value ...]

where the accepted thermostat keywords are the same as for :doc:`fix nvt <fix_nh>`,
and the classical constant-potential keywords are::

   keyword = mu or ne or dedn or ne_velocity

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all uvt temp 300.0 300.0 100.0 mu 0.0 0.0 100.0 ne 1.0 dedn v_dEdN

Description
"""""""""""

``fix uvt`` is a thin wrapper around the coupled classical uVT
implementation in :doc:`fix nvt <fix_nh>`.  It uses the same
Nose-Hoover chain machinery as ``fix nvt`` but requires the chemical
potential control keywords ``mu``, ``ne``, and ``dedn`` so that the
electron-number degree of freedom is propagated together with the
particle thermostat.

The ``dedn`` value may be supplied as an equal-style variable
(``v_name``), a global compute (``c_ID``), or a global fix
(``f_ID``).  If the source provides a global vector, the first
component can be selected with the usual ``[index]`` syntax, e.g.
``c_model[1]`` or ``f_state[2]``.

The thermostat keywords follow the same meaning as in
:doc:`fix nvt <fix_nh>`.  The ``mu``, ``ne``, ``dedn``, and
``ne_velocity`` keywords are specific to ``fix uvt``.

Restrictions
""""""""""""

``fix uvt`` requires temperature control and cannot be combined with
pressure control.
