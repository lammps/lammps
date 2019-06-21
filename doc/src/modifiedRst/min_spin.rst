.. index:: min\_style spin

min\_style spin command
=======================

Syntax
""""""


.. parsed-literal::

   min_style spin

Examples
""""""""


.. parsed-literal::

   min_style 	spin

Description
"""""""""""

Apply a minimization algorithm to use when a :doc:`minimize <minimize>`
command is performed.

Style *spin* defines a damped spin dynamics with an adaptive
timestep, according to:

.. math::

   :align: center

with lambda a damping coefficient (similar to a Gilbert
damping).
Lambda can be defined by setting the *alpha\_damp* keyword with the 
:doc:`min\_modify <min_modify>` command.

The minimization procedure solves this equation using an
adaptive timestep. The value of this timestep is defined 
by the largest precession frequency that has to be solved in the 
system:

.. math::

   :align: center

with *\|omega\|\_\ *max*\ * the norm of the largest precession frequency
in the system (across all processes, and across all replicas if a
spin/neb calculation is performed).

Kappa defines a discretization factor *discrete\_factor* for the 
definition of this timestep. 
*discrete\_factor* can be defined with the :doc:`min\_modify <min_modify>`
command.

.. note::

   The *spin* style replaces the force tolerance by a torque
   tolerance. See :doc:`minimize <minimize>` for more explanation.

Restrictions
""""""""""""


This minimization procedure is only applied to spin degrees of
freedom for a frozen lattice configuration.

Related commands
""""""""""""""""

:doc:`min\_style <min_style>`, :doc:`minimize <minimize>`, 
:doc:`min\_modify <min_modify>`

Default
"""""""

The option defaults are *alpha\_damp* = 1.0 and *discrete\_factor* =
10.0.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
