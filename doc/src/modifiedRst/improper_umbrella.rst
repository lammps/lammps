.. index:: improper\_style umbrella

improper\_style umbrella command
================================

improper\_style umbrella/omp command
====================================

Syntax
""""""


.. parsed-literal::

   improper_style umbrella

Examples
""""""""


.. parsed-literal::

   improper_style umbrella
   improper_coeff 1 100.0 180.0

Description
"""""""""""

The *umbrella* improper style uses the following potential, which is
commonly referred to as a classic inversion and used in the
:doc:`DREIDING <Howto_bioFF>` force field:

.. math source doc: src/Eqs/improper_umbrella.tex
.. math::

   E=\frac{1}{2}K\left( \frac{1}{\sin\omega_0}\right) ^2 \left( \cos\omega - \cos\omega_0\right) ^2 \qquad \omega_0 \neq 0^o

>>>image was here
   E=K\left( 1-cos\omega\right)  \qquad \omega_0 = 0^o


where K is the force constant and omega is the angle between the IL
axis and the IJK plane:

.. image:: Eqs/umbrella.jpg
   :align: center

If omega0 = 0 the potential term has a minimum for the planar
structure.  Otherwise it has two minima at +/- omega0, with a barrier
in between.

See :ref:`(Mayo) <umbrella-Mayo>` for a description of the DREIDING force field.

The following coefficients must be defined for each improper type via
the :doc:`improper\_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* K (energy)
* omega0 (degrees)


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`improper\_coeff <improper_coeff>`

**Default:** none


----------


.. _umbrella-Mayo:



**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990),


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
