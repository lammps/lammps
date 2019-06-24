.. index:: pair\_style sph/taitwater/morris

pair\_style sph/taitwater/morris command
========================================

Syntax
""""""


.. parsed-literal::

   pair_style sph/taitwater/morris

Examples
""""""""


.. parsed-literal::

   pair_style sph/taitwater/morris
   pair_coeff \* \* 1000.0 1430.0 1.0 2.4

Description
"""""""""""

The sph/taitwater/morris style computes pressure forces between SPH
particles according to Tait's equation of state:

.. math source doc: src/Eqs/pair_sph_tait.tex
.. math::

   p = B [(\frac{\rho}{\rho_0})^{\gamma} - 1]


where gamma = 7 and B = c\_0\^2 rho\_0 / gamma, with rho\_0 being the
reference density and c\_0 the reference speed of sound.

This pair style also computes laminar viscosity :ref:`(Morris) <Morris>`.

See `this PDF guide <USER/sph/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above.

* rho0 reference density (mass/volume units)
* c0 reference soundspeed (distance/time units)
* nu dynamic viscosity (mass\*distance/time units)
* h kernel function cutoff (distance units)


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This style does not support the :doc:`pair\_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.  Thus, you need to re-specify the pair\_style and
pair\_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run\_style respa <run_style>` command.  It does not support the *inner*\ ,
*middle*\ , *outer* keywords.

Restrictions
""""""""""""


This pair style is part of the USER-SPH package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, pair\_sph/rhosum

**Default:** none


----------


.. _Morris:



**(Morris)** Morris, Fox, Zhu, J Comp Physics, 136, 214-226 (1997).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
