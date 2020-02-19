.. index:: pair\_style born/coul/dsf/cs

pair\_style born/coul/dsf/cs command
====================================

pair\_style born/coul/long/cs command
=====================================

pair\_style born/coul/long/cs/gpu command
=========================================

pair\_style born/coul/wolf/cs command
=====================================

pair\_style born/coul/wolf/cs/gpu command
=========================================

pair\_style buck/coul/long/cs command
=====================================

pair\_style coul/long/cs command
================================

pair\_style coul/long/cs/gpu command
====================================

pair\_style coul/wolf/cs command
================================

pair\_style lj/cut/coul/long/cs command
=======================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *born/coul/dsf/cs* or *born/coul/long/cs* or *born/coul/wolf/cs* or *buck/coul/long/cs* or *coul/long/cs* or *coul/wolf/cs* or *lj/cut/coul/long/cs*
* args = list of arguments for a particular style


.. parsed-literal::

     *born/coul/dsf/cs* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (distance units)
     *born/coul/long/cs* args = cutoff (cutoff2)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *born/coul/wolf/cs* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for Buckingham (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *buck/coul/long/cs* args = cutoff (cutoff2)
       cutoff = global cutoff for Buckingham (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *coul/long* args = cutoff
       cutoff = global cutoff for Coulombic (distance units)
     *coul/wolf* args = alpha cutoff
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for Coulombic (distance units)
     *lj/cut/coul/long/cs* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style born/coul/dsf/cs 0.1 10.0 12.0
   pair_coeff \* \*   0.0 1.00 0.00 0.00 0.00
   pair_coeff 1 1 480.0 0.25 0.00 1.05 0.50

   pair_style born/coul/long/cs 10.0 8.0
   pair_coeff 1 1 6.08 0.317 2.340 24.18 11.51

   pair_style born/coul/wolf/cs 0.25 10.0 12.0
   pair_coeff \* \*   0.0 1.00 0.00 0.00 0.00
   pair_coeff 1 1 480.0 0.25 0.00 1.05 0.50

   pair_style buck/coul/long/cs 10.0
   pair_style buck/coul/long/cs 10.0 8.0
   pair_coeff \* \* 100.0 1.5 200.0
   pair_coeff 1 1 100.0 1.5 200.0 9.0

   pair_style coul/long/cs 10.0
   pair_coeff \* \*

   pair_style coul/wolf/cs 0.2 9.0
   pair_coeff \* \*

   pair_style lj/cut/coul/long/cs 10.0
   pair_style lj/cut/coul/long/cs 10.0 8.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

Description
"""""""""""

These pair styles are designed to be used with the adiabatic
core/shell model of :ref:`(Mitchell and Finchham) <MitchellFinchham2>`.  See
the :doc:`Howto coreshell <Howto_coreshell>` doc page for an overview of
the model as implemented in LAMMPS.

All the styles are identical to the corresponding pair style without
the "/cs" in the name:

* :doc:`pair_style born/coul/dsf <pair_born>`
* :doc:`pair_style born/coul/long <pair_born>`
* :doc:`pair_style born/coul/wolf <pair_born>`
* :doc:`pair_style buck/coul/long <pair_buck>`
* :doc:`pair_style coul/long <pair_coul>`
* :doc:`pair_style coul/wolf <pair_coul>`
* :doc:`pair_style lj/cut/coul/long <pair_lj>`

except that they correctly treat the special case where the distance
between two charged core and shell atoms in the same core/shell pair
approach r = 0.0.

Styles with a "/long" in the name are used with a long-range solver
for Coulombic interactions via the :doc:`kspace_style <kspace_style>`
command.  They require special treatment of the short-range Coulombic
interactions within the cor/shell model.

Specifically, the short-range Coulomb interaction between a core and
its shell should be turned off using the
:doc:`special_bonds <special_bonds>` command by setting the 1-2 weight
to 0.0, which works because the core and shell atoms are bonded to
each other.  This induces a long-range correction approximation which
fails at small distances (~< 10e-8). Therefore, the Coulomb term which
is used to calculate the correction factor is extended by a minimal
distance (r\_min = 1.0-6) when the interaction between a core/shell
pair is treated, as follows

.. image:: Eqs/pair_cs.jpg
   :align: center

where C is an energy-conversion constant, Qi and Qj are the charges on
the core and shell, epsilon is the dielectric constant and r\_min is the
minimal distance.

For styles that are not used with a long-range solver, i.e. those with
"/dsf" or "/wolf" in the name, the only correction is the addition of
a minimal distance to avoid the possible r = 0.0 case for a core/shell
pair.


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


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

See the corresponding doc pages for pair styles without the "cs"
suffix to see how mixing, shifting, tabulation, tail correction,
restarting, and rRESPA are handled by theses pair styles.


----------


Restrictions
""""""""""""


These pair styles are part of the CORESHELL package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style born <pair_born>`,
:doc:`pair_style buck <pair_buck>`

**Default:** none


----------


.. _MitchellFinchham2:



**(Mitchell and Finchham)** Mitchell, Finchham, J Phys Condensed Matter,
5, 1031-1038 (1993).
