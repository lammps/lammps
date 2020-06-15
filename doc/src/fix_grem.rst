.. index:: fix grem

fix grem command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID grem lambda eta H0 thermostat-ID

* ID, group-ID are documented in :doc:`fix <fix>` command
* grem = style name of this fix command
* lambda = intercept parameter of linear effective temperature function
* eta = slope parameter of linear effective temperature function
* H0 = shift parameter of linear effective temperature function
* thermostat-ID = ID of Nose-Hoover thermostat or barostat used in simulation

Examples
""""""""

.. code-block:: LAMMPS

   fix             fxgREM all grem 400 -0.01 -30000 fxnpt
   thermo_modify   press fxgREM_press

   fix             fxgREM all grem 502 -0.15 -80000 fxnvt

Description
"""""""""""

This fix implements the molecular dynamics version of the generalized
replica exchange method (gREM) originally developed by :ref:`(Kim) <Kim2010>`,
which uses non-Boltzmann ensembles to sample over first order phase
transitions. The is done by defining replicas with an enthalpy
dependent effective temperature

.. math::

  T_{eff} = \lambda + \eta (H - H_0)

with :math:`\eta` negative and steep enough to only intersect the
characteristic microcanonical temperature (Ts) of the system once,
ensuring a unimodal enthalpy distribution in that replica.
:math:`\lambda` is the intercept and effects the generalized ensemble
similar to how temperature effects a Boltzmann ensemble. :math:`H_0`
is a reference enthalpy, and is typically set as the lowest desired
sampled enthalpy.  Further explanation can be found in our recent
papers :ref:`(Malolepsza) <Malolepsza>`.

This fix requires a Nose-Hoover thermostat fix reference passed to the
grem as *thermostat-ID*\ . Two distinct temperatures exist in this
generalized ensemble, the effective temperature defined above, and a
kinetic temperature that controls the velocity distribution of
particles as usual. Either constant volume or constant pressure
algorithms can be used.

The fix enforces a generalized ensemble in a single replica
only. Typically, this ideology is combined with replica exchange with
replicas differing by :math:`\lambda` only for simplicity, but this is not
required. A multi-replica simulation can be run within the LAMMPS
environment using the :doc:`temper/grem <temper_grem>` command. This
utilizes LAMMPS partition mode and requires the number of available
processors be on the order of the number of desired replicas. A
100-replica simulation would require at least 100 processors (1 per
world at minimum). If many replicas are needed on a small number of
processors, multi-replica runs can be run outside of LAMMPS.  An
example of this can be found in examples/USER/misc/grem and has no
limit on the number of replicas per processor. However, this is very
inefficient and error prone and should be avoided if possible.

In general, defining the generalized ensembles is unique for every
system. When starting a many-replica simulation without any knowledge
of the underlying microcanonical temperature, there are several tricks
we have utilized to optimize the process.  Choosing a less-steep
:math:`\eta` yields broader distributions, requiring fewer replicas to
map the microcanonical temperature.  While this likely struggles from
the same sampling problems gREM was built to avoid, it provides quick
insight to Ts.  Initially using an evenly-spaced :math:`\lambda`
distribution identifies regions where small changes in enthalpy lead
to large temperature changes. Replicas are easily added where needed.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`thermo_modify <thermo_modify>` *press* option is supported
by this fix to add the rescaled kinetic pressure as part of
:doc:`thermodynamic output <thermo_style>`.

Restrictions
""""""""""""

This fix is part of the USER-MISC package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`temper/grem <temper_grem>`, :doc:`fix nvt <fix_nh>`, :doc:`fix npt <fix_nh>`, :doc:`thermo_modify <thermo_modify>`

**Default:** none

----------

.. _Kim2010:

**(Kim)** Kim, Keyes, Straub, J Chem. Phys, 132, 224107 (2010).

.. _Malolepsza:

**(Malolepsza)** Malolepsza, Secor, Keyes, J Phys Chem B 119 (42),
13379-13384 (2015).
