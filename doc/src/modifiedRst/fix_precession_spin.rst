.. index:: fix precession/spin

fix precession/spin command
===========================

Syntax
""""""


.. parsed-literal::

   fix ID group precession/spin style args

* ID, group are documented in :doc:`fix <fix>` command
* precession/spin = style name of this fix command
* style = *zeeman* or *anisotropy* or *cubic*
  
  .. parsed-literal::
  
       *zeeman* args = H x y z
         H = intensity of the magnetic field (in Tesla)
         x y z = vector direction of the field
       *anisotropy* args = K x y z
         K = intensity of the magnetic anisotropy (in eV)
         x y z = vector direction of the anisotropy

  
  .. parsed-literal::
  
       *cubic* args = K1 K2c n1x n1y n1x n2x n2y n2z n3x n3y n3z 
         K1 and K2c = intensity of the magnetic anisotropy (in eV)
         n1x to n3z = three direction vectors of the cubic anisotropy



Examples
""""""""


.. parsed-literal::

   fix 1 all precession/spin zeeman 0.1 0.0 0.0 1.0
   fix 1 3 precession/spin anisotropy 0.001 0.0 0.0 1.0
   fix 1 iron precession/spin cubic 0.001 0.0005 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
   fix 1 all precession/spin zeeman 0.1 0.0 0.0 1.0 anisotropy 0.001 0.0 0.0 1.0

Description
"""""""""""

This fix applies a precession torque to each magnetic spin in the group.

Style *zeeman* is used for the simulation of the interaction
between the magnetic spins in the defined group and an external
magnetic field:

.. math source doc: src/Eqs/force_spin_zeeman.tex
.. math::

   :align: center

with mu0 the vacuum permeability, muB the Bohr magneton (muB = 5.788 eV/T
in metal units).

Style *anisotropy* is used to simulate an easy axis or an easy plane
for the magnetic spins in the defined group:

.. math source doc: src/Eqs/force_spin_aniso.tex
.. math::

   :align: center

with n defining the direction of the anisotropy, and K (in eV) its intensity.
If K>0, an easy axis is defined, and if K<0, an easy plane is defined.

Style *cubic* is used to simulate a cubic anisotropy, with three
possible easy axis for the magnetic spins in the defined group:

.. math source doc: src/Eqs/fix_spin_cubic.tex
.. math::

   :align: center

with K1 and K2c (in eV) the intensity coefficients and 
n1, n2 and n3 defining the three anisotropic directions
defined by the command (from n1x to n3z). 
For n1 = (100), n2 = (010), and n3 = (001), K1 < 0 defines an 
iron type anisotropy (easy axis along the (001)-type cube
edges), and K1 > 0 defines a nickel type anisotropy (easy axis
along the (111)-type cube diagonals). 
K2\^c > 0 also defines easy axis along the (111)-type cube
diagonals.
See chapter 2 of :ref:`(Skomski) <Skomski1>` for more details on cubic
anisotropies.

In all cases, the choice of (x y z) only imposes the vector
directions for the forces. Only the direction of the vector is 
important; it's length is ignored (the entered vectors are
normalized).

Those styles can be combined within one single command line.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

By default, the energy associated to this fix is not added to the potential
energy of the system.
The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this fix
to add this magnetic potential energy to the potential energy of the system,


.. parsed-literal::

   fix             1 all precession/spin zeeman 1.0 0.0 0.0 1.0
   fix_modify      1 energy yes

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.

No information about this fix is written to :doc:`binary restart files <restart>`.

Restrictions
""""""""""""


The *precession/spin* style is part of the SPIN package.  This style
is only enabled if LAMMPS was built with this package, and if the
atom\_style "spin" was declared.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`atom\_style spin <atom_style>`

**Default:** none


----------


.. _Skomski1:



**(Skomski)** Skomski, R. (2008). Simple models of magnetism.
Oxford University Press.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
