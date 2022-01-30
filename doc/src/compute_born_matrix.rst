.. index:: compute born/matrix

compute born/matrix command
===========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID born/matrix

* ID, group-ID are documented in :doc:`compute <compute>` command
* born/matrix = style name of this compute command
* the keyword *numdiff* may be appended

    .. parsed-literal::

       *numdiff* values = delta virial-ID
  
Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all born/matrix
   compute 1 all born/matrix numdiff 1.0e-4 myvirial

Description
"""""""""""

Define a compute that calculates
:math:`\frac{\partial{}^2U}{\partial\varepsilon_{i}\partial\varepsilon_{j}}` the
second derivatives of the potential energy :math:`U` with regard to strain
tensor :math:`\varepsilon` elements. These values are related to:

.. math::

    C^{B}_{i,j}=\frac{1}{V}\frac{\partial{}^2U}{\partial{}\varepsilon_{i}\partial\varepsilon_{j}}

also called the Born term of elastic constants in the stress-stress fluctuation
formalism.  This quantity can be used to compute the elastic constant tensor.
Using the symmetric Voigt notation, the elastic constant tensor can be written
as a 6x6 symmetric matrix:

.. math::

    C_{i,j} = \langle{}C^{B}_{i,j}\rangle
             + \frac{V}{k_{B}T}\left(\langle\sigma_{i}\sigma_{j}\rangle\right.
             \left.- \langle\sigma_{i}\rangle\langle\sigma_{j}\rangle\right)
             + \frac{Nk_{B}T}{V}
               \left(\delta_{i,j}+(\delta_{1,i}+\delta_{2,i}+\delta_{3,i})\right.
               \left.*(\delta_{1,j}+\delta_{2,j}+\delta_{3,j})\right)

In the above expression, :math:`\sigma` stands for the virial stress
tensor, :math:`\delta` is the Kronecker delta and the usual notation apply for
the number of particle, the temperature and volume respectively :math:`N`,
:math:`T` and :math:`V`. :math:`k_{B}` is the Boltzmann constant.

The Born term is a symmetric 6x6 matrix by construction and as such can be
expressed as 21 independent terms. The terms are ordered corresponding to the
following matrix element:

.. math::

    \begin{matrix}
       C_{1}  & C_{7}   & C_{8}  & C_{9}  & C_{10} & C_{11} \\
       C_{7}  & C_{2}   & C_{12} & C_{13} & C_{14} & C_{15} \\
       \vdots & C_{12}  & C_{3}  & C_{16} & C_{17} & C_{18} \\
       \vdots & C_{13}  & C_{16} & C_{4}  & C_{19} & C_{20} \\
       \vdots & \vdots  & \vdots & C_{19} & C_{5}  & C_{21} \\
       \vdots & \vdots  & \vdots & \vdots & C_{21} & C_{6}
    \end{matrix}

in this matrix the indices of :math:`C_{k}` value are the corresponding index
:math:`k` in the compute output. Each term comes from the sum of every
interactions derivatives in the system as explained in :ref:`(VanWorkum)
<VanWorkum>` or :ref:`(Voyiatzis) <Voyiatzis>`.

The output can be accessed using usual Lammps routines:

.. code-block:: LAMMPS

   compute 1 all born/matrix
   compute 2 all pressure NULL virial
   variable S1 equal -c_2[1]
   variable S2 equal -c_2[2]
   variable S3 equal -c_2[3]
   variable S4 equal -c_2[4]
   variable S5 equal -c_2[5]
   variable S6 equal -c_2[6]
   fix 1 all ave/time 1 1 1 v_S1 v_S2 v_S3 v_S4 v_S5 v_S6 c_1[*] file born.out

In this example, the file *born.out* will contain the information needed to
compute the first and second terms of the elastic constant matrix in a post
processing procedure. The other required quantities can be accessed using any
other *LAMMPS* usual method.

NOTE: In the above :math:`C_{i,j}` computation, the term involving the virial
stress tensor :math:`\sigma` is the covariance between each elements. In a
solid the virial stress can have large variations between timesteps and average
values can be slow to converge. This term is better computed using
instantaneous values.

The *numdiff* keyword uses finite differences of energy to numerically
approximate the derivative. This is useful when using interaction styles
for which the analytical derivatives have not been implemented.
The keyword requirs the additional values *delta* and *virial-ID*
giving the size of the applied strain and the ID of the pressure compute
that provides the virial tensor, requiring that it use the virial
keyword e.g.

.. code-block:: LAMMPS

   compute myvirial all pressure NULL virial
   compute 1 all born/matrix numdiff 1.0e-4 myvirial


**Output info:**

This compute calculates a global array with the number of rows=21.
The values are ordered as explained above. These values can be used
by any command that uses a global values from a compute as input. See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "extensive".

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

The Born term can be decomposed as a product of two terms. The first one
is a general term which depends on the configuration. The second one is
specific to every interaction composing your force field (non-bonded,
bonds, angle...). Currently not all LAMMPS interaction styles
implement the *born_matrix*
method giving first and second order derivatives and LAMMPS will
exit with an error if this compute is used with such interactions,
unless the *numdiff* option is also used.

Default
"""""""

none

----------

.. _VanWorkum:

**(Van Workum)** K. Van Workum et al., J. Chem. Phys. 125 144506 (2006)

.. _Voyiatzis:

**(Voyiatzis)** E. Voyiatzis, Computer Physics Communications 184(2013)27-33
