.. index:: pair_style smatb
.. index:: pair_style smatb/single

pair_style smatb command
=========================

pair_style smatb/single command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *smatb* or *smatb/single*
* args = none

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smatb
   pair_coeff 1 1 2.88 10.35 4.178 0.210 1.818 4.07293506 4.9883063257983666

   pair_style smatb/single
   pair_coeff 1 1 2.88 10.35 4.178 0.210 1.818 4.07293506 4.9883063257983666


Description
"""""""""""

.. versionadded:: 4May2022

The *smatb* and *smatb/single* styles compute the Second Moment
Approximation to the Tight Binding :ref:`(Cyrot) <Cyrot>`,
:ref:`(Gupta) <Gupta>`, :ref:`(Rosato) <Rosato>`, given by

.. math::
      E_{i}  = \sum_{j,R_{ij}\leq R_{c}} \alpha(R_{ij}) - \sqrt{\sum_{j,R_{ij}\leq R_{c}}\Xi^2(R_{ij})}

:math:`R_{ij}` is the distance between the atom :math:`i` and :math:`j`.
And the two functions :math:`\alpha\left(r\right)` and :math:`\Xi\left(r\right)` are:

.. math::
   \alpha\left(r\right)=\left\lbrace\begin{array}{ll}
      A e^{-p \left(\frac{r}{R_{0}}-1\right)} & r < R_{sc}\\
      a_3\left(r-R_{c}\right)^3+a_4\left(r-R_{c}\right)^4
      +a_5\left(r-R_{c}\right)^5& R_{sc} < r < R_{c}
      \end{array}
      \right.

.. math::
      \Xi\left(r\right)=\left\lbrace\begin{array}{ll}
      \xi e^{-q \left(\frac{r}{R_{0}}-1\right)} & r < R_{sc}\\
      x_3\left(r-R_{c}\right)^3+x_4\left(r-R_{c}\right)^4
      +x_5\left(r-R_{c}\right)^5& R_{sc} < r < R_{c}
      \end{array}
      \right.


The polynomial coefficients :math:`a_3`, :math:`a_4`, :math:`a_5`,
:math:`x_3`, :math:`x_4`, :math:`x_5` are computed by LAMMPS: the two
exponential terms and their first and second derivatives are smoothly
reduced to zero, from the inner cutoff :math:`R_{sc}` to the outer
cutoff :math:`R_{c}`.

The *smatb/single* style is an optimization when using only a single atom type.

Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`R_{0}` (distance units)
* :math:`p` (dimensionless)
* :math:`q` (dimensionless)
* :math:`A` (energy units)
* :math:`\xi` (energy units)
* :math:`R_{cs}` (distance units)
* :math:`R_{c}` (distance units)


Note that: :math:`R_{0}` is the nearest neighbor distance, usually coincides
with the diameter of the atoms

See the :doc:`run_style <run_style>` command for details.

----------

Mixing info
"""""""""""

For atom type pairs I,J and I != J the coefficients are not automatically mixed.

----------

Restrictions
""""""""""""

These pair styles are part of the SMTBQ package and are only enabled
if LAMMPS is built with that package.  See the :doc:`Build package <Build_package>` page for more info.

These pair styles require the :doc:`newton <newton>` setting to be "on" for pair interactions.

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Cyrot:

**(Cyrot)**  Cyrot-Lackmann and Ducastelle, Phys Rev. B, 4, 2406-2412 (1971).

.. _Gupta:

**(Gupta)** Gupta ,Phys Rev. B, 23, 6265-6270 (1981).

.. _Rosato:

**(Rosato)** Rosato and Guillope  and Legrand, Philosophical Magazine A, 59.2, 321-336 (1989).

