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

* style = *smatb*
* args = none

.. parsed-literal::

     *smatb*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smatb
   pair_coeff 1 1 2.88 10.35 4.178 0.210 1.818 4.07293506 4.9883063257983666


Description
"""""""""""

The *smatb* styles compute the Second Moment Approximation to the Tight Binding
:ref:`(Cyrot) <Cyrot>`, :ref:`(Gupta) <Gupta>`, :ref:`(Rosato) <Rosato>`,
given by

.. math::

      E_{ik}  =\left\lbrace\begin{array}{ll}
       A e^{-p \left(\frac{R_{ik}}{R_{0}}-1\right)}
      -\sqrt{\sum_{j,R_{ij}\leq R_c}\Xi^2e^{-2q\left(\frac{R_{ij}}{R_{0}}-1\right)}}& R_{ij} < R_{sc}\\
      {\left(a_3\left(R_{ik}-R_c\right)^3+a_4\left(R_{ik}-R_c\right)^4
      +a_5\left(R_{ik}-R_c\right)^5\right)
      -\sqrt{\sum_{j,R_{ij}\leq R_c}\left(x_3\left(R_{ij}-R_c\right)^3
      +x_4\left(R_{ij}-R_c\right)^4+x_5\left(R_{ij}-R_c\right)^5\right)^2}} & R_{sc} < R_{ij} < r_c
      \end{array}
      \right.

The polynomial coefficients :math:`a_3`, :math:`a_4`, :math:`a_5`, :math:`x_3`,
:math:`x_4`, :math:`x_5` are computed by LAMMPS smoothly
link the two exponentials and their first and second order derivatives to zero 
from the inner cutoff :math:`R_{sc}` to the outer cutoff :math:`R_c`.


Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`R_{0}` (distance units)
* :math:`p` 
* :math:`q`
* :math:`A` (energy units)
* :math:`\Xi` (energy units)
* :math:`R_{cs}` (distance units) 
* :math:`R_c` (distance units)


Note that: :math:`R_{0}` is the nearest neighbour distance, usually coincides
with the diameter of the atoms

See the :doc:`run_style <run_style>` command for details.

----------

Mixing info
"""""""""""

For atom type pairs I,J and I != J the coefficients are not automatically mixed.

----------

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

