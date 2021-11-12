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
   pair_coeff 1 1 2.88 10.35	4.178	0.210	1.818	4.07293506	4.9883063257983666

Description
"""""""""""

The *smatb* styles compute the Second Moment Aproximation to the Tight Binding
:ref:`(Cyrot) <Cyrot>`, :ref:`(Gupta) <Gupta>`, :ref:`(Rosato) <Rosato>`,
given by

.. math::

      E_{ik} &= A e^{-p \left(\frac{R_{ik}}{R_{0}}-1\right)} -\sqrt{
            \sum_{j,R_{ij}\leq R_c}\Xi^2
            e^{-2q\left(\frac{R_{ij}}{R_{0}}-1\right)}}.\qquad R_{ij} < R_{sc}\\
      E_{ik} &=\left(a_3\left(R_{ik}-R_c\right)^3
            +a_4\left(R_{ik}-R_c\right)^4
            +a_5\left(R_{ik}-R_c\right)^5\right) -\sqrt{
            \sum_{j,R_{ij}\left(
            x_3\left(R_{ij}-R_c\right)^3
            +x_4\left(R_{ij}-R_c\right)^4
            +x_5\left(R_{ij}-R_c\right)^5
            \right)^2}. \qquad R_{sc} < R_{ij} < r_c

The polynomial coefficients a3, a4, a5, x3, x4, x5 are computed by LAMMPS smoothly
link the two exponential to zero from the inner cutoff :math:`R_{sc}` to the
outer cutoff :math:`R_c`.


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


Note that : when :math:`R_{cs} < R_{ij} < R_c` the exponentials are substituted
by a polynomial in the form:
:math:`P(R_{ij}) = P_3\left(R_{ij}-R_c\right)^3+P_4\left(R_{ij}-R_c\right)^4+P_5\left(R_{ij}-R_c\right)^5`,
where :math:`P_3`, :math:`P_4` nd :math:`P_5` are calculated at the
initialization of the potential in order to guarantee the continuity of the 
function, of the first and of the second derivateves in :math:`R_{cs}`,
moreover the vaule of :math:`P(R_{ij})` and of its first and second derivative
is automatically zero in the end of the cutoff :math:`R_c` due to its polynomial
form.

See the :doc:`run_style <run_style>` command for details.

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

