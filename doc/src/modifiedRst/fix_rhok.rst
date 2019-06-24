.. index:: fix rhok

fix rhok command
================


.. parsed-literal::

   fix ID group-ID rhok nx ny nz K a

* ID, group-ID are documented in :doc:`fix <fix>` command
* nx, ny, nz = k-vector of collective density field
* K = spring constant of bias potential
* a = anchor point of bias potential

Examples
""""""""


.. parsed-literal::

   fix bias all rhok 16 0 0 4.0 16.0
   fix 1 all npt temp 0.8 0.8 4.0 z 2.2 2.2 8.0
   # output of 4 values from fix rhok: U_bias rho_k_RE  rho_k_IM  \|rho_k\|
   thermo_style custom step temp pzz lz f_bias f_bias[1] f_bias[2] f_bias[3]

Description
"""""""""""

The fix applies a force to atoms given by the potential

.. math source doc: src/Eqs/fix_rhok.tex
.. math::

   U &=&  \frac{1}{2} K (|\rho_{\vec{k}}| - a)^2 \\
   \rho_{\vec{k}} &=& \sum_j^N \exp(-i\vec{k} \cdot \vec{r}_j )/\sqrt{N} \\
   \vec{k} &=& (2\pi n_x /L_x , 2\pi n_y  /L_y , 2\pi n_z/L_z ) 


as described in :ref:`(Pedersen) <Pedersen>`.

This field, which biases configurations with long-range order, can be
used to study crystal-liquid interfaces and determine melting
temperatures :ref:`(Pedersen) <Pedersen>`.

An example of using the interface pinning method is located in the
*examples/USER/misc/rhok* directory.

Restrictions
""""""""""""


This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`thermo\_style <thermo_style>`

**Default:** none


----------


.. _Pedersen:



**(Pedersen)** Pedersen, J. Chem. Phys., 139, 104102 (2013).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
