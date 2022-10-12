.. index:: fix drude

fix drude command
=================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID drude flag1 flag2 ... flagN

* ID, group-ID are documented in :doc:`fix <fix>` command
* drude = style name of this fix command
* flag1 flag2 ... flagN = Drude flag for each atom type (1 to N) in the system

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all drude 1 1 0 1 0 2 2 2
   fix 1 all drude C C N C N D D D

Example input scripts available: examples/PACKAGES/drude

Description
"""""""""""

Assign each atom type in the system to be one of 3 kinds of atoms
within the Drude polarization model. This fix is designed to be used
with the :doc:`thermalized Drude oscillator model <Howto_drude>`.
Polarizable models in LAMMPS are described on the :doc:`Howto polarizable <Howto_polarizable>` doc page.

The three possible types can be designated with an integer (0,1,2)
or capital letter (N,C,D):

* 0 or N = non-polarizable atom (not part of Drude model)
* 1 or C = Drude core
* 2 or D = Drude electron

Restrictions
""""""""""""

This fix should be invoked before any other commands that implement
the Drude oscillator model, such as :doc:`fix langevin/drude <fix_langevin_drude>`, :doc:`fix tgnvt/drude <fix_tgnh_drude>`, :doc:`fix drude/transform <fix_drude_transform>`, :doc:`compute temp/drude <compute_temp_drude>`, :doc:`pair_style thole <pair_thole>`.

Related commands
""""""""""""""""

:doc:`fix langevin/drude <fix_langevin_drude>`, :doc:`fix tgnvt/drude <fix_tgnh_drude>`, :doc:`fix drude/transform <fix_drude_transform>`, :doc:`compute temp/drude <compute_temp_drude>`, :doc:`pair_style thole <pair_thole>`

Default
"""""""

none
