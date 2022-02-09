.. index:: fix electrode/thermo

fix electrode/thermo command
============================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID charge_udpate potential eta temp T_v tau_v keyword values ...

* ID = name of fix
* group-ID = name of group fix is applied to
* potential = electric potential in Volts
* eta = reciprocal width of electrode charge densities
* T_v temperature parameter of thermo-potentiostat
* tau_v time parameter of thermo-potentiostat

.. parsed-literal::

    *symm(etry) on/off*
        turn on/off charge neutrality constraint
    *couple <group-ID> <potential>*
        group-ID = add another group of electrode atoms
        potential = electric potential in volts applied to this electrode
    *etypes <types>*
        specify atom types exclusive to the electrode for optimized neighbor
        lists
    *ffield on/off*
        turn on/off finite-field implementation
    *write_mat <filename>*
        write elastance matrix to file
    *write_inv <filename>*
        write inverted matrix to file
    *read_mat <filename>*
        read elastance matrix from file
    *read_inv <filename>*
        read inverted matrix from file


Examples
""""""""

.. code-block:: LAMMPS

   fix fxthermo electrodes electrode/thermo 0.0 1.805 symm on
   fix fxthermo bot electrode/thermo -1.0 1.805 couple top 1.0 write_inv inv.csv symm on

Description
"""""""""""

Lorem ipsum

.. code-block:: LAMMPS

   kspace_modify slab <slab_factor>
   kspace_modify wire <wire_factor>
   kspace_modify slab ew2d

.. warning::

   Atom positions of electrode particles have to be fixed at all times.

----------
