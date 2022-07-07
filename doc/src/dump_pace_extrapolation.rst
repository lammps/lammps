.. index:: dump pace/extrapolation

dump pace/extrapolation  command
================================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID pace/extrapolation N file args


* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* N = dump every this many timesteps (use 1 to potentially dump on every step, see below)
* file = name of file to write dump info to
* args = list of arguments (see :doc:`dump custom <dump>`)

Examples
""""""""

.. code-block:: LAMMPS

    dump pace all pace/extrapolation 1 extrapolation.dat id type mass x y z c_pace_gamma
    dump pace all pace/extrapolation 1 extrapolation.lammpsbin id type mass x y z c_pace_gamma

Description
"""""""""""

Dump a snapshot of atom coordinates if extrapolation grade, computed by
:doc:`compute pace/extrapolation <compute_pace_extrapolation>`, exceeds *gamma_lower_bound* threshold,
provided in :doc:`pair_style pace/extrapolation <pair_pace>`.

.. note::

   To be able to use this dump, you need to setup  :doc:`pair_style pace/extrapolation <pair_pace>`
   and  :doc:`compute pace/extrapolation <compute_pace_extrapolation>` beforehand

----------

Related commands
""""""""""""""""

:doc:`pair_style pace/extrapolation <pair_pace>`, :doc:`compute pace/extrapolation <compute_pace_extrapolation>`
