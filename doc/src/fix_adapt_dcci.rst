.. index:: fix adapt/dcci

fix adapt/dcci command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID adapt/dcci lambda attribute args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* adapt/dcci = style name of this fix command
* lambda = intercept parameter of linear effective temperature function
* one or more attribute/arg pairs may be appended
* attribute = *pair*

 .. parsed-literal::

       *pair* args = pstyle pparam I J v_name
         pstyle = pair style name (e.g., lj/cut)
         fscale = scaling force parameter to adapt over time
         I,J = type pair(s) to set parameter for        

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all adapt/dcci ${lambda} pair lj/cut fscale * * 
   fix 1 all adapt/dcci ${lambda} pair meam fscale * * 

Description
"""""""""""

This fix is modified version of the :doc:`fix adapt <fix_adapt>` command, where the forces
on the atoms are dynamically scaled during the simulation for each step by the *lambda* parameter, which is controlled by the :doc:`dcci <dcci>` command in order to calculate the coxistence conditions while run a non equilibrium dynamics :ref:`(deKoning) <deKoning>`. Further explanation can be found in our recent paper :ref:`(Cajahuaringa) <Cajahuaringa>`.

----------

The *pair* keyword enables various parameters of potentials defined by
the :doc:`pair_style <pair_style>` command to be changed, if the pair
style supports it.  Note that the :doc:`pair_style <pair_style>` and
:doc:`pair_coeff <pair_coeff>` commands must be used in the usual manner
to specify these parameters initially; the fix adapt command simply
overrides the parameters.

The *pstyle* argument is the name of the pair style.  If
:doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` is used,
*pstyle* should be a sub-style name. If there are multiple
sub-styles using the same pair style, then *pstyle* should be specified
as "style:N", where *N* is which instance of the pair style you wish to
adapt (e.g., the first or second).  For example, *pstyle* could be
specified as "soft" or "lubricate" or "lj/cut:1" or "lj/cut:2".  The
*fscale* argument is the name of the parameter to indicate the scaling force. See the doc pages for individual pair styles and their energy formulas for the meaning of these parameters:

+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut <pair_lj>`                                                      | fscale                                           | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`meam <pair_meam>`                                                      | fscale                                           | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`deepmd <pair_meam>`                                                    | fscale                                           | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+

.. note::

   It is easy to add new pairwise potentials and their parameters
   to this list.  All it typically takes is adding an extract() method to
   the pair\_\*.cpp file associated with the potential.


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix works in combination with the :doc:`dcci <dcci>` command.

No information about this fix is written to :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.  

This fix computes a global scalar quantity which can be accessed by various :doc:`output commands <Howto_output>`. The scalar is an \lambda parameter which is controlled by :doc:`dcci <dcci>` 

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the REPLICA package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`dcci <dcci>`

Default
"""""""

none

----------

.. _Cajahuaringa:

**(Cajahuaringa)** Cajahuaringa and Antonelli, Comput. Mater. Sci., 111275, 207 (2022).

.. _deKoning:

**(deKoning)** de Koning, Antonelli and Sidney, J Chem Phys 115, 11025 (2001).
