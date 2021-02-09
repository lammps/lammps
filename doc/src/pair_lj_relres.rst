.. index:: pair_style lj/relres

pair_style lj/relres command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/relres Rsi Rso Rci Rco

* Rsi = inner switching distance - boundary up to which LJ potential of fine-grained model is applied (distance units)
* Rso = outer switching distance - boundary beyond which LJ potential of coarse-grained model is applied (distance units)  
* Rci = inner cutoff beyond which force smoothing is applied (distance units)
* Rco = outer cutoff for lj/relres interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/relres 4.0 5.0 8.0 10.0
   pair_coeff * * 0.5 1.0 0.0 1.0
   pair_coeff 1 1 1.5 1.2 3.0 1.2 3.0 3.5 6.0 7.0

Description
"""""""""""

Style *lj/relres* computes a LJ interaction with RelRes methodology developed by :ref:`Chaimovich at al.<Chaimovich1>`
This methodology applies fine-grained model between near neighbors (up to :math:`r_{si}` boundary) and a simplified coarse-grained model 
for far neighbors (beyond :math:`r_{so}` boundary) allowing significant improvement in computational efficiency while preserving correctness 
of simulation results.

.. math::

   E = \left\{\begin{array}{lr}
        4 \epsilon^{FG} \left[ \left(\frac{\sigma^{FG}}{r}\right)^{12} - \left(\frac{\sigma^{FG}}{r}\right)^6 \right]-G_{si}, & r< r_{si} \\
        \sum_{m=0}^{4} C_{sm}\left(r-r_{si}\right)^m-G_{so} ,   &   r_{si}\leq r< r_{so} \\
        4 \epsilon^{CG} \left[ \left(\frac{\sigma^{CG}}{r}\right)^{12} -     \left(\frac{\sigma^{CG}}{r}\right)^6 \right]-G_c, &    r_{so}\leq r<r_{ci} \\
        \sum_{m=0}^{4} C_{cm}\left(r-r_{ci}\right)^m -G_c, &  r_{ci}\leq r< r_{co} \\
        0, &  r\geq r_{co}\end{array}\right.

Between :math:`r_{si}` and :math:`r_{so}` the polynomial smoothing is applied in a way that the force and its 1st derivative are not discontinued 
at switching between fine- and coarse-grained potentials (between :math:`r_{si}` and :math:`r_{so}`) and at cutoff (between :math:`r_{ci}` and :math:`r_{co}`). 
The corresponding polynomial coefficients :math:`C_{sm}` and :math:`C_{cm}` and shifting constants :math:`G_{si}`, :math:`G_{so}` and :math:`G_{c}` are computed by LAMMPS accordingly. 
To avoid smoothing, the inner switching distance :math:`r_{si}` parameter should be set equal to the outer switching distance :math:`r_{so}` parameter 
(:math:`r_{si}=r_{so}`). Similarly, to avoid smoothing at cutoff, inner and outer cutoff parameters should be set equal (:math:`r_{ci}=r_{co}`).
Details can be found in :ref:`(Chaimovich) <Chaimovich2>`.

.. note::

   Energy and force resulting from this methodology can be plotted via the
   :doc:`pair_write <pair_write>` command to see the effect.

In implementation of *lj/relres* style, atoms are grouped in the way that one of the atoms in the group plays the role of a coarse-grained site for the calculation 
of interactions beyond :math:`r_{so}` distance while continuing to play the role of a fine-grained site for shorter distances. 
This atom must be defined as a different atom type. Other atoms in the group participate in the fine-grained interactions only.

The following coefficients must be defined for each pair of atom
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as will be described below:

* :math:`\epsilon^{FG}` (energy units)
* :math:`\sigma^{FG}` (distance units)
* :math:`\epsilon^{CG}` (energy units)
* :math:`\sigma^{CG}` (distance units)

For atom types that are used as fine-grained sites only, :math:`\epsilon^{CG}` must be set to 0 (zero). 
For atom types that are used as coarse-grained sites only (if any), :math:`\epsilon^{FG}` must be set to 0 (zero).

Additional parameters can be defined to specify different :math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`, :math:`r_{co}` for a particular set of atom types:

* :math:`r_{si}` (distance units)
* :math:`r_{so}` (distance units)
* :math:`r_{ci}` (distance units)
* :math:`r_{co}` (distance units)

These parameters are optional and they are used to override global values defined in the pair_style command. 
If this override option is used, all four values must be specified.  If not specified, the global values for :math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`, and :math:`r_{co}` are used.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J with I != J, the :math:`\epsilon^{FG}`, :math:`\sigma^{FG}`, :math:`\epsilon^{CG}`, :math:`\sigma^{CG}`, :math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`, and :math:`r_{co}`
parameters for this pair style can be mixed, if not defined explicitly.
All parameters are mixed according to the pair_modify mix option.  The
default mix value is *geometric*\ , and it is recommended to use with this *lj/relres* style.  
See the "pair_modify" command for details. 

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction. It is recommended to set this option to *yes*\ . 
Otherwise, the shifting constant :math:`G_{c}` is set to zero. Constants :math:`G_{si}` and :math:`G_{so}` are not impacted by this option.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since the energy of the pair interaction is smoothed to 0.0
at the cutoff.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Chaimovich1:

**(Chaimovich at al.)** A.Chaimovich, C. Peter and K. Kremer, J. Chem. Phys. 143, 243107
(2015).

.. _Chaimovich2:

**(Chaimovich)** M.Chaimovich and A. Chaimovich, J. Chem. Theory Comput. 17, 1045-1059
(2021).

