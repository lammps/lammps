.. index:: pair_style ylz
.. index:: pair_style ylz/gpu
.. index:: pair_style ylz/intel
.. index:: pair_style ylz/omp
 
pair_style ylc command
===========================

Accelerator Variants: *ylz/gpu*, *ylz/intel*, *ylz/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style ylz cutoff


* cutoff = global cutoff for interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style   ylz  2.6
   pair_coeff   *  *  1.0  1.0  4  3  0.0  2.6  
   

Description
"""""""""""

The *ylz* (Yuan-Li-Zhang) 
:ref:`(Yuan) <Yuan>` style computes anisotropic interactions between pairs of particles considering the relative particle orientations via the formulas 

.. math::

   U ( \mathbf{r}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) =\left\{\begin{matrix} \mathbf{u}_R(r)+\left [ 1-\phi (\mathbf{r\hat{}}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) \right ]\epsilon, ~~ r<\mathbf{r}_{min} \\ \mathbf{u}_A(r)\phi (\mathbf{r\hat{}}_{ij}, \mathbf{n}_i, \mathbf{n}_j ),~~  \mathbf{r}_{min}<r<\mathbf{r}_{c} \\ \end{matrix}\right.\\\\ \phi (\mathbf{r\hat{}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )=1+\mu (a(\mathbf{r\hat{}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )-1) \\\\ a(\mathbf{r\hat{}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )=(\mathbf{n}_i\times\mathbf{r\hat{}}_{ij} )\cdot (\mathbf{n}_j\times\mathbf{r\hat{}}_{ij} )+sin\mathbf{\theta}_0(\mathbf{n}_i-\mathbf{n}_j)\cdot \mathbf{r\hat{}}_{ij}\\\\ \mathbf{u}_R(r)=\epsilon \left [ \left ( \frac{{r}_{min}}{r} \right )^{4}-2\left ( \frac{{r}_{min}}{r}\right )^{2} \right ] \\\\ \mathbf{u}_A(r)=-\epsilon\;cos^{2\zeta }\left ( \frac{\pi}{2}\frac{\left ( {r}-{r}_{min} \right )}{\left ( {r}_{c}-{r}_{min} \right )} \right ) \\ 

where :math:`r_{i}` and :math:`r_{j}` are the center position vectors of particles i and j, respectively, :math:`r_{ij}=r_{i}-r_{j}` is the inter-particle distance vector, :math:`r=\left|r_{ij} \right|` and :math:`{r\hat{}}_{ij}=\mathbf{r}_{ij}/r`. The unit vectors :math:`n_{i}` and :math:`n_{j}` represent the axes of symmetry of particles i and j, respectively. :math:`u_R` and :math:`u_A` are the repulsive and attractive potentials. :math:`\phi` is an angular function which depends on the relative orientation between pair particles. :math:`\mu` is the parameter related to bending rigidity, :math:`\theta_{0}` is the parameter related to the spontaneous curvature, and :math:`\epsilon` is the energy unit, respectively. The :math:`\zeta` controls the slope of the attractive branch and :math:`{r}_{c}`is the cutoff radius. :math:`r_{min}` is the distance which minimizes the potential energy :math:`u_{A}(r)`and :math:`r_{min}=2^{1/6}\sigma`, where :math:`\sigma` is the length unit.

Use of this pair style requires the NVE, NVT, or NPT fixes with the *asphere* extension (e.g. :doc:`fix nve/asphere <fix_nve_asphere>`) in order to integrate particle rotation.  Additionally, :doc:`atom_style ellipsoid <atom_style>` should be used since it defines the rotational state of each particle.

The following coefficients must be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` = well depth (energy units)
* :math:`\sigma` = minimum effective particle radii (distance units)
* :math:`\zeta` = tune parameter for the slope of the attractive branch
* :math:`\mu` = parameter related to bending rigidity
* :math:`\sin(\theta _{0})` is the sine of parameter related to the spontaneous curvature
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global
cutoff specified in the pair_style command is used.






----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients and cutoff distance for this pair style can be mixed.  The default mix value is *geometric*\ .  See the "pair_modify" command for details.

The :doc:`pair_modify <pair_modify>` table option is not relevant for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>` tail option for adding long-range tail corrections to energy and pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *ylz* style is part of the ASPHERE package.  It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

These pair styles require that atoms store torque and a quaternion to represent their orientation, as defined by the :doc:`atom_style <atom_style>`.  It also require they store a per-type :doc:`shape <set>`.  The particles cannot store a per-particle diameter.

This pair style requires that atoms be ellipsoids as defined by the
:doc:`atom_style ellipsoid <atom_style>` command.

Particles acted on by the potential can be finite-size aspherical or
spherical particles, or point particles.  Spherical particles have all
3 of their shape parameters equal to each other.  Point particles have
all 3 of their shape parameters equal to 0.0.


Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix nve/asphere <fix_nve_asphere>`,
:doc:`compute temp/asphere <compute_temp_asphere>`, :doc:`pair_style resquared <pair_resquared>`, :doc:`pair_style gayberne <pair_resquared>`

Default
"""""""

none

----------

.. _Yuan:

**(Yuan)** Yuan, Huang, Li, Lykotrafitis, Zhang, Phys. Rev. E, 82, 011905(2010).
