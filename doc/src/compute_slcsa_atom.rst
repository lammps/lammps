.. index:: compute slcsa/atom

compute slcsa/atom command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID slcsa/atom twojmax nclasses db_mean_descriptor_file lda_file lr_decision_file lr_bias_file maha_file value

* ID, group-ID are documented in :doc:`compute <compute>` command
* slcsa/atom = style name of this compute command
* twojmax = band limit for bispectrum components (non-negative integer)
* nclasses = number of crystal structures used in the database for the classifier SL-CSA
* db_mean_descriptor_file = file name of file containing the database mean descriptor
* lda_file = file name of file containing the linear discriminant analysis matrix for dimension reduction
* lr_decision_file = file name of file containing the scaling matrix for logistic regression classification
* lr_bias_file = file name of file containing the bias vector for logistic regression classification
* maha_file = file name of file containing for each crystal structure: the Mahalanobis distance threshold for sanity check purposes, the average reduced descriptor and the inverse of the corresponding covariance matrix
* c_ID[*] = compute ID of previously required *compute sna/atom* command

Examples
""""""""

.. code-block:: LAMMPS

   compute b1 all sna/atom 9.0 0.99363 8 0.5 1.0 rmin0 0.0 nnn 24 wmode 1 delta 0.3
   compute b2 all slcsa/atom 8 4 mean_descriptors.dat lda_scalings.dat lr_decision.dat lr_bias.dat maha_thresholds.dat c_b1[*]

Description
"""""""""""

.. versionadded:: TBD

Define a computation that performs the Supervised Learning Crystal
Structure Analysis (SL-CSA) from :ref:`(Lafourcade) <Lafourcade2023_1>`
for each atom in the group. The SL-CSA tool takes as an input a per-atom
descriptor (bispectrum) that is computed through the *compute sna/atom*
command and then proceeds to a dimension reduction step followed by a
logistic regression in order to assign a probable crystal structure to
each atom in the group. The SL-CSA tool is pre-trained on a database
containing :math:`C` distinct crystal structures from which a crystal
structure classifier is derived and a tutorial to build such a tool is
available at `SL-CSA <https://github.com/lafourcadep/SL-CSA>`_.

The first step of the SL-CSA tool consists in performing a dimension
reduction of the per-atom descriptor :math:`\mathbf{B}^i \in
\mathbb{R}^{D}` through the Linear Discriminant Analysis (LDA) method,
leading to a new projected descriptor
:math:`\mathbf{x}^i=\mathrm{P}_\mathrm{LDA}(\mathbf{B}^i):\mathbb{R}^D
\rightarrow \mathbb{R}^{d=C-1}`:

.. math::

   \mathbf{x}^i = \mathbf{C}^T_\mathrm{LDA} \cdot (\mathbf{B}^i - \mu^\mathbf{B}_\mathrm{db})

where :math:`\mathbf{C}^T_\mathrm{LDA} \in \mathbb{R}^{D \times d}` is
the reduction coefficients matrix of the LDA model read in file
*lda_file*, :math:`\mathbf{B}^i \in \mathbb{R}^{D}` is the bispectrum of
atom :math:`i` and :math:`\mu^\mathbf{B}_\mathrm{db} \in \mathbb{R}^{D}`
is the average descriptor of the entire database. The latter is computed
from the average descriptors of each crystal structure read from the
file *mean_descriptors_file*.

The new projected descriptor with dimension :math:`d=C-1` allows for a
good separation of different crystal structures fingerprints in the
latent space.

Once the dimension reduction step is performed by means of LDA, the new
descriptor :math:`\mathbf{x}^i \in \mathbb{R}^{d=C-1}` is taken as an
input for performing a multinomial logistic regression (LR) which
provides a score vector
:math:`\mathbf{s}^i=\mathrm{P}_\mathrm{LR}(\mathbf{x}^i):\mathbb{R}^d
\rightarrow \mathbb{R}^C` defined as:

.. math::

   \mathbf{s}^i = \mathbf{b}_\mathrm{LR} + \mathbf{D}_\mathrm{LR} \cdot {\mathbf{x}^i}^T

with :math:`\mathbf{b}_\mathrm{LR} \in \mathbb{R}^C` and
:math:`\mathbf{D}_\mathrm{LR} \in \mathbb{R}^{C \times d}` the bias
vector and decision matrix of the LR model after training both read in
files *lr_fil1* and *lr_file2* respectively.

Finally, a probability vector
:math:`\mathbf{p}^i=\mathrm{P}_\mathrm{LR}(\mathbf{x}^i):\mathbb{R}^d
\rightarrow \mathbb{R}^C` is defined as:

.. math::

   \mathbf{p}^i = \frac{\mathrm{exp}(\mathbf{s}^i)}{\sum\limits_{j} \mathrm{exp}(s^i_j) }

from which the crystal structure assigned to each atom with descriptor
:math:`\mathbf{B}^i` and projected descriptor :math:`\mathbf{x}^i` is
computed as the *argmax* of the probability vector
:math:`\mathbf{p}^i`. Since the logistic regression step systematically
attributes a crystal structure to each atom, a sanity check is needed to
avoid misclassification. To this end, a per-atom Mahalanobis distance to
each crystal structure *CS* present in the database is computed:

.. math::

   d_\mathrm{Mahalanobis}^{i \rightarrow \mathrm{CS}} = \sqrt{(\mathbf{x}^i - \mathbf{\mu}^\mathbf{x}_\mathrm{CS})^\mathrm{T} \cdot \mathbf{\Sigma}^{-1}_\mathrm{CS} \cdot (\mathbf{x}^i - \mathbf{\mu}^\mathbf{x}_\mathrm{CS}) }

where :math:`\mathbf{\mu}^\mathbf{x}_\mathrm{CS} \in \mathbb{R}^{d}` is
the average projected descriptor of crystal structure *CS* in the
database and where :math:`\mathbf{\Sigma}_\mathrm{CS} \in \mathbb{R}^{d
\times d}` is the corresponding covariance matrix. Finally, if the
Mahalanobis distance to crystal structure *CS* for atom *i* is greater
than the pre-determined threshold, no crystal structure is assigned to
atom *i*. The Mahalanobis distance thresholds are read in file
*maha_file* while the covariance matrices are read in file
*covmat_file*.

The `SL-CSA <https://github.com/lafourcadep/SL-CSA>`_ framework provides
an automatic computation of the different matrices and thresholds
required for a proper classification and writes down all the required
files for calling the *compute slcsa/atom* command.

The *compute slcsa/atom* command requires that the :doc:`compute
sna/atom <compute_sna_atom>` command is called before as it takes the
resulting per-atom bispectrum as an input. In addition, it is crucial
that the value *twojmax* is set to the same value of the value *twojmax*
used in the *compute sna/atom* command, as well as that the value
*nclasses* is set to the number of crystal structures used in the
database to train the SL-CSA tool.

Output info
"""""""""""

By default, this compute computes the Mahalanobis distances to the
different crystal structures present in the database in addition to
assigning a crystal structure for each atom as a per-atom vector, which
can be accessed by any command that uses per-atom values from a compute
as input.  See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute sna/atom <compute_sna_atom>`

Default
"""""""

none

----------

.. _Lafourcade2023_1:

**(Lafourcade)** Lafourcade, Maillet, Denoual, Duval, Allera, Goryaeva, and Marinica,
`Comp. Mat. Science, 230, 112534 (2023) <https://doi.org/10.1016/j.commatsci.2023.112534>`_
