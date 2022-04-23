.. index:: pair_style rann

pair_style rann command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rann
   pair_coeff file Type1_element Type2_element Type3_element...

Examples
""""""""

.. code-block:: LAMMPS

    pair_style rann
    pair_coeff ** Mg.rann Mg
    pair_coeff ** MgAlalloy.rann Mg Mg Al Mg

Description
"""""""""""

Pair style *rann* computes pairwise interactions for a variety of
materials using rapid atomistic neural network (RANN) potentials
(:ref:`Dickel <Dickel>` , :ref:`Nitol <Nitol>`).  Neural network
potentials work by first generating a series of symmetry functions
i.e. structural fingerprints from the neighbor list and then using these
values as the input layer of a neural network.  There is a single output
neuron in the final layer which is the energy.  Atomic forces are found
by analytical derivatives computed via back-propagation.  For alloy
systems, each element has a unique network.

Potential file syntax
"""""""""""""""""""""

The RANN potential is defined by a single text file which contains all
the fitting parameters for the alloy system.  The potential file also
defines the active fingerprints, network architecture, activation
functions, etc.  The potential file is divided into several sections
which are identified by one of the following keywords:

* atomtypes
* mass
* fingerprintsperelement
* fingerprints
* fingerprintconstants
* screening (optional)
* networklayers
* layersize
* weight
* bias
* activationfunctions
* calibrationparameters (ignored)

The '#' character is treated as a comment marker, similar to LAMMPS
input scripts.  Sections are not required to follow a rigid ordering,
but do require previous definition of prerequisite information.  E.g.,
fingerprintconstants for a particular fingerprint must follow the
fingerprints definition; layersize for a particular layer must follow
the declaration of network layers.

*atomtypes* are defined as follows using element keywords separated by spaces.

.. code-block::

    atomtypes:
    Fe Mg Al etc.

*mass* must be specified for each element keyword as follows:

.. code-block::

    mass:Mg:
    24.305
    mass:Fe:
    55.847
    mass:Al:
    26.982

*fingerprintsperelement* specifies how many fingerprints are active for
computing the energy of a given atom.  This number must be specified for
each element keyword.  Active elements for each fingerprint depend upon
the type of the central atom and the neighboring atoms.  Pairwise
fingerprints may be defined for a Mg atom based exclusively on its Al
neighbors, for example.  Bond fingerprints may use two neighbor lists of
different element types.  In computing fingerprintsperelement from all
defined fingerprints, only the fingerprints defined for atoms of a
particular element should be considered, regardless of the elements used
in its neighbor list.  In the following code, for example, some
fingerprints may compute pairwise fingerprints summing contributions
about Fe atoms based on a neighbor list of exclusively Al atoms, but if
there are no fingerprints summing contributions of all neighbors about a
central Al atom, then fingerprintsperelement of Al is zero:

.. code-block::

    fingerprintsperelement:Mg:
    5
    fingerprintsperelement:Fe:
    2
    fingerprintsperelement:Al:
    0

*fingerprints* specifies the active fingerprints for a certain element
combination.  Pair fingerprints are specified for two elements, while
bond fingerprints are specified for three elements.  Only one
fingerprints header should be used for an individual combination of
elements.  The ordering of the fingerprints in the network input layer
is determined by the order of element combinations specified by
subsequent *fingerprints* lines, and the order of the fingerprints
defined for each element combination.  Multiple fingerprints of the same
style or different ones may be specified.  If the same style and element
combination is used for multiple fingerprints, they should have
different id numbers.  The first element specifies the atoms for which
this fingerprint is computed while the other(s) specify which atoms to
use in the neighbor lists for the computation.  Switching the second and
third element type in bond fingerprints has no effect on the
computation:

.. code-block::

    fingerprints:Mg_Mg:
    radial_0 radialscreened_0 radial_1
    fingerprints:Mg_Al_Fe:
    bond_0 bondspin_0
    fingerprints:Mg_Al:
    radial_0 radialscreened_0

The following fingerprint styles are currently defined.  See the
:ref:`formulation section <fingerprints>` below for their definitions:

* radial
* radialscreened
* radialspin
* radialscreenedspin
* bond
* bondscreened
* bondspin
* bondscreenedspin

*fingerprintconstants*  specifies the meta-parameters for a defined fingerprint.  For all radial styles, re, rc,
alpha, dr, o, and n must be specified.  re should usually be the stable interatomic distance, rc is the cutoff
radius, dr is the cutoff smoothing distance, o is the lowest radial power term (which may be negative), and n
is the highest power term.  The total length of the fingerprint vector is (n-o+1).  alpha is a list of decay parameters
used for exponential decay of radial contributions.  It may be set proportionally to the bulk modulus similarly
to MEAM potentials, but other values may provided better fitting in special cases.  Bond style fingerprints require
specification of re, rc, alphak, dr, k, and m.  Here m is the power of the bond cosines and k is the number of
decay parameters.  Cosine powers go from 0 to m-1 and are each computed for all values of alphak.  Thus the total
length of the fingerprint vector is m*k.

.. code-block::

   fingerprintconstants:Mg_Mg:radialscreened_0:re:
   3.193592
   fingerprintconstants:Mg_Mg:radialscreened_0:rc:
   6.000000
   fingerprintconstants:Mg_Mg:radialscreened_0:alpha:
   5.520000 5.520000 5.520000 5.520000 5.520000
   fingerprintconstants:Mg_Mg:radialscreened_0:dr:
   2.806408
   fingerprintconstants:Mg_Mg:radialscreened_0:o:
   -1
   fingerprintconstants:Mg_Mg:radialscreened_0:n:
   3

*screening* specifies the Cmax and Cmin values used in the screening
fingerprints. Contributions form neighbors to the fingerprint are
omitted if they are blocked by a closer neighbor, and reduced if they
are partially blocked.  Larger values of Cmin correspond to neighbors
being blocked more easily.  Cmax cannot be greater than 3, and Cmin
cannot be greater than Cmax or less than zero.  Screening may be omitted
in which case the default values Cmax = 2.8, Cmin = 0.8 are used.  Since
screening is a bond computation, it is specified separately for each
combination of three elements in which the latter two may be
interchanged with no effect.

.. code-block::

    screening:Mg_Mg_Mg:Cmax:
    2.700000
    screening:Mg_Mg_Mg:Cmin:
    0.400000

*networklayers* species the size of the neural network for each atom.
It counts both the input and output layer and so is 2 + \<hidden layers\>.

.. code-block::

   networklayers:Mg:
   3

*layersize* specifies the length of each layer, including the input
layer and output layer.  The input layer is layer 0.  The size of the
input layer size must match the summed length of all the fingerprints
for that element, and the output layer size must be 1:

.. code-block::

    layersize:Mg:0:
    14
    layersize:Mg:1:
    20
    layersize:Mg:2:
    1


*weight* specifies the weight for a given element and layer.  Weight
cannot be specified for the output layer.  The weight of layer i is a
*m* x *n* matrix where *m* is the layer size of *i* and *n* is the layer size of
*i*\ +1:

.. code-block::

   weight:Mg:0:
   w11 w12 w13 ...
   w21 w22 w23 ...
   ...

*bias* specifies the bias for a given element and layer.  Bias cannot be
specified for the output layer.  The bias of layer i is a nx1 vector
where n is the layer size of i+1:

.. code-block::

   bias:Mg:0:
   b1
   b2
   b3
   ...

*activationfunctions* specifies the activation function for a given
element and layer.  Activation functions cannot be specified for the
output layer:

.. code-block::

    activationfunctions:Mg:0:
    sigI
    activationfunctions:Mg:1:
    linear

The following activation styles are currently specified.  See the
:ref:`formulation section <activations>` below for their definitions.

* sigI

* linear

*calibrationparameters* specifies a number of parameters used to calibrate the potential. These are ignored
by LAMMPS.

Formulation
"""""""""""

In the RANN formulation, the total energy of a system of atoms
is given by:

.. math::

    E = \sum_{\alpha} E^{\alpha}\\\\
    E^{\alpha} = {}^{N}\!A^{\alpha}\\\\
    {}^{n+1}\!A_i^{\alpha} = {}^{n}\!F\left({}^{n}\!W_{ij}{\;}^{n}\!A_j^{\alpha}+{}^{n}\!B_i\right)\\\\
    {}^{0}\!A_i^{\alpha} = \left[\begin{array}{c} {}^1\!S\!f^\alpha\\ {}^2\!S\!f^\alpha \\...\\\end{array}\right]

Here :math:`E^\alpha` is the energy of atom :math:`\alpha`,
:math:`{}^n\!F()`, :math:`{}^n\!W_{ij}` and :math:`{}^n\!B_i` are the
activation function, weight matrix and bias vector of the n-th layer
respectively.  The inputs to the first layer are a collection of
structural fingerprints which are collected and reshaped into a single
long vector.  The individual fingerprints may be defined in any order
and have various shapes and sizes.  Multiple fingerprints of the same
type and varying parameters may also be defined in the input layer.

Eight types of structural fingerprints are currently defined.  In the
following, :math:`\beta` and :math:`\gamma` span the full neighbor list
of atom :math:`\alpha`.  :math:`\delta_i` are decay meta-parameters, and
:math:`r_e` is a meta-parameter roughly proportional to the first
neighbor distance.  :math:`r_c` and :math:`dr` are the neighbor cutoff
distance and cutoff smoothing distance respectively.
:math:`S^{\alpha\beta}` is the MEAM screening function :ref:`(Baskes)
<Baskes97>`, :math:`s_i^\alpha` and :math:`s_i^\beta` are the atom spin
vectors :ref:`(Tranchida) <Tranchida7>`.  :math:`r^{\alpha\beta}` is the
distance from atom :math:`\alpha` to atom :math:`\beta`, and
:math:`\theta^{\alpha\beta\gamma}` is the bond angle:

.. math ::

    cos\left(\theta^{\alpha\beta\gamma}\right)=\frac{\mathbf{r}^{\alpha\beta} \cdot \mathbf{r}^{\alpha\gamma}}{r^{\alpha\beta}r^{\alpha\gamma}}

:math:`S^{\alpha\beta}` is defined as :ref:`(Baskes) <Baskes97>`:

.. math::

    X^{\gamma\beta} = \left(\frac{r^{\gamma\beta}}{r^{\alpha\beta}}\right)^2\\
    \\
    X^{\alpha\gamma} = \left(\frac{r^{\alpha\gamma}}{r^{\alpha\beta}}\right)^2\\
    \\
    C = \frac{2\left(X^{\alpha\gamma}+X^{\gamma\beta}\right)-\left(X^{\alpha\gamma}-X^{\gamma\beta}\right)^2-1}{1-\left(X^{\alpha\gamma}-X^{\gamma\beta}\right)^2}\\
    \\
    f_c(x) = \left[\begin{array}{l}  1 \; \: x \geq 1\\ \left(1-\left(1-x\right)^4\right)^2 \; \: 0<x<1\\0\; \; x\leq0\end{array}\right.\\
    \\
    S^{\alpha\beta\gamma} = f_c\left(\frac{C-C_{min}}{C_{max}-C_{min}}\right)\\
    \\
    S^{\alpha\beta} = \prod_\gamma S^{\alpha\beta\gamma}\\


The structural fingerprints are computed as follows:

.. _fingerprints:

* **radial**

.. math::

    {}^r\!S\!f_i^\alpha = \sum_{\beta} \left(\frac{r^{\alpha\beta}}{r_e}\right)^ie^{-\delta_i \frac{r^{\alpha\beta}}{r_e}}f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)

* **bond**

.. math::

    {}^b\!S\!f_{ij}^\alpha = \sum_{\beta}\sum_{\gamma} \left(cos(\theta_{\alpha\beta\gamma})\right)^ie^{-\delta_j \frac{r^{\alpha\beta}}{r_e}}e^{-\delta_j \frac{r^{\alpha\gamma}}{r_e}}f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)f_c\left(\frac{r_c-r^{\alpha\gamma}}{dr}\right)

* **radialscreened**

.. math::

    {}^{rsc}\!S\!f_i^\alpha = \sum_{\beta} \left(\frac{r^{\alpha\beta}}{r_e}\right)^ie^{-\delta_i \frac{r^{\alpha\beta}}{r_e}}S^{\alpha\beta}f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)

* **bondscreened**

.. math::

    {}^{bsc}\!S\!f_{ij}^\alpha = \sum_{\beta}\sum_{\gamma} \left(cos(\theta_{\alpha\beta\gamma})\right)^ie^{-\delta_j \frac{r^{\alpha\beta}}{r_e}}e^{-\delta_j \frac{r^{\alpha\gamma}}{r_e}}S^{\alpha\beta}S^{\alpha\gamma}f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)f_c\left(\frac{r_c-r^{\alpha\gamma}}{dr}\right)

* **radialspin**

.. math::

    {}^{rsp}\!S\!f_i^\alpha = \sum_{\beta} \left(\frac{r^{\alpha\beta}}{r_e}\right)^ie^{-\delta_i \frac{r^{\alpha\beta}}{r_e}}\left(\mathbf{s^\alpha \cdot s^\beta}\right)f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)

*  **bondspin**

.. math::

   {}^{bsp}\!S\!f_{ij}^\alpha = \sum_{\beta}\sum_{\gamma} \left(cos(\theta_{\alpha\beta\gamma})\right)^ie^{-\delta_j \frac{r^{\alpha\beta}}{r_e}}e^{-\delta_j \frac{r^{\alpha\gamma}}{r_e}}\left(\mathbf{s^\alpha \cdot s^\beta}\right)\left(\mathbf{s^\alpha \cdot s^\gamma}\right)f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)f_c\left(\frac{r_c-r^{\alpha\gamma}}{dr}\right)

* **radialscreenedspin**

.. math::

    {}^{rscsp}\!S\!f_i^\alpha = \sum_{\beta} \left(\frac{r^{\alpha\beta}}{r_e}\right)^ie^{-\delta_i \frac{r^{\alpha\beta}}{r_e}}S^{\alpha\beta}\left(\mathbf{s^\alpha \cdot s^\beta}\right)f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)

* **bondscreenedspin**

.. math::

    {}^{bscsp}\!S\!f_{ij}^\alpha = \sum_{\beta}\sum_{\gamma} \left(cos(\theta_{\alpha\beta\gamma})\right)^ie^{-\delta_j \frac{r^{\alpha\beta}}{r_e}}e^{-\delta_j \frac{r^{\alpha\gamma}}{r_e}}S^{\alpha\beta}S^{\alpha\gamma}\left(\mathbf{s^\alpha \cdot s^\beta}\right)\left(\mathbf{s^\alpha \cdot s^\gamma}\right)f_c\left(\frac{r_c-r^{\alpha\beta}}{dr}\right)f_c\left(\frac{r_c-r^{\alpha\gamma}}{dr}\right)

The activation functions are computed as follows:

.. _activations:

* **sigI**

.. math::

    F^{sigI}(x) = 0.1x+0.9ln\left(e^x+1\right)

* **linear**

.. math::

   F^{linear}(x) = x

Restrictions
""""""""""""

Pair style *rann* is part of the ML-RANN package.  It is only enabled if LAMMPS was built with that
package.  Additionally, if any spin fingerprint styles are used LAMMPS must be built with the SPIN
package as well.

Pair style *rann* does not support computing per-atom stress or using :doc:`pair_modify nofdotr <pair_modify>`.

Defaults
""""""""""""

Cmin = 0.8, Cmax = 2.8.

----------


.. _Baskes97:

**(Baskes)** Baskes,
Materials Chemistry and Physics, 50(2), 152-158, (1997).

.. _Dickel:

**(Dickel)** Dickel, Francis, and Barrett,
Computational Materials Science 171 (2020): 109157.

.. _Nitol:

**(Nitol)** Nitol, Dickel, and Barrett,
Computational Materials Science 188 (2021): 110207.

.. _Tranchida7:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).

