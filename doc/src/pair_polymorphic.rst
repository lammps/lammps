.. index:: pair_style polymorphic

pair_style polymorphic command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style polymorphic

style = *polymorphic*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style polymorphic
   pair_coeff * * TlBr_msw.polymorphic Tl Br
   pair_coeff * * AlCu_eam.polymorphic Al Cu
   pair_coeff * * GaN_tersoff.polymorphic Ga N
   pair_coeff * * GaN_sw.polymorphic GaN

Description
"""""""""""

The *polymorphic* pair style computes a 3-body free-form potential
(:ref:`Zhou <Zhou3>`) for the energy E of a system of atoms as

.. math::

   E & = \frac{1}{2}\sum_{i=1}^{i=N}\sum_{j=1}^{j=N}\left[\left(1-\delta_{ij}\right)\cdot U_{IJ}\left(r_{ij}\right)-\left(1-\eta_{ij}\right)\cdot F_{IJ}\left(r_{ij}\right)\cdot V_{IJ}\left(r_{ij}\right)\right] \\
   X_{ij} & = \sum_{k=i_1,k\neq i,j}^{i_N}W_{IK}\left(r_{ik}\right)\cdot G_{JIK}\left(\theta_{jik}\right)\cdot P_{IK}\left(\Delta r_{jik}\right) \\
   \Delta r_{jik} & = r_{ij}-\xi_{IJ}\cdot r_{ik}

where I, J, K represent species of atoms i, j, and k, :math:`i_1, ...,
i_N` represents a list of *i*\ 's neighbors, :math:`\delta_{ij}` is a
Dirac constant (i.e., :math:`\delta_{ij} = 1` when :math:`i = j`, and
:math:`\delta_{ij} = 0` otherwise), :math:`\eta_{ij}` is similar
constant that can be set either to :math:`\eta_{ij} = \delta_{ij}` or
:math:`\eta_{ij} = 1 - \delta_{ij}` depending on the potential type,
:math:`U_{IJ}(r_{ij})`, :math:`V_{IJ}(r_{ij})`, :math:`W_{IK}(r_{ik})`
are pair functions, :math:`G_{JIK}(\cos(\theta))` is an angular
function, :math:`P_{IK}(\Delta r_{jik})` is a function of atomic spacing
differential :math:`\Delta r_{jik} = r_{ij} - \xi_{IJ} \cdot r_{ik}`
with :math:`\xi_{IJ}` being a pair-dependent parameter, and
:math:`F_{IJ}(X_{ij})` is a function of the local environment variable
:math:`X_{ij}`. This generic potential is fully defined once the
constants :math:`\eta_{ij}` and :math:`\xi_{IJ}`, and the six functions
:math:`U_{IJ}(r_{ij})`, :math:`V_{IJ}(r_{ij})`, :math:`W_{IK}(r_{ik})`,
:math:`G_{JIK}(\cos(\theta))`, :math:`P_{IK}(\Delta r_{jik})`, and
:math:`F_{IJ}(X_{ij})` are given. Note that these six functions are all
one dimensional, and hence can be provided in an analytic or tabular
form. This allows users to design different potentials solely based on a
manipulation of these functions. For instance, the potential reduces to
Stillinger-Weber potential (:ref:`SW <SW>`) if we set

.. math::

   \left\{\begin{array}{l}
   \eta_{ij} = \delta_{ij},\xi_{IJ}=0 \\
   U_{IJ}\left(r\right)=A_{IJ}\cdot\epsilon_{IJ}\cdot \left(\frac{\sigma_{IJ}}{r}\right)^q\cdot \left[B_{IJ}\cdot \left(\frac{\sigma_{IJ}}{r}\right)^{p-q}-1\right]\cdot exp\left(\frac{\sigma_{IJ}}{r-a_{IJ}\cdot \sigma_{IJ}}\right) \\
   V_{IJ}\left(r\right)=\sqrt{\lambda_{IJ}\cdot \epsilon_{IJ}}\cdot exp\left(\frac{\gamma_{IJ}\cdot \sigma_{IJ}}{r-a_{IJ}\cdot \sigma_{IJ}}\right) \\
   F_{IJ}\left(X\right)=-X \\
   P_{IJ}\left(\Delta r\right)=1 \\
   W_{IJ}\left(r\right)=\sqrt{\lambda_{IJ}\cdot \epsilon_{IJ}}\cdot exp\left(\frac{\gamma_{IJ}\cdot \sigma_{IJ}}{r-a_{IJ}\cdot \sigma_{IJ}}\right) \\
   G_{JIK}\left(\theta\right)=\left(cos\theta+\frac{1}{3}\right)^2
   \end{array}\right.

The potential reduces to Tersoff types of potential
(:ref:`Tersoff <Tersoff>` or :ref:`Albe <poly-Albe>`) if we set

.. math::

   \left\{\begin{array}{l}
   \eta_{ij}=\delta_{ij},\xi_{IJ}=1 \\
   U_{IJ}\left(r\right)=\frac{D_{e,IJ}}{S_{IJ}-1}\cdot exp\left[-\beta_{IJ}\sqrt{2S_{IJ}\left(r-r_{e,IJ}\right)}\right]\cdot f_{c,IJ}\left(r\right) \\
   V_{IJ}\left(r\right)=\frac{S_{IJ}\cdot D_{e,IJ}}{S_{IJ}-1}\cdot exp\left[-\beta_{IJ}\sqrt{\frac{2}{S_{IJ}}\left(r-r_{e,IJ}\right)}\right]\cdot f_{c,IJ}\left(r\right) \\
   F_{IJ}\left(X\right)=\left(1+X\right)^{-\frac{1}{2}} \\
   P_{IJ}\left(\Delta r\right)=exp\left(2\mu_{IK}\cdot \Delta r\right) \\
   W_{IJ}\left(r\right)=f_{c,IK}\left(r\right) \\
   G_{JIK}\left(\theta\right)=\gamma_{IK}\left[1+\frac{c_{IK}^2}{d_{IK}^2}-\frac{c_{IK}^2}{d_{IK}^2+\left(h_{IK}+cos\theta\right)^2}\right]
   \end{array}\right.

.. math::

   f_{c,IJ}=\left\{\begin{array}{lr}
   1, & r\leq r_{s,IJ} \\
   \frac{1}{2}+\frac{1}{2} cos \left[\frac{\pi \left(r-r_{s,IJ}\right)}{r_{c,IJ}-r_{s,IJ}}\right], & r_{s,IJ}<r<r_{c,IJ} \\
   0, & r \geq r_{c,IJ} \\
   \end{array}\right.

The potential reduces to Rockett-Tersoff (:ref:`Wang <Wang3>`) type if we set

.. math::

   \left\{\begin{array}{l}
   \eta_{ij}=\delta_{ij},\xi_{IJ}=1 \\
   U_{IJ}\left(r\right)=\left\{\begin{array}{lr}
   A_{IJ}\cdot exp\left(-\lambda_{1,IJ}\cdot r\right)\cdot f_{c,IJ}\left(r\right), & r\leq r_{s,1,IJ} \\
   A_{IJ}\cdot exp\left(-\lambda_{1,IJ}\cdot r\right)\cdot f_{c,IJ}\left(r\right)\cdot f_{c,1,IJ}\left(r\right), & r_{s,1,IJ}<r<r_{c,1,IJ} \\
   0, & r\ge r_{c,1,IJ}
   \end{array}\right. \\
   V_{IJ}\left(r\right)=\left\{\begin{array}{lr}
   B_{IJ} \cdot exp\left(-\lambda_{2,IJ}\cdot r\right)\cdot f_{c,IJ}\left(r\right), & r\le r_{s,1,IJ} \\
   B_{IJ} \cdot exp\left(-\lambda_{2,IJ}\cdot r\right)\cdot f_{c,IJ}\left(r\right)+A_{IJ}\cdot exp\left(-\lambda_{1,IJ}\cdot r\right)\cdot & \\ ~~~~~~ f_{c,IJ}\left(r\right)\cdot \left[1-f_{c,1,IJ}\left(r\right)\right], & r_{s,1,IJ}<r<r_{c,1,IJ} \\
   B_{IJ} \cdot exp\left(-\lambda_{2,IJ}\cdot r\right)\cdot f_{c,IJ}\left(r\right)+A_{IJ}\cdot exp\left(-\lambda_{1,IJ}\cdot r\right)\cdot & \\ ~~~~~~ f_{c,IJ}\left(r\right) & r \ge r_{c,1,IJ}
   \end{array}\right. \\
   F_{IJ}\left(X\right)=\left[1+\left(\beta_{IJ}\cdot X\right)^{n_{IJ}}\right]^{-\frac{1}{2n_{IJ}}} \\
   P_{IJ}\left(\Delta r\right)=exp\left(\lambda_{3,IK}\cdot \Delta r^3\right) \\
   W_{IJ}\left(r\right)=f_{c,IK}\left(r\right) \\
   G_{JIK}\left(\theta\right)=1+\frac{c_{IK}^2}{d_{IK}^2}-\frac{c_{IK}^2}{d_{IK}^2+\left(h_{IK}+cos\theta\right)^2}
   \end{array}\right.

.. math::

   f_{c,IJ}=\left\{\begin{array}{lr}
   1, & r\leq r_{s,IJ} \\
   \frac{1}{2}+\frac{1}{2} cos \left[\frac{\pi \left(r-r_{s,IJ}\right)}{r_{c,IJ}-r_{s,IJ}}\right], & r_{s,IJ}<r<r_{c,IJ} \\
   0, & r \geq r_{c,IJ} \\
   \end{array}\right.

.. math::

   f_{c,1,IJ}=\left\{\begin{array}{lr}
   1, & r\leq r_{s,1,IJ} \\
   \frac{1}{2}+\frac{1}{2} cos \left[\frac{\pi \left(r-r_{s,1,IJ}\right)}{r_{c,1,IJ}-r_{s,1,IJ}}\right], & r_{s,1,IJ}<r<r_{c,1,IJ} \\
   0, & r \geq r_{c,1,IJ} \\
   \end{array}\right.

The potential becomes embedded atom method (:ref:`Daw <poly-Daw>`) if we set

.. math::

   \left\{\begin{array}{l}
   \eta_{ij}=1-\delta_{ij},\xi_{IJ}=0 \\
   U_{IJ}\left(r\right)=\phi_{IJ}\left(r\right) \\
   V_{IJ}\left(r\right)=1 \\
   F_{II}\left(X\right)=-2F_I\left(X\right) \\
   P_{IJ}\left(\Delta r\right)=1 \\
   W_{IJ}\left(r\right)=f_{K}\left(r\right) \\
   G_{JIK}\left(\theta\right)=1
   \end{array}\right.

In the embedded atom method case, :math:`\phi_{IJ}(r_{ij})` is the pair
energy, :math:`F_I(X)` is the embedding energy, *X* is the local
electron density, and :math:`f_K(r)` is the atomic electron density function.

If the tabulated functions are created using the parameters of sw,
tersoff, and eam potentials, the polymorphic pair style will produce
the same global properties (energies and stresses) and the same forces
as the sw, tersoff, and eam pair styles. The polymorphic pair style
also produces the same atom properties (energies and stresses) as the
corresponding tersoff and eam pair styles. However, due to a different
partition of global properties to atom properties, the polymorphic
pair style will produce different atom properties (energies and
stresses) as the sw pair style. This does not mean that polymorphic
pair style is different from the sw pair style in this case. It just
means that the definitions of the atom energies and atom stresses are
different.

Only a single pair_coeff command is used with the polymorphic style
which specifies an potential file for all needed elements. These are
mapped to LAMMPS atom types by specifying N additional arguments after
the filename in the pair_coeff command, where N is the number of
LAMMPS atom types:

* filename
* N element names = mapping of Tersoff elements to atom types

See the pair_coeff doc page for alternate ways to specify the path for
the potential file.  Several files for polymorphic potentials are
included in the potentials directory of the LAMMPS distribution.  They
have a "poly" suffix.

As an example, imagine the SiC_tersoff.poly file has tabulated
functions for Si-C tersoff potential. If your LAMMPS simulation has 4
atoms types and you want the 1st 3 to be Si, and the 4th to be C, you
would use the following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * SiC_tersoff.poly Si Si Si C

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom
types. The first three Si arguments map LAMMPS atom types 1,2,3 to the
Si element in the polymorphic file. The final C argument maps LAMMPS
atom type 4 to the C element in the polymorphic file. If a mapping
value is specified as NULL, the mapping is not performed. This can be
used when an polymorphic potential is used as part of the hybrid pair
style. The NULL values are placeholders for atom types that will be
used with other potentials.

Potential files in the potentials directory of the LAMMPS distribution
have a ".poly" suffix. At the beginning of the files, an unlimited
number of lines starting with '#' are used to describe the potential
and are ignored by LAMMPS. The next line lists two numbers:

.. parsed-literal::

   ntypes :math:`\eta`

Here ntypes represent total number of species defined in the potential
file, and :math:`\eta = 0` or 1. The number ntypes must equal the total
number of different species defined in the pair_coeff command. When
:math:`\eta = 1`, :math:\eta_{ij}` defined in the potential functions
above is set to :math:`1 - \delta_{ij}`, otherwise :math:`\eta_{ij}` is
set to :math:`\delta_{ij}`. The next ntypes lines each lists two numbers
and a character string representing atomic number, atomic mass, and name
of the species of the ntypes elements:

.. parsed-literal::

   atomic_number atomic-mass element (1)
   atomic_number atomic-mass element (2)
   ...
   atomic_number atomic-mass element (ntypes)

The next ntypes\*(ntypes+1)/2 lines contain two numbers:

.. parsed-literal::

   cut :math:`xi` (1)
   cut :math:`xi` (2)
   ...
   cut :math:`xi` (ntypes\*(ntypes+1)/2)

Here cut means the cutoff distance of the pair functions, :math:`\xi` is
the same as defined in the potential functions above. The
ntypes\*(ntypes+1)/2 lines are related to the pairs according to the
sequence of first ii (self) pairs, i = 1, 2, ..., ntypes, and then then
ij (cross) pairs, i = 1, 2, ..., ntypes-1, and j = i+1, i+2, ..., ntypes
(i.e., the sequence of the ij pairs follows 11, 22, ..., 12, 13, 14,
..., 23, 24, ...).

The final blocks of the potential file are the U, V, W, P, G, and F
functions are listed sequentially. First, U functions are given for
each of the ntypes\*(ntypes+1)/2 pairs according to the sequence
described above. For each of the pairs, nr values are listed. Next,
similar arrays are given for V, W, and P functions. Then G functions
are given for all the ntypes\*ntypes\*ntypes ijk triplets in a natural
sequence i from 1 to ntypes, j from 1 to ntypes, and k from 1 to
ntypes (i.e., ijk = 111, 112, 113, ..., 121, 122, 123 ..., 211, 212,
...). Each of the ijk functions contains ng values. Finally, the F
functions are listed for all ntypes\*(ntypes+1)/2 pairs, each
containing nx values. Either analytic or tabulated functions can be
specified. Currently, constant, exponential, sine and cosine analytic
functions are available which are specified with: constant c1 , where
f(x) = c1 exponential c1 c2 , where f(x) = c1 exp(c2\*x) sine c1 c2 ,
where f(x) = c1 sin(c2\*x) cos c1 c2 , where f(x) = c1 cos(c2\*x)
Tabulated functions are specified by spline n x1 x2, where n=number of
point, (x1,x2)=range and then followed by n values evaluated uniformly
over these argument ranges.  The valid argument ranges of the
functions are between 0 <= r <= cut for the U(r), V(r), W(r)
functions, -cutmax <= delta_r <= cutmax for the P(delta_r) functions,
-1 <= :math:`\cos\theta` <= 1 for the G(:math:`\cos\theta`) functions,
and 0 <= X <= maxX for the F(X) functions.

**Mixing, shift, table tail correction, restart**\ :

This pair styles does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write their information to :doc:`binary restart files <restart>`, since it is stored in potential files. Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

----------

Restrictions
""""""""""""

If using create_atoms command, atomic masses must be defined in the
input script. If using read_data, atomic masses must be defined in the
atomic structure data file.

This pair style is part of the MANYBODY package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

This pair potential requires the :doc:`newtion <newton>` setting to be
"on" for pair interactions.

The potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`. You can use
any LAMMPS units, but you would need to create your own potential
files.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

----------

.. _Zhou3:

**(Zhou)** X. W. Zhou, M. E. Foster, R. E. Jones, P. Yang, H. Fan, and
F. P. Doty, J. Mater. Sci. Res., 4, 15 (2015).

.. _SW:

**(SW)** F. H. Stillinger-Weber, and T. A. Weber, Phys. Rev. B, 31, 5262 (1985).

.. _Tersoff:

**(Tersoff)** J. Tersoff, Phys. Rev. B, 39, 5566 (1989).

.. _poly-Albe:

**(Albe)** K. Albe, K. Nordlund, J. Nord, and A. Kuronen, Phys. Rev. B,
66, 035205 (2002).

.. _Wang3:

**(Wang)** J. Wang, and A. Rockett, Phys. Rev. B, 43, 12571 (1991).

.. _poly-Daw:

**(Daw)** M. S. Daw, and M. I. Baskes, Phys. Rev. B, 29, 6443 (1984).
