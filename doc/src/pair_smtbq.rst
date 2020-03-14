.. index:: pair_style smtbq

pair_style smtbq command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style smtbq

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smtbq
   pair_coeff * * ffield.smtbq.Al2O3 O Al

Description
"""""""""""

This pair style computes a variable charge SMTB-Q (Second-Moment
tight-Binding QEq) potential as described in :ref:`SMTB-Q_1 <SMTB-Q_1>` and
:ref:`SMTB-Q_2 <SMTB-Q_2>`. Briefly, the energy of metallic-oxygen systems
is given by three contributions:

.. math::

   E_{tot} & =  E_{ES} + E_{OO} + E_{MO} \\
   E_{ES}  & =  \sum_i{\biggl[ \chi_{i}^{0}Q_i + \frac{1}{2}J_{i}^{0}Q_{i}^{2} +
   \frac{1}{2} \sum_{j\neq i}{ J_{ij}(r_{ij})f_{cut}^{R_{coul}}(r_{ij})Q_i Q_j } \biggr] } \\
   E_{OO}  & =  \sum_{i,j}^{i,j = O}{\biggl[Cexp( -\frac{r_{ij}}{\rho} ) - Df_{cut}^{r_1^{OO}r_2^{OO}}(r_{ij}) exp(Br_{ij})\biggr]}  \\
   E_{MO}  & =  \sum_i{E_{cov}^{i} + \sum_{j\neq i}{ Af_{cut}^{r_{c1}r_{c2}}(r_{ij})exp\bigl[-p(\frac{r_{ij}}{r_0} -1) \bigr] } }

where :math:`E_{tot}` is the total potential energy of the system,
:math:`E_{ES}` is the electrostatic part of the total energy,
:math:`E_{OO}` is the interaction between oxygen atoms and
:math:`E_{MO}` is a short-range interaction between metal and oxygen
atoms. This interactions depend on interatomic distance :math:`r_{ij}`
and/or the charge :math:`Q_{i}` of atoms *i*\ . Cut-off function enables
smooth convergence to zero interaction.

The parameters appearing in the upper expressions are set in the
ffield.SMTBQ.Syst file where Syst corresponds to the selected system
(e.g. field.SMTBQ.Al2O3). Examples for :math:`\mathrm{TiO_2}`,
:math:`\mathrm{Al_2O_3}` are provided.  A single pair_coeff command
is used with the SMTBQ styles which provides the path to the potential
file with parameters for needed elements. These are mapped to LAMMPS
atom types by specifying additional arguments after the potential
filename in the pair_coeff command. Note that atom type 1 must always
correspond to oxygen atoms. As an example, to simulate a :math:`\mathrm{TiO_2}` system,
atom type 1 has to be oxygen and atom type 2 Ti. The following
pair_coeff command should then be used:

.. code-block:: LAMMPS

   pair_coeff * * PathToLammps/potentials/ffield.smtbq.TiO2 O Ti

The electrostatic part of the energy consists of two components

self-energy of atom *i* in the form of a second order charge dependent
polynomial and a long-range Coulombic electrostatic interaction. The
latter uses the wolf summation method described in :ref:`Wolf <Wolf2>`,
spherically truncated at a longer cutoff, :math:`R_{coul}`. The
charge of each ion is modeled by an orbital Slater which depends on
the principal quantum number (\ *n*\ ) of the outer orbital shared by the
ion.

Interaction between oxygen, :math:`E_{OO}`, consists of two parts,
an attractive and a repulsive part. The attractive part is effective
only at short range (< :math:`r_2^{OO}`). The attractive
contribution was optimized to study surfaces reconstruction
(e.g. :ref:`SMTB-Q_2 <SMTB-Q_2>` in :math:`\mathrm{TiO_2}`) and is not necessary
for oxide bulk modeling. The repulsive part is the Pauli interaction
between the electron clouds of oxygen. The Pauli repulsion and the
coulombic electrostatic interaction have same cut off value. In the
ffield.SMTBQ.Syst, the keyword *'buck'* allows to consider only the
repulsive O-O interactions. The keyword *'buckPlusAttr'* allows to
consider the repulsive and the attractive O-O interactions.

The short-range interaction between metal-oxygen, :math:`E_{MO}` is
based on the second moment approximation of the density of states with
a N-body potential for the band energy term,
:math:`E^i_{cov}`, and a Born-Mayer type repulsive terms
as indicated by the keyword *'second_moment'* in the
ffield.SMTBQ.Syst. The energy band term is given by:

.. math::

   E_{cov}^{i(i=M,O)} & = - \biggl\{\eta_i(\mu \xi^{0})^2 f_{cut}^{r_{c1}r_{c2}}(r_{ij})
   \biggl( \sum_{j(j=O,M)}{ exp[ -2q(\frac{r_{ij}}{r_0} - 1)] } \biggr)
   \delta Q_i \bigl( 2\frac{n_0}{\eta_i} - \delta Q_i \bigr) \biggr\}^{1/2} \\
   \delta Q_i & =  | Q_i^{F} | - | Q_i |

where :math:`\eta_i` is the stoichiometry of atom *i*\ ,
:math:`\delta Q_i` is the charge delocalization of atom *i*\ ,
compared to its formal charge
:math:`Q^F_i`. :math:`n_0`, the number of hybridized
orbitals, is calculated with to the atomic orbitals shared
:math:`d_i` and the stoichiometry
:math:`\eta_i`. :math:`r_{c1}` and :math:`r_{c2}` are the two
cutoff radius around the fourth neighbors in the cutoff function.

In the formalism used here, :math:`\xi^0` is the energy
parameter. :math:`\xi^0` is in tight-binding approximation the
hopping integral between the hybridized orbitals of the cation and the
anion. In the literature we find many ways to write the hopping
integral depending on whether one takes the point of view of the anion
or cation. These are equivalent vision. The correspondence between the
two visions is explained in appendix A of the article in the
SrTiO3 :ref:`SMTB-Q_3 <SMTB-Q_3>` (parameter :math:`\beta` shown in
this article is in fact the :math:`\beta_O`). To summarize the
relationship between the hopping integral :math:`\xi^O`  and the
others, we have in an oxide :math:`\mathrm{C_n O_m}` the following
relationship:

.. math::

   \xi^0 & = \frac{\xi_O}{m} = \frac{\xi_C}{n} \\
   \frac{\beta_O}{\sqrt{m}} & = \frac{\beta_C}{\sqrt{n}} = \xi^0 \frac{\sqrt{m}+\sqrt{n}}{2}

Thus parameter :math:`\mu`, indicated above, is given by :math:`\mu = \frac{1}{2}(\sqrt{n}+\sqrt{m})`

The potential offers the possibility to consider the polarizability of
the electron clouds of oxygen by changing the slater radius of the
charge density around the oxygen atoms through the parameters *rBB, rB and
rS* in the ffield.SMTBQ.Syst. This change in radius is performed
according to the method developed by E. Maras
:ref:`SMTB-Q_2 <SMTB-Q_2>`. This method needs to determine the number of
nearest neighbors around the oxygen. This calculation is based on
first (:math:`r_{1n}`) and second (:math:`r_{2n}`) distances
neighbors.

The SMTB-Q potential is a variable charge potential. The equilibrium
charge on each atom is calculated by the electronegativity
equalization (QEq) method. See :ref:`Rick <Rick3>` for further detail. One
can adjust the frequency, the maximum number of iterative loop and the
convergence of the equilibrium charge calculation. To obtain the
energy conservation in NVE thermodynamic ensemble, we recommend to use
a convergence parameter in the interval 10e-5 -
10e-6 eV.

The ffield.SMTBQ.Syst files are provided for few systems. They consist
of nine parts and the lines beginning with '#' are comments (note that
the number of comment lines matter). The first sections are on the
potential parameters and others are on the simulation options and
might be modified. Keywords are character type and must be enclosed in
quotation marks ('').

1) Number of different element in the oxide:

* N_elem= 2 or 3
* Divider line

2) Atomic parameters

For the anion (oxygen)

* Name of element (char) and stoichiometry in oxide
* Formal charge and mass of element
* Principal quantum number of outer orbital n), electronegativity (:math:`\chi^0_i`) and hardness (:math:`J^0_i`)
* Ionic radius parameters  : max coordination number (\ *coordBB* = 6 by default), bulk coordination number *(coordB)*\ , surface coordination number  *(coordS)* and *rBB, rB and rS*  the slater radius for each coordination number. (**note : If you don't want to change the slater radius, use three identical radius values**)
* Number of orbital shared by the element in the oxide (:math:`d_i`)
* Divider line

For each cations (metal):

* Name of element (char) and stoichiometry in oxide
* Formal charge and mass of element
* Number of electron in outer orbital *(ne)*\ , electronegativity (:math:`\chi^0_i`), hardness (:math:`J^0_i`) and :math:`r_{Slater}` the slater radius for the cation.
* Number of orbitals shared by the elements in the oxide (:math:`d_i`)
* Divider line

3) Potential parameters:

* Keyword for element1, element2 and interaction potential
  ('second_moment' or 'buck' or 'buckPlusAttr') between element 1
  and 2.  If the potential is 'second_moment', specify 'oxide' or
  'metal' for metal-oxygen or metal-metal interactions respectively.
* Potential parameter:

  - If type of potential is 'second_moment' : A (eV), *p*,
    :math:`\zeta^0` (eV) and *q*, :math:`r_{c1} (\mathrm{\mathring{A}})`, :math:`r_{c2}
    (\mathrm{\mathring{A}})` and :math:`r_0 (\mathrm{\mathring{A}})`
  - If type of potential is 'buck' : *C* (eV) and :math:`\rho (\mathrm{\mathring{A}})`
  - If type of potential is 'buckPlusAttr' : *C* (eV) and :math:`\rho
    (\mathrm{\mathring{A}})` *D* (eV), *B* :math:`(\mathrm{\mathring{A}}^{-1})`, :math:`r^{OO}_1 (\mathrm{\mathring{A}})` and
    :math:`r^{OO}_2 (\mathrm{\mathring{A}})`
* Divider line

4) Tables parameters:

* Cutoff radius for the Coulomb interaction (:math:`R_{coul}`)
* Starting radius (:math:`r_{min} = 1,18845 \mathrm{\mathring{A}}`) and increments
  (:math:`dr = 0.001 \mathrm{\mathring{A}}`) for creating the potential table.
* Divider line

5) Rick model parameter:

* *Nevery* : parameter to set the frequency of the charge
  resolution. The charges are evaluated each *Nevery* time steps.
* Max number of iterative loop (\ *loopmax*\ ) and convergence criterion
  (\ *prec*\ ) in eV of the charge resolution
* Divider line

6) Coordination parameter:

* First (:math:`r_{1n}`) and second (:math:`r_{2n}`) neighbor distances
  in angstrom
* Divider line

7) Charge initialization mode:

* Keyword (\ *QInitMode*\ ) and initial oxygen charge
  (:math:`Q_{init}`). If keyword = 'true', all oxygen charges are
  initially set equal to :math:`Q_{init}`. The charges on the cations
  are initially set in order to respect the neutrality of the box. If
  keyword = 'false', all atom charges are initially set equal to 0 if
  you use the :doc:`create_atoms <create_atoms>` command or the charge
  specified in the file structure using :doc:`read_data <read_data>`
  command.
* Divider line

8) Mode for the electronegativity equalization (Qeq)

* Keyword (\ *mode*\ ) followed by:

  - QEqAll  (one QEq group) \|   no parameters
  - QEqAllParallel (several QEq groups) \|   no parameters
  - Surface \|   zlim   (QEq only for z>zlim)

* Parameter if necessary
* Divider line

9) Verbose

* If you want the code to work in verbose mode or not : 'true' or 'false'
* If you want to print or not in the file 'Energy_component.txt' the
  three main contributions to the energy of the system according to the
  description presented above : 'true' or 'false' and
  :math:`N_{Energy}`. This option writes to the file every
  :math:`N_{Energy}` time steps. If the value is 'false' then
  :math:`N_{Energy} = 0`. The file takes into account the possibility to
  have several QEq groups *g* then it writes: time step, number of atoms
  in group *g*\ , electrostatic part of energy, :math:`E_{ES}`, the
  interaction between oxygen, :math:`E_{OO}`, and short range
  metal-oxygen interaction, :math:`E_{MO}`.
* If you want to print to the file 'Electroneg_component.txt' the
  electronegativity component (:math:`\frac{\partial E_{tot}}{\partial
  Q_i}`) or not: 'true' or 'false' and :math:`N_{Electroneg}`. This
  option writes to the file every :math:`N_{Electroneg}` time steps. If
  the value is 'false' then :math:`N_{Electroneg} = 0`.  The file
  consist of atom number *i*\ , atom type (1 for oxygen and # higher
  than 1 for metal), atom position: *x*\ , *y* and *z*\ , atomic charge
  of atom *i*\ , electrostatic part of atom *i* electronegativity,
  covalent part of atom *i* electronegativity, the hopping integral of
  atom *i* :math:`(Z\beta^2)_i` and box electronegativity.

.. note::

   This last option slows down the calculation dramatically.  Use
   only with a single processor simulation.

----------

**Mixing, shift, table, tail correction, restart, rRESPA info:**

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
needs to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

**Restriction:**

This pair style is part of the USER-SMTBQ package and is only enabled
if LAMMPS is built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This potential requires using atom type 1 for oxygen and atom type
higher than 1 for metal atoms.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The SMTB-Q potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.

----------

**Citing this work:**

Please cite related publication: N. Salles, O. Politano, E. Amzallag
and R. Tetot, Comput. Mater. Sci. 111 (2016) 181-189

----------

.. _SMTB-Q_1:

**(SMTB-Q_1)** N. Salles, O. Politano, E. Amzallag, R. Tetot,
Comput. Mater. Sci. 111 (2016) 181-189

.. _SMTB-Q_2:

**(SMTB-Q_2)** E. Maras, N. Salles, R. Tetot, T. Ala-Nissila,
H. Jonsson, J. Phys. Chem. C 2015, 119, 10391-10399

.. _SMTB-Q_3:

**(SMTB-Q_3)** R. Tetot, N. Salles, S. Landron, E. Amzallag, Surface
Science 616, 19-8722 28 (2013)

.. _Wolf2:

**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).

.. _Rick3:

**(Rick)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys 101, 6141
(1994).
