.. index:: pair_style mesocnt
.. index:: pair_style mesocnt/viscous

pair_style mesocnt command
==========================

pair_style mesocnt/viscous command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style neigh_cutoff mode neigh_mode

* style = *mesocnt* or *mesocnt/viscous*
* neigh_cutoff = neighbor list cutoff (distance units)
* mode = *chain* or *segment* (optional)
* neigh_mode = *id* or *topology* (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style mesocnt 30.0
   pair_coeff * * C_10_10.mesocnt 2

   pair_style mesocnt/viscous 60.0 chain topology
   pair_coeff * * C_10_10.mesocnt 0.001 20.0 0.2 2 4

Description
"""""""""""

Style *mesocnt* implements a mesoscopic potential for the interaction
of carbon nanotubes (CNTs), or other quasi-1D objects such as other
kinds of nanotubes or nanowires. In this potential, CNTs are modelled
as chains of cylindrical segments in which each infinitesimal surface
element interacts with all other CNT surface elements with the
Lennard-Jones (LJ) term adopted from the :doc:`airebo <pair_airebo>`
style. The interaction energy is then computed by integrating over the
surfaces of all interacting CNTs.

In LAMMPS, cylindrical segments are represented by bonds. Each segment
is defined by its two end points ("nodes") which correspond to atoms
in LAMMPS. For the exact functional form of the potential and
implementation details, the reader is referred to the original papers
:ref:`(Volkov1) <Volkov1>` and :ref:`(Volkov2) <Volkov2>`.

.. versionchanged:: 15Sep2022

The potential supports two modes, *segment* and *chain*. By default,
*chain* mode is enabled.  In *segment* mode, interactions are
pair-wise between all neighboring segments based on a segment-segment
approach (keyword *segment* in pair_style command).  In *chain* mode,
interactions are calculated between each segment and infinitely or
semi-infinitely long CNTs as described in :ref:`(Volkov1) <Volkov1>`.
Chains of segments are converted to these (semi-)infinite CNTs bases
on an approximate chain approach outlined in :ref:`(Volkov2)
<Volkov2>`. Hence, interactions are calculated on a segment-chain
basis (keyword *chain* in the pair_style command).  Using *chain* mode
allows to simplify the computation of the interactions significantly
and reduces the computational times to the same order of magnitude as
for regular bead spring models where beads interact with the standard
:doc:`pair_lj/cut <pair_lj>` potential. However, this method is only
valid when the curvature of the CNTs in the system is small.  When
CNTs are buckled (see :doc:`angle_mesocnt <angle_mesocnt>`), local
curvature can be very high and the pair_style automatically switches
to *segment* mode for interactions involving buckled CNTs.

The potential further implements two different neighbor list
construction modes. Mode *id* uses atom and mol IDs to construct
neighbor lists while *topology* modes uses only the bond topology of
the system. While *id* mode requires bonded atoms to have consecutive
LAMMPS atom IDs and atoms in different CNTs to have different LAMMPS
molecule IDs, *topology* mode has no such requirement. Using *id* mode
is faster and is enabled by default.

.. note::

  Neighbor *id* mode requires all CNTs in the system to have distinct
  LAMMPS molecule IDs and bonded atoms to have consecutive LAMMPS atom
  IDs. If this is not possible (e.g. in simulations of CNT rings),
  *topology* mode needs to be enabled in the pair_style command.

.. versionadded:: 15Sep2022

In addition to the LJ interactions described above, style
*mesocnt/viscous* explicitly models friction between neighboring
segments. Friction forces are a function of the relative velocity
between a segment and its neighboring approximate chain (even in
*segment* mode) and only act along the axes of the interacting segment
and chain. In this potential, friction forces acting per unit length
of a nanotube segment are modelled as a shifted logistic function:

.. math::

   F^{\text{FRICTION}}(v) / L = \frac{F^{\text{max}}}{1 +
   \exp(-k(v-v_0))} - \frac{F^{\text{max}}}{1 + \exp(k v_0)}

----------

In the pair_style command, the modes described above can be toggled
using the *segment* or *chain* keywords.  The neighbor list cutoff
defines the cutoff within which atoms are included in the neighbor
list for constructing neighboring CNT chains.  This is different from
the potential cutoff, which is directly calculated from parameters
specified in the potential file. We recommend using a neighbor list
cutoff of at least 3 times the maximum segment length used in the
simulation to ensure proper neighbor chain construction.

.. note::

   CNT ends are treated differently by all *mesocnt* styles. Atoms on
   CNT ends need to be assigned different LAMMPS atom types than atoms
   not on CNT ends.

Style *mesocnt* requires tabulated data provided in a single ASCII
text file, as well as a list of integers corresponding to all LAMMPS
atom types representing CNT ends:

* filename
* :math:`N` CNT end atom types

For example, if your LAMMPS simulation of (10, 10) nanotubes has 4
atom types where atom types 1 and 3 are assigned to 'inner' nodes and
atom types 2 and 4 are assigned to CNT end nodes, the pair_coeff
command would be:

.. code-block:: LAMMPS

   pair_coeff * * C_10_10.mesocnt 2 4

Likewise, style *mesocnt/viscous* also requires the same information
as style *mesocnt*, with the addition of 3 parameters for the viscous
friction forces as listed above:

* filename
* :math:`F^{\text{max}}`
* :math:`k`
* :math:`v_0`
* :math:`N` CNT end atom types

Using the same example system as with style *mesocnt* with the
addition of friction, the pair_coeff command is:

.. code-block:: LAMMPS

   pair_coeff * * C_10_10.mesocnt 0.03 20.0 0.20 2 4

Potential files for CNTs can be readily generated using the freely
available code provided on

.. parsed-literal::

   https://github.com/phankl/cntpot

Using the same approach, it should also be possible to generate
potential files for other 1D systems mentioned above.

.. note::

   Because of their size, *mesocnt* style potential files are not
   bundled with LAMMPS.  When compiling LAMMPS from source code, the
   file ``C_10_10.mesocnt`` should be downloaded separately from
   `https://download.lammps.org/potentials/C_10_10.mesocnt
   <https://download.lammps.org/potentials/C_10_10.mesocnt>`_

   The first line of the potential file provides a time stamp and
   general information. The second line lists four integers giving the
   number of data points provided in the subsequent four data
   tables. The third line lists four floating point numbers: the CNT
   radius R, the LJ parameter sigma and two numerical parameters
   delta1 and delta2. These four parameters are given in
   Angstroms. This is followed by four data tables each separated by a
   single empty line. The first two tables have two columns and list
   the parameters uInfParallel and Gamma respectively.  The last two
   tables have three columns giving data on a quadratic array and list
   the parameters Phi and uSemiParallel respectively.  uInfParallel
   and uSemiParallel are given in eV/Angstrom, Phi is given in eV and
   Gamma is unitless.

   If a simulation produces many warnings about segment-chain
   interactions falling outside the interpolation range, we recommend
   generating a potential file with lower values of delta1 and delta2.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These pair styles does not support mixing.

These pair styles does not support the :doc:`pair_modify
<pair_modify>` shift, table, and tail options.

These pair styles do not write their information to :doc:`binary
restart files <restart>`, since it is stored in tabulated potential
files.  Thus, you need to re-specify the pair_style and pair_coeff
commands in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

These styles are part of the MESONT package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

These pair styles require the :doc:`newton <newton>` setting to be
"on" for pair interactions.

These pair styles require all 3 :doc:`special_bonds lj <special_bonds>`
settings to be non-zero for proper neighbor list construction.

Pair style *mesocnt/viscous* requires you to use the :doc:`comm_modify
vel yes <comm_modify>` command so that velocities are stored by ghost
atoms.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`bond_style mesocnt <bond_mesocnt>`,
:doc:`angle_style mesocnt <angle_mesocnt>`

Default
"""""""

mode = chain, neigh_mode = id

----------

.. _Volkov1:

**(Volkov1)** Volkov and Zhigilei, J Phys Chem C, 114, 5513 (2010).

.. _Volkov2:

**(Volkov2)** Volkov, Simov and Zhigilei, APS Meeting Abstracts,
Q31.013 (2008).
