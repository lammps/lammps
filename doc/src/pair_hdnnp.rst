.. index:: pair_style hdnnp

pair_style hdnnp command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style hdnnp keyword value ...

* zero or more keyword/value pairs may be appended
* keyword = *dir* or *showew* or *showewsum* or *maxew* or *resetew* or *cflength* or *cfenergy*
* value depends on the preceding keyword:

.. parsed-literal::
   
   *emap* value = mapping
        mapping = Element mapping from LAMMPS atom types to n2p2 elements
   *dir* value = directory
       directory = Path to HDNNP configuration files
   *showew* value = *yes* or *no*
   *showewsum* value = summary
             summary = Write EW summary every this many timesteps (*0* turns summary off)
   *maxew* value = threshold
         threshold = Maximum number of EWs allowed
   *resetew* value = *yes* or *no*
   *cflength* value = length
            length = Length unit conversion factor
   *cfenergy* value = energy
            energy = Energy unit conversion factor

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hdnnp showew yes showewsum 100 maxew 1000 resetew yes cflength 1.8897261328 cfenergy 0.0367493254 emap "1:H,2:O"
   pair_coeff * * 6.35

   pair_style hdnnp dir "./" showewsum 10000
   pair_coeff * * 6.01

Description
"""""""""""

This pair style adds an interaction based on the high-dimensional neural network
potential (HDNNP) method as presented in :ref:`(Behler and Parrinello 2007)
<Behler_Parrinello_2007>`. HDNNPs are machine learning potentials which require
careful training of neural networks prior to application in MD simulations. The
pair style uses an interface to the *n2p2* library :ref:`(Singraber, Behler and
Dellago 2019) <Singraber_Behler_Dellago_2019>` which is available on Github
`here <https://github.com/CompPhysVienna/n2p2>`__. Please see the *n2p2*
`documentation <https://compphysvienna.github.io/n2p2/>`__ for further details.
*n2p2* (and hence this pair style) is compatible with neural network potentials
trained with its own tools :ref:`(Singraber et al 2019) <Singraber_et_al_2019>`
and with `RuNNer <https://www.uni-goettingen.de/de/560580.html>`__.

The maximum cutoff radius of all symmetry functions (the atomic environment
descriptors of HDNNPs) is the only argument of the *pair_coeff* command which
should be invoked with asterisk wild-cards only:

.. code-block:: LAMMPS

   pair_coeff * * cutoff

.. note::

   The cutoff must be given in LAMMPS length units, even if the neural network
   potential has been trained using a different unit system (see remarks about the
   *cflength* and *cfenergy* keywords below for details).

The numeric value may be slightly larger than the actual maximum symmetry
function cutoff radius (to account for rounding errors when converting units),
but must not be smaller.

----

If provided, the keyword *emap* determines the mapping from LAMMPS atom types to
n2p2 elements. The format is a comma-separated list of ``atom type:element``
pairs, e.g. ``"1:Cu,2:Zn"`` will map atom types 1 and 2 to elements Cu and Zn,
respectively. Atom types not present in the list will be completely ignored by
the HDNNP. The keyword *emap* is mandatory in a "hybrid" setup (see
:doc:`pair_hybrid <pair_hybrid>`) with "extra" atom types in the simulation
which are not handled by the HDNNP.

.. warning::

   Without an explicit mapping it is by default assumed that the atom type
   specifications in LAMMPS configuration files are consistent with the ordering
   of elements in the *n2p2* library. Thus, without the *emap* keyword present
   atom types must be sorted in order of ascending atomic number, e.g. the only
   correct mapping for a configuration containing hydrogen, oxygen and zinc
   atoms would be:
   
   +---------+---------------+-------------+
   | Element | Atomic number | LAMMPS type |
   +=========+===============+=============+
   |       H |             1 |           1 |
   +---------+---------------+-------------+
   |       O |             8 |           2 |
   +---------+---------------+-------------+
   |      Zn |            30 |           3 |
   +---------+---------------+-------------+

Use the *dir* keyword to specify the directory containing the HDNNP configuration
files. The directory must contain ``input.nn`` with neural network and symmetry
function setup, ``scaling.data`` with symmetry function scaling data and
``weights.???.data`` with weight parameters for each element.

The keyword *showew* can be used to turn on/off the display of extrapolation
warnings (EWs) which are issued whenever a symmetry function value is out of
bounds defined by minimum/maximum values in ``scaling.data``. An extrapolation
warning may look like this:

.. code-block:: LAMMPS

   ### NNP EXTRAPOLATION WARNING ### STRUCTURE:      0 ATOM:       119 ELEMENT: Cu SYMFUNC:   32 TYPE:  3 VALUE:  2.166E-02 MIN:  2.003E-05 MAX:  1.756E-02

stating that the value 2.166E-02 of symmetry function 32 of type 3 (angular
narrow), element Cu (see the log file for a symmetry function listing) was out
of bounds (maximum in ``scaling.data`` is 1.756E-02) for atom 119. Here, the
atom index refers to the LAMMPS tag (global index) and the structure index is
used to print out the MPI rank the atom belongs to.

.. note::

   The *showew* keyword should only be set to *yes* for debugging purposes.
   Extrapolation warnings may add lots of overhead as they are communicated each
   timestep. Also, if the simulation is run in a region where the HDNNP was not
   correctly trained, lots of extrapolation warnings may clog log files and the
   console. In a production run use *showewsum* instead.

The keyword *showewsum* can be used to get an overview of extrapolation warnings
occurring during an MD simulation. The argument specifies the interval at which
extrapolation warning summaries are displayed and logged. An EW summary may look
like this:

.. code-block:: LAMMPS

   ### NNP EW SUMMARY ### TS:        100 EW         11 EWPERSTEP  1.100E-01

Here, at timestep 100 the occurrence of 11 extrapolation warnings since the last
summary is reported, which corresponds to an EW rate of 0.11 per timestep.
Setting *showewsum* to 0 deactivates the EW summaries.

A maximum number of allowed extrapolation warnings can be specified with the
*maxew* keyword. If the number of EWs exceeds the *maxew* argument the
simulation is stopped. Note however that this is merely an approximate threshold
since the check is only performed at the end of each timestep and each MPI
process counts individually to minimize communication overhead.

The keyword *resetew* alters the behavior of the above mentioned *maxew*
threshold. If *resetew* is set to *yes* the threshold is applied on a
per-timestep basis and the internal EW counters are reset at the beginning of
each timestep. With *resetew* set to *no* the counters accumulate EWs along the
whole trajectory.

If the training of a neural network potential has been performed with different
physical units for length and energy than those set in LAMMPS, it is still
possible to use the potential when the unit conversion factors are provided via
the *cflength* and *cfenergy* keywords. If for example, the HDNNP was
parameterized with Bohr and Hartree training data and symmetry function
parameters (i.e. distances and energies in "input.nn" are given in Bohr and
Hartree) but LAMMPS is set to use *metal* units (Angstrom and eV) the correct
conversion factors are:

.. code-block:: LAMMPS

   cflength 1.8897261328

   cfenergy 0.0367493254

Thus, arguments of *cflength* and *cfenergy* are the multiplicative factors
required to convert lengths and energies given in LAMMPS units to respective
quantities in native HDNNP units (1 Angstrom = 1.8897261328 Bohr, 1 eV =
0.0367493254 Hartree).

----

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support mixing. The :doc:`pair_coeff <pair_coeff>` command
should only be invoked with asterisk wild cards (see above).

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.
Thus, you need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa
<run_style>` command.  It does not support the *inner*\ , *middle*\ , *outer*
keywords.

Restrictions
""""""""""""

This pair style is part of the USER-HDNNP package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Please report bugs and feature requests to the `n2p2 GitHub issue page
<https://github.com/CompPhysVienna/n2p2/issues>`__.

Related commands
^^^^^^^^^^^^^^^^

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_hybrid <pair_hybrid>`, :doc:`units <units>`

Default
^^^^^^^

The default options are *dir* = "hdnnp/", *showew* = yes, *showewsum* = 0, *maxew*
= 0, *resetew* = no, *cflength* = 1.0, *cfenergy* = 1.0. The default atom type
mapping is determined automatically according to ascending atomic number of
present elements (see above).

----

.. _Behler_Parrinello_2007:

**(Behler and Parrinello 2007)** `Behler, J.; Parrinello, M. Generalized
Neural-Network Representation of High-Dimensional Potential-Energy Surfaces.
Phys. Rev. Lett.  2007, 98 (14), 146401.
<https://doi.org/10.1103/PhysRevLett.98.146401>`__

.. _Singraber_Behler_Dellago_2019:

**(Singraber, Behler and Dellago 2019)** `Singraber, A.; Behler, J.; Dellago, C.
Library-Based LAMMPS Implementation of High-Dimensional Neural Network
Potentials. J. Chem.  Theory Comput. 2019, 15 (3), 1827-1840
<https://doi.org/10.1021/acs.jctc.8b00770>`__

.. _Singraber_et_al_2019:

**(Singraber et al 2019)** `Singraber, A.; Morawietz, T.; Behler, J.; Dellago,
C. Parallel Multistream Training of High-Dimensional Neural Network Potentials.
J. Chem. Theory Comput.  2019, 15 (5), 3075-3092.
<https://doi.org/10.1021/acs.jctc.8b01092>`__
