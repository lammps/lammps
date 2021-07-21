.. index:: pair_style mesont/tpm

pair_style mesont/tpm command
=============================

Syntax
""""""


.. parsed-literal::

   pair_style mesont/tpm cut table_path BendingMode TPMType

* cut = the cutoff distance
* table_path = the path to the potential table
* BendingMode = the parameter defining the type of the bending potential for nanotubes: 0 - harmonic bending :ref:`(Srivastava) <Srivastava>`, 1 - anharmonic potential of bending and bending-buckling :ref:`(Zhigilei1) <Zhigilei1>`
* TPMType = the parameter determining the type of the inter-tube interaction term: 0 - segment-segment approach, 1 - segment-chain approach :ref:`(Zhigilei2 <Zhigilei2>`, :ref:`Zhigilei3) <Zhigilei3>`

The segment-segment approach is approximately 5 times slower than segment-chain approximation.
The parameter BendingMode also affects the calculation of the inter-tube interaction term when TPMType = 1. In this case, when BendingMode = 1, each continuous chain of segments is additionally replaced by a number of sub-chains divided by bending buckling kinks.

Examples
""""""""


.. parsed-literal::

   pair_style mesont/tpm 30.0 MESONT-TABTP_10_10.xrs 0 0

Description
"""""""""""

The tubular potential model (TPM) force field is designed for mesoscopic
simulations of interacting flexible nanotubes. The force field is based on the
mesoscopic computational model suggested in Ref. :ref:`(Srivastava) <Srivastava>`.
In this model, each nanotube is represented by a chain of mesoscopic elements
in the form of stretchable cylindrical segments, where each segment consists
of multiple atoms. Each nanotube is divided into segments by a sequence of
nodes placed on the nanotube centerline. This sequence of nodes determines the
spatial position of the cylindrical segments and defines the configuration of
the entire tube.

The potential force field that controls the dynamic behavior of a system of
interacting nanotubes is given by the following equation defining the potential
energy of the system:

.. math::

   U = U_{str} + U_{bnd} + U_{vdW}

where :math:`U_{str}`  is the harmonic potential describing the stretching of nanotube
:ref:`(Srivastava) <Srivastava>`, :math:`U_{bnd}`  is the potential for nanotube bending
:ref:`(Srivastava) <Srivastava>` and bending-buckling :ref:`(Zhigilei1) <Zhigilei1>`, and
:math:`U_{vdW}`  is the potential describing van-der Waals interaction between nanotubes
:ref:`(Zhigilei2 <Zhigilei2>`, :ref:`Zhigilei3) <Zhigilei3>`. The stretching energy, :math:`U_{str}` ,
is given by the sum of stretching energies of individual nanotube segments.
The bending energy, :math:`U_{bnd}` , is given by the sum of bending energies in all
internal nanotube nodes. The tube-tube interaction energy, :math:`U_{vdW}` , is calculated
based on the tubular potential method suggested in Ref. :ref:`(Zhigilei2) <Zhigilei2>`.
The tubular potential method is briefly described below.

The interaction between two straight nanotubes of arbitrary length and
orientation is described by the approximate tubular potential developed in
:ref:`(Zhigilei3) <Zhigilei3>`. This potential approximates the results of direct
integration of carbon-carbon interatomic potential over the surfaces of the
interacting nanotubes, with the force sources homogeneously distributed over
the nanotube surfaces. The input data for calculation of tubular potentials
are partially tabulated. For single-walled CNTs of arbitrary chirality, the
tabulated potential data can be generated in the form of ASCII files
TPMSSTP.xrs and TPMA.xrs by the stand-alone code TMDPotGen included in the
tool directory of LAMMPS release. The potential provided with LAMMPS release,
MESONT-TABTP_10_10.xrs, is tabulated for (10,10) nanotubes.

Calculations of the interaction between curved or bent nanotubes are performed
on either segment-segment or segment-chain basis. In the first case, activated
when parameter TPMType is equal to 0, the tubular potential is calculated for
each pair of interacting mesoscopic segments. In this case, however, small
potential barriers for inter-tube sliding are introduced. While relatively
small, these barriers are still larger than the ones that originate from the
atomic-scale corrugation in atomistic modeling of inter-tube interaction. The
latter are too weak to prevent room-temperature rearrangements of defect-free
CNT, while the artificial mesoscopic barriers due to the segment-segment
interaction can impede sliding of nanotubes with respect to each other and
affect the kinetics of structural rearrangements in a system of nanotubes at
moderate mesoscopic temperatures. In the second case, activated when parameter
TPMType is equal to 1, the inter-tube interaction term is calculated based on
the segment-chain approach. In this case, for each NT segment, the list of its
neighboring segments is divided into short continuous chains of segments
belonging to individual nanotubes. For each pair of a segment and a chain, the
curved chain is approximated by a straight equivalent nanotube based on the
weighted approach suggested in Ref. :ref:`(Zhigilei2) <Zhigilei2>`. Finally, the
interaction between the segment and straight equivalent chain is calculated
based on the tubular potential. In this case, and in the absence of bending
buckling (i.e., when parameter BendingMode is equal to 0), the tubular
potential method ensures the absence of corrugation of the effective inter-tube
interaction potential for curved nanotubes and eliminates any barriers for the
inter-tube sliding. As a result, the tubular potential method can describe the
spontaneous self-assembly of nanotubes into continuous networks of bundles
:ref:`(Zhigilei1 <Zhigilei1>`, :ref:`Zhigilei3) <Zhigilei3>`.


----------


The TMD force field has been used for generation of nanotube films, fibers,
and vertically aligned forests of nanotubes. Mesoscopic dynamic simulations
were used to prepare realistic structures of continuous networks of nanotube
bundles and to study their structural and mechanical properties
:ref:`(Zhigilei1 <Zhigilei1>`, :ref:`Zhigilei3 <Zhigilei3>`, :ref:`Zhigilei4 <Zhigilei4>`,
:ref:`Zhigilei5 <Zhigilei5>`, :ref:`Zhigilei6) <Zhigilei6>`. With
additional models for heat transfer, this force filed was also used to
study the thermal transport properties of carbon nanotube films
:ref:`(Zhigilei7 <Zhigilei7>`, :ref:`Zhigilei8 <Zhigilei8>`, :ref:`Zhigilei8) <Zhigilei8>`.
The methods for modeling of
the mechanical energy dissipation into heat (energy exchange between the
dynamic degrees of freedom of the mesoscopic model and the energy of atomic
vibrations that are not explicitly represented in the model)
:ref:`(Zhigilei10) <Zhigilei10>` and mesoscopic description of covalent cross-links
between nanotubes :ref:`(Banna) <Banna>` have also been developed but are not
included in this first release of the LAMMPS implementation of the force field.
Further details can be found in references provided below.

The MESONT package also provides TMDGen code designed to generate initial samples
composed of straight and dispersed nanotubes of given chirality and length at a
given material density, which is available in tools directory. In the generated
samples, nanotubes are distributed with random positions and orientations. Both
periodic and free boundary conditions are available along each axis of the
system of coordinates. All parameters in the sample files generated with TMDGen
are given in metal :doc:`units <units>`.

Restrictions
""""""""""""


This pair style is a part of the MSEONT package, and it is only enabled if
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page for more information.

This pair potential requires use of :doc:`mesont atomic style <atom_style>`.

This pair potential requires the :doc:`newton <newton>` setting to be "on" for
pair interactions.

The cutoff distance should be set to be at least :math:`max\left[2L,\sqrt{L^2/2+(2R+T_{cut})^2}\right]` ,
where L is the maximum segment length, R is the maximum tube radius, and
:math:`T_{cut}` = 10.2 A is the maximum distance between the surfaces of interacting
segments. Because of the use of extended chain concept at CNT ends, the recommended
cutoff is 3L.

.. note::

   Because of their size, *mesont* style potential files
   are not bundled with LAMMPS.   When compiling LAMMPS from
   source code, the file ``TABTP_10_10.mesont`` should be downloaded
   transparently from `https://download.lammps.org/potentials/TABTP_10_10.mesont <https://download.lammps.org/potentials/TABTP_10_10.mesont>`_

The ``TABTP_10_10.mesont`` potential file is parameterized for metal :doc:`units <units>`.
You can use the carbon nanotube mesoscopic force field with any LAMMPS units,
but you would need to create your own potential files with coefficients listed in
appropriate units, if your simulation does not use "metal" units.

The chirality parameters set during system generation must match the values
specified during generation of the potential tables.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

----------

.. _Srivastava:

**(Srivastava)** Zhigilei, Wei, Srivastava, Phys. Rev. B 71, 165417 (2005).

.. _Zhigilei1:

**(Zhigilei1)** Volkov and Zhigilei, ACS Nano 4, 6187 (2010).

.. _Zhigilei2:

**(Zhigilei2)** Volkov, Simov, Zhigilei, ASME paper IMECE2008, 68021 (2008).

.. _Zhigilei3:

**(Zhigilei3)** Volkov, Zhigilei, J. Phys. Chem. C 114, 5513 (2010).

.. _Zhigilei4:

**(Zhigilei4)** Wittmaack, Banna, Volkov, Zhigilei, Carbon 130, 69 (2018).

.. _Zhigilei5:

**(Zhigilei5)** Wittmaack, Volkov, Zhigilei, Compos. Sci. Technol. 166, 66 (2018).

.. _Zhigilei6:

**(Zhigilei6)** Wittmaack, Volkov, Zhigilei, Carbon 143, 587 (2019).

.. _Zhigilei7:

**(Zhigilei7)** Volkov, Zhigilei, Phys. Rev. Lett. 104, 215902 (2010).

.. _Zhigilei8:

**(Zhigilei8)** Volkov, Shiga, Nicholson, Shiomi, Zhigilei, J. Appl. Phys. 111, 053501 (2012).

.. _Zhigilei9:

**(Zhigilei9)** Volkov, Zhigilei, Appl. Phys. Lett. 101, 043113 (2012).

.. _Zhigilei10:

**(Zhigilei10)** Jacobs, Nicholson, Zemer, Volkov, Zhigilei, Phys. Rev. B 86, 165414 (2012).

.. _Banna:

**(Banna)** Volkov, Banna, Comp. Mater. Sci. 176, 109410 (2020).

