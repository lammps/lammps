.. index:: pair\_style cnt/tpm

pair\_style cnt/tpm command
============================

Syntax
""""""


.. parsed-literal::

   pair_style cnt/tpm cut table\_path BendingMode TPMType 

* cut = the cutoff distance
* table\_path = the path to the potential table, the default value is ./
* BendingMode = the parameter defining the type of the bending potential for nanotubes: 0 - harmonic bending :ref:`[1] <Srivastava>`, 1 - anharmonic potential of bending and bending-buckling :ref:`[2] <Zhigilei1>`
* TPMType = the parameter determining the type of the inter-tube interaction term: 0 - segment-segment approach, 1 - segment-chain approach :ref:`[3 <Zhigilei2>`, :ref:`4] <Zhigilei3>`

The parameter BendingMode also affects the calculation of the inter-tube interaction term when TPMType = 1. In this case, when BendingMode = 1, each continuous chain of segments is additionally replaced by a number of sub-chains divided by bending buckling kinks.

Examples
""""""""


.. parsed-literal::

   pair_style cnt/tpm 25.0 ./ 0 0

Description
"""""""""""

The tubular potential model (TPM) force field is designed for mesoscopic
simulations of interacting flexible nanotubes. The force field is based on the
mesoscopic computational model suggested in Ref. :ref:`[1] <Srivastava>`.
In this model, each nanotube is represented by a chain of mesoscopic elements
in the form of stretchable cylindrical segments, where each segment consists
of multiple atoms. Each nanotube is divided into segments by a sequence of
nodes placed on the nanotube centerline. This sequence of nodes determines the
spatial position of the cylindrical segments and defines the configuration of
the entire tube.

The potential force field that controls the dynamic behavior of a system of
interacting nanotubes is given by the following equation defining the potential
energy of the system:

.. image:: Eqs/cnt_eq.jpg
   :align: center

where U\_str is the harmonic potential describing the stretching of CNTs
:ref:`[1] <Srivastava>`, U\_bnd is the potential for nanotube bending
:ref:`[1] <Srivastava>` and bending-buckling :ref:`[2] <Zhigilei1>`, and
U\_vdW is the potential describing van-der Waals interaction between nanotubes
:ref:`[3 <Zhigilei2>`, :ref:`4] <Zhigilei3>`. The stretching energy, U\_str,
is given by the sum of stretching energies of individual nanotube segments.
The bending energy, U\_bnd, is given by the sum of bending energies in all
internal nanotube nodes. The tube-tube interaction energy, U\_vdW, is calculated
based on the tubular potential method suggested in Ref. :ref:`[3] <Zhigilei2>`.
The tubular potential method is briefly described below.

The interaction between two straight nanotubes of arbitrary length and
orientation is described by the approximate tubular potential developed in
:ref:`[4] <Zhigilei3>`. This potential approximates the results of direct
integration of carbon-carbon interatomic potential over the surfaces of the
interacting nanotubes, with the force sources homogeneously distributed over
the nanotube surfaces. The input data for calculation of tubular potentials
are partially tabulated. For single-walled CNTs of arbitrary chirality, the
tabulated potential data can be generated in the form of ASCII files
TPMSSTP.xrs and TPMA.xrs by the stand-alone code TMDPotGen included in the
tool directory of LAMMPS release. The potential provided with LAMMPS release,
CNT\_10\_10, is tabulated for (10,10) nanotubes.

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
the segment-chain approach. In this case, for each CNT segment, the list of its
neighboring segments is divided into short continuous chains of segments
belonging to individual nanotubes. For each pair of a segment and a chain, the
curved chain is approximated by a straight equivalent nanotube based on the
weighted approach suggested in Ref. :ref:`[3] <Zhigilei2>`. Finally, the
interaction between the segment and straight equivalent chain is calculated
based on the tubular potential. In this case, and in the absence of bending
buckling (i.e., when parameter BendingMode is equal to 0), the tubular
potential method ensures the absence of corrugation of the effective inter-tube
interaction potential for curved nanotubes and eliminates any barriers for the
inter-tube sliding. As a result, the tubular potential method can describe the
spontaneous self-assembly of nanotubes into continuous networks of bundles
:ref:`[2 <Zhigilei1>`, :ref:`4] <Zhigilei3>`.


----------


The TMD force field has been used for generation of nanotube films, fibers,
and vertically aligned forests of nanotubes. Mesoscopic dynamic simulations
were used to prepare realistic structures of continuous networks of nanotube
bundles and to study their structural and mechanical properties
:ref:`[2 <Zhigilei1>`, :ref:`4 <Zhigilei3>` - :ref:`7] <Zhigilei6>`. With
additional models for heat transfer, this force filed was also used to
study the thermal transport properties of carbon nanotube films
:ref:`[8 <Zhigilei7>` - :ref:`10] <Zhigilei9>`. The methods for modeling of
the mechanical energy dissipation into heat (energy exchange between the
dynamic degrees of freedom of the mesoscopic model and the energy of atomic
vibrations that are not explicitly represented in the model) 
:ref:`[11] <Zhigilei10>` and mesoscopic description of covalent cross-links
between nanotubes :ref:`[12] <Banna>` have also been developed but are not
included in this first release of the LAMMPS implementation of the force field.
Further details can be found in references provided below.

The CNT package also provides TMDGen code designed to generate initial samples
composed of straight and dispersed nanotubes of given chirality and length at a
given material density, which is availible in tools directory. In the generated
samples, nanotubes are distributed with random positions and orientations. Both
periodic and free boundary conditions are available along each axis of the
system of coordinates. All parameters in the sample files generated with TMDGen
are given in metal :doc:`units <units>`.

Restrictions
""""""""""""


This pair style is a part of the USER-CNT package, and it is only enabled if
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page for more information.

This pair potential requires use of :doc:`cnt atomic style <atom_style>`.

This pair potential requires the :doc:`newton <newton>` setting to be "on" for
pair interactions.

The cutoff distance should be set to be at least:

.. image:: Eqs/cnt_cut.jpg
   :align: center

where L is the maximum segment length, R is the maximum tube radius, and
T_cut = 10.2 A is the maximum distance between the surfaces of interacting
segments.

The TPMSSTP.xrs and TPMA.xrs potential files provided with LAMMPS (see the
potentials directory) are parameterized for metal :doc:`units <units>`.
You can use the carbon nanotube mesoscopic force field with any LAMMPS units,
but you would need to create your own TPMSSTP.xrs and TPMA.xrs potential files
with coefficients listed in appropriate units, if your simulation
does not use "metal" units.

The chirality parameters set during system generation must match the values
specified during generation of the potential tables.

This pair style has not been developed to support :doc:`hybrid <pair_hybrid>`
pair style and has never been tested for this style.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

----------


.. _Srivastava:



**[1]** Zhigilei, Wei, Srivastava, Phys. Rev. B 71, 165417 (2005).

.. _Zhigilei1:



**[2]** Volkov and Zhigilei, ACS Nano 4, 6187 (2010).

.. _Zhigilei2:



**[3]** Volkov, Simov, Zhigilei, ASME paper IMECE2008, 68021 (2008).

.. _Zhigilei3:



**[4]** Volkov, Zhigilei, J. Phys. Chem. C 114, 5513 (2010).

.. _Zhigilei4:



**[5]** Wittmaack, Banna, Volkov, Zhigilei, Carbon 130, 69 (2018).

.. _Zhigilei5:



**[6]** Wittmaack, Volkov, Zhigilei, Compos. Sci. Technol. 166, 66 (2018).

.. _Zhigilei6:



**[7]** Wittmaack, Volkov, Zhigilei, Carbon 143, 587 (2019).

.. _Zhigilei7:



**[8]** Volkov, Zhigilei, Phys. Rev. Lett. 104, 215902 (2010).

.. _Zhigilei8:



**[9]** Volkov, Shiga, Nicholson, Shiomi, Zhigilei, J. Appl. Phys. 111, 053501 (2012).

.. _Zhigilei9:



**[10]** Volkov, Zhigilei, Appl. Phys. Lett. 101, 043113 (2012).

.. _Zhigilei10:



**[11]** Jacobs, Nicholson, Zemer, Volkov, Zhigilei, Phys. Rev. B 86, 165414 (2012).

.. _Banna:



**[12]** Volkov, Banna, Comp. Mater. Sci. 176, 109410 (2020).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
