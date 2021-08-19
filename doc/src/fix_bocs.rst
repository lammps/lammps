.. index:: fix bocs

fix bocs command
================

Syntax
""""""

.. parsed-literal::

   fix  ID group-ID bocs keyword values ...

   keyword = *temp* or *cgiso* or *analytic* or *linear_spline* or *cubic_spline*
     *temp* values = Tstart Tstop Tdamp
     *cgiso* values = Pstart Pstop Pdamp
     *basis set*
       *analytic* values = V_avg N_particles N_coeff Coeff_1 Coeff_2 ... Coeff_N
       *linear_spline* values = input_filename
       *cubic_spline* values = input_filename

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all bocs temp 300.0 300.0 100.0 cgiso 0.986 0.986 1000.0 analytic 66476.015 968 2 245030.10 8962.20

   fix 1 all bocs temp 300.0 300.0 100.0 cgiso 0.986 0.986 1000.0 cubic_spline input_Fv.dat

   thermo_modify press 1_press

Description
"""""""""""

These commands incorporate a pressure correction as described by
Dunn and Noid in :ref:`(Dunn1) <bocs-Dunn1>` to the standard MTTK
barostat by Martyna et. al. in :ref:`(Martyna) <bocs-Martyna>` .
The first half of the command mimics a standard fix npt command:

.. code-block:: LAMMPS

   fix 1 all bocs temp Tstart Tstop Tcoupl cgiso Pstart Pstop Pdamp

The two differences are replacing *npt* with *bocs*, and replacing
*iso*\ /\ *aniso*\ /\ *etc* with *cgiso*\ .
The rest of the command details what form you would like to use for
the pressure correction equation. The choices are: *analytic*, *linear_spline*,
or *cubic_spline*.

With either spline method, the only argument that needs to follow it
is the name of a file that contains the desired pressure correction
as a function of volume. The file must be formatted so each line has:

.. parsed-literal::

   Volume_i, PressureCorrection_i

Note both the COMMA and the SPACE separating the volume's
value and its corresponding pressure correction. The volumes in the file
must be uniformly spaced. Both the volumes and the pressure corrections
should be provided in the proper units, e.g. if you are using *units real*,
the volumes should all be in cubic angstroms, and the pressure corrections
should all be in atmospheres. Furthermore, the table should start/end at a
volume considerably smaller/larger than you expect your system to sample
during the simulation. If the system ever reaches a volume outside of the
range provided, the simulation will stop.

With the *analytic* option, the arguments are as follows:

.. parsed-literal::

   ... analytic V_avg N_particles N_coeff Coeff_1 Coeff_2 ... Coeff_N

Note that *V_avg* and *Coeff_i* should all be in the proper units, e.g. if you
are using *units real*, *V_avg* should be in cubic angstroms, and the
coefficients should all be in atmospheres \* cubic angstroms.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the cumulative global energy change to :doc:`binary
restart files <restart>`.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the fix continues in an uninterrupted
fashion.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

The cumulative energy change in the system imposed by this fix is
included in the :doc:`thermodynamic output <thermo_style>` keywords
*ecouple* and *econserve*.  See the :doc:`thermo_style <thermo_style>`
doc page for details.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the same
cumulative energy change due to this fix described in the previous
paragraph.  The scalar value calculated by this fix is "extensive".

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

As this is computing a (modified) pressure, group-ID should be *all*\ .

The pressure correction has only been tested for use with an isotropic
pressure coupling in 3 dimensions.

By default, LAMMPS will still report the normal value for the pressure
if the pressure is printed via a *thermo* command, or if the pressures
are written to a file every so often. In order to have LAMMPS report the
modified pressure, you must include the *thermo_modify* command given in
the examples. For the last argument in the command, you should put
XXXX_press, where XXXX is the ID given to the fix bocs command (in the
example, the ID of the fix bocs command is 1 ).

This fix is part of the BOCS package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Further information
"""""""""""""""""""

For more details about the pressure correction and the entire BOCS software
package, visit the `BOCS package on GitHub <bocsgithub_>`_ and read the release
paper by Dunn et. al. :ref:`(Dunn2) <bocs-Dunn2>` .

.. _bocsgithub: https://github.com/noid-group/BOCS

----------

.. _bocs-Dunn1:

**(Dunn1)** Dunn and Noid, J Chem Phys, 143, 243148 (2015).

.. _bocs-Martyna:

**(Martyna)** Martyna, Tobias, and Klein, J Chem Phys, 101, 4177 (1994).

.. _bocs-Dunn2:

**(Dunn2)** Dunn, Lebold, DeLyser, Rudzinski, and Noid, J. Phys. Chem. B, 122, 3363 (2018).
