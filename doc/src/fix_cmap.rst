.. index:: fix cmap

fix cmap command
================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID cmap filename

* ID, group-ID are documented in :doc:`fix <fix>` command
* cmap = style name of this fix command
* filename = force-field file with CMAP coefficients

Examples
""""""""

.. code-block:: LAMMPS

   fix            myCMAP all cmap ../potentials/cmap36.data
   read_data      proteinX.data fix myCMAP crossterm CMAP
   fix_modify     myCMAP energy yes

Description
"""""""""""

This command enables CMAP 5-body interactions to be added to
simulations which use the CHARMM force field.  These are relevant for
any CHARMM model of a peptide or protein sequences that is 3 or more
amino-acid residues long; see :ref:`(Buck) <Buck>` and :ref:`(Brooks)
<Brooks2>` for details, including the analytic energy expressions for
CMAP interactions.  The CMAP 5-body terms add additional potential
energy contributions to pairs of overlapping phi-psi dihedrals of
amino-acids, which are important to properly represent their
conformational behavior.

The examples/cmap directory has a sample input script and data file
for a small peptide, that illustrates use of the fix cmap command.

As in the example above, this fix should be used before reading a data
file that contains a listing of CMAP interactions.  The *filename*
specified should contain the CMAP parameters for a particular version
of the CHARMM force field.  Two such files are including in the
lammps/potentials directory: charmm22.cmap and charmm36.cmap.

The data file read by the "read_data" must contain the topology of all
the CMAP interactions, similar to the topology data for bonds, angles,
dihedrals, etc.  Specially it should have a line like this in its
header section:

.. parsed-literal::

   N crossterms

where N is the number of CMAP 5-body interactions.  It should also
have a section in the body of the data file like this with N lines:

.. parsed-literal::

   CMAP

          1       1       8      10      12      18      20
          2       5      18      20      22      25      27
          [...]
          N       3     314     315     317      318    330

The first column is an index from 1 to N to enumerate the CMAP 5-atom
tuples; it is ignored by LAMMPS.  The second column is the "type" of
the interaction; it is an index into the CMAP force field file.  The
remaining 5 columns are the atom IDs of the atoms in the two 4-atom
dihedrals that overlap to create the CMAP interaction.  Note that the
"crossterm" and "CMAP" keywords for the header and body sections match
those specified in the read_data command following the data file name;
see the :doc:`read_data <read_data>` page for more details.

A data file containing CMAP 5-body interactions can be generated from
a PDB file using the charmm2lammps.pl script in the tools/ch2lmp
directory of the LAMMPS distribution.  The script must be invoked with
the optional "-cmap" flag to do this; see the tools/ch2lmp/README file
for more information.  The same conversion script also creates the
file of CMAP coefficient data which is read by this command.

The potential energy associated with CMAP interactions can be output
as described below.  It can also be included in the total potential
energy of the system, as output by the :doc:`thermo_style
<thermo_style>` command, if the :doc:`fix_modify energy <fix_modify>`
command is used, as in the example above.  See the note below about
how to include the CMAP energy when performing an :doc:`energy
minimization <minimize>`.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the list of CMAP cross-terms to :doc:`binary restart
files <restart>`.  See the :doc:`read_restart <read_restart>` command
for info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the potential energy of the CMAP interactions to both
the global potential energy and peratom potential energies of the
system as part of :doc:`thermodynamic output <thermo_style>` or
output by the :doc:`compute pe/atom <compute_pe_atom>` command.  The
default setting for this fix is :doc:`fix_modify energy yes
<fix_modify>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by
this fix to add the contribution due to the CMAP interactions to both
the global pressure and per-atom stress of the system via the
:doc:`compute pressure <compute_pressure>` and :doc:`compute
stress/atom <compute_stress_atom>` commands.  The former can be
accessed by :doc:`thermodynamic output <thermo_style>`.  The default
setting for this fix is :doc:`fix_modify virial yes <fix_modify>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy discussed above.  The scalar value calculated by this fix is
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its forces. Default is the
outermost level.

.. note::

   For energy minimization, if you want the potential energy
   associated with the CMAP terms forces to be included in the total
   potential energy of the system (the quantity being minimized), you
   MUST not disable the :doc:`fix_modify <fix_modify>` *energy* option
   for this fix.

Restrictions
""""""""""""

To function as expected this fix command must be issued *before* a
:doc:`read_data <read_data>` command but *after* a
:doc:`read_restart <read_restart>` command.

This fix can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`fix_modify <fix_modify>`, :doc:`read_data <read_data>`

Default
"""""""

none

----------

.. _Buck:

**(Buck)** Buck, Bouguet-Bonnet, Pastor, MacKerell Jr., Biophys J, 90, L36
(2006).

.. _Brooks2:

**(Brooks)** Brooks, Brooks, MacKerell Jr., J Comput Chem, 30, 1545 (2009). https://doi.org/10.1002/jcc.21287
