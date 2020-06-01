.. index:: fix qmmm

fix qmmm command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID qmmm

* ID, group-ID are documented in :doc:`fix <fix>` command
* qmmm = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 qmol qmmm

Description
"""""""""""

This fix provides functionality to enable a quantum
mechanics/molecular mechanics (QM/MM) coupling of LAMMPS to a quantum
mechanical code.  The current implementation only supports an ONIOM
style mechanical coupling to the `Quantum ESPRESSO <espresso_>`_ plane
wave DFT package.  Electrostatic coupling is in preparation and the
interface has been written in a manner that coupling to other QM codes
should be possible without changes to LAMMPS itself.

.. _espresso: http://www.quantum-espresso.org

The interface code for this is in the lib/qmmm directory of the LAMMPS
distribution and is being made available at this early stage of
development in order to encourage contributions for interfaces to
other QM codes.  This will allow the LAMMPS side of the implementation
to be adapted if necessary before being finalized.

Details about how to use this fix are currently documented in the
description of the QM/MM interface code itself in lib/qmmm/README.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global scalar or vector or per-atom
quantities are stored by this fix for access by various :doc:`output commands <Howto_output>`.  No parameter of this fix can be used
with the *start/stop* keywords of the :doc:`run <run>` command.  This
fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the USER-QMMM package.  It is only enabled if
LAMMPS was built with that package. It also requires building a
library provided with LAMMPS.  See the :doc:`Build package <Build_package>` doc page for more info.

The fix is only functional when LAMMPS is built as a library and
linked with a compatible QM program and a QM/MM front end into a QM/MM
executable.  See the lib/qmmm/README file for details.

**Related commands:** none

**Default:** none
