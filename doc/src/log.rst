.. index:: log

log command
===========

Syntax
""""""

.. parsed-literal::

   log file keyword

* file = name of new logfile
* keyword = *append* if output should be appended to logfile (optional)

Examples
""""""""

.. code-block:: LAMMPS

   log log.equil
   log log.equil append

Description
"""""""""""

This command closes the current LAMMPS log file, opens a new file with
the specified name, and begins logging information to it.  If the
specified file name is *none*, then no new log file is opened.  If the
optional keyword *append* is specified, then output will be appended
to an existing log file, instead of overwriting it.

If multiple processor partitions are being used, the file name should
be a variable, so that different processors do not attempt to write to
the same log file.

The file "log.lammps" is the default log file for a LAMMPS run.  The
name of the initial log file can also be set by the :doc:`-log command-line switch <Run_options>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

none


Default
"""""""

The default LAMMPS log file is named log.lammps
