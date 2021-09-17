LAMMPS open-source license
--------------------------

GPL version of LAMMPS
^^^^^^^^^^^^^^^^^^^^^

LAMMPS is an open-source code, available free-of-charge, and distributed
under the terms of the `GNU Public License Version 2 <gpl_>`_ (GPLv2),
which means you can use or modify the code however you wish for your own
purposes, but have to adhere to certain rules when redistributing it -
specifically in binary form - or are distributing software derived
from it or that includes parts of it.

LAMMPS comes with no warranty of any kind.

As each source file states in its header, it is a copyrighted code, and
thus not in the public domain. For more information about open-source
software and open-source distribution, see `www.gnu.org <gnuorg_>`_
or `www.opensource.org <opensource_>`_.  The legal text of the GPL as it
applies to LAMMPS is in the LICENSE file included in the LAMMPS distribution.

.. _gpl: https://github.com/lammps/lammps/blob/master/LICENSE

.. _lgpl: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html

.. _gnuorg: http://www.gnu.org

.. _opensource: http://www.opensource.org

Here is a more specific summary of what the GPL means for LAMMPS users:

(1) Anyone is free to use, copy, modify, or extend LAMMPS in any way they
choose, including for commercial purposes.

(2) If you **distribute** a modified version of LAMMPS, it must remain
open-source, meaning you are required to distribute **all** of it under
the terms of the GPL.  You should clearly annotate such a modified code
as a derivative version of LAMMPS.

(3) If you release any code that includes or uses LAMMPS source code,
then it must also be open-sourced, meaning you distribute it under
the terms of the GPL.  You may write code that interfaces LAMMPS to
a differently licensed library.  In that case the code that provides
the interface must be licensed GPL, but not necessarily that library
unless you are distributing binaries that require the library to run.

(4) If you give LAMMPS files to someone else, the GPL LICENSE file and
source file headers (including the copyright and GPL notices) should
remain part of the code.


LGPL version of LAMMPS
^^^^^^^^^^^^^^^^^^^^^^

We occasionally make stable LAMMPS releases available under the `GNU
Lesser Public License v2.1 <lgpl_>`_.  This is on request only and with
non-LGPL compliant files removed.  This allows uses linking non-GPL
compatible software with the (otherwise unmodified) LAMMPS library
or loading it dynamically at runtime.  Any **modifications** to
the LAMMPS code however, even with the LGPL licensed version, must still
be made available under the same open source terms as LAMMPS itself.
