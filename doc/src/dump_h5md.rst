.. index:: dump h5md

dump h5md command
=================

Syntax
""""""


.. parsed-literal::

   dump ID group-ID h5md N file.h5 args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be imaged
* h5md = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file.h5 = name of file to write to

.. parsed-literal::

   args = list of data elements to dump, with their dump "sub-intervals"
     position options
     image
     velocity options
     force options
     species options
     file_from ID: do not open a new file, re-use the already opened file from dump ID
     box value = *yes* or *no*
     create_group value = *yes* or *no*
     author value = quoted string



Note that at least one element must be specified and image may only be
present if position is specified first.

For the elements *position*\ , *velocity*\ , *force* and *species*\ , a
sub-interval may be specified to write the data only every N\_element
iterations of the dump (i.e. every N\*N\_element time steps). This is
specified by this option directly following the element declaration:


.. parsed-literal::

   every N_element



Examples
""""""""


.. parsed-literal::

   dump h5md1 all h5md 100 dump_h5md.h5 position image
   dump h5md1 all h5md 100 dump_h5md.h5 position velocity every 10
   dump h5md1 all h5md 100 dump_h5md.h5 velocity author "John Doe"

Description
"""""""""""

Dump a snapshot of atom coordinates every N timesteps in the
`HDF5 <HDF5-ws_>`_ based `H5MD <h5md_>`_ file format :ref:`(de Buyl) <h5md_cpc>`.
HDF5 files are binary, portable and self-describing.  This dump style
will write only one file, on the root node.

Several dumps may write to the same file, by using file\_from and
referring to a previously defined dump.  Several groups may also be
stored within the same file by defining several dumps.  A dump that
refers (via *file\_from*) to an already open dump ID and that concerns
another particle group must specify *create\_group yes*.

.. _h5md: http://nongnu.org/h5md/



Each data element is written every N\*N\_element steps. For *image*\ , no
sub-interval is needed as it must be present at the same interval as
*position*\ .  *image* must be given after *position* in any case.  The
box information (edges in each dimension) is stored at the same
interval than the *position* element, if present. Else it is stored
every N steps.

.. note::

   Because periodic boundary conditions are enforced only on
   timesteps when neighbor lists are rebuilt, the coordinates of an atom
   written to a dump file may be slightly outside the simulation box.

**Use from write\_dump:**

It is possible to use this dump style with the
:doc:`write_dump <write_dump>` command.  In this case, the sub-intervals
must not be set at all.  The write\_dump command can be used either to
create a new file or to add current data to an existing dump file by
using the *file\_from* keyword.

Typically, the *species* data is fixed. The following two commands
store the position data every 100 timesteps, with the image data, and
store once the species data in the same file.


.. parsed-literal::

   dump h5md1 all h5md 100 dump.h5 position image
   write_dump all h5md dump.h5 file_from h5md1 species


----------


Restrictions
""""""""""""


The number of atoms per snapshot cannot change with the h5md style.
The position data is stored wrapped (box boundaries not enforced, see
note above).  Only orthogonal domains are currently supported. This is
a limitation of the present dump h5md command and not of H5MD itself.

The *h5md* dump style is part of the USER-H5MD package. It is only
enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info. It also requires
(i) building the ch5md library provided with LAMMPS (See the :doc:`Build package <Build_package>` doc page for more info.) and (ii) having
the `HDF5 <HDF5-ws_>`_ library installed (C bindings are sufficient) on
your system.  The library ch5md is compiled with the h5cc wrapper
provided by the HDF5 library.

.. _HDF5-ws: http://www.hdfgroup.org/HDF5/




----------


Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`


----------


.. _h5md\_cpc:



**(de Buyl)** de Buyl, Colberg and Hofling, H5MD: A structured,
efficient, and portable file format for molecular data,
Comp. Phys. Comm. 185(6), 1546-1553 (2014) -
`[arXiv:1308.6382] <http://arxiv.org/abs/1308.6382/>`_.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
