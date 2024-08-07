.. index:: geturl

geturl command
==============

Syntax
""""""

.. code-block:: LAMMPS

   geturl url keyword args ...

* url = URL of the file to download
* zero or more keyword argument pairs may be provided
* keyword = *output* or *verify* or *overwrite*

  .. parsed-literal::

     *output* filename = write to *filename* instead of inferring the name from the URL
     *verify* yes/no = verify SSL certificate and hostname if *yes*, do not if *no*
     *overwrite* yes/no = if *yes* overwrite the output file in case it exists, do not if *no*

Examples
""""""""

.. code-block:: LAMMPS

   geturl https://www.ctcms.nist.gov/potentials/Download/1990--Ackland-G-J-Vitek-V--Cu/2/Cu2.eam.fs
   geturl https://github.com/lammps/lammps/blob/develop/bench/in.lj output in.bench-lj

Description
"""""""""""

Download a file from an URL to the local disk. This is implemented with
the `libcurl library <https:://curl.se/libcurl/>`_ which supports a
large variety of protocols including "http", "https", "ftp", "scp",
"sftp", "file".  The transfer will only be performed on MPI rank 0.

The *output* keyword can be used to set the filename. By default, the last part
of the URL is used.

The *verify* keyword allows to turn on or off whether ``libcurl`` will validate
the SSL certificate and hostname for encrypted connections.  Turning this off
may be required when using a proxy or connecting to a server with a self-signed
SSL certificate.

The *overwrite* keyword determines whether a file should be overwritten if it
already exists.  If the argument is *no*, then the download will be skipped
if the file exists.

----------

Restrictions
""""""""""""

This command is part of the EXTRA-COMMAND package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.  It also requires that LAMMPS was
built with support for `the libcurl library
<https://curl.se/libcurl/>`_.  See the page about :ref:`Compiling LAMMPS
with libcurl support <libcurl>` for further info.  If support for
libcurl is not included, using *geturl* will trigger an error.

Related commands
""""""""""""""""

:doc:`shell <shell>`

Default
"""""""

*verify* = yes, *overwrite* = yes
