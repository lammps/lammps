.. index:: geturl

geturl command
==============

Syntax
""""""

.. code-block:: LAMMPS

   geturl url keyword args ...

* url = URL of the file to download
* zero or more keyword argument pairs may be provided
* keyword = *output* or *verify* or *overwrite* or *verbose*

  .. parsed-literal::

     *output* filename = write to *filename* instead of inferring the name from the URL
     *verify* yes/no = verify SSL certificate and hostname if *yes*, do not if *no*
     *overwrite* yes/no = if *yes* overwrite the output file in case it exists, do not if *no*
     *verbose* yes/no = if *yes* write verbose debug output from libcurl to screen, do not if *no*

Examples
""""""""

.. code-block:: LAMMPS

   geturl https://www.ctcms.nist.gov/potentials/Download/1990--Ackland-G-J-Vitek-V--Cu/2/Cu2.eam.fs
   geturl https://github.com/lammps/lammps/blob/develop/bench/in.lj output in.bench-lj

Description
"""""""""""

.. versionadded:: 29Aug2024

Download a file from an URL to the local disk. This is implemented with
the `libcurl library <https:://curl.se/libcurl/>`_ which supports a
large variety of protocols including "http", "https", "ftp", "scp",
"sftp", "file".  The transfer will only be performed on MPI rank 0.

The *output* keyword can be used to set the filename. By default, the last part
of the URL is used.

The *verify* keyword determines whether ``libcurl`` will validate the
SSL certificate and hostname for encrypted connections.  Turning this
off may be required when using a proxy or connecting to a server with a
self-signed SSL certificate.

The *overwrite* keyword determines whether a file should be overwritten if it
already exists.  If the argument is *no*, then the download will be skipped
if the file exists.

The *verbose* keyword determines whether a detailed protocol of the steps
performed by libcurl is written to the screen.  Using the argument *yes*
can be used to debug connection issues when the *geturl* command does not
behave as expected.  If the argument is *no*, geturl will operate silently
and only report the error status number provided by libcurl, in case of a
failure.

.. _geturl_proxy:

.. admonition:: Using *geturl* with proxies for http or https
   :class: note

   The `libcurl library <https:://curl.se/libcurl/>`_ supports `routing
   traffic through proxies
   <https://everything.curl.dev/usingcurl/proxies/env.html>`_ by setting
   suitable environment variables (e.g. ``http_proxy`` or
   ``https_proxy``) as required by some institutional or corporate
   security protocols.  In that case you probably also want to use the
   *verify* *no* setting.

   Using a proxy may also be needed if you are running on an HPC cluster
   where only the login or head nodes have access to the internet, but
   not the compute nodes.  In this case the following input can be adapted
   and used for your local HPC environment:

   .. code-block:: LAMMPS

      variable headnode getenv PBS_O_HOST     # use SLURM_SUBMIT_HOST when using SLURM instead of Torque/PBS
      shell ssh -N -f -D 8001 ${headnode}     # start SOCKS5 proxy with backgrounded ssh connection to cluster head node
      shell putenv http_proxy=socks5://localhost:8001 https_proxy=socks5://localhost:8001
      geturl https://download.lammps.org/tars/SHA256SUMS # download a file using proxy
      shell head SHA256SUMS                # check if the download was successful

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
