.. index:: kim\_query

kim\_query command
==================

Syntax
""""""


.. parsed-literal::

   kim_query variable query_function web_query_flags

* variable = name of a (string style) variable where the result of the query is stored
* query\_function = name of the OpenKIM web API query function to be used
* web\_query\_flags = a series of keyword=value pairs that represent the web query; supported keywords depend on query function

Examples
""""""""


.. parsed-literal::

   kim_query latconst get_test_result test=TE_156715955670 model=MO_800509458712 &
     prop=structure-cubic-crystal-npt species=["Al"] keys=["a"] units=["angstrom"]

Description
"""""""""""

The kim\_query command allows to retrieve properties from the OpenKIM
through a web query. The result is stored in a string style
:doc:`variable <variable>`, the name of which must be given as the first
argument of the kim\_query command.  The second required argument is the
name of the actual query function (e.g. *get\_test\_result*).  All following
arguments are parameters handed over to the web query in the format
*keyword=value*\ .  The list of supported keywords and the type of how
the value has to be encoded depends on the query function used.  This
mirrors the functionality available on the OpenKIM webpage at
`https://query.openkim.org <https://query.openkim.org/>`_

Restrictions
""""""""""""


This command is part of the KIM package.  It is only enabled if
LAMMPS was built with that package.  Furthermore, its correct
functioning depends on compiling LAMMPS with libcurl support.
See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_style kim <pair_kim>`, :doc:`variable <variable>`


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
