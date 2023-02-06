.. index:: shell

shell command
=============

Syntax
""""""

.. parsed-literal::

   shell command args

* command = *cd* or *mkdir* or *mv* or *rm* or *rmdir* or *putenv* or arbitrary command

  .. parsed-literal::

       *cd* arg = dir
         dir = directory to change to
       *mkdir* args = dir1 dir2 ...
         dir1,dir2 = one or more directories to create
       *mv* args = old new
         old = old filename
         new = new filename or destination folder
       *rm* args = [-f] file1 file2 ...
         -f = turn off warnings (optional)
         file1,file2 = one or more filenames to delete
       *rmdir* args = dir1 dir2 ...
         dir1,dir2 = one or more directories to delete
       *putenv* args = var1=value1 var2=value2
         var=value = one of more definitions of environment variables
       anything else is passed as a command to the shell for direct execution

Examples
""""""""

.. code-block:: LAMMPS

   shell cd sub1
   shell cd ..
   shell mkdir tmp1 tmp2/tmp3
   shell rmdir tmp1 tmp2
   shell mv log.lammps hold/log.1
   shell rm TMP/file1 TMP/file2
   shell putenv LAMMPS_POTENTIALS=../../potentials
   shell my_setup file1 10 file2
   shell my_post_process 100 dump.out

Description
"""""""""""

Execute a shell command.  A few simple file-based shell commands are
supported directly, in Unix-style syntax.  Any command not listed
above is passed as-is to the C-library system() call, which invokes
the command in a shell.  To use the external executable instead of
the built-in version one needs to use a full path, for example
*/bin/rm* instead of *rm*.  The built-in commands will also work
on operating systems, that do not - by default - provide the
corresponding external executables (like *mkdir* on Windows).

This command provides a ways to invoke custom commands or executables
from your input script.  For example, you can move files around in
preparation for the next section of the input script.  Or you can run a
program that pre-processes data for input into LAMMPS.  Or you can run a
program that post-processes LAMMPS output data.

With the exception of *cd*, all commands, including ones invoked via a
system() call, are executed by only a single processor, so that
files/directories are not being manipulated by multiple processors
concurrently which may result in unexpected errors or corrupted files.

The *cd* command changes the current working directory similar to
the ``cd`` command.  All subsequent LAMMPS commands that read/write files
will use the new directory.  All processors execute this command.

The *mkdir* command creates directories similar to the Unix ``mkdir -p``
command.  That is, it will attempt to create the entire path of
subdirectories if they do not exist yet.

The *mv* command renames a file and/or moves it to a new directory.
It cannot rename files across filesystem boundaries or between drives.

The *rm* command deletes file similar to the Unix ``rm`` command.

The *rmdir* command deletes directories similar to Unix ``rmdir`` command.
If a directory is not empty, its contents are also removed recursively
similar to the Unix ``rm -r`` command.

The *putenv* command defines or updates an environment variable directly.
Since this command does not pass through the shell, no shell variable
expansion or globbing is performed, only the usual substitution for
LAMMPS variables defined with the :doc:`variable <variable>` command is
performed.  The resulting string is then used literally.

Any other command is passed as-is to the shell along with its arguments as
one string, invoked by the C-library system() call.  For example,
these lines in your input script:

.. code-block:: LAMMPS

   variable n equal 10
   variable foo string file2
   shell my_setup file1 $n ${foo}

would be the same as invoking

.. code-block:: bash

   my_setup file1 10 file2

from a command-line prompt.  The executable program "my_setup" is run
with 3 arguments: file1 10 file2.

Restrictions
""""""""""""

LAMMPS will do a best effort to detect errors and print suitable
warnings, but due to the nature of delegating commands to the C-library
system() call, this is not always reliable.

Related commands
""""""""""""""""

none


Default
"""""""

none
