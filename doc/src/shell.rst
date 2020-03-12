.. index:: shell

shell command
=============

Syntax
""""""

.. parsed-literal::

   shell cmd args

* cmd = *cd* or *mkdir* or *mv* or *rm* or *rmdir* or *putenv* or arbitrary command

  .. parsed-literal::

       *cd* arg = dir
         dir = directory to change to
       *mkdir* args = dir1 dir2 ...
         dir1,dir2 = one or more directories to create
       *mv* args = old new
         old = old filename
         new = new filename
       *rm* args = file1 file2 ...
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
   shell mkdir tmp1 tmp2 tmp3
   shell rmdir tmp1
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
the command in a shell.

This is means to invoke other commands from your input script.  For
example, you can move files around in preparation for the next section
of the input script.  Or you can run a program that pre-processes data
for input into LAMMPS.  Or you can run a program that post-processes
LAMMPS output data.

With the exception of *cd*\ , all commands, including ones invoked via a
system() call, are executed by only a single processor, so that
files/directories are not being manipulated by multiple processors.

The *cd* cmd executes the Unix "cd" command to change the working
directory.  All subsequent LAMMPS commands that read/write files will
use the new directory.  All processors execute this command.

The *mkdir* cmd executes the Unix "mkdir" command to create one or
more directories.

The *mv* cmd executes the Unix "mv" command to rename a file and/or
move it to a new directory.

The *rm* cmd executes the Unix "rm" command to remove one or more
files.

The *rmdir* cmd executes the Unix "rmdir" command to remove one or
more directories.  A directory must be empty to be successfully
removed.

The *putenv* cmd defines or updates an environment variable directly.
Since this command does not pass through the shell, no shell variable
expansion or globbing is performed, only the usual substitution for
LAMMPS variables defined with the :doc:`variable <variable>` command is
performed.  The resulting string is then used literally.

Any other cmd is passed as-is to the shell along with its arguments as
one string, invoked by the C-library system() call.  For example,
these lines in your input script:

.. code-block:: LAMMPS

   variable n equal 10
   variable foo string file2
   shell my_setup file1 $n ${foo}

would be the same as invoking

.. code-block:: bash

   % my_setup file1 10 file2

from a command-line prompt.  The executable program "my\_setup" is run
with 3 arguments: file1 10 file2.

Restrictions
""""""""""""

LAMMPS does not detect errors or print warnings when any of these
commands execute.  E.g. if the specified directory does not exist,
executing the *cd* command will silently do nothing.

**Related commands:** none

**Default:** none
