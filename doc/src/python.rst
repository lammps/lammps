.. index:: python

python command
==============

Syntax
""""""


.. parsed-literal::

   python func keyword args ...

* func = name of Python function
* one or more keyword/args pairs must be appended
  
  .. parsed-literal::
  
     keyword = *invoke* or *input* or *return* or *format* or *length* or *file* or *here* or *exists* or *source*
       *invoke* arg = none = invoke the previously defined Python function
       *input* args = N i1 i2 ... iN
         N = # of inputs to function
         i1,...,iN = value, SELF, or LAMMPS variable name
           value = integer number, floating point number, or string
           SELF = reference to LAMMPS itself which can be accessed by Python function
           variable = v_name, where name = name of LAMMPS variable, e.g. v_abc
       *return* arg = varReturn
         varReturn = v_name  = LAMMPS variable name which return value of function will be assigned to
       *format* arg = fstring with M characters
         M = N if no return value, where N = # of inputs
         M = N+1 if there is a return value
         fstring = each character (i,f,s,p) corresponds in order to an input or return value
         'i' = integer, 'f' = floating point, 's' = string, 'p' = SELF
       *length* arg = Nlen
         Nlen = max length of string returned from Python function
       *file* arg = filename
         filename = file of Python code, which defines func
       *here* arg = inline
         inline = one or more lines of Python code which defines func
                  must be a single argument, typically enclosed between triple quotes
       *exists* arg = none = Python code has been loaded by previous python command
       *source* arg = *filename* or *inline*
         filename = file of Python code which will be executed immediately
         inline = one or more lines of Python code which will be executed immediately
                  must be a single argument, typically enclosed between triple quotes



Examples
""""""""


.. parsed-literal::

   python pForce input 2 v_x 20.0 return v_f format fff file force.py
   python pForce invoke

   python factorial input 1 myN return v_fac format ii here """
   def factorial(n):
     if n == 1: return n
     return n \* factorial(n-1)
    """

   python loop input 1 SELF return v_value format pf here """
   def loop(lmpptr,N,cut0):
     from lammps import lammps
     lmp = lammps(ptr=lmpptr)

     # loop N times, increasing cutoff each time

     for i in range(N):
       cut = cut0 + i\*0.1
       lmp.set_variable("cut",cut)                 # set a variable in LAMMPS
       lmp.command("pair_style lj/cut ${cut}")   # LAMMPS commands
       lmp.command("pair_coeff \* \* 1.0 1.0")
       lmp.command("run 100")
    """

Description
"""""""""""

Define a Python function or execute a previously defined function or
execute some arbitrary python code.
Arguments, including LAMMPS variables, can be passed to the function
from the LAMMPS input script and a value returned by the Python
function to a LAMMPS variable.  The Python code for the function can
be included directly in the input script or in a separate Python file.
The function can be standard Python code or it can make "callbacks" to
LAMMPS through its library interface to query or set internal values
within LAMMPS.  This is a powerful mechanism for performing complex
operations in a LAMMPS input script that are not possible with the
simple input script and variable syntax which LAMMPS defines.  Thus
your input script can operate more like a true programming language.

Use of this command requires building LAMMPS with the PYTHON package
which links to the Python library so that the Python interpreter is
embedded in LAMMPS.  More details about this process are given below.

There are two ways to invoke a Python function once it has been
defined.  One is using the *invoke* keyword.  The other is to assign
the function to a :doc:`python-style variable <variable>` defined in
your input script.  Whenever the variable is evaluated, it will
execute the Python function to assign a value to the variable.  Note
that variables can be evaluated in many different ways within LAMMPS.
They can be substituted for directly in an input script.  Or they can
be passed to various commands as arguments, so that the variable is
evaluated during a simulation run.

A broader overview of how Python can be used with LAMMPS is given on
the :doc:`Python <Python_head>` doc page.  There is an examples/python
directory which illustrates use of the python command.


----------


The *func* setting specifies the name of the Python function.  The
code for the function is defined using the *file* or *here* keywords
as explained below. In case of the *source* keyword, the name of
the function is ignored.

If the *invoke* keyword is used, no other keywords can be used, and a
previous python command must have defined the Python function
referenced by this command.  This invokes the Python function with the
previously defined arguments and return value processed as explained
below.  You can invoke the function as many times as you wish in your
input script.

If the *source* keyword is used, no other keywords can be used.
The argument can be a filename or a string with python commands,
either on a single line enclosed in quotes, or as multiple lines
enclosed in triple quotes. These python commands will be passed
to the python interpreter and executed immediately without registering
a python function for future execution.

The *input* keyword defines how many arguments *N* the Python function
expects.  If it takes no arguments, then the *input* keyword should
not be used.  Each argument can be specified directly as a value,
e.g. 6 or 3.14159 or abc (a string of characters).  The type of each
argument is specified by the *format* keyword as explained below, so
that Python will know how to interpret the value.  If the word SELF is
used for an argument it has a special meaning.  A pointer is passed to
the Python function which it converts into a reference to LAMMPS
itself.  This enables the function to call back to LAMMPS through its
library interface as explained below.  This allows the Python function
to query or set values internal to LAMMPS which can affect the
subsequent execution of the input script.  A LAMMPS variable can also
be used as an argument, specified as v\_name, where "name" is the name
of the variable.  Any style of LAMMPS variable can be used, as defined
by the :doc:`variable <variable>` command.  Each time the Python
function is invoked, the LAMMPS variable is evaluated and its value is
passed to the Python function.

The *return* keyword is only needed if the Python function returns a
value.  The specified *varReturn* must be of the form v\_name, where
"name" is the name of a python-style LAMMPS variable, defined by the
:doc:`variable <variable>` command.  The Python function can return a
numeric or string value, as specified by the *format* keyword.

As explained on the :doc:`variable <variable>` doc page, the definition
of a python-style variable associates a Python function name with the
variable.  This must match the *func* setting for this command.  For
example these two commands would be self-consistent:


.. parsed-literal::

   variable foo python myMultiply
   python myMultiply return v_foo format f file funcs.py

The two commands can appear in either order in the input script so
long as both are specified before the Python function is invoked for
the first time.

The *format* keyword must be used if the *input* or *return* keyword
is used.  It defines an *fstring* with M characters, where M = sum of
number of inputs and outputs.  The order of characters corresponds to
the N inputs, followed by the return value (if it exists).  Each
character must be one of the following: "i" for integer, "f" for
floating point, "s" for string, or "p" for SELF.  Each character
defines the type of the corresponding input or output value of the
Python function and affects the type conversion that is performed
internally as data is passed back and forth between LAMMPS and Python.
Note that it is permissible to use a :doc:`python-style variable <variable>` in a LAMMPS command that allows for an
equal-style variable as an argument, but only if the output of the
Python function is flagged as a numeric value ("i" or "f") via the
*format* keyword.

If the *return* keyword is used and the *format* keyword specifies the
output as a string, then the default maximum length of that string is
63 characters (64-1 for the string terminator).  If you want to return
a longer string, the *length* keyword can be specified with its *Nlen*
value set to a larger number (the code allocates space for Nlen+1 to
include the string terminator).  If the Python function generates a
string longer than the default 63 or the specified *Nlen*\ , it will be
truncated.


----------


Either the *file*\ , *here*\ , or *exists* keyword must be used, but only
one of them.  These keywords specify what Python code to load into the
Python interpreter.  The *file* keyword gives the name of a file,
which should end with a ".py" suffix, which contains Python code.  The
code will be immediately loaded into and run in the "main" module of
the Python interpreter.  Note that Python code which contains a
function definition does not "execute" the function when it is run; it
simply defines the function so that it can be invoked later.

The *here* keyword does the same thing, except that the Python code
follows as a single argument to the *here* keyword.  This can be done
using triple quotes as delimiters, as in the examples above.  This
allows Python code to be listed verbatim in your input script, with
proper indentation, blank lines, and comments, as desired.  See the
:doc:`Commands parse <Commands_parse>` doc page, for an explanation of
how triple quotes can be used as part of input script syntax.

The *exists* keyword takes no argument.  It means that Python code
containing the required Python function defined by the *func* setting,
is assumed to have been previously loaded by another python command.

Note that the Python code that is loaded and run must contain a
function with the specified *func* name.  To operate properly when
later invoked, the function code must match the *input* and
*return* and *format* keywords specified by the python command.
Otherwise Python will generate an error.


----------


This section describes how Python code can be written to work with
LAMMPS.

Whether you load Python code from a file or directly from your input
script, via the *file* and *here* keywords, the code can be identical.
It must be indented properly as Python requires.  It can contain
comments or blank lines.  If the code is in your input script, it
cannot however contain triple-quoted Python strings, since that will
conflict with the triple-quote parsing that the LAMMPS input script
performs.

All the Python code you specify via one or more python commands is
loaded into the Python "main" module, i.e. \__main\__.  The code can
define global variables or statements that are outside of function
definitions.  It can contain multiple functions, only one of which
matches the *func* setting in the python command.  This means you can
use the *file* keyword once to load several functions, and the
*exists* keyword thereafter in subsequent python commands to access
the other functions previously loaded.

A Python function you define (or more generally, the code you load)
can import other Python modules or classes, it can make calls to other
system functions or functions you define, and it can access or modify
global variables (in the "main" module) which will persist between
successive function calls.  The latter can be useful, for example, to
prevent a function from being invoke multiple times per timestep by
different commands in a LAMMPS input script that access the returned
python-style variable associated with the function.  For example,
consider this function loaded with two global variables defined
outside the function:


.. parsed-literal::

   nsteplast = -1
   nvaluelast = 0

   def expensive(nstep):
     global nsteplast,nvaluelast
     if nstep == nsteplast: return nvaluelast
     nsteplast = nstep
     # perform complicated calculation
     nvalue = ...
     nvaluelast = nvalue
     return nvalue

Nsteplast stores the previous timestep the function was invoked
(passed as an argument to the function).  Nvaluelast stores the return
value computed on the last function invocation.  If the function is
invoked again on the same timestep, the previous value is simply
returned, without re-computing it.  The "global" statement inside the
Python function allows it to overwrite the global variables.

Note that if you load Python code multiple times (via multiple python
commands), you can overwrite previously loaded variables and functions
if you are not careful.  E.g. if the code above were loaded twice, the
global variables would be re-initialized, which might not be what you
want.  Likewise, if a function with the same name exists in two chunks
of Python code you load, the function loaded second will override the
function loaded first.

It's important to realize that if you are running LAMMPS in parallel,
each MPI task will load the Python interpreter and execute a local
copy of the Python function(s) you define.  There is no connection
between the Python interpreters running on different processors.
This implies three important things.

First, if you put a print statement in your Python function, you will
see P copies of the output, when running on P processors.  If the
prints occur at (nearly) the same time, the P copies of the output may
be mixed together.  Welcome to the world of parallel programming and
debugging.

Second, if your Python code loads modules that are not pre-loaded by
the Python library, then it will load the module from disk.  This may
be a bottleneck if 1000s of processors try to load a module at the
same time.  On some large supercomputers, loading of modules from disk
by Python may be disabled.  In this case you would need to pre-build a
Python library that has the required modules pre-loaded and link
LAMMPS with that library.

Third, if your Python code calls back to LAMMPS (discussed in the
next section) and causes LAMMPS to perform an MPI operation requires
global communication (e.g. via MPI\_Allreduce), such as computing the
global temperature of the system, then you must insure all your Python
functions (running independently on different processors) call back to
LAMMPS.  Otherwise the code may hang.


----------


Your Python function can "call back" to LAMMPS through its
library interface, if you use the SELF input to pass Python
a pointer to LAMMPS.  The mechanism for doing this in your
Python function is as follows:


.. parsed-literal::

   def foo(lmpptr,...):
     from lammps import lammps
     lmp = lammps(ptr=lmpptr)
     lmp.command('print "Hello from inside Python"')
     ...

The function definition must include a variable (lmpptr in this case)
which corresponds to SELF in the python command.  The first line of
the function imports the Python module lammps.py in the python dir of
the distribution.  The second line creates a Python object "lmp" which
wraps the instance of LAMMPS that called the function.  The
"ptr=lmpptr" argument is what makes that happen.  The third line
invokes the command() function in the LAMMPS library interface.  It
takes a single string argument which is a LAMMPS input script command
for LAMMPS to execute, the same as if it appeared in your input
script.  In this case, LAMMPS should output


.. parsed-literal::

   Hello from inside Python

to the screen and log file.  Note that since the LAMMPS print command
itself takes a string in quotes as its argument, the Python string
must be delimited with a different style of quotes.

The :doc:`Pytnon library <Python_library>` doc page describes the syntax
for how Python wraps the various functions included in the LAMMPS
library interface.

A more interesting example is in the examples/python/in.python script
which loads and runs the following function from examples/python/funcs.py:


.. parsed-literal::

   def loop(N,cut0,thresh,lmpptr):
     print "LOOP ARGS",N,cut0,thresh,lmpptr
     from lammps import lammps
     lmp = lammps(ptr=lmpptr)
     natoms = lmp.get_natoms()

     for i in range(N):
       cut = cut0 + i\*0.1

       lmp.set_variable("cut",cut)                 # set a variable in LAMMPS
       lmp.command("pair_style lj/cut ${cut}")     # LAMMPS command
       #lmp.command("pair_style lj/cut %d" % cut)  # LAMMPS command option

       lmp.command("pair_coeff \* \* 1.0 1.0")       # ditto
       lmp.command("run 10")                       # ditto
       pe = lmp.extract_compute("thermo_pe",0,0)   # extract total PE from LAMMPS
       print "PE",pe/natoms,thresh
       if pe/natoms < thresh: return

with these input script commands:


.. parsed-literal::

   python          loop input 4 10 1.0 -4.0 SELF format iffp file funcs.py
   python          loop invoke

This has the effect of looping over a series of 10 short runs (10
timesteps each) where the pair style cutoff is increased from a value
of 1.0 in distance units, in increments of 0.1.  The looping stops
when the per-atom potential energy falls below a threshold of -4.0 in
energy units.  More generally, Python can be used to implement a loop
with complex logic, much more so than can be created using the LAMMPS
:doc:`jump <jump>` and :doc:`if <if>` commands.

Several LAMMPS library functions are called from the loop function.
Get\_natoms() returns the number of atoms in the simulation, so that it
can be used to normalize the potential energy that is returned by
extract\_compute() for the "thermo\_pe" compute that is defined by
default for LAMMPS thermodynamic output.  Set\_variable() sets the
value of a string variable defined in LAMMPS.  This library function
is a useful way for a Python function to return multiple values to
LAMMPS, more than the single value that can be passed back via a
return statement.  This cutoff value in the "cut" variable is then
substituted (by LAMMPS) in the pair\_style command that is executed
next.  Alternatively, the "LAMMPS command option" line could be used
in place of the 2 preceding lines, to have Python insert the value
into the LAMMPS command string.

.. note::

   When using the callback mechanism just described, recognize that
   there are some operations you should not attempt because LAMMPS cannot
   execute them correctly.  If the Python function is invoked between
   runs in the LAMMPS input script, then it should be OK to invoke any
   LAMMPS input script command via the library interface command() or
   file() functions, so long as the command would work if it were
   executed in the LAMMPS input script directly at the same point.

However, a Python function can also be invoked during a run, whenever
an associated LAMMPS variable it is assigned to is evaluated.  If the
variable is an input argument to another LAMMPS command (e.g. :doc:`fix setforce <fix_setforce>`), then the Python function will be invoked
inside the class for that command, in one of its methods that is
invoked in the middle of a timestep.  You cannot execute arbitrary
input script commands from the Python function (again, via the
command() or file() functions) at that point in the run and expect it
to work.  Other library functions such as those that invoke computes
or other variables may have hidden side effects as well.  In these
cases, LAMMPS has no simple way to check that something illogical is
being attempted.

The same applies to Python functions called during a simulation run at
each time step using :doc:`fix python/invoke <fix_python_invoke>`.


----------


If you run Python code directly on your workstation, either
interactively or by using Python to launch a Python script stored in a
file, and your code has an error, you will typically see informative
error messages.  That is not the case when you run Python code from
LAMMPS using an embedded Python interpreter.  The code will typically
fail silently.  LAMMPS will catch some errors but cannot tell you
where in the Python code the problem occurred.  For example, if the
Python code cannot be loaded and run because it has syntax or other
logic errors, you may get an error from Python pointing to the
offending line, or you may get one of these generic errors from
LAMMPS:


.. parsed-literal::

   Could not process Python file
   Could not process Python string

When the Python function is invoked, if it does not return properly,
you will typically get this generic error from LAMMPS:


.. parsed-literal::

   Python function evaluation failed

Here are three suggestions for debugging your Python code while
running it under LAMMPS.

First, don't run it under LAMMPS, at least to start with!  Debug it
using plain Python.  Load and invoke your function, pass it arguments,
check return values, etc.

Second, add Python print statements to the function to check how far
it gets and intermediate values it calculates.  See the discussion
above about printing from Python when running in parallel.

Third, use Python exception handling.  For example, say this statement
in your Python function is failing, because you have not initialized the
variable foo:


.. parsed-literal::

   foo += 1

If you put one (or more) statements inside a "try" statement,
like this:


.. parsed-literal::

   import exceptions
   print "Inside simple function"
   try:
     foo += 1      # one or more statements here
   except Exception, e:
     print "FOO error:",e

then you will get this message printed to the screen:


.. parsed-literal::

   FOO error: local variable 'foo' referenced before assignment

If there is no error in the try statements, then nothing is printed.
Either way the function continues on (unless you put a return or
sys.exit() in the except clause).


----------


Restrictions
""""""""""""


This command is part of the PYTHON package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Building LAMMPS with the PYTHON package will link LAMMPS with the
Python library on your system.  Settings to enable this are in the
lib/python/Makefile.lammps file.  See the lib/python/README file for
information on those settings.

If you use Python code which calls back to LAMMPS, via the SELF input
argument explained above, there is an extra step required when
building LAMMPS.  LAMMPS must also be built as a shared library and
your Python function must be able to load the Python module in
python/lammps.py that wraps the LAMMPS library interface.  These are
the same steps required to use Python by itself to wrap LAMMPS.
Details on these steps are explained on the :doc:`Python <Python_head>`
doc page.  Note that it is important that the stand-alone LAMMPS
executable and the LAMMPS shared library be consistent (built from the
same source code files) in order for this to work.  If the two have
been built at different times using different source files, problems
may occur.

Related commands
""""""""""""""""

:doc:`shell <shell>`, :doc:`variable <variable>`, :doc:`fix python/invoke <fix_python_invoke>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
