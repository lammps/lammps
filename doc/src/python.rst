.. index:: python

python command
==============

Syntax
""""""

.. code-block:: LAMMPS

   python mode keyword args ...

* mode = *source* or *name* of Python function

  if mode is *source*:

  .. parsed-literal::

     keyword = *here* or name of a *Python file*
       *here* arg = inline
         inline = one or more lines of Python code which will be executed immediately
                  must be a single argument, typically enclosed between triple quotes
       *Python file* = name of a file with Python code which will be executed immediately

* if *mode* is *name* of a Python function:

  .. parsed-literal::

     one or more keywords with/without arguments must be appended
     keyword = *invoke* or *input* or *return* or *format* or *length* or *file* or *here* or *exists*
       *invoke* arg = logreturn (optional)
         invoke the previously-defined Python function
         if logreturn is specified, print the return value of the invoked function to the screen and logfile
       *input* args = N i1 i2 ... iN
         N = # of inputs to function
         i1,...,iN = value, SELF, or LAMMPS variable name
           value = integer number, floating point number, or string
           SELF = reference to LAMMPS itself which can then be accessed by Python function
           variable = v_name, where name = name of a LAMMPS variable, e.g. v_abc
           internal variable = iv_name, where name = name of a LAMMPS internal-style variable, e.g. iv_xyz
       *return* arg = varReturn
         varReturn = v_name  = LAMMPS variable name which the return value of the Python function will be assigned to
       *format* arg = fstring with M characters
         M = N if no return value, where N = # of inputs
         M = N+1 if there is a return value
         fstring = each character (i,f,s,p) corresponds (in order) to an input or return value
           'i' = integer, 'f' = floating point, 's' = string, 'p' = SELF
       *length* arg = Nlen
         Nlen = max length of string returned from Python function
       *file* arg = filename
         filename = file of Python code, which defines the Python function
       *here* arg = inline
         inline = one or more lines of Python code which defines the Python function
                  must be a single argument, typically enclosed between triple quotes
       *exists* arg = none = Python code has been loaded by previous python command

Examples
""""""""

.. code-block:: LAMMPS

   python pForce input 2 v_x 20.0 return v_f format fff file force.py
   python pForce invoke logreturn

   python factorial input 1 myN return v_fac format ii here """
   def factorial(n):
     if n == 1: return n
     return n * factorial(n-1)
    """

   python loop input 1 SELF return v_value format pf here """
   def loop(lmpptr,N,cut0):
     from lammps import lammps
     lmp = lammps(ptr=lmpptr)

     # loop N times, increasing cutoff each time

     for i in range(N):
       cut = cut0 + i*0.1
       lmp.set_variable("cut",cut)               # set a variable in LAMMPS
       lmp.command("pair_style lj/cut ${cut}")   # LAMMPS commands
       lmp.command("pair_coeff * * 1.0 1.0")
       lmp.command("run 100")
   """

   python source funcdef.py

   python source here "from lammps import lammps"


Description
"""""""""""

The *python* command interfaces LAMMPS with an embedded Python
interpreter and enables executing arbitrary python code in that
interpreter.  This can be done immediately, by using *mode* = *source*.
Or execution can be deferred, by registering a Python function for later
execution, by using *mode* = *name* of a Python function.

Later execution can be triggered in one of two ways.  One is to use the
python command again with its *invoke* keyword.  The other is to trigger
the evaluation of a python-style, equal-style, vector-style, or
atom-style variable.  A python-style variable invokes its associated
Python function; its return value becomes the value of the python-style
variable.  Equal-, vector-, and atom-style variables can use a Python
function wrapper in their formulas which encodes the python-style
variable name, and specifies arguments (which themselves can be numeric
formulas) to pass to the Python function associated with the
python-style variable.

As explained on the :doc:`variable <variable>` doc page, the definition
of a python-style variable associates a Python function name with the
variable.  Its specification must match the *mode* argument of the
*python* command for the Python function name.  For example these two
commands would be consistent:

.. code-block:: LAMMPS

   variable foo python myMultiply
   python myMultiply return v_foo format f file funcs.py

The two commands can appear in either order in the input script so long
as both are specified before the Python function is invoked for the
first time.

Note that python-style, equal-style, vector-style, and atom-style
variables can be used in many different ways within LAMMPS.  They can be
evaluated directly in an input script, effectively replacing the
variable with its value.  Or they can be passed to various commands as
arguments, so that the variable is evaluated multiple times during a
simulation run.  See the :doc:`variable <variable>` command doc page for
more details on variable styles which enable Python function evaluation.

The Python code for a Python function can be included directly in the
input script or in a separate Python file.  The function can be standard
Python code or it can make "callbacks" to LAMMPS through its library
interface to query or set internal values within LAMMPS.  This is a
powerful mechanism for performing complex operations in a LAMMPS input
script that are not possible with the simple input script and variable
syntax which LAMMPS defines.  Thus your input script can operate more
like a true programming language.

Use of this command requires building LAMMPS with the PYTHON package
which links to the Python library so that the Python interpreter is
embedded in LAMMPS.  More details about this process are given below.

A broader overview of how Python can be used with LAMMPS is given in the
:doc:`Use Python with LAMMPS <Python_head>` section of the
documentation.  There is also an ``examples/python`` directory which
illustrates use of the python command.

----------

The first argument to the *python* command is the *mode* setting, which
is either *source* or the *name* of a Python function.

.. versionchanged:: 22Dec2022

If *source* is used, it is followed by either the *here* keyword or a
file name containing Python code.  The *here* keyword is followed by a
single *inline* argument which is a string containing one or more python
commands.  The string can either be on the same line as the *python*
command, enclosed in quotes, or it can be multiple lines enclosed in
triple quotes.

In either case, the in-line code or the file contents are passed to the
python interpreter and executed immediately.  The code will be loaded
into and run in the "main" module of the Python interpreter.  This
allows running arbitrary Python code at any time while processing the
LAMMPS input file.  This can be used to pre-load Python modules,
initialize global variables, define functions or classes, or perform
operations using the Python programming language.  The Python code will
be executed in parallel on all the MPI processes being used to run
LAMMPS.  Note that no arguments can be passed to the executed Python
code.

If the *mode* setting is the *name* of a Python function, then it will
be registered with LAMMPS for future execution (or can already be
defined, see the *exists* keyword).  One or more keywords must follow
the *mode* function name.  One of the keywords must be *invoke*, *file*,
*here*, or *exists*, which specifies what Python code to load into the
Python interpreter.  Note that only one of those 4 keywords is allowed
since their operations are mutually exclusive.

----------

If the *invoke* keyword is used, no other keywords can be used.  A
previous *python* command must have registered the Python function
referenced by this command, which can then be invoked multiple times in
an input script via the *invoke* keyword.  Each invocation passes
current values for arguments to the Python function.  A return value of
the Python function will be ignored unless the Python function is linked
to a :doc:`python style variable <variable>` with the *return* keyword.
This return value can be logged to the screen and logfile by adding the
optional *logreturn* argument to the *invoke* keyword.  In that case a
message with the name of the python command and the return value is
printed.  Note that return values of python functions are otherwise
*only* accessible when the function is invoked indirectly by evaluating
its associated :doc:`python style variable <variable>`, as described
below.

The *file* keyword gives the name of a file containing Python code,
which should end with a ".py" suffix.  The code will be immediately
loaded into and run in the "main" module of the Python interpreter.  The
Python code will be executed in parallel on all MPI processes.  Note
that Python code which contains a function definition does NOT "execute"
the function when it is run; it simply defines the function so that it
can be invoked later.

The *here* keyword does the same thing, except that the Python code
follows as a single argument to the *here* keyword.  This can be done
using triple quotes as delimiters, as in the examples above and below.
This allows Python code to be listed verbatim in your input script, with
proper indentation, blank lines, and comments, as desired.  See the
:doc:`Commands parse <Commands_parse>` doc page, for an explanation of
how triple quotes can be used as part of input script syntax.

The *exists* keyword takes no argument.  It simply means that Python
code containing the needed Python function has already been loaded into
the LAMMPS Python interpreter, for example by previous *python source*
command or in a file that was loaded previously with the *file*
keyword. This allows use of a single file of Python code which contains
multiple functions, any of which can be used in the same (or different)
input scripts (see below).

Note that the Python code that is loaded and run by the *file* or *here*
keyword must contain a function with the specified function *name*.  To
operate properly when the function is later invoked, the code for the
function must match the *input* and *return* and *format* keywords
specified by the python command.  Otherwise Python will generate an
error.

----------

The other keywords which can be used with the *python* command are
*input*, *return*, *format*, and *length*.

The *input* keyword defines how many arguments *N* the Python function
expects.  If it takes no arguments, then the *input* keyword should not
be used.  Each argument can be specified directly as a value, e.g. '6'
or '3.14159' or 'abc' (a string of characters).  The type of each
argument is specified by the *format* keyword as explained below, so
that Python will know how to interpret the value.  If the word SELF is
used for an argument it has a special meaning.  A pointer is passed to
the Python function which it can convert into a reference to LAMMPS
itself using the :doc:`LAMMPS Python module <Python_module>`.  This
enables the function to call back to LAMMPS through its library
interface as explained below.  This allows the Python function to query
or set values internal to LAMMPS which can affect the subsequent
execution of the input script.

A LAMMPS variable can also be used as an *input* argument, specified as
v_name, where "name" is the name of the variable defined in the input
script.  Any style of LAMMPS variable returning a scalar or a string can
be used, as defined by the :doc:`variable <variable>` command.  The
style of variable must be consistent with the *format* keyword
specification for the type of data that is passed to Python.  Each time
the Python function is invoked, the LAMMPS variable is evaluated and its
value is passed as an argument to the Python function.  Note that a
python-style variable can be used as an argument, which means that the a
Python function can use arguments which invoke other Python functions.

A LAMMPS internal-style variable can also be used as an *input*
argument, specified as iv_name, where "name" is the name of the
internal-style variable.  The internal-style variable does not have to
be defined in the input script (though it can be); if it is not defined,
this command creates an :doc:`internal-style variable <variable>` with
the specified name.

An internal-style variable must be used when an equal-style,
vector-style, or atom-style variable triggers the invocation of the
Python function defined by this command, by including a Python function
wrapper with arguments in its formula.  Each of the arguments must be
specified as an internal-style variable via the *input* keyword.

In brief, the syntax for a Python function wrapper in a variable formula
is ``py_varname(arg1,arg2,...argN)``, where "varname" is the name of a
python-style variable associated with a Python function defined by this
command.  One or more arguments to the function wrapper can themselves
be sub-formulas which the variable command will evaluate and pass as
arguments to the Python function.  This is done by assigning the numeric
result for each argument to an internal-style variable; thus the *input*
keyword must specify the arguments as internal-style variables and their
format (see below) as "f" for floating point.  This is because LAMMPS
variable formulas are calculated with floating point arithmetic (any
integer values are converted to floating point).  Note that the Python
function can also have additional inputs, also specified by the *input*
keyword, which are NOT arguments in the Python function wrapper.  See
the example below for the ``mixedargs`` Python function.

See the :doc:`variable <variable>` command doc page for full details on
formula syntax including for Python function wrappers.  Examples using
Python function wrappers are shown below.  Note that as explained above
with python-style variables, Python function wrappers can be nested; a
sub-formula for an argument can contain its own Python function wrapper
which invokes another Python function.

The *return* keyword is only needed if the Python function returns a
value.  The specified *varReturn* is of the form v_name, where "name" is
the name of a python-style LAMMPS variable, defined by the
:doc:`variable <variable>` command.  The Python function can return a
numeric or string value, as specified by the *format* keyword.  This
return value is *only* accessible when its associated python-style
variable is evaluated.  When the *invoke* keyword is used, the return
value of the python function is ignored unless the optional *logreturn*
argument is specified.

The *format* keyword must be used if the *input* or *return* keywords
are used.  It defines an *fstring* with M characters, where M = sum of
number of inputs and outputs.  The order of characters corresponds to
the N inputs, followed by the return value (if it exists).  Each
character must be one of the following: "i" for integer, "f" for
floating point, "s" for string, or "p" for SELF.  Each character defines
the type of the corresponding input or output value of the Python
function and affects the type conversion that is performed internally as
data is passed back and forth between LAMMPS and Python.  Note that it
is permissible to use a :doc:`python-style variable <variable>` in a
LAMMPS command that allows for an equal-style variable as an argument,
but only if the output of the Python function is flagged as a numeric
value ("i" or "f") via the *format* keyword.

If the *return* keyword is used and the *format* keyword specifies the
output as a string, then the default maximum length of that string is 63
characters (64-1 for the string terminator).  If you want to return a
longer string, the *length* keyword can be specified with its *Nlen*
value set to a larger number.  LAMMPS will then allocate Nlen+1 space to
include the string terminator.  If the Python function generates a
string longer than the default 63 or the specified *Nlen*, it will be
truncated.

----------

This section describes how Python code can be written to work with
LAMMPS.

Whether you load Python code from a file or directly from your input
script, via the *file* and *here* keywords, the code can be identical.
It must be indented properly as Python requires.  It can contain
comments or blank lines.  If the code is in your input script, it cannot
however contain triple-quoted Python strings, since that will conflict
with the triple-quote parsing that the LAMMPS input script performs.

All the Python code you specify via one or more python commands is
loaded into the Python "main" module, i.e. ``__name__ == '__main__'``.
The code can define global variables, define global functions, define
classes or execute statements that are outside of function definitions.
It can contain multiple functions, only one of which matches the *func*
setting in the python command.  This means you can use the *file*
keyword once to load several functions, and the *exists* keyword
thereafter in subsequent python commands to register the other functions
that were previously loaded with LAMMPS.

A Python function you define (or more generally, the code you load) can
import other Python modules or classes, it can make calls to other
system functions or functions you define, and it can access or modify
global variables (in the "main" module) which will persist between
successive function calls.  The latter can be useful, for example, to
prevent a function from being invoke multiple times per timestep by
different commands in a LAMMPS input script that access the returned
python-style variable associated with the function.  For example,
consider this function loaded with two global variables defined outside
the function:

.. code-block:: python

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

The variable 'nsteplast' stores the previous timestep the function was
invoked (passed as an argument to the function).  The variable
'nvaluelast' stores the return value computed on the last function
invocation.  If the function is invoked again on the same timestep, the
previous value is simply returned, without re-computing it.  The
"global" statement inside the Python function allows it to overwrite the
global variables from within the local context of the function.

Also note that if you load Python code multiple times (via multiple
python commands), you can overwrite previously loaded variables and
functions if you are not careful.  E.g. if the code above were loaded
twice, the global variables would be re-initialized, which might not be
what you want.  Likewise, if a function with the same name exists in two
chunks of Python code you load, the function loaded second will override
the function loaded first.

It's important to realize that if you are running LAMMPS in parallel,
each MPI task will load the Python interpreter and execute a local copy
of the Python function(s) you define.  There is no connection between
the Python interpreters running on different processors.  This implies
three important things.

First, if you put a print or other statement creating output to the
screen in your Python function, you will see P copies of the output,
when running on P processors.  If the prints occur at (nearly) the same
time, the P copies of the output may be mixed together.

It is possible to avoid this issue, by passing the pointer to the
current LAMMPS class instance to the Python function via the {input}
SELF argument described above.  The Python function can then use the
Python interface to the LAMMPS library interface, and determine the MPI
rank of the current process.  The Python code can then ensure output
will only appear on MPI rank 0.  The following LAMMPS input demonstrates
how this could be done. The text 'Hello, LAMPS!' should be printed only
once, even when running LAMMPS in parallel.

.. code-block:: LAMMPS

   python python_hello input 1 SELF format p here """
   def python_hello(handle):
       from lammps import lammps
       lmp = lammps(ptr=handle)
       me = lmp.extract_setting('world_rank')
       if me == 0:
           print('Hello, LAMMPS!')
   """

   python python_hello invoke

Second, if your Python code loads Python modules that are not pre-loaded
by the Python library, then it will load the module from disk.  This may
be a bottleneck if 1000s of processors try to load a module at the same
time.  On some large supercomputers, loading of modules from disk by
Python may be disabled.  In this case you would need to pre-build a
Python library that has the required modules pre-loaded and link LAMMPS
with that library.

Third, if your Python code calls back to LAMMPS (discussed in the next
section) and causes LAMMPS to perform an MPI operation requires global
communication (e.g. via MPI_Allreduce), such as computing the global
temperature of the system, then you must ensure all your Python
functions (running independently on different processors) call back to
LAMMPS.  Otherwise the code may hang.

----------

As mentioned above, a Python function can "call back" to LAMMPS through
its library interface, if the SELF input is used to pass Python a
pointer to LAMMPS.  The mechanism for doing this is as follows:

.. code-block:: python

   def foo(handle,...):
     from lammps import lammps
     lmp = lammps(ptr=handle)
     lmp.command('print "Hello from inside Python"')
     ...

The function definition must include a variable ('handle' in this case)
which corresponds to SELF in the *python* command.  The first line of
the function imports the :doc:`"lammps" Python module <Python_module>`.
The second line creates a Python object ``lmp`` which wraps the instance
of LAMMPS that called the function.  The 'ptr=handle' argument is what
makes that happen.  The third line invokes the command() function in the
LAMMPS library interface.  It takes a single string argument which is a
LAMMPS input script command for LAMMPS to execute, the same as if it
appeared in your input script.  In this case, LAMMPS should output

.. parsed-literal::

   Hello from inside Python

to the screen and log file.  Note that since the LAMMPS print command
itself takes a string in quotes as its argument, the Python string must
be delimited with a different style of quotes.

The :doc:`Python_head` page describes the syntax for how Python wraps
the various functions included in the LAMMPS library interface.

A more interesting example is in the ``examples/python/in.python``
script which loads and runs the following function from
``examples/python/funcs.py``:

.. code-block:: python

   def loop(N,cut0,thresh,lmpptr):
     print("LOOP ARGS", N, cut0, thresh, lmpptr)
     from lammps import lammps
     lmp = lammps(ptr=lmpptr)
     natoms = lmp.get_natoms()

     for i in range(N):
       cut = cut0 + i*0.1

       lmp.set_variable("cut",cut)                 # set a variable in LAMMPS
       lmp.command("pair_style lj/cut ${cut}")     # LAMMPS command
       #lmp.command("pair_style lj/cut %d" % cut)  # alternate form of LAMMPS command

       lmp.command("pair_coeff * * 1.0 1.0")       # ditto
       lmp.command("run 10")                       # ditto
       pe = lmp.extract_compute("thermo_pe",0,0)   # extract total PE from LAMMPS
       print("PE", pe/natoms, thresh)
       if pe/natoms < thresh: return

with these input script commands:

.. code-block:: LAMMPS

   python          loop input 4 10 1.0 -4.0 SELF format iffp file funcs.py
   python          loop invoke

This has the effect of looping over a series of 10 short runs (10
timesteps each) where the pair style cutoff is increased from a value of
1.0 in distance units, in increments of 0.1.  The looping stops when the
per-atom potential energy falls below a threshold of -4.0 in energy
units.  More generally, Python can be used to implement a loop with
complex logic, much more so than can be created using the LAMMPS
:doc:`jump <jump>` and :doc:`if <if>` commands.

Several LAMMPS library functions are called from the loop function.
Get_natoms() returns the number of atoms in the simulation, so that it
can be used to normalize the potential energy that is returned by
extract_compute() for the "thermo_pe" compute that is defined by default
for LAMMPS thermodynamic output.  Set_variable() sets the value of a
string variable defined in LAMMPS.  This library function is a useful
way for a Python function to return multiple values to LAMMPS, more than
the single value that can be passed back via a return statement.  This
cutoff value in the "cut" variable is then substituted (by LAMMPS) in
the pair_style command that is executed next.  Alternatively, the
"alternate form of LAMMPS command" line could be used in place of the 2
preceding lines, to have Python insert the value into the LAMMPS command
string.

.. note::

   When using the callback mechanism just described, recognize that
   there are some operations you should not attempt because LAMMPS
   cannot execute them correctly.  If the Python function is invoked
   between runs in the LAMMPS input script, then it should be OK to
   invoke any LAMMPS input script command via the library interface
   command() or file() functions, so long as the command would work if
   it were executed in the LAMMPS input script directly at the same
   point.


----------

As noted above, a Python function can be invoked during a run, whenever
an associated python-style variable it is assigned to is evaluated.

If the variable is an input argument to another LAMMPS command
(e.g. :doc:`fix setforce <fix_setforce>`), then the Python function will
be invoked inside the class for that command, possibly in one of its
methods that is invoked in the middle of a timestep.  You cannot execute
arbitrary input script commands from the Python function (again, via the
command() or file() functions) at that point in the run and expect it to
work.  Other library functions such as those that invoke computes or
other variables may have hidden side effects as well.  In these cases,
LAMMPS has no simple way to check that something illogical is being
attempted.

The same constraints apply to Python functions called during a
simulation run at each time step using the :doc:`fix python/invoke
<fix_python_invoke>` command.

----------

As noted above, a Python function can also be invoked within the formula
for an equal-style, vector-style, or atom-style variable.  This means
the Python function will be invoked whenever that variable is invoked.
In the case of a vector-style variable, the Python function can be
invoked once per element of the global vector.  In the case of an
atom-style variable, the Python function can be invoked once per atom.

Here are three simple examples using equal-, vector-, and atom-style
variables to trigger execution of a Python function:

.. code-block:: LAMMPS

   variable        foo python truncate
   python          truncate return v_foo input 1 iv_arg format fi here """
   def truncate(x):
     return int(x)
   """
   variable        ptrunc equal py_foo(press)
   print           "TRUNCATED pressure = ${ptrunc}"

The Python ``truncate`` function simply converts a floating-point value
to an integer value.  When the LAMMPS print command evaluates the
equal-style ``ptrunc`` variable, the current thermodynamic pressure is
passed to the Python function.  The truncated value is output to the
screen and logfile by the print command.  Note that the *input* keyword
for the *python* command, specifies an internal-style variable named
"arg" as iv_arg which is required to invoke the Python function from a
Python function wrapper.

The last 2 lines can be replaced by these to define a vector-style
variable which invokes the same Python ``truncate`` function:

.. code-block:: LAMMPS

  compute         ke all temp
  variable        ke vector c_ke
  variable        ketrunc vector py_foo(v_ke)
  thermo_style    custom step temp epair v_ketrunc[*6]

The vector-style variable ``ketrunc`` invokes the Python ``truncate``
function on each of the 6 components of the global kinetic energy tensor
calculated by the :doc:`compute ke <compute_ke>` command.  The 6
truncated values will be printed with thermo output to the screen and
log file.

Or the last 2 lines of the equal-style variable example can be replaced
by these to define atom-style variables which invoke the same Python
``truncate`` function:

.. code-block:: LAMMPS

   variable        xtrunc atom py_foo(x)
   variable        ytrunc atom py_foo(y)
   variable        ztrunc atom py_foo(z)
   dump            1 all custom 100 tmp.dump id x y z v_xtrunc v_ytrunc v_ztrunc

When the dump command invokes the 3 atom-style variables, their
arguments x,y,z to the Python function wrapper are the current per-atom
coordinates of each atom.  The Python ``truncate`` function is thus
invoked 3 times for each atom, and the truncated coordinate values for
each atom are written to the dump file.

Note that when using a Python function wrapper in a variable, arguments
can be passed to the Python function either from the variable formula or
by *input* keyword to the :doc:`python command <python>`.  For example,
consider these (made up) commands:

.. code-block:: LAMMPS

   variable        foo python mixedargs
   python          mixedargs return v_foo input 6 7.5 v_myValue iv_arg1 iv_argy iv_argz v_flag &
                   format fffffsf here """
   def mixedargs(a,b,x,y,z,flag):
     ...
     return result
   """
   variable        flag string optionABC
   variable        myValue equal "2.0*temp*c_pe"
   compute         pe all pe
   compute         peatom all pe/atom
   variable        field atom py_foo(x+3.0,sqrt(y),(z-zlo)*c_peatom)

They define a Python ``mixedargs`` function with 6 arguments.  Three of
them are internal-style variables, which the variable formula calculates
as numeric values for each atom and passes to the function.  In this
example, these arguments are themselves small formulas containing the
x,y,z coordinates of each atom as well as a per-atom compute (c_peratom)
and thermodynamic keyword (zlo).

The other three arguments ``(7.5,v_myValue,v_flag)`` are defined by the
*python* command.  The first and last are constant values ("7.5" and the
``optionABC`` string).  The second argument (``myValue``) is the result
of an equal-style variable formula which accesses the system temperature
and potential energy.

The "result" returned by each invocation of the Python ``mixedargs``
function becomes the per-atom value in the atom-style "field" variable,
which could be output to a dump file or used elsewhere in the input
script.

----------

If you run Python code directly on your workstation, either
interactively or by using Python to launch a Python script stored in a
file, and your code has an error, you will typically see informative
error messages.  That is not the case when you run Python code from
LAMMPS using an embedded Python interpreter.  The code will typically
fail silently.  LAMMPS will catch some errors but cannot tell you where
in the Python code the problem occurred.  For example, if the Python
code cannot be loaded and run because it has syntax or other logic
errors, you may get an error from Python pointing to the offending line,
or you may get one of these generic errors from LAMMPS:

.. parsed-literal::

   Could not process Python file
   Could not process Python string

When the Python function is invoked, if it does not return properly,
you will typically get this generic error from LAMMPS:

.. parsed-literal::

   Python function evaluation failed

Here are three suggestions for debugging your Python code while running
it under LAMMPS.

First, don't run it under LAMMPS, at least to start with!  Debug it
using plain Python.  Load and invoke your function, pass it arguments,
check return values, etc.

Second, add Python print statements to the function to check how far it
gets and intermediate values it calculates.  See the discussion above
about printing from Python when running in parallel.

Third, use Python exception handling.  For example, say this statement
in your Python function is failing, because you have not initialized the
variable foo:

.. code-block:: python

   foo += 1

If you put one (or more) statements inside a "try" statement, like this:

.. code-block:: python

   import exceptions
   print("Inside simple function")
   try:
     foo += 1      # one or more statements here
   except Exception as e:
     print("FOO error:", e)

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
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Building LAMMPS with the PYTHON package will link LAMMPS with the Python
library on your system.  Settings to enable this are in the
lib/python/Makefile.lammps file.  See the lib/python/README file for
information on those settings.

If you use Python code which calls back to LAMMPS, via the SELF input
argument explained above, there is an extra step required when building
LAMMPS.  LAMMPS must also be built as a shared library and your Python
function must be able to load the :doc:`"lammps" Python module
<Python_module>` that wraps the LAMMPS library interface.

These are the same steps required to use Python by itself to wrap
LAMMPS.  Details on these steps are explained on the :doc:`Python
<Python_head>` doc page.  Note that it is important that the stand-alone
LAMMPS executable and the LAMMPS shared library be consistent (built
from the same source code files) in order for this to work.  If the two
have been built at different times using different source files,
problems may occur.

Another limitation of calling back to Python from the LAMMPS module
using the *python* command in a LAMMPS input is that both, the Python
interpreter and LAMMPS, must be linked to the same Python runtime as a
shared library.  If the Python interpreter is linked to Python
statically (which seems to happen with Conda) then loading the shared
LAMMPS library will create a second python "main" module that hides the
one from the Python interpreter and all previous defined function and
global variables will become invisible.

Related commands
""""""""""""""""

:doc:`shell <shell>`, :doc:`variable <variable>`,
:doc:`fix python/invoke <fix_python_invoke>`

Default
"""""""

none
