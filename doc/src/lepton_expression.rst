
Lepton expression syntax and features
"""""""""""""""""""""""""""""""""""""

Lepton supports the following operators in expressions:

.. table_from_list::
   :columns: 14

   * \+
   * Add
   *
   * \-
   * Subtract
   *
   * \*
   * Multiply
   *
   * \/
   * Divide
   *
   * \^
   * Power

The following mathematical functions are available:

.. table_from_list::
   :columns: 4

   * sqrt(x)
   * Square root
   * exp(x)
   * Exponential
   * log(x)
   * Natural logarithm
   * sin(x)
   * Sine (angle in radians)
   * cos(x)
   * Cosine (angle in radians)
   * sec(x)
   * Secant (angle in radians)
   * csc(x)
   * Cosecant (angle in radians)
   * tan(x)
   * Tangent (angle in radians)
   * cot(x)
   * Cotangent (angle in radians)
   * asin(x)
   * Inverse sine (in radians)
   * acos(x)
   * Inverse cosine (in radians)
   * atan(x)
   * Inverse tangent (in radians)
   * sinh(x)
   * Hyperbolic sine
   * cosh(x)
   * Hyperbolic cosine
   * tanh(x)
   * Hyperbolic tangent
   * erf(x)
   * Error function
   * erfc(x)
   * Complementary Error function
   * abs(x)
   * Absolute value
   * min(x,y)
   * Minimum of two values
   * max(x,y)
   * Maximum of two values
   * delta(x)
   * delta(x) is 1 for `x = 0`, otherwise 0
   * step(x)
   * step(x) is 0 for `x < 0`, otherwise 1

Numbers may be given in either decimal or exponential form.  All of the
following are valid numbers: `5`, `-3.1`, `1e6`, and `3.12e-2`.

As an extension to the standard Lepton syntax, it is also possible to
use LAMMPS :doc:`variables <variable>` in the format "v_name".  Before
evaluating the expression, "v_name" will be replaced with the value of
the variable "name".  This is compatible with all kinds of scalar
variables, but not with vectors, arrays, local, or per-atom
variables.  If necessary, a custom scalar variable needs to be defined
that can access the desired (single) item from a non-scalar variable.
As an example, the following lines will instruct LAMMPS to ramp
the force constant for a harmonic bond from 100.0 to 200.0 during the
next run:

.. code-block:: LAMMPS

   variable fconst equal ramp(100.0, 200)
   bond_style lepton
   bond_coeff 1 1.5 "v_fconst * (r^2)"

An expression may be followed by definitions for intermediate values that appear in the
expression. A semicolon ";" is used as a delimiter between value definitions. For example,
the expression:

.. code-block:: C

   a^2+a*b+b^2; a=a1+a2; b=b1+b2

is exactly equivalent to

.. code-block:: C

   (a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2

The definition of an intermediate value may itself involve other
intermediate values. Whitespace and quotation characters ('\'' and '"')
are ignored.  All uses of a value must appear *before* that value's
definition.  For efficiency reasons, the expression string is parsed,
optimized, and then stored in an internal, pre-parsed representation for
evaluation.

Evaluating a Lepton expression is typically between 2.5 and 5 times
slower than the corresponding compiled and optimized C++ code.  If
additional speed or GPU acceleration (via GPU or KOKKOS) is required,
the interaction can be represented as a table.  Suitable table files
can be created either internally using the :doc:`pair_write <pair_write>`
or :doc:`bond_write <bond_write>` command or through the Python scripts
in the :ref:`tools/tabulate <tabulate>` folder.
