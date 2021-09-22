
Utility functions
-----------------

The ``utils`` sub-namespace inside the ``LAMMPS_NS`` namespace provides
a collection of convenience functions and utilities that perform common
tasks that are required repeatedly throughout the LAMMPS code like
reading or writing to files with error checking or translation of
strings into specific types of numbers with checking for validity.  This
reduces redundant implementations and encourages consistent behavior.

I/O with status check and similar functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The the first two functions are wrappers around the corresponding C
library calls ``fgets()`` or ``fread()``.  They will check if there
were errors on reading or an unexpected end-of-file state was reached.
In that case, the functions will stop with an error message, indicating
the name of the problematic file, if possible unless the *error* argument
is a NULL pointer.

The :cpp:func:`fgets_trunc` function will work similar for ``fgets()``
but it will read in a whole line (i.e. until the end of line or end
of file), but store only as many characters as will fit into the buffer
including a final newline character and the terminating NULL byte.
If the line in the file is longer it will thus be truncated in the buffer.
This function is used by :cpp:func:`read_lines_from_file` to read individual
lines but make certain they follow the size constraints.

The :cpp:func:`read_lines_from_file` function will read the requested
number of lines of a maximum length into a buffer and will return 0
if successful or 1 if not. It also guarantees that all lines are
terminated with a newline character and the entire buffer with a
NULL character.

----------

.. doxygenfunction:: sfgets
   :project: progguide

.. doxygenfunction:: sfread
   :project: progguide

.. doxygenfunction:: fgets_trunc
   :project: progguide

.. doxygenfunction:: read_lines_from_file
   :project: progguide

----------

String to number conversions with validity check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These functions should be used to convert strings to numbers. They are
are strongly preferred over C library calls like ``atoi()`` or
``atof()`` since they check if the **entire** provided string is a valid
(floating-point or integer) number, and will error out instead of
silently returning the result of a partial conversion or zero in cases
where the string is not a valid number.  This behavior allows to more
easily detect typos or issues when processing input files.

The *do_abort* flag should be set to ``true`` in case  this function
is called only on a single MPI rank, as that will then trigger the
a call to ``Error::one()`` for errors instead of ``Error::all()``
and avoids a "hanging" calculation when run in parallel.

Please also see :cpp:func:`is_integer() <LAMMPS_NS::utils::is_integer>`
and :cpp:func:`is_double() <LAMMPS_NS::utils::is_double>` for testing
strings for compliance without conversion.

----------

.. doxygenfunction:: numeric
   :project: progguide

.. doxygenfunction:: inumeric
   :project: progguide

.. doxygenfunction:: bnumeric
   :project: progguide

.. doxygenfunction:: tnumeric
   :project: progguide


String processing
^^^^^^^^^^^^^^^^^

The following are functions to help with processing strings
and parsing files or arguments.

----------

.. doxygenfunction:: strdup
   :project: progguide

.. doxygenfunction:: trim
   :project: progguide

.. doxygenfunction:: trim_comment
   :project: progguide

.. doxygenfunction:: has_utf8
   :project: progguide

.. doxygenfunction:: utf8_subst
   :project: progguide

.. doxygenfunction:: count_words(const char *text)
   :project: progguide

.. doxygenfunction:: count_words(const std::string &text)
   :project: progguide

.. doxygenfunction:: count_words(const std::string &text, const std::string &separators)
   :project: progguide

.. doxygenfunction:: trim_and_count_words
   :project: progguide

.. doxygenfunction:: split_words
   :project: progguide

.. doxygenfunction:: split_lines
   :project: progguide

.. doxygenfunction:: strmatch
   :project: progguide

.. doxygenfunction:: strfind
   :project: progguide

.. doxygenfunction:: is_integer
   :project: progguide

.. doxygenfunction:: is_double
   :project: progguide

File and path functions
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: guesspath
   :project: progguide

.. doxygenfunction:: path_basename
   :project: progguide

.. doxygenfunction:: path_join
   :project: progguide

.. doxygenfunction:: file_is_readable
   :project: progguide

Potential file functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: get_potential_file_path
   :project: progguide

.. doxygenfunction:: get_potential_date
   :project: progguide

.. doxygenfunction:: get_potential_units
   :project: progguide

.. doxygenfunction:: get_supported_conversions
   :project: progguide

.. doxygenfunction:: get_conversion_factor
   :project: progguide

.. doxygenfunction:: open_potential(const std::string &name, LAMMPS *lmp, int *auto_convert)
   :project: progguide

Argument processing
^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: bounds
   :project: progguide

.. doxygenfunction:: expand_args
   :project: progguide

Convenience functions
^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: logmesg(LAMMPS *lmp, const S &format, Args&&... args)
   :project: progguide

.. doxygenfunction:: logmesg(LAMMPS *lmp, const std::string &mesg)
   :project: progguide

.. doxygenfunction:: getsyserror
   :project: progguide

.. doxygenfunction:: check_packages_for_style
   :project: progguide

.. doxygenfunction:: timespec2seconds
   :project: progguide

.. doxygenfunction:: date2num
   :project: progguide

.. doxygenfunction:: current_date
   :project: progguide

Customized standard functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: binary_search
   :project: progguide

.. doxygenfunction:: merge_sort
   :project: progguide

---------------------------

Tokenizer classes
-----------------

The purpose of the tokenizer classes is to simplify the recurring task
of breaking lines of text down into words and/or numbers.
Traditionally, LAMMPS code would be using the ``strtok()`` function from
the C library for that purpose, but that function has two significant
disadvantages: 1) it cannot be used concurrently from different LAMMPS
instances since it stores its status in a global variable and 2) it
modifies the string that it is processing.  These classes were
implemented to avoid both of these issues and also to reduce the amount
of code that needs to be written.

The basic procedure is to create an instance of the tokenizer class with
the string to be processed as an argument and then do a loop until all
available tokens are read.  The constructor has a default set of
separator characters, but that can be overridden. The default separators
are all "whitespace" characters, i.e. the space character, the tabulator
character, the carriage return character, the linefeed character, and
the form feed character.

.. code-block:: C++
   :caption: Tokenizer class example listing entries of the PATH environment variable

   #include "tokenizer.h"
   #include <cstdlib>
   #include <string>
   #include <iostream>

   using namespace LAMMPS_NS;

   int main(int, char **)
   {
       const char *path = getenv("PATH");

       if (path != nullptr) {
           Tokenizer p(path,":");
           while (p.has_next())
               std::cout << "Entry: " << p.next() << "\n";
       }
       return 0;
   }

Most tokenizer operations cannot fail except for
:cpp:func:`LAMMPS_NS::Tokenizer::next` (when used without first
checking with :cpp:func:`LAMMPS_NS::Tokenizer::has_next`) and
:cpp:func:`LAMMPS_NS::Tokenizer::skip`.  In case of failure, the class
will throw an exception, so you may need to wrap the code using the
tokenizer into a ``try`` / ``catch`` block to handle errors.  The
:cpp:class:`LAMMPS_NS::ValueTokenizer` class may also throw an exception
when a (type of) number is requested as next token that is not
compatible with the string representing the next word.

.. code-block:: C++
   :caption: ValueTokenizer class example with exception handling

   #include "tokenizer.h"
   #include <cstdlib>
   #include <string>
   #include <iostream>

   using namespace LAMMPS_NS;

   int main(int, char **)
   {
       const char *text = "1 2 3 4 5 20.0 21 twentytwo 2.3";
       double num1(0),num2(0),num3(0),num4(0);

       ValueTokenizer t(text);
       // read 4 doubles after skipping over 5 numbers
       try {
           t.skip(5);
           num1 = t.next_double();
           num2 = t.next_double();
           num3 = t.next_double();
           num4 = t.next_double();
       } catch (TokenizerException &e) {
           std::cout << "Reading numbers failed: " << e.what() << "\n";
       }
       std::cout << "Values: " << num1 << " " << num2 << " " << num3 << " " << num4 << "\n";
       return 0;
   }

This code example should produce the following output:

.. code-block::

   Reading numbers failed: Not a valid floating-point number: 'twentytwo'
   Values: 20 21 0 0

----------

.. doxygenclass:: LAMMPS_NS::Tokenizer
   :project: progguide
   :members:

.. doxygenclass:: LAMMPS_NS::TokenizerException
   :project: progguide
   :members:

.. doxygenclass:: LAMMPS_NS::ValueTokenizer
   :project: progguide
   :members:

.. doxygenclass:: LAMMPS_NS::InvalidIntegerException
   :project: progguide
   :members: what

.. doxygenclass:: LAMMPS_NS::InvalidFloatException
   :project: progguide
   :members: what

----------


Argument parsing classes
---------------------------

The purpose of argument parsing classes it to simplify and unify how
arguments of commands in LAMMPS are parsed and to make abstractions of
repetitive tasks.

The :cpp:class:`LAMMPS_NS::ArgInfo` class provides an abstraction
for parsing references to compute or fix styles, variables or custom
integer or double properties handled by :doc:`fix property/atom <fix_property_atom>`.
These would start with a "c\_", "f\_", "v\_", "d\_", "d2\_", "i\_", or "i2\_"
followed by the ID or name of than instance and may be postfixed with
one or two array indices "[<number>]" with numbers > 0.

A typical code segment would look like this:

.. code-block:: C++
   :caption: Usage example for ArgInfo class

   int nvalues = 0;
   for (iarg = 0; iarg < nargnew; iarg++) {
     ArgInfo argi(arg[iarg]);

     which[nvalues] = argi.get_type();
     argindex[nvalues] = argi.get_index1();
     ids[nvalues] = argi.copy_name();

     if ((which[nvalues] == ArgInfo::UNKNOWN)
          || (which[nvalues] == ArgInfo::NONE)
          || (argi.get_dim() > 1))
       error->all(FLERR,"Illegal compute XXX command");

     nvalues++;
   }

----------

.. doxygenclass:: LAMMPS_NS::ArgInfo
   :project: progguide
   :members:


----------

File reader classes
-------------------

The purpose of the file reader classes is to simplify the recurring task
of reading and parsing files. They can use the
:cpp:class:`LAMMPS_NS::ValueTokenizer` class to process the read in
text.  The :cpp:class:`LAMMPS_NS::TextFileReader` is a more general
version while :cpp:class:`LAMMPS_NS::PotentialFileReader` is specialized
to implement the behavior expected for looking up and reading/parsing
files with potential parameters in LAMMPS.  The potential file reader
class requires a LAMMPS instance, requires to be run on MPI rank 0 only,
will use the :cpp:func:`LAMMPS_NS::utils::get_potential_file_path`
function to look up and open the file, and will call the
:cpp:class:`LAMMPS_NS::Error` class in case of failures to read or to
convert numbers, so that LAMMPS will be aborted.

.. code-block:: C++
   :caption: Use of PotentialFileReader class in pair style coul/streitz

    PotentialFileReader reader(lmp, file, "coul/streitz");
    char * line;

    while((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);
        std::string iname = values.next_string();

        int ielement;
        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");
        }

        params[nparams].ielement = ielement;
        params[nparams].chi = values.next_double();
        params[nparams].eta = values.next_double();
        params[nparams].gamma = values.next_double();
        params[nparams].zeta = values.next_double();
        params[nparams].zcore = values.next_double();

      } catch (TokenizerException & e) {
        error->one(FLERR, e.what());
      }
      nparams++;
    }

A file that would be parsed by the reader code fragment looks like this:

.. parsed-literal::

   # DATE: 2015-02-19 UNITS: metal CONTRIBUTOR: Ray Shan CITATION: Streitz and Mintmire, Phys Rev B, 50, 11996-12003 (1994)
   #
   # X (eV)                J (eV)          gamma (1/\AA)   zeta (1/\AA)    Z (e)

   Al      0.000000        10.328655       0.000000        0.968438        0.763905
   O       5.484763        14.035715       0.000000        2.143957        0.000000


----------

.. doxygenclass:: LAMMPS_NS::TextFileReader
   :project: progguide
   :members:

.. doxygenclass:: LAMMPS_NS::PotentialFileReader
   :project: progguide
   :members:

----------

Memory pool classes
-------------------

The memory pool classes are used for cases where otherwise many
small memory allocations would be needed and where the data would
be either all used or all freed.  One example for that is the
storage of neighbor lists.  The memory management strategy is
based on the assumption that allocations will be in chunks of similar
sizes.  The allocation is then not done per individual call for a
reserved chunk of memory, but for a "page" that can hold multiple
chunks of data.  A parameter for the maximum chunk size must be
provided, as that is used to determine whether a new page of memory
must be used.

The :cpp:class:`MyPage <LAMMPS_NS::MyPage>` class offers two ways to
reserve a chunk: 1) with :cpp:func:`get() <LAMMPS_NS::MyPage::get>` the
chunk size needs to be known in advance, 2) with :cpp:func:`vget()
<LAMMPS_NS::MyPage::vget>` a pointer to the next chunk is returned, but
its size is registered later with :cpp:func:`vgot()
<LAMMPS_NS::MyPage::vgot>`.

.. code-block:: C++
   :caption: Example of using :cpp:class:`MyPage <LAMMPS_NS::MyPage>`

      #include "my_page.h"
      using namespace LAMMPS_NS;

      MyPage<double> *dpage = new MyPage<double>;
      // max size of chunk: 256, size of page: 10240 doubles (=81920 bytes)
      dpage->init(256,10240);

      double **build_some_lists(int num)
      {
          dpage->reset();
          double **dlist = new double*[num];
          for (int i=0; i < num; ++i) {
              double *dptr = dpage.vget();
              int jnum = 0;
              for (int j=0; j < jmax; ++j) {
                  // compute some dvalue for eligible loop index j
                  dptr[j] = dvalue;
                  ++jnum;
              }
              if (dpage.status() != 0) {
                  // handle out of memory or jnum too large errors
              }
              dpage.vgot(jnum);
              dlist[i] = dptr;
          }
          return dlist;
      }

----------

.. doxygenclass:: LAMMPS_NS::MyPage
   :project: progguide
   :members:

.. doxygenclass:: LAMMPS_NS::MyPoolChunk
   :project: progguide
   :members:

----------

Eigensolver functions
---------------------

The ``MathEigen`` sub-namespace of the ``LAMMPS_NS`` namespace contains
functions and classes for eigensolvers. Currently only the
:cpp:func:`jacobi3 function <MathEigen::jacobi3>` is used in various
places in LAMMPS.  That function is built on top of a group of more
generic eigensolvers that are maintained in the ``math_eigen_impl.h``
header file.  This header contains the implementation of three template
classes:

#. "Jacobi" calculates all of the eigenvalues and eigenvectors
   of a dense, symmetric, real matrix.

#. The "PEigenDense" class only calculates the principal eigenvalue
   (ie. the largest or smallest eigenvalue), and its corresponding
   eigenvector.  However it is much more efficient than "Jacobi" when
   applied to large matrices (larger than 13x13).  PEigenDense also can
   understand complex-valued Hermitian matrices.

#. The "LambdaLanczos" class is a generalization of "PEigenDense" which can be
   applied to arbitrary sparse matrices.

The "math_eigen_impl.h" code is an amalgamation of `jacobi_pd
<https://github.com/jewettaij/jacobi_pd>`_ by Andrew Jewett at Scripps
Research (under CC0-1.0 license) and `Lambda Lanczos
<https://github.com/mrcdr/lambda-lanczos>`_ by Yuya Kurebayashi at
Tohoku University (under MIT license)

----------

.. doxygenfunction:: MathEigen::jacobi3(double const *const *mat, double *eval, double **evec)
   :project: progguide

.. doxygenfunction:: MathEigen::jacobi3(double const mat[3][3], double *eval, double evec[3][3])
   :project: progguide

---------------------------

Communication buffer coding with *ubuf*
---------------------------------------

LAMMPS uses communication buffers where it collects data from various
class instances and then exchanges the data with neighboring sub-domains.
For simplicity those buffers are defined as ``double`` buffers and
used for doubles and integer numbers. This presents a unique problem
when 64-bit integers are used.  While the storage needed for a ``double``
is also 64-bit, it cannot be used by a simple assignment.  To get around
that limitation, LAMMPS uses the :cpp:union:`ubuf <LAMMPS_NS::ubuf>`
union.  It is used in the various "pack" and "unpack" functions in the
LAMMPS classes to store and retrieve integers that may be 64-bit from
the communication buffers.

---------------------------

.. doxygenunion:: LAMMPS_NS::ubuf
   :project: progguide

