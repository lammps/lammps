# Outline of include file conventions in LAMMPS

This purpose of this document is to provide a point of reference
for LAMMPS developers and contributors as to what include files
and definitions to put where into LAMMPS source.
Last change 2019-07-05

## Table of Contents

  * [Motivation](#motivation)
  * [Rules](#rules)
  * [Tools](#tools)
  * [Legacy Code](#legacy-code)

## Motivation

The conventions outlined in this document are supposed to help make
maintenance of the LAMMPS software easier.  By trying to achieve
consistency across files contributed by different developers, it will
become easier for the code maintainers to modify and adjust files and,
overall, the chance for errors or portability issues will be reduced.
The rules employed are supposed to minimize naming conflicts and
simplify dependencies between files and thus speed up compilation. They
may, as well, make otherwise hidden dependencies visible.

## Rules

Below are the various rules that are applied.  Not all are enforced
strictly and automatically.  If there are no significant side effects,
exceptions may be possible for cases where a full compliance to the
rules may require a large effort compared to the benefit.

### Core Files Versus Package Files

All rules listed below are most strictly observed for core LAMMPS files,
which are the files that are not part of a package, and the files of the
packages MOLECULE, MANYBODY, KSPACE, and RIGID.  On the other end of
the spectrum are USER packages and legacy packages that predate these
rules and thus may not be fully compliant.  Also, new contributions
will be checked more closely, while existing code will be incrementally
adapted to the rules as time and required effort permits.

### System Versus Local Header Files

All system- or library-provided include files are included with angular
brackets (examples: `#include <cstring>` or `#include <mpi.h>`) while
include files provided with LAMMPS are included with double quotes
(examples: `#include "pointers.h"` or `#include "compute_temp.h"`).

For headers declaring functions of the C-library, the corresponding
C++ versions should be included (examples: `#include <cstdlib>` or
`#include <cctypes>` instead of `#include <stdlib.h>` or
`#include<ctypes.h>` ).

### C++ Standard Compliance

LAMMPS core files use standard conforming C++ compatible with the
C++11 standard, unless explicitly noted.  Also, LAMMPS uses the C-style
stdio library for I/O instead of iostreams.  Since using both at the
same time can cause problems, iostreams should be avoided where possible.

### Lean Header Files

Header files will typically contain the definition of a (single) class.
These header files should have as few include statements as possible.
This is particularly important for classes that implement a "style" and
thus use a macro of the kind `SomeStyle(some/name,SomeName)`. These will
all be included in the auto-generated `"some_style.h"` files which 
results in a high potential for direct or indirect symbol name clashes.

In the ideal case, the header would only include one file defining the
parent class. That would typically be either `#include "pointers.h"` for
the `Pointers` class, or a header of a class derived from it like
`#include "pair.h"` for the `Pair` class and so on.  References to other
classes inside the class should be make through pointers, for which forward
declarations (inside the `LAMMPS_NS` or the new class' namespace) can
be employed.  The full definition will then be included into the corresponding
implementation file.  In the given example from above, the header file
would be called `some_name.h` and the implementation `some_name.cpp` (all
lower case with underscores, while the class itself would be in camel case
and no underscores `SomeName`, and the style name with lower case names separated by
a forward slash).

### Implementation Files

In the implementation files (typically, those would have the same base name
as the corresponding header with a .cpp extension instead of .h) include
statements should follow the "include what you use" principle.

### Order of Include Statements

Include files should be included in this order:
* the header matching the implementation (`some_class.h` for file `some_class.cpp`)
* mpi.h
* system and library headers (anything that is using angular brackets; C-library headers first, then C++)
* LAMMPS local headers (preferably in alphabetical order)

### Special Cases and Exceptions

#### pointers.h

The `pointer.h` header file also includes `cstdio` and `lmptype.h`
(and through it `stdint.h`, `intttypes.h`, cstdlib, and `climits`).
This means any header including `pointers.h` can assume that `FILE`,
`NULL`, `INT_MAX` are defined.

## Tools

The [Include What You Use tool](https://include-what-you-use.org/)
can be used to provide supporting information about compliance with
the rules listed here.  There are some limitations and the IWYU tool
may give incorrect advice.  The tools is activated by setting the
CMake variable `CMAKE_CXX_INCLUDE_WHAT_YOU_USE` variable to the
path of the `include-what-you-use` command.  When activated, the
tool will be run after each compilation and provide suggestions for
which include files should be added or removed.

## Legacy Code

A lot of code predates the application of the rules in this document
and the rules themselves are a moving target.  So there are going to be
significant chunks of code that do not fully comply.  This applies
for example to the USER-REAXC, or the USER-ATC package.  The LAMMPS
developers are dedicated to make an effort to improve the compliance
and welcome volunteers wanting to help with the process.

