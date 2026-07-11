/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Coverity Scan modeling file.

   This file is NOT part of the LAMMPS build.  It is compiled on its own
   with the Coverity analysis toolkit and uploaded to the Coverity Scan
   project so that the server-side analysis applies the models below:

       cov-make-library -of user_models.xmldb coverity_model.cpp

   See tools/coverity/README.md for the build and upload workflow.

   A model describes, to the Coverity engine, only the behavior of a
   function that matters for defect detection -- using the special
   __coverity_*__ modeling intrinsics -- not its real implementation.
   Coverity matches a model to the actual function by signature
   (namespace + class + name + parameter types), so each stub below must
   declare exactly the same signature as the corresponding LAMMPS source.

   Keep these signatures in sync with the headers cited in each comment.
   A drifted signature stops matching silently (no error) and the model
   goes inert.
------------------------------------------------------------------------- */

#include <cstdint>

namespace LAMMPS_NS {

// Must match the typedef in src/lmptype.h for the analyzed build.
// The Coverity Scan build uses -D LAMMPS_SIZES=SMALLBIG (see
// .github/workflows/coverity.yml), where bigint == int64_t.
typedef int64_t bigint;

/* ----------------------------------------------------------------------
   Custom memory allocator -- src/memory.h

   The templated Memory::create / grow / destroy helpers all funnel
   through these three primitives, so modeling the primitives teaches
   Coverity the allocate/free pairing for the entire custom allocator and
   removes spurious "resource leak" and "use after free" reports that come
   from the analyzer not recognizing it as an allocator.
------------------------------------------------------------------------- */

class Memory {
 public:
  void *smalloc(bigint n, const char *);
  void *srealloc(void *, bigint n, const char *);
  void sfree(void *);
};

void *Memory::smalloc(bigint n, const char *)
{
  // returns a fresh allocation of n bytes; tying the size to n also
  // informs OVERRUN/buffer-size analysis of the allocation extent.
  return __coverity_alloc__(n);
}

void *Memory::srealloc(void *ptr, bigint n, const char *)
{
  // standard realloc model: the old block is released and a new block of
  // n bytes is returned.
  __coverity_free__(ptr);
  return __coverity_alloc__(n);
}

void Memory::sfree(void *ptr)
{
  __coverity_free__(ptr);
}

/* ----------------------------------------------------------------------
   NOT modeled here on purpose:

   Error::all / Error::one / Error::done are already declared [[noreturn]]
   in src/error.h, an attribute cov-analyze honors directly.  An in-source
   attribute is preferable to a model -- it also informs the compiler,
   clang-tidy, and every other tool, and it is version controlled and
   reviewed.  Rule of thumb: when we own the code, annotate it in the
   source ([[noreturn]], [[nodiscard]], assert, explicit checks) rather
   than adding a model.  Reserve this file for code we cannot annotate or
   where no suitable attribute exists.
------------------------------------------------------------------------- */

}    // namespace LAMMPS_NS
