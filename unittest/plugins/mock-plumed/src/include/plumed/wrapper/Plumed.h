/* Stand-in for the wrapper header of PLUMED (src/wrapper/Plumed.h), only for testing that
   the PLUMED plugin of LAMMPS configures, compiles, links and loads.  It declares only what
   fix plumed uses, and the class does nothing.  The real cmd() is type checked and also
   takes some arguments by value; this one takes any pointer. */

#ifndef __PLUMED_wrapper_Plumed_h
#define __PLUMED_wrapper_Plumed_h

#ifdef __cplusplus
namespace PLMD {
class Plumed {
 public:
  Plumed();
  ~Plumed();
  void cmd(const char *key, const void *val = nullptr);
};
}    // namespace PLMD
#endif

#endif
