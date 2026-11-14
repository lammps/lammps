// Stand-in for the PLUMED library: does nothing, except reporting an API version that
// fix plumed accepts, so that the fix can also be created.

#include "plumed/wrapper/Plumed.h"

#include <cstring>

namespace PLMD {
Plumed::Plumed() = default;
Plumed::~Plumed() = default;

void Plumed::cmd(const char *key, const void *val)
{
  if (val && (strcmp(key, "getApiVersion") == 0)) *(int *) const_cast<void *>(val) = 10;
}
}    // namespace PLMD
