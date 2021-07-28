// Hash function hashlittle()
// from lookup3.c, by Bob Jenkins, May 2006, Public Domain
// bob_jenkins@burtleburtle.net

#ifndef LMP_HASHLITTLE_H
#define LMP_HASHLITTLE_H

#include <cstddef>
#include <cstdint>

namespace LAMMPS_NS {
uint32_t hashlittle(const void *key, size_t length, uint32_t);
}
#endif
