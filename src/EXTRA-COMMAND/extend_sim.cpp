#include "extend_sim.h"
#include "domain.h"
#include "error.h"
using namespace LAMMPS_NS;

void ExtendSim::command(int narg, char **arg)
{
  if (domain->box_exist)
    error->all(FLERR, "extend_sim must be run before the simulation box exists");
  error->all(FLERR, "extend_sim: not implemented yet");
}
