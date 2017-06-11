#include "meam.h"

using namespace LAMMPS_NS;

void
MEAM::meam_cleanup(void)
{

  memory->destroy(this->phirar6);
  memory->destroy(this->phirar5);
  memory->destroy(this->phirar4);
  memory->destroy(this->phirar3);
  memory->destroy(this->phirar2);
  memory->destroy(this->phirar1);
  memory->destroy(this->phirar);
  memory->destroy(this->phir);
}