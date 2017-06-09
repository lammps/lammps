#include "meam.h"

void
MEAM::meam_cleanup_(void)
{

  meam_deallocate(this->phirar6);
  meam_deallocate(this->phirar5);
  meam_deallocate(this->phirar4);
  meam_deallocate(this->phirar3);
  meam_deallocate(this->phirar2);
  meam_deallocate(this->phirar1);
  meam_deallocate(this->phirar);
  meam_deallocate(this->phir);
}