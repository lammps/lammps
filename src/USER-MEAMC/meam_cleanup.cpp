extern "C" {
#include "meam.h"

void
meam_cleanup_(void)
{

  deallocate(meam_data.phirar6);
  deallocate(meam_data.phirar5);
  deallocate(meam_data.phirar4);
  deallocate(meam_data.phirar3);
  deallocate(meam_data.phirar2);
  deallocate(meam_data.phirar1);
  deallocate(meam_data.phirar);
  deallocate(meam_data.phir);
}
}