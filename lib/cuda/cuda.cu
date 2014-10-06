#include "cuda_precision.h"
#include "cuda_shared.h"
#include "cuda_cu.h"

void Cuda_Cuda_GetCompileSettings(cuda_shared_data* sdata)
{
  sdata->compile_settings.prec_glob = sizeof(CUDA_CFLOAT) / 4;
  sdata->compile_settings.prec_x = sizeof(X_CFLOAT) / 4;
  sdata->compile_settings.prec_v = sizeof(V_CFLOAT) / 4;
  sdata->compile_settings.prec_f = sizeof(F_CFLOAT) / 4;
  sdata->compile_settings.prec_pppm = sizeof(PPPM_CFLOAT) / 4;
  sdata->compile_settings.prec_fft = sizeof(FFT_CFLOAT) / 4;

#ifdef FFT_CUFFT
  sdata->compile_settings.cufft = 1;
#else
  sdata->compile_settings.cufft = 0;
#endif

  sdata->compile_settings.arch = CUDA_ARCH;

}
