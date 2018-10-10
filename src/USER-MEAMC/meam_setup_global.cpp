#include "meam.h"
#include <cmath>
using namespace LAMMPS_NS;

void
MEAM::meam_setup_global(int nelt, lattice_t* lat, double* z, int* ielement, double* /*atwt*/, double* alpha,
                        double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                        double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                        int* ibar)
{

  int i;
  double tmplat[maxelt];

  this->neltypes = nelt;

  for (i = 0; i < nelt; i++) {
    this->lattce_meam[i][i] = lat[i];

    this->Z_meam[i] = z[i];
    this->ielt_meam[i] = ielement[i];
    this->alpha_meam[i][i] = alpha[i];
    this->beta0_meam[i] = b0[i];
    this->beta1_meam[i] = b1[i];
    this->beta2_meam[i] = b2[i];
    this->beta3_meam[i] = b3[i];
    tmplat[i] = alat[i];
    this->Ec_meam[i][i] = esub[i];
    this->A_meam[i] = asub[i];
    this->t0_meam[i] = t0[i];
    this->t1_meam[i] = t1[i];
    this->t2_meam[i] = t2[i];
    this->t3_meam[i] = t3[i];
    this->rho0_meam[i] = rozero[i];
    this->ibar_meam[i] = ibar[i];

    if (this->lattce_meam[i][i] == FCC)
      this->re_meam[i][i] = tmplat[i] / sqrt(2.0);
    else if (this->lattce_meam[i][i] == BCC)
      this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 2.0;
    else if (this->lattce_meam[i][i] == HCP)
      this->re_meam[i][i] = tmplat[i];
    else if (this->lattce_meam[i][i] == DIM)
      this->re_meam[i][i] = tmplat[i];
    else if (this->lattce_meam[i][i] == DIA)
      this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 4.0;
    else {
      //           error
    }
  }

  // Set some defaults
  this->rc_meam = 4.0;
  this->delr_meam = 0.1;
  setall2d(this->attrac_meam, 0.0);
  setall2d(this->repuls_meam, 0.0);
  setall3d(this->Cmax_meam, 2.8);
  setall3d(this->Cmin_meam, 2.0);
  setall2d(this->ebound_meam, pow(2.8, 2) / (4.0 * (2.8 - 1.0)));
  setall2d(this->delta_meam, 0.0);
  setall2d(this->nn2_meam, 0);
  setall2d(this->zbl_meam, 1);
  this->gsmooth_factor = 99.0;
  this->augt1 = 1;
  this->ialloy = 0;
  this->mix_ref_t = 0;
  this->emb_lin_neg = 0;
  this->bkgd_dyn = 0;
  this->erose_form = 0;
}
