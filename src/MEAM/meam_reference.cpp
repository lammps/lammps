// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "meam.h"

#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MEAM_NS;

void hcp_shpfcn(const double sthe, const double cthe, double (&s)[3])
{
  s[2] = 1.0 / 3.0;
}

void dia_shpfcn(const double sthe, const double cthe, double (&s)[3])
{
  s[2] = 32.0 / 9.0;
}

void dim_shpfcn(const double sthe, const double cthe, double (&s)[3])
{
  s[0] = 1.0;
  s[1] = 2.0 / 3.0;
  //        s(3) = 1.d0 // this should be 0.4 unless (1-legendre) is multiplied in the density calc.
  s[2] = 0.40;    // this is (1-legendre) where legendre = 0.6 in dynamo is accounted.
}

void lin_shpfcn(const double sthe, const double cthe, double (&s)[3])
{
  s[0] = 0.0;
  // FIXME: comment and value are in conflict
  s[1] = 8.0 / 3.0;    // 4*(co**4 + si**4 - 1.0/3.0) in zig become 4*(1-1/3)
  s[2] = 0.0;
}

void zig_shpfcn(const double sthe, const double cthe, double (&s)[3])
{
  s[0] = 4.0 * pow(cthe, 2);
  s[1] = 4.0 * (pow(cthe, 4) + pow(sthe, 4) - 1.0 / 3.0);
  s[2] = 4.0 * (pow(cthe, 2) * (3 * pow(sthe, 4) + pow(cthe, 4)));
  s[2] = s[2] - 0.6 * s[0];    //legend in dyn, 0.6 is default value.
}

// Ensure these are in the same order as lattice_t (for now)!
const reference_lattice_t MEAM_NS::lattice_defs[MAXLAT] = {
  {
    .name                    = "fcc",
    .single                  = true,
    .re                      = 1.0 / sqrt(2.0),
    .Zij                     = 12,
    .Zij2                    = 6,
    .Nscr2                   = 4,
    .ratio_2nn               = sqrt(2.0),
  },
  {
    .name                    = "bcc",
    .single                  = true,
    .re                      = sqrt(3.0) / 2.0,
    .Zij                     = 8,
    .Zij2                    = 6,
    .Nscr2                   = 4,
    .ratio_2nn               = 2.0 / sqrt(3.0),
  },
  {
    .name                    = "hcp",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 12,
    .Zij2                    = 6,
    .Nscr2                   = 4,
    .ratio_2nn               = sqrt(2.0),
    .shpfcn                  = hcp_shpfcn,
  },
  {
    .name                    = "dim",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 1,
    .shpfcn                  = dim_shpfcn,
  },
  {
    .name                    = "dia",
    .single                  = true,
    .re                      = sqrt(3.0) / 4.0,
    .Zij                     = 4,
    .Zij2                    = 12,
    .Nscr2                   = 1,
    .ratio_2nn               = sqrt(8.0 / 3.0),
    .shpfcn                  = dia_shpfcn,
  },
  {
    .name                    = "dia3",
    .single                  = true,
    .re                      = sqrt(3.0) / 4.0,
    .Zij                     = 4,
    .Zij2                    = 12,
    .Nscr2                   = 4,
    .ratio_2nn               = sqrt(11.0 / 3.0),
    .shpfcn                  = dia_shpfcn,
    .nn3                     = true,
  },
  {
    .name                    = "b1",
    .Zij                     = 6,
    .Zij2                    = 12,
    .Nscr2                   = 2,
    .ratio_2nn               = sqrt(2.0),
  },
  {
    .name                    = "c11",
    .Zij                     = 10,
  },
  {
    .name                    = "l12",
    .Zij                     = 12,
    .Zij2                    = 6,
    .Nscr2                   = 4,
    .ratio_2nn               = sqrt(2.0),
  },
  {
    .name                    = "b2",
    .Zij                     = 8,
    .Zij2                    = 6,
    .Nscr2                   = 4,
    .ratio_2nn               = 2.0 / sqrt(3.0),
  },
  {
    .name                    = "ch4",
    .re                      = 1.0,
    .Zij                     = 4,
    // CH4 actually needs shape factor for diamond for C, dimer for H
    .shpfcn                  = dia_shpfcn,
  },
  {
    .name                    = "lin",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 2,
    .shpfcn                  = lin_shpfcn,
  },
  {
    .name                    = "zig",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 2,
    .Zij2                    = 2,
    .Nscr2                   = 1,
    .ratio_2nn               = 2.0,
    .ratio_2nn_angular       = true,
    .shpfcn                  = zig_shpfcn,
  },
  {
    .name                    = "tri",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 2,
    .Zij2                    = 1,
    .Nscr2                   = 2,
    .ratio_2nn               = 2.0,
    .ratio_2nn_angular       = true,
    .shpfcn                  = zig_shpfcn,
  },
  {
    .name                    = "sc",
    .single                  = true,
    .re                      = 1.0,
    .Zij                     = 6,
    .Zij2                    = 12,
    .Nscr2                   = 2,
    .ratio_2nn               = sqrt(2.0),
  },
};

static_assert(MAXLAT == sizeof(lattice_defs) / sizeof(reference_lattice_t));

/* ----------------------------------------------------------------------
   Convert lattice name from input to to lattice_t
   return false on failure
   return true and set lat on success
------------------------------------------------------------------------- */

bool MEAM_NS::str_to_lat(const std::string & str, bool single, lattice_t& lat)
{
  for (int i=0; i<MAXLAT; i++) {
    const reference_lattice_t& def = lattice_defs[i];

    if (single && !def.single)
      continue;

    if (str == def.name) {
      lat = (lattice_t)i;
      return true;
    }
  }
  return false;
}


/* ----------------------------------------------------------------------
   Number of first neighbors for reference structure
------------------------------------------------------------------------- */

int MEAM_NS::get_Zij(const lattice_t latt)
{
  const int lidx = (int)latt;
  if (latt < 0 || latt >= MAXLAT)
    return -1;

  return lattice_defs[lidx].Zij;
}

/* ----------------------------------------------------------------------
   Number of second neighbors for the reference structure
     a = distance ratio R1/R2 (a2nn in dynamo)
     numscr = number of atoms that screen the 2NN bond
     S = second neighbor screening function (xfac, a part of b2nn in dynamo)
------------------------------------------------------------------------- */

int MEAM_NS::get_Zij2(const lattice_t latt, const double cmin, const double cmax, const double stheta,
                      double &arat, double &S)
{
  double C, sijk;
  int Zij2, numscr;

  const int lidx = (int)latt;
  if (latt < 0 || latt >= MAXLAT)
    return -1;
  const reference_lattice_t& def = lattice_defs[lidx];

  Zij2 = def.Zij2;
  numscr = def.Nscr2;
  arat = def.ratio_2nn;
  if (def.ratio_2nn_angular) {
    arat *= stheta;
  }

  if (iszero(arat)) {
    // error
    return -1;
  }

  // Compute screening for each first neighbor
  if (def.nn3) {
    C = 1.0;
  } else {
    C = 4.0 / (arat * arat) - 1.0;
  }
  // There are numscr first neighbors screening the second neighbors
  sijk = Csijk(C, cmin, cmax);
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}

/* ----------------------------------------------------------------------
   Number of second neighbors for the reference structure
   for the case of different elements on alternating positions
------------------------------------------------------------------------- */

int MEAM_NS::get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S)
{

  double sijk;
  int numscr = 0, Zij2 = 0;

  switch (latt) {
    case ZIG:    //zig-zag for b11s and b22s term
    case TRI:    //trimer for b11s
      Zij2 = 2;
      numscr = 1;
      break;
    default:
      // unknown lattic flag in get Zij
      //        call error('Lattice not defined in get_Zij.')
      break;
  }
  sijk = Csijk(1.0, cmin, cmax);
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}

/* ----------------------------------------------------------------------
   Shape factors for various configurations
------------------------------------------------------------------------- */

void MEAM_NS::get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3])
{
  const int lidx = (int)latt;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 0.0;

  if (latt < 0 || latt >= MAXLAT) {
    return;
  }
  const reference_lattice_t& def = lattice_defs[lidx];

  if (def.shpfcn == nullptr) {
    return;
  }

  def.shpfcn(sthe, cthe, s);
}

/* ----------------------------------------------------------------------
   Average weighting factors for the reference structure
------------------------------------------------------------------------- */

void MEAM::get_tavref(double* t11av, double* t21av, double* t31av, double* t12av, double* t22av, double* t32av,
                      double t11, double t21, double t31, double t12, double t22, double t32, double r, int a,
                      int b, lattice_t latt)
{
  double rhoa01, rhoa02, a1, a2;

  //     For ialloy = 2, no averaging is done
  if (ialloy == 2) {
    *t11av = t11;
    *t21av = t21;
    *t31av = t31;
    *t12av = t12;
    *t22av = t22;
    *t32av = t32;
  } else switch (latt)  {
    case FCC:
    case BCC:
    case DIA:
    case DIA3:
    case HCP:
    case B1:
    case DIM:
    case B2:
    case CH4:
    case LIN:
    case ZIG:
    case TRI:
    case SC:
      //     all neighbors are of the opposite type
      *t11av = t12;
      *t21av = t22;
      *t31av = t32;
      *t12av = t11;
      *t22av = t21;
      *t32av = t31;
      break;
    default:
      a1 = r / re_meam[a][a] - 1.0;
      a2 = r / re_meam[b][b] - 1.0;
      rhoa01 = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
      rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
      if (latt == L12) {
        double rho01 = 8 * rhoa01 + 4 * rhoa02;
        *t12av = t11;
        *t22av = t21;
        *t32av = t31;
        if (rho01 > 0.0) {
          *t11av = (8 * t11 * rhoa01 + 4 * t12 * rhoa02) / rho01;
          *t21av = (8 * t21 * rhoa01 + 4 * t22 * rhoa02) / rho01;
          *t31av = (8 * t31 * rhoa01 + 4 * t32 * rhoa02) / rho01;
        } else {
          // limit for rhoa01 and rhoa02 -> 0.0. Should not happen.
          *t11av = (2.0 * t11 + t12) / 3.0;
          *t21av = (2.0 * t21 + t22) / 3.0;
          *t31av = (2.0 * t31 + t32) / 3.0;
        }
      } else {
        //      call error('Lattice not defined in get_tavref.')
      }
  }
}

/* ----------------------------------------------------------------------
   Calculate density functions, assuming reference configuration
------------------------------------------------------------------------- */

void MEAM::get_densref(double r, int a, int b, double* rho01, double* rho11, double* rho21, double* rho31,
                       double* rho02, double* rho12, double* rho22, double* rho32,
                       double* rho1m1, double* rho2m1, double* rho3m1,
                       double* rho1m2, double* rho2m2, double* rho3m2)
{
  double a1, a2;
  double t1ma, t2ma, t3ma, t1mb, t2mb, t3mb;
  double s[3];
  lattice_t lat;
  int Zij,Zij2nn;
  double rhoa01nn, rhoa02nn;
  double rhoa01, rhoa11, rhoa21, rhoa31;
  double rhoa02, rhoa12, rhoa22, rhoa32;
  double arat, scrn, denom;
  double C, s111, s112, s221, S11, S22;
  // msmeam
  double rhoa1m1, rhoa2m1, rhoa3m1, rhoa1m2, rhoa2m2, rhoa3m2;

  t1ma = t2ma = t3ma = t1mb = t2mb = t3mb = 1.0;
  a1 = r / re_meam[a][a] - 1.0;
  a2 = r / re_meam[b][b] - 1.0;

  if (msmeamflag) {
    // the rho variables are multiplied by t here since ialloy not needed in msmeam
    t1ma = t1_meam[a];
    t2ma = t2_meam[a];
    t3ma = t3_meam[a];
    t1mb = t1_meam[b];
    t2mb = t2_meam[b];
    t3mb = t3_meam[b];
    // msmeam specific rho vars
    rhoa1m1 = rho0_meam[a] * t1m_meam[a] * MathSpecial::fm_exp(-beta1m_meam[a] * a1);
    rhoa2m1 = rho0_meam[a] * t2m_meam[a] * MathSpecial::fm_exp(-beta2m_meam[a] * a1);
    rhoa3m1 = rho0_meam[a] * t3m_meam[a] * MathSpecial::fm_exp(-beta3m_meam[a] * a1);
    rhoa1m2 = rho0_meam[b] * t1m_meam[b] * MathSpecial::fm_exp(-beta1m_meam[b] * a2);
    rhoa2m2 = rho0_meam[b] * t2m_meam[b] * MathSpecial::fm_exp(-beta2m_meam[b] * a2);
    rhoa3m2 = rho0_meam[b] * t3m_meam[b] * MathSpecial::fm_exp(-beta3m_meam[b] * a2);
  }
  rhoa01 = rho0_meam[a]        * MathSpecial::fm_exp(-beta0_meam[a] * a1);
  rhoa11 = rho0_meam[a] * t1ma * MathSpecial::fm_exp(-beta1_meam[a] * a1);
  rhoa21 = rho0_meam[a] * t2ma * MathSpecial::fm_exp(-beta2_meam[a] * a1);
  rhoa31 = rho0_meam[a] * t3ma * MathSpecial::fm_exp(-beta3_meam[a] * a1);
  rhoa02 = rho0_meam[b]        * MathSpecial::fm_exp(-beta0_meam[b] * a2);
  rhoa12 = rho0_meam[b] * t1mb * MathSpecial::fm_exp(-beta1_meam[b] * a2);
  rhoa22 = rho0_meam[b] * t2mb * MathSpecial::fm_exp(-beta2_meam[b] * a2);
  rhoa32 = rho0_meam[b] * t3mb * MathSpecial::fm_exp(-beta3_meam[b] * a2);

  lat = lattce_meam[a][b];

  Zij = get_Zij(lat);

  *rho11 = 0.0;
  *rho21 = 0.0;
  *rho31 = 0.0;
  *rho12 = 0.0;
  *rho22 = 0.0;
  *rho32 = 0.0;
  if (msmeamflag) {
    *rho1m1 = 0.0;
    *rho2m1 = 0.0;
    *rho3m1 = 0.0;
    *rho1m2 = 0.0;
    *rho2m2 = 0.0;
    *rho3m2 = 0.0;
  }

  // keep track of density components separately; combine in the calling subroutine
  switch (lat) {
    case FCC:
      *rho01 = 12.0 * rhoa02;
      *rho02 = 12.0 * rhoa01;
      break;
    case BCC:
      *rho01 = 8.0 * rhoa02;
      *rho02 = 8.0 * rhoa01;
      break;
    case B1:
    case SC:
      *rho01 = 6.0 * rhoa02;
      *rho02 = 6.0 * rhoa01;
      break;
    case DIA:
    case DIA3:
      *rho01 = 4.0 * rhoa02;
      *rho02 = 4.0 * rhoa01;
      *rho31 = 32.0 / 9.0 * rhoa32 * rhoa32;
      *rho32 = 32.0 / 9.0 * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho3m1 = 32.0 / 9.0 * rhoa3m2 * rhoa3m2;
        *rho3m2 = 32.0 / 9.0 * rhoa3m1 * rhoa3m1;
      }
      break;
    case HCP:
      *rho01 = 12 * rhoa02;
      *rho02 = 12 * rhoa01;
      *rho31 = 1.0 / 3.0 * rhoa32 * rhoa32;
      *rho32 = 1.0 / 3.0 * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho3m1 = 1.0 / 3.0 * rhoa3m2 * rhoa3m2;
        *rho3m2 = 1.0 / 3.0 * rhoa3m1 * rhoa3m1;
      }
      break;
    case DIM:
      get_shpfcn(DIM, 0, 0, s);
      *rho01 = rhoa02;
      *rho02 = rhoa01;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho31 = s[2] * rhoa32 * rhoa32;
      *rho32 = s[2] * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho1m1 = s[0] * rhoa1m2 * rhoa1m2;
        *rho1m2 = s[0] * rhoa1m1 * rhoa1m1;
        *rho2m1 = s[1] * rhoa2m2 * rhoa2m2;
        *rho2m2 = s[1] * rhoa2m1 * rhoa2m1;
        *rho3m1 = s[2] * rhoa3m2 * rhoa3m2;
        *rho3m2 = s[2] * rhoa3m1 * rhoa3m1;
      }
      break;
    case C11:
      *rho01 = rhoa01;
      *rho02 = rhoa02;
      *rho11 = rhoa11;
      *rho12 = rhoa12;
      *rho21 = rhoa21;
      *rho22 = rhoa22;
      *rho31 = rhoa31;
      *rho32 = rhoa32;
      if (msmeamflag) {
        *rho1m1 = rhoa1m1;
        *rho1m2 = rhoa1m2;
        *rho2m1 = rhoa2m1;
        *rho2m2 = rhoa2m2;
        *rho3m1 = rhoa3m1;
        *rho3m2 = rhoa3m2;
      }
      break;
    case L12:
      *rho01 = 8 * rhoa01 + 4 * rhoa02;
      *rho02 = 12 * rhoa01;
      if (!msmeamflag && ialloy == 1) {
        *rho21 = 8. / 3. * MathSpecial::square(rhoa21 * t2_meam[a] - rhoa22 * t2_meam[b]);
        denom = 8 * rhoa01 * MathSpecial::square(t2_meam[a]) + 4 * rhoa02 * MathSpecial::square(t2_meam[b]);
        if (denom > 0.)
          *rho21 = *rho21 / denom * *rho01;
      } else {
        *rho21 = 8. / 3. * (rhoa21 - rhoa22) * (rhoa21 - rhoa22);
      }
      if (msmeamflag) {
        *rho2m1 = 8. / 3. * (rhoa2m1 - rhoa2m2) * (rhoa2m1 - rhoa2m2);
      }
      break;
    case B2:
      *rho01 = 8.0 * rhoa02;
      *rho02 = 8.0 * rhoa01;
      break;
    case CH4:
      *rho01 = 4.0 * rhoa02; //in assumption that 'a' represent carbon
      *rho02 = rhoa01;       //in assumption that 'b' represent hydrogen

      get_shpfcn(DIM, 0, 0, s); //H
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;

      get_shpfcn(CH4, 0, 0, s); //C
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case LIN:
      *rho01 = rhoa02*Zij;
      *rho02 = rhoa01*Zij;

      get_shpfcn(LIN, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case ZIG:
      *rho01 = rhoa02*Zij;
      *rho02 = rhoa01*Zij;

      get_shpfcn(ZIG, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case TRI:
      *rho01 = rhoa02;
      *rho02 = rhoa01*Zij;

      get_shpfcn(TRI, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      s[0] = 1.0;
      s[1] = 2.0/3.0;
      s[2] = 1.0 - 0.6*s[0];

      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;


    // default:
    //        call error('Lattice not defined in get_densref.')
  }

  if (nn2_meam[a][b] == 1) {


    Zij2nn = get_Zij2(lat, Cmin_meam[a][a][b], Cmax_meam[a][a][b],
                      stheta_meam[a][b], arat, scrn);

    a1 = arat * r / re_meam[a][a] - 1.0;
    a2 = arat * r / re_meam[b][b] - 1.0;

    rhoa01nn = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
    rhoa02nn = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);

    if (lat == L12) {
      //     As usual, L12 thinks it's special; we need to be careful computing
      //     the screening functions
      C = 1.0;
      s111 = get_sijk(C, a, a, a);
      s112 = get_sijk(C, a, a, b);
      s221 = get_sijk(C, b, b, a);
      S11 = s111 * s111 * s112 * s112;
      S22 = s221 * s221 * s221 * s221;
      *rho01 = *rho01 + 6 * S11 * rhoa01nn;
      *rho02 = *rho02 + 6 * S22 * rhoa02nn;

    } else {
      //     For other cases, assume that second neighbor is of same type,
      //     first neighbor may be of different type

      *rho01 = *rho01 + Zij2nn * scrn * rhoa01nn;

      //     Assume Zij2nn and arat don't depend on order, but scrn might
      Zij2nn = get_Zij2(lat, Cmin_meam[b][b][a], Cmax_meam[b][b][a],
                        stheta_meam[a][b], arat, scrn);
      *rho02 = *rho02 + Zij2nn * scrn * rhoa02nn;
    }
  }
}
