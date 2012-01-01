# ifndef _TCPDEFS_H
# define _TCPDEFS_H

//Constants
const double sqr_2_over_3 = 0.816496580927726;
const double two_over_sqr_pi = 1.12837916709551;

# ifndef PLASMA_UNITS // using GridMD units: A, eV, m_proton=1
const double h_plank=0.06442021524615668;
const double h_sq=h_plank*h_plank;
const double m_electron=1./1836.1527556560675;
const double m_proton=1.;
const double coul_pref=14.39965172693122; // e^2/(4*pi*eps0), eV*A
const double eV_to_K=11604.447517053462; // Temperature in K correspondent to 1eV
# endif

# endif