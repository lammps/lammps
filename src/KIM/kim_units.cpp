/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ellad B. Tadmor (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.0.2 (and newer) package
------------------------------------------------------------------------- */

#include <iostream>
#include <math.h>
#include <string>
#include <map>
using namespace std;

namespace
{

// Constants of nature and basic conversion factors
// Source: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
//         Working with NIST values even when there are newer values for
//         compatibility with LAMMPS

/*----------------------
   Fundamental constants
------------------------ */
double const boltz_si = 1.38064852e-23;  // [J K^-1] Boltzmann's factor
                                         //   (NIST value)
double const Nav = 6.022140857e23;       // [unitless] Avogadro's number
                                         //   (NIST value)
// double const Nav = 6.02214076e23;     // [unitless] Avogadro's number
                                         //   (official value May 2019)
double const me_si = 9.10938356e-31;     // [kg] electron rest mass
                                         //   (NIST value)
// double me_si = 9.10938291e-31;        // [kg] electron rest mass
double const e_si = 1.6021766208e-19;    // [C] elementary charge
                                         //   (charge of an electron/proton)
                                         //   (NIST value)

/*----------------------
   Distance units
------------------------ */
double const bohr_si = 5.2917721067e-11;  // [m] Bohr unit (distance between
                                          //   nucleus and electron in H)
                                          //   (NIST value)
double const angstrom_si = 1e-10;         // [m] Angstrom
double const centimeter_si = 1e-2;        // [m] centimeter
double const micrometer_si = 1e-6;        // [m] micrometer (micron)
double const nanometer_si = 1e-9;         // [m] nanometer

/*----------------------
   Mass units
------------------------ */
double const gram_per_mole_si = 1e-3/Nav;  // [kg] gram per mole
double const amu_si = 1e-3/Nav;            // [kg] atomic mass unit (molecular
                                           //   weight) For example, the mean
                                           //   molecular weight of water
                                           //   is 18.015 atomic mass units
                                           //   (amu), so one mole of water
                                           //   weight 18.015 grams.
double const gram_si = 1e-3;               // [kg] gram
double const picogram_si = 1e-15;          // [kg] picogram
double const attogram_si = 1e-21;          // [kg[ attogram

/*----------------------
   Time units
------------------------ */
double const atu_si = 2.418884326509e-17;  // [s] atomic time unit
                                           //   ( = hbar/E_h where E_h is the
                                           //   Hartree energy) (NIST value)
double const atu_electron_si = atu_si*sqrt(amu_si/me_si);
                                          // [s] atomic time unit
                                          //   used in electron system (see https://sourceforge.net/p/lammps/mailman/lammps-users/thread/BCA2BDB2-BA03-4280-896F-1E6120EF47B2%40caltech.edu/)
double const microsecond_si = 1e-6;        // [s] microsecond
double const nanosecond_si = 1e-9;         // [s] nanosecond
double const picosecond_si = 1e-12;        // [s] picosecond
double const femtosecond_si = 1e-15;       // [s] femtosecond

/*----------------------
   Density units
------------------------ */
double const gram_per_centimetercu_si =
             gram_si/pow(centimeter_si,3);  // [kg/m^3] gram/centimeter^3
double const amu_per_bohrcu_si = amu_si/pow(bohr_si,3);  // [kg/m^3] amu/bohr^3
double const picogram_per_micrometercu_si =
             picogram_si/pow(micrometer_si,3); // [kg/m^3] picogram/micrometer^3
double const attogram_per_nanometercu_si =
             attogram_si/pow(nanometer_si,3);  // [kg/m^3] attogram/nanometer^3

/*----------------------
   Energy/torque units
------------------------ */
double const kcal_si = 4184.0;              // [J] kilocalorie (heat energy
                                            //   involved in warming up one
                                            //   kilogram of water by one
                                            //   degree Kelvin)
double const ev_si = 1.6021766208e-19;      // [J] electon volt (amount of
                                            //   energy gained or lost by the
                                            //   charge of a single electron
                                            //   moving across an electric
                                            //   potential difference of one
                                            //   volt.) (NIST value)
double const hartree_si = 4.359744650e-18;  // [J] Hartree (approximately the
                                            //   electric potential energy of
                                            //   the hydrogen atom in its
                                            //   ground state) (NIST value)
double const kcal_per_mole_si = kcal_si/Nav;// [J] kcal/mole
double const erg_si = 1e-7;                 // [J] erg
double const dyne_centimeter_si = 1e-7;     // [J[ dyne*centimeter
double const picogram_micrometersq_per_microsecondsq_si =
             picogram_si*pow(micrometer_si,2)/pow(microsecond_si,2);
                                            // [J] picogram*micrometer^2/
                                            //   microsecond^2
double const attogram_nanometersq_per_nanosecondsq_si =
             attogram_si*pow(nanometer_si,2)/pow(nanosecond_si,2);
                                            // [J] attogram*nanometer^2/
                                            //   nanosecond^2

/*----------------------
   Velocity units
------------------------ */
double const angstrom_per_femtosecond_si =
             angstrom_si/femtosecond_si;      // [m/s] Angstrom/femtosecond
double const angstrom_per_picosecond_si =
             angstrom_si/picosecond_si;       // [m/s] Angstrom/picosecond
double const micrometer_per_microsecond_si =
             micrometer_si/microsecond_si;    // [m/s] micrometer/microsecond
double const nanometer_per_nanosecond_si =
             nanometer_si/nanosecond_si;      // [m/s] nanometer/nanosecond
double const centimeter_per_second_si =
             centimeter_si;                  // [m/s] centimeter/second
double const bohr_per_atu_electron_si =
             bohr_si/atu_electron_si;        // [m/s] bohr/atu

/*----------------------
   Force units
------------------------ */
double const kcal_per_mole_angstrom_si =
             kcal_per_mole_si/angstrom_si;       // [N] kcal/(mole*Angstrom)
double const ev_per_angstrom_si =
             ev_si/angstrom_si;                  // [N] eV/Angstrom
double const dyne_si =
             dyne_centimeter_si/centimeter_si;   // [N] dyne
double const hartree_per_bohr_si =
             hartree_si/bohr_si;                 // [N] hartree/bohr
double const picogram_micrometer_per_microsecondsq_si =
             picogram_si*micrometer_si/pow(microsecond_si,2);
                                                 // [N] picogram*micrometer/
                                                 //   microsecond^2
double const attogram_nanometer_per_nanosecondsq_si =
             attogram_si*nanometer_si/pow(nanosecond_si,2);
                                                 // [N] attogram*nanometer/
                                                 //   nanosecond^2

/*----------------------
   Pressure units
------------------------ */
double const atmosphere_si = 101325.0; // [Pa] standard atmosphere (NIST value)
double const bar_si = 1e5;             // [Pa] bar
double const dyne_per_centimetersq_si =
             dyne_centimeter_si/pow(centimeter_si,3);
                                       // [Pa] dyne/centimeter^2
double const picogram_per_micrometer_microsecondsq_si =
             picogram_si/(micrometer_si*pow(microsecond_si,2));
                                       // [Pa] picogram/(micrometer*
                                       //   microsecond^2)
double const attogram_per_nanometer_nanosecondsq_si =
             attogram_si/(nanometer_si*pow(nanosecond_si,2));
                                       // [Pa] attogram/(nanometer*nanosecond^2)

/*----------------------
   Viscosity units
------------------------ */
double const poise_si = 0.1;                  // [Pa*s] Poise
double const amu_per_bohr_femtosecond_si =
             amu_si/(bohr_si*femtosecond_si); // [Pa*s] amu/(bohr*femtosecond)
double const picogram_per_micrometer_microsecond_si =
             picogram_si/(micrometer_si*microsecond_si);
                                              // [Pa*s] picogram/(micrometer*
                                              //   microsecond)
double const attogram_per_nanometer_nanosecond_si =
             attogram_si/(nanometer_si*nanosecond_si);
                                              // [Pa*s] attogram/(nanometer*
                                              //   nanosecond)

/*----------------------
   Charge units
------------------------ */
double const echarge_si = e_si;                   // [C] electron charge unit
double const statcoulomb_si = e_si/4.8032044e-10; // [C] Statcoulomb or esu
                                                  //   (value from LAMMPS units
                                                  //   documentation)
double const picocoulomb_si = 1e-12;              // [C] picocoulomb

/*----------------------
   Dipole units
------------------------ */
double const electron_angstrom_si = echarge_si*angstrom_si;
                                                 // [C*m] electron*angstrom
double const statcoulomb_centimeter_si = statcoulomb_si*centimeter_si;
                                                 // [C*m] statcoulomb*centimeter
double const debye_si = 1e-18*statcoulomb_centimeter_si;
                                                 // [C*m] Debye
double const picocoulomb_micrometer_si = picocoulomb_si*micrometer_si;
                                                 // [C*m] picocoulomb*micrometer
double const electron_nanometer_si = echarge_si*nanometer_si;
                                                 // [C*m] electron*nanometer

/*----------------------
   Electric field units
------------------------ */
double const volt_per_angstrom_si = 1.0/angstrom_si;// [V/m] volt/angstrom
double const statvolt_per_centimeter_si =
             erg_si/(statcoulomb_si*centimeter_si); // [V/m] statvolt/centimeter
double const volt_per_centimeter_si =
             1.0/centimeter_si;                     // [V/m] volt/centimeter
double const volt_per_micrometer_si =
             1.0/micrometer_si;                     // [V/m] volt/micrometer
double const volt_per_nanometer_si =
             1.0/nanometer_si;                      // [V/m] volt/nanometer

// Define enumerations
enum sys_type
{
   real      = 1,
   metal     = 2,
   si        = 3,
   cgs       = 4,
   electron  = 5,
   micro     = 6,
   nano      = 7
};

enum unit_type
{
    mass         = 1,
    distance     = 2,
    time         = 3,
    energy       = 4,
    velocity     = 5,
    force        = 6,
    torque       = 7,
    temperature  = 8,
    pressure     = 9,
    viscosity    = 10,
    charge       = 11,
    dipole       = 12,
    efield       = 13,
    density      = 14
};

enum units
{
    // mass
    gram_per_mole = 101,
    kilogram      = 102,
    gram          = 103,
    amu           = 104,
    picogram      = 105,
    attogram      = 106,
    // distance
    angstrom    = 201,
    meter       = 202,
    centimeter  = 203,
    bohr        = 204,
    micrometer  = 205,
    nanometer   = 206,
    // time
    femtosecond  = 301,
    picosecond   = 302,
    second       = 303,
    microsecond  = 304,
    nanosecond   = 305,
    // energy
    kcal_per_mole = 401,
    ev            = 402,
    joule         = 403,
    erg           = 404,
    hartree       = 405,
    picogram_micrometersq_per_microsecondsq = 406,
    attogram_nanometersq_per_nanosecondsq   = 407,
    // velocity
    angstrom_per_femtosecond   = 501,
    angstrom_per_picosecond    = 502,
    meter_per_second           = 503,
    centimeter_per_second      = 504,
    bohr_per_atu               = 505,
    micrometer_per_microsecond = 506,
    nanometer_per_nanosecond   = 507,
    // force
    kcal_per_mole_angstrom = 601,
    ev_per_angstrom        = 602,
    newton                 = 603,
    dyne                   = 604,
    hartree_per_bohr       = 605,
    picogram_micrometer_per_microsecondsq = 606,
    attogram_nanometer_per_nanosecondsq   = 607,
    // torque
    newton_meter    = 701,
    dyne_centimeter = 702,
    // temperature
    kelvin         = 801,
    // pressure
    atmosphere = 901,
    bar        = 902,
    pascal     = 903,
    dyne_per_centimetersq                 = 904,
    picogram_per_micrometer_microsecondsq = 905,
    attogram_per_nanometer_nanosecondsq   = 906,
    // viscosity
    poise                    = 1001,
    pascal_second            = 1002,
    amu_per_bohr_femtosecond = 1003, // electron system, not in docs, GUESSED
    picogram_per_micrometer_microsecond = 1004,
    attogram_per_nanometer_nanosecond   = 1005,
    // charge
    echarge     = 1101,
    coulomb     = 1102,
    statcoulomb = 1103,
    picocoulomb = 1104,
    // dipole
    electron_angstrom      = 1201,
    coulomb_meter          = 1202,
    statcoulomb_centimeter = 1203,
    debye                  = 1204,
    picocoulomb_micrometer = 1205,
    electron_nanometer     = 1206,
    // electric field
    volt_per_angstrom       = 1301,
    volt_per_meter          = 1302,
    statvolt_per_centimeter = 1303,
    volt_per_centimeter     = 1304,
    volt_per_micrometer     = 1305,
    volt_per_nanometer      = 1306,
    // density
    gram_per_centimetercu     = 1401,
    kilogram_per_metercu      = 1402,
    amu_per_bohrcu            = 1403, // electron system, not in docs, GUESSED
    picogram_per_micrometercu = 1404,
    attogram_per_nanometercu  = 1405
};

// Define dictionaries
map<string, sys_type> system_dic;
map<string, unit_type> unit_dic;
map<unit_type, units> units_real_dic;
map<unit_type, units> units_metal_dic;
map<unit_type, units> units_si_dic;
map<unit_type, units> units_cgs_dic;
map<unit_type, units> units_electron_dic;
map<unit_type, units> units_micro_dic;
map<unit_type, units> units_nano_dic;

/* ---------------------------------------------------------------------- */

void initialize_dictionaries()
{
  system_dic["real"]     = real;
  system_dic["metal"]    = metal;
  system_dic["si"]       = si;
  system_dic["cgs"]      = cgs;
  system_dic["electron"] = electron;
  system_dic["micro"]    = micro;
  system_dic["nano"]     = nano;

  unit_dic["mass"]       = mass;
  unit_dic["distance"]   = distance;
  unit_dic["time"]       = time;
  unit_dic["energy"]     = energy;
  unit_dic["velocity"]   = velocity;
  unit_dic["force"]      = force;
  unit_dic["torque"]     = torque;
  unit_dic["temperature"]= temperature;
  unit_dic["pressure"]   = pressure;
  unit_dic["viscosity"]  = viscosity;
  unit_dic["charge"]     = charge;
  unit_dic["dipole"]     = dipole;
  unit_dic["efield"]     = efield;
  unit_dic["density"]    = density;

  units_real_dic[mass]       = gram_per_mole;
  units_real_dic[distance]   = angstrom;
  units_real_dic[time]       = femtosecond;
  units_real_dic[energy]     = kcal_per_mole;
  units_real_dic[velocity]   = angstrom_per_femtosecond;
  units_real_dic[force]      = kcal_per_mole_angstrom;
  units_real_dic[torque]     = kcal_per_mole;
  units_real_dic[temperature]= kelvin;
  units_real_dic[pressure]   = atmosphere;
  units_real_dic[viscosity]  = poise;
  units_real_dic[charge]     = echarge;
  units_real_dic[dipole]     = electron_angstrom;
  units_real_dic[efield]     = volt_per_angstrom;
  units_real_dic[density]    = gram_per_centimetercu;

  units_metal_dic[mass]       = gram_per_mole;
  units_metal_dic[distance]   = angstrom;
  units_metal_dic[time]       = picosecond;
  units_metal_dic[energy]     = ev;
  units_metal_dic[velocity]   = angstrom_per_picosecond;
  units_metal_dic[force]      = ev_per_angstrom;
  units_metal_dic[torque]     = ev;
  units_metal_dic[temperature]= kelvin;
  units_metal_dic[pressure]   = bar;
  units_metal_dic[viscosity]  = poise;
  units_metal_dic[charge]     = echarge;
  units_metal_dic[dipole]     = electron_angstrom;
  units_metal_dic[efield]     = volt_per_angstrom;
  units_metal_dic[density]    = gram_per_centimetercu;

  units_si_dic[mass]       = kilogram;
  units_si_dic[distance]   = meter;
  units_si_dic[time]       = second;
  units_si_dic[energy]     = joule;
  units_si_dic[velocity]   = meter_per_second;
  units_si_dic[force]      = newton;
  units_si_dic[torque]     = newton_meter;
  units_si_dic[temperature]= kelvin;
  units_si_dic[pressure]   = pascal;
  units_si_dic[viscosity]  = pascal_second;
  units_si_dic[charge]     = coulomb;
  units_si_dic[dipole]     = coulomb_meter;
  units_si_dic[efield]     = volt_per_meter;
  units_si_dic[density]    = kilogram_per_metercu;

  units_cgs_dic[mass]       = gram;
  units_cgs_dic[distance]   = centimeter;
  units_cgs_dic[time]       = second;
  units_cgs_dic[energy]     = erg;
  units_cgs_dic[velocity]   = centimeter_per_second;
  units_cgs_dic[force]      = dyne;
  units_cgs_dic[torque]     = dyne_centimeter;
  units_cgs_dic[temperature]= kelvin;
  units_cgs_dic[pressure]   = dyne_per_centimetersq;
  units_cgs_dic[viscosity]  = poise;
  units_cgs_dic[charge]     = statcoulomb;
  units_cgs_dic[dipole]     = statcoulomb_centimeter;
  units_cgs_dic[efield]     = statvolt_per_centimeter;
  units_cgs_dic[density]    = gram_per_centimetercu;

  units_electron_dic[mass]       = amu;
  units_electron_dic[distance]   = bohr;
  units_electron_dic[time]       = femtosecond;
  units_electron_dic[energy]     = hartree;
  units_electron_dic[velocity]   = bohr_per_atu;
  units_electron_dic[force]      = hartree_per_bohr;
  units_electron_dic[torque]     = hartree;              // unknown, GUESSED
  units_electron_dic[temperature]= kelvin;
  units_electron_dic[pressure]   = pascal;
  units_electron_dic[viscosity]  = pascal_second;        // unknown, GUESSED
  units_electron_dic[charge]     = echarge;
  units_electron_dic[dipole]     = debye;
  units_electron_dic[efield]     = volt_per_centimeter;
  units_electron_dic[density]    = amu_per_bohrcu;       // unknown, GUESSED

  units_micro_dic[mass]       = picogram;
  units_micro_dic[distance]   = micrometer;
  units_micro_dic[time]       = microsecond;
  units_micro_dic[energy]     = picogram_micrometersq_per_microsecondsq;
  units_micro_dic[velocity]   = micrometer_per_microsecond;
  units_micro_dic[force]      = picogram_micrometer_per_microsecondsq;
  units_micro_dic[torque]     = picogram_micrometersq_per_microsecondsq;
  units_micro_dic[temperature]= kelvin;
  units_micro_dic[pressure]   = picogram_per_micrometer_microsecondsq;
  units_micro_dic[viscosity]  = picogram_per_micrometer_microsecond;
  units_micro_dic[charge]     = picocoulomb;
  units_micro_dic[dipole]     = picocoulomb_micrometer;
  units_micro_dic[efield]     = volt_per_micrometer;
  units_micro_dic[density]    = picogram_per_micrometercu;

  units_nano_dic[mass]       = attogram;
  units_nano_dic[distance]   = nanometer;
  units_nano_dic[time]       = nanosecond;
  units_nano_dic[energy]     = attogram_nanometersq_per_nanosecondsq;
  units_nano_dic[velocity]   = nanometer_per_nanosecond;
  units_nano_dic[force]      = attogram_nanometer_per_nanosecondsq;
  units_nano_dic[torque]     = attogram_nanometersq_per_nanosecondsq;
  units_nano_dic[temperature]= kelvin;
  units_nano_dic[pressure]   = attogram_per_nanometer_nanosecondsq;
  units_nano_dic[viscosity]  = attogram_per_nanometer_nanosecond;
  units_nano_dic[charge]     = echarge;
  units_nano_dic[dipole]     = electron_nanometer;
  units_nano_dic[efield]     = volt_per_nanometer;
  units_nano_dic[density]    = attogram_per_nanometercu;

}

/* ---------------------------------------------------------------------- */

// Get the enumeration for the unit of type `unit_type_enum`
// for LAMMPS system `system_enum`.
units get_lammps_system_unit(sys_type system_enum, unit_type unit_type_enum)
{
  switch(system_enum) {
    case real :
      return units_real_dic[unit_type_enum];
    case metal :
      return units_metal_dic[unit_type_enum];
    case si :
      return units_si_dic[unit_type_enum];
    case cgs :
      return units_cgs_dic[unit_type_enum];
    case electron :
      return units_electron_dic[unit_type_enum];
    case micro :
      return units_micro_dic[unit_type_enum];
    case nano :
    default :    // This is here to a prevent a compiler warning
      return units_nano_dic[unit_type_enum];
  }
}

/* ---------------------------------------------------------------------- */

// Mass conversion
double get_mass_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[kilogram][kilogram]           = 1.0;
  conv[kilogram][gram_per_mole]      = 1.0/gram_per_mole_si;
  conv[kilogram][gram]               = 1.0/gram_si;
  conv[kilogram][amu]                = 1.0/amu_si;
  conv[kilogram][picogram]           = 1.0/picogram_si;
  conv[kilogram][attogram]           = 1.0/attogram_si;

  to_si = 1.0/conv[kilogram][gram_per_mole];
  conv[gram_per_mole][kilogram]      = to_si*conv[kilogram][kilogram];
  conv[gram_per_mole][gram_per_mole] = 1.0;
  conv[gram_per_mole][gram]          = to_si*conv[kilogram][gram];
  conv[gram_per_mole][amu]           = to_si*conv[kilogram][amu];
  conv[gram_per_mole][picogram]      = to_si*conv[kilogram][picogram];
  conv[gram_per_mole][attogram]      = to_si*conv[kilogram][attogram];

  to_si = 1.0/conv[kilogram][gram];
  conv[gram][kilogram]               = to_si*conv[kilogram][kilogram];
  conv[gram][gram_per_mole]          = to_si*conv[kilogram][gram_per_mole];
  conv[gram][gram]                   = 1.0;
  conv[gram][amu]                    = to_si*conv[kilogram][amu];
  conv[gram][picogram]               = to_si*conv[kilogram][picogram];
  conv[gram][attogram]               = to_si*conv[kilogram][attogram];

  to_si = 1.0/conv[kilogram][amu];
  conv[amu][kilogram]                = to_si*conv[kilogram][kilogram];
  conv[amu][gram_per_mole]           = to_si*conv[kilogram][gram_per_mole];
  conv[amu][gram]                    = to_si*conv[kilogram][gram];
  conv[amu][amu]                     = 1.0;
  conv[amu][picogram]                = to_si*conv[kilogram][picogram];
  conv[amu][attogram]                = to_si*conv[kilogram][attogram];

  to_si = 1.0/conv[kilogram][picogram];
  conv[picogram][kilogram]           = to_si*conv[kilogram][kilogram];
  conv[picogram][gram_per_mole]      = to_si*conv[kilogram][gram_per_mole];
  conv[picogram][gram]               = to_si*conv[kilogram][gram];
  conv[picogram][amu]                = to_si*conv[kilogram][amu];
  conv[picogram][picogram]           = 1.0;
  conv[picogram][attogram]           = to_si*conv[kilogram][attogram];

  to_si = 1.0/conv[kilogram][attogram];
  conv[attogram][kilogram]           = to_si*conv[kilogram][kilogram];
  conv[attogram][gram_per_mole]      = to_si*conv[kilogram][gram_per_mole];
  conv[attogram][gram]               = to_si*conv[kilogram][gram];
  conv[attogram][amu]                = to_si*conv[kilogram][amu];
  conv[attogram][picogram]           = to_si*conv[kilogram][picogram];
  conv[attogram][attogram]           = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Distance conversion
double get_distance_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[meter][meter]           = 1.0;
  conv[meter][angstrom]        = 1.0/angstrom_si;
  conv[meter][centimeter]      = 1.0/centimeter_si;
  conv[meter][bohr]            = 1.0/bohr_si;
  conv[meter][micrometer]      = 1.0/micrometer_si;
  conv[meter][nanometer]       = 1.0/nanometer_si;

  to_si = 1.0/conv[meter][angstrom];
  conv[angstrom][meter]        = to_si*conv[meter][meter];
  conv[angstrom][angstrom]     = 1.0;
  conv[angstrom][centimeter]   = to_si*conv[meter][centimeter];
  conv[angstrom][bohr]         = to_si*conv[meter][bohr];
  conv[angstrom][micrometer]   = to_si*conv[meter][micrometer];
  conv[angstrom][nanometer]    = to_si*conv[meter][nanometer];

  to_si = 1.0/conv[meter][centimeter];
  conv[centimeter][meter]      = to_si*conv[meter][meter];
  conv[centimeter][angstrom]   = to_si*conv[meter][angstrom];
  conv[centimeter][centimeter] = 1.0;
  conv[centimeter][bohr]       = to_si*conv[meter][bohr];
  conv[centimeter][micrometer] = to_si*conv[meter][micrometer];
  conv[centimeter][nanometer]  = to_si*conv[meter][nanometer];

  to_si = 1.0/conv[meter][bohr];
  conv[bohr][meter]            = to_si*conv[meter][meter];
  conv[bohr][angstrom]         = to_si*conv[meter][angstrom];
  conv[bohr][centimeter]       = to_si*conv[meter][centimeter];
  conv[bohr][bohr]             = 1.0;
  conv[bohr][micrometer]       = to_si*conv[meter][micrometer];
  conv[bohr][nanometer]        = to_si*conv[meter][nanometer];

  to_si = 1.0/conv[meter][micrometer];
  conv[micrometer][meter]      = to_si*conv[meter][meter];
  conv[micrometer][angstrom]   = to_si*conv[meter][angstrom];
  conv[micrometer][centimeter] = to_si*conv[meter][centimeter];
  conv[micrometer][bohr]       = to_si*conv[meter][bohr];
  conv[micrometer][micrometer] = 1.0;
  conv[micrometer][nanometer]  = to_si*conv[meter][nanometer];

  to_si = 1.0/conv[meter][nanometer];
  conv[nanometer][meter]       = to_si*conv[meter][meter];
  conv[nanometer][angstrom]    = to_si*conv[meter][angstrom];
  conv[nanometer][centimeter]  = to_si*conv[meter][centimeter];
  conv[nanometer][bohr]        = to_si*conv[meter][bohr];
  conv[nanometer][micrometer]  = to_si*conv[meter][micrometer];
  conv[nanometer][nanometer]   = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Time conversion
double get_time_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[second][second]           = 1.0;
  conv[second][femtosecond]      = 1.0/femtosecond_si;
  conv[second][picosecond]       = 1.0/picosecond_si;
  conv[second][microsecond]      = 1.0/microsecond_si;
  conv[second][nanosecond]       = 1.0/nanosecond_si;

  to_si = 1.0/conv[second][femtosecond];
  conv[femtosecond][second]      = to_si*conv[second][second];
  conv[femtosecond][femtosecond] = 1.0;
  conv[femtosecond][picosecond]  = to_si*conv[second][picosecond];
  conv[femtosecond][microsecond] = to_si*conv[second][microsecond];
  conv[femtosecond][nanosecond]  = to_si*conv[second][nanosecond];

  to_si = 1.0/conv[second][picosecond];
  conv[picosecond][second]       = to_si*conv[second][second];
  conv[picosecond][femtosecond]  = to_si*conv[second][femtosecond];
  conv[picosecond][picosecond]   = 1.0;
  conv[picosecond][microsecond]  = to_si*conv[second][microsecond];
  conv[picosecond][nanosecond]   = to_si*conv[second][nanosecond];

  to_si = 1.0/conv[second][microsecond];
  conv[microsecond][second]      = to_si*conv[second][second];
  conv[microsecond][femtosecond] = to_si*conv[second][femtosecond];
  conv[microsecond][picosecond]  = to_si*conv[second][picosecond];
  conv[microsecond][microsecond] = 1.0;
  conv[microsecond][nanosecond]  = to_si*conv[second][nanosecond];

  to_si = 1.0/conv[second][nanosecond];
  conv[nanosecond][second]       = to_si*conv[second][second];
  conv[nanosecond][femtosecond]  = to_si*conv[second][femtosecond];
  conv[nanosecond][picosecond]   = to_si*conv[second][picosecond];
  conv[nanosecond][microsecond]  = to_si*conv[second][microsecond];
  conv[nanosecond][nanosecond]   = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Energy conversion
double get_energy_conversion_factor(units from_unit_enum, units to_unit_enum)
{

  map<units, map<units, double> > conv;
  double to_si;

  conv[joule][joule]                                                                     = 1.0;
  conv[joule][kcal_per_mole]                                                             = 1.0/kcal_per_mole_si;
  conv[joule][ev]                                                                        = 1.0/ev_si;
  conv[joule][erg]                                                                       = 1.0/erg_si;
  conv[joule][hartree]                                                                   = 1.0/hartree_si;
  conv[joule][picogram_micrometersq_per_microsecondsq]                                   = 1.0/picogram_micrometersq_per_microsecondsq_si;
  conv[joule][attogram_nanometersq_per_nanosecondsq]                                     = 1.0/attogram_nanometersq_per_nanosecondsq_si;

  to_si = 1.0/conv[joule][kcal_per_mole];
  conv[kcal_per_mole][joule]                                                             = to_si*conv[joule][joule];
  conv[kcal_per_mole][kcal_per_mole]                                                     = 1.0;
  conv[kcal_per_mole][ev]                                                                = to_si*conv[joule][ev];
  conv[kcal_per_mole][erg]                                                               = to_si*conv[joule][erg];
  conv[kcal_per_mole][hartree]                                                           = to_si*conv[joule][hartree];
  conv[kcal_per_mole][picogram_micrometersq_per_microsecondsq]                           = to_si*conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[kcal_per_mole][attogram_nanometersq_per_nanosecondsq]                             = to_si*conv[joule][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[joule][ev];
  conv[ev][joule]                                                                        = to_si*conv[joule][joule];
  conv[ev][kcal_per_mole]                                                                = to_si*conv[joule][kcal_per_mole];
  conv[ev][ev]                                                                           = 1.0;
  conv[ev][erg]                                                                          = to_si*conv[joule][erg];
  conv[ev][hartree]                                                                      = to_si*conv[joule][hartree];
  conv[ev][picogram_micrometersq_per_microsecondsq]                                      = to_si*conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[ev][attogram_nanometersq_per_nanosecondsq]                                        = to_si*conv[joule][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[joule][erg];
  conv[erg][joule]                                                                       = to_si*conv[joule][joule];
  conv[erg][kcal_per_mole]                                                               = to_si*conv[joule][kcal_per_mole];
  conv[erg][ev]                                                                          = to_si*conv[joule][ev];
  conv[erg][erg]                                                                         = 1.0;
  conv[erg][hartree]                                                                     = to_si*conv[joule][hartree];
  conv[erg][picogram_micrometersq_per_microsecondsq]                                     = to_si*conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[erg][attogram_nanometersq_per_nanosecondsq]                                       = to_si*conv[joule][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[joule][hartree];
  conv[hartree][joule]                                                                   = to_si*conv[joule][joule];
  conv[hartree][kcal_per_mole]                                                           = to_si*conv[joule][kcal_per_mole];
  conv[hartree][ev]                                                                      = to_si*conv[joule][ev];
  conv[hartree][erg]                                                                     = to_si*conv[joule][erg];
  conv[hartree][hartree]                                                                 = 1.0;
  conv[hartree][picogram_micrometersq_per_microsecondsq]                                 = to_si*conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[hartree][attogram_nanometersq_per_nanosecondsq]                                   = to_si*conv[joule][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[picogram_micrometersq_per_microsecondsq][joule]                                   = to_si*conv[joule][joule];
  conv[picogram_micrometersq_per_microsecondsq][kcal_per_mole]                           = to_si*conv[joule][kcal_per_mole];
  conv[picogram_micrometersq_per_microsecondsq][ev]                                      = to_si*conv[joule][ev];
  conv[picogram_micrometersq_per_microsecondsq][erg]                                     = to_si*conv[joule][erg];
  conv[picogram_micrometersq_per_microsecondsq][hartree]                                 = to_si*conv[joule][hartree];
  conv[picogram_micrometersq_per_microsecondsq][picogram_micrometersq_per_microsecondsq] = 1.0;
  conv[picogram_micrometersq_per_microsecondsq][attogram_nanometersq_per_nanosecondsq]   = to_si*conv[joule][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[joule][attogram_nanometersq_per_nanosecondsq];
  conv[attogram_nanometersq_per_nanosecondsq][joule]                                     = to_si*conv[joule][joule];
  conv[attogram_nanometersq_per_nanosecondsq][kcal_per_mole]                             = to_si*conv[joule][kcal_per_mole];
  conv[attogram_nanometersq_per_nanosecondsq][ev]                                        = to_si*conv[joule][ev];
  conv[attogram_nanometersq_per_nanosecondsq][erg]                                       = to_si*conv[joule][erg];
  conv[attogram_nanometersq_per_nanosecondsq][hartree]                                   = to_si*conv[joule][hartree];
  conv[attogram_nanometersq_per_nanosecondsq][picogram_micrometersq_per_microsecondsq]   = to_si*conv[joule][picogram_micrometersq_per_microsecondsq];
  conv[attogram_nanometersq_per_nanosecondsq][attogram_nanometersq_per_nanosecondsq]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Velocity conversion
double get_velocity_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[meter_per_second][meter_per_second]                     = 1.0;
  conv[meter_per_second][angstrom_per_femtosecond]             = 1.0/angstrom_per_femtosecond_si;
  conv[meter_per_second][angstrom_per_picosecond]              = 1.0/angstrom_per_picosecond_si;
  conv[meter_per_second][centimeter_per_second]                = 1.0/centimeter_per_second_si;
  conv[meter_per_second][bohr_per_atu]                         = 1.0/bohr_per_atu_electron_si;
  conv[meter_per_second][micrometer_per_microsecond]           = 1.0/micrometer_per_microsecond_si;
  conv[meter_per_second][nanometer_per_nanosecond]             = 1.0/nanometer_per_nanosecond_si;

  to_si = 1.0/conv[meter_per_second][angstrom_per_femtosecond];
  conv[angstrom_per_femtosecond][meter_per_second]             = to_si*conv[meter_per_second][meter_per_second];
  conv[angstrom_per_femtosecond][angstrom_per_femtosecond]     = 1.0;
  conv[angstrom_per_femtosecond][angstrom_per_picosecond]      = to_si*conv[meter_per_second][angstrom_per_picosecond];
  conv[angstrom_per_femtosecond][centimeter_per_second]        = to_si*conv[meter_per_second][centimeter_per_second];
  conv[angstrom_per_femtosecond][bohr_per_atu]                 = to_si*conv[meter_per_second][bohr_per_atu];
  conv[angstrom_per_femtosecond][micrometer_per_microsecond]   = to_si*conv[meter_per_second][micrometer_per_microsecond];
  conv[angstrom_per_femtosecond][nanometer_per_nanosecond]     = to_si*conv[meter_per_second][nanometer_per_nanosecond];

  to_si = 1.0/conv[meter_per_second][angstrom_per_picosecond];
  conv[angstrom_per_picosecond][meter_per_second]              = to_si*conv[meter_per_second][meter_per_second];
  conv[angstrom_per_picosecond][angstrom_per_femtosecond]      = to_si*conv[meter_per_second][angstrom_per_femtosecond];
  conv[angstrom_per_picosecond][angstrom_per_picosecond]       = 1.0;
  conv[angstrom_per_picosecond][centimeter_per_second]         = to_si*conv[meter_per_second][centimeter_per_second];
  conv[angstrom_per_picosecond][bohr_per_atu]                  = to_si*conv[meter_per_second][bohr_per_atu];
  conv[angstrom_per_picosecond][micrometer_per_microsecond]    = to_si*conv[meter_per_second][micrometer_per_microsecond];
  conv[angstrom_per_picosecond][nanometer_per_nanosecond]      = to_si*conv[meter_per_second][nanometer_per_nanosecond];

  to_si = 1.0/conv[meter_per_second][centimeter_per_second];
  conv[centimeter_per_second][meter_per_second]                = to_si*conv[meter_per_second][meter_per_second];
  conv[centimeter_per_second][angstrom_per_femtosecond]        = to_si*conv[meter_per_second][angstrom_per_femtosecond];
  conv[centimeter_per_second][angstrom_per_picosecond]         = to_si*conv[meter_per_second][angstrom_per_picosecond];
  conv[centimeter_per_second][centimeter_per_second]           = 1.0;
  conv[centimeter_per_second][bohr_per_atu]                    = to_si*conv[meter_per_second][bohr_per_atu];
  conv[centimeter_per_second][micrometer_per_microsecond]      = to_si*conv[meter_per_second][micrometer_per_microsecond];
  conv[centimeter_per_second][nanometer_per_nanosecond]        = to_si*conv[meter_per_second][nanometer_per_nanosecond];

  to_si = 1.0/conv[meter_per_second][bohr_per_atu];
  conv[bohr_per_atu][meter_per_second]                         = to_si*conv[meter_per_second][meter_per_second];
  conv[bohr_per_atu][angstrom_per_femtosecond]                 = to_si*conv[meter_per_second][angstrom_per_femtosecond];
  conv[bohr_per_atu][angstrom_per_picosecond]                  = to_si*conv[meter_per_second][angstrom_per_picosecond];
  conv[bohr_per_atu][centimeter_per_second]                    = to_si*conv[meter_per_second][centimeter_per_second];
  conv[bohr_per_atu][bohr_per_atu]                             = 1.0;
  conv[bohr_per_atu][micrometer_per_microsecond]               = to_si*conv[meter_per_second][micrometer_per_microsecond];
  conv[bohr_per_atu][nanometer_per_nanosecond]                 = to_si*conv[meter_per_second][nanometer_per_nanosecond];

  to_si = 1.0/conv[meter_per_second][micrometer_per_microsecond];
  conv[micrometer_per_microsecond][meter_per_second]           = to_si*conv[meter_per_second][meter_per_second];
  conv[micrometer_per_microsecond][angstrom_per_femtosecond]   = to_si*conv[meter_per_second][angstrom_per_femtosecond];
  conv[micrometer_per_microsecond][angstrom_per_picosecond]    = to_si*conv[meter_per_second][angstrom_per_picosecond];
  conv[micrometer_per_microsecond][centimeter_per_second]      = to_si*conv[meter_per_second][centimeter_per_second];
  conv[micrometer_per_microsecond][bohr_per_atu]               = to_si*conv[meter_per_second][bohr_per_atu];
  conv[micrometer_per_microsecond][micrometer_per_microsecond] = 1.0;
  conv[micrometer_per_microsecond][nanometer_per_nanosecond]   = to_si*conv[meter_per_second][nanometer_per_nanosecond];

  to_si = 1.0/conv[meter_per_second][nanometer_per_nanosecond];
  conv[nanometer_per_nanosecond][meter_per_second]             = to_si*conv[meter_per_second][meter_per_second];
  conv[nanometer_per_nanosecond][angstrom_per_femtosecond]     = to_si*conv[meter_per_second][angstrom_per_femtosecond];
  conv[nanometer_per_nanosecond][angstrom_per_picosecond]      = to_si*conv[meter_per_second][angstrom_per_picosecond];
  conv[nanometer_per_nanosecond][centimeter_per_second]        = to_si*conv[meter_per_second][centimeter_per_second];
  conv[nanometer_per_nanosecond][bohr_per_atu]                 = to_si*conv[meter_per_second][bohr_per_atu];
  conv[nanometer_per_nanosecond][micrometer_per_microsecond]   = to_si*conv[meter_per_second][micrometer_per_microsecond];
  conv[nanometer_per_nanosecond][nanometer_per_nanosecond]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Force conversion
double get_force_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[newton][newton]                                                               = 1.0;
  conv[newton][kcal_per_mole_angstrom]                                               = 1.0/kcal_per_mole_angstrom_si;
  conv[newton][ev_per_angstrom]                                                      = 1.0/ev_per_angstrom_si;
  conv[newton][dyne]                                                                 = 1.0/dyne_si;
  conv[newton][hartree_per_bohr]                                                     = 1.0/hartree_per_bohr_si;
  conv[newton][picogram_micrometer_per_microsecondsq]                                = 1.0/picogram_micrometer_per_microsecondsq_si;
  conv[newton][attogram_nanometer_per_nanosecondsq]                                  = 1.0/attogram_nanometer_per_nanosecondsq_si;

  to_si = 1.0/conv[newton][kcal_per_mole_angstrom];
  conv[kcal_per_mole_angstrom][newton]                                               = to_si*conv[newton][newton];
  conv[kcal_per_mole_angstrom][kcal_per_mole_angstrom]                               = 1.0;
  conv[kcal_per_mole_angstrom][ev_per_angstrom]                                      = to_si*conv[newton][ev_per_angstrom];
  conv[kcal_per_mole_angstrom][dyne]                                                 = to_si*conv[newton][dyne];
  conv[kcal_per_mole_angstrom][hartree_per_bohr]                                     = to_si*conv[newton][hartree_per_bohr];
  conv[kcal_per_mole_angstrom][picogram_micrometer_per_microsecondsq]                = to_si*conv[newton][picogram_micrometer_per_microsecondsq];
  conv[kcal_per_mole_angstrom][attogram_nanometer_per_nanosecondsq]                  = to_si*conv[newton][attogram_nanometer_per_nanosecondsq];

  to_si = 1.0/conv[newton][ev_per_angstrom];
  conv[ev_per_angstrom][newton]                                                      = to_si*conv[newton][newton];
  conv[ev_per_angstrom][kcal_per_mole_angstrom]                                      = to_si*conv[newton][kcal_per_mole_angstrom];
  conv[ev_per_angstrom][ev_per_angstrom]                                             = 1.0;
  conv[ev_per_angstrom][dyne]                                                        = to_si*conv[newton][dyne];
  conv[ev_per_angstrom][hartree_per_bohr]                                            = to_si*conv[newton][hartree_per_bohr];
  conv[ev_per_angstrom][picogram_micrometer_per_microsecondsq]                       = to_si*conv[newton][picogram_micrometer_per_microsecondsq];
  conv[ev_per_angstrom][attogram_nanometer_per_nanosecondsq]                         = to_si*conv[newton][attogram_nanometer_per_nanosecondsq];

  to_si = 1.0/conv[newton][dyne];
  conv[dyne][newton]                                                                 = to_si*conv[newton][newton];
  conv[dyne][kcal_per_mole_angstrom]                                                 = to_si*conv[newton][kcal_per_mole_angstrom];
  conv[dyne][ev_per_angstrom]                                                        = to_si*conv[newton][ev_per_angstrom];
  conv[dyne][dyne]                                                                   = 1.0;
  conv[dyne][hartree_per_bohr]                                                       = to_si*conv[newton][hartree_per_bohr];
  conv[dyne][picogram_micrometer_per_microsecondsq]                                  = to_si*conv[newton][picogram_micrometer_per_microsecondsq];
  conv[dyne][attogram_nanometer_per_nanosecondsq]                                    = to_si*conv[newton][attogram_nanometer_per_nanosecondsq];

  to_si = 1.0/conv[newton][hartree_per_bohr];
  conv[hartree_per_bohr][newton]                                                     = to_si*conv[newton][newton];
  conv[hartree_per_bohr][kcal_per_mole_angstrom]                                     = to_si*conv[newton][kcal_per_mole_angstrom];
  conv[hartree_per_bohr][ev_per_angstrom]                                            = to_si*conv[newton][ev_per_angstrom];
  conv[hartree_per_bohr][dyne]                                                       = to_si*conv[newton][dyne];
  conv[hartree_per_bohr][hartree_per_bohr]                                           = 1.0;
  conv[hartree_per_bohr][picogram_micrometer_per_microsecondsq]                      = to_si*conv[newton][picogram_micrometer_per_microsecondsq];
  conv[hartree_per_bohr][attogram_nanometer_per_nanosecondsq]                        = to_si*conv[newton][attogram_nanometer_per_nanosecondsq];

  to_si = 1.0/conv[newton][picogram_micrometer_per_microsecondsq];
  conv[picogram_micrometer_per_microsecondsq][newton]                                = to_si*conv[newton][newton];
  conv[picogram_micrometer_per_microsecondsq][kcal_per_mole_angstrom]                = to_si*conv[newton][kcal_per_mole_angstrom];
  conv[picogram_micrometer_per_microsecondsq][ev_per_angstrom]                       = to_si*conv[newton][ev_per_angstrom];
  conv[picogram_micrometer_per_microsecondsq][dyne]                                  = to_si*conv[newton][dyne];
  conv[picogram_micrometer_per_microsecondsq][hartree_per_bohr]                      = to_si*conv[newton][hartree_per_bohr];
  conv[picogram_micrometer_per_microsecondsq][picogram_micrometer_per_microsecondsq] = 1.0;
  conv[picogram_micrometer_per_microsecondsq][attogram_nanometer_per_nanosecondsq]   = to_si*conv[newton][attogram_nanometer_per_nanosecondsq];

  to_si = 1.0/conv[newton][attogram_nanometer_per_nanosecondsq];
  conv[attogram_nanometer_per_nanosecondsq][newton]                                  = to_si*conv[newton][newton];
  conv[attogram_nanometer_per_nanosecondsq][kcal_per_mole_angstrom]                  = to_si*conv[newton][kcal_per_mole_angstrom];
  conv[attogram_nanometer_per_nanosecondsq][ev_per_angstrom]                         = to_si*conv[newton][ev_per_angstrom];
  conv[attogram_nanometer_per_nanosecondsq][dyne]                                    = to_si*conv[newton][dyne];
  conv[attogram_nanometer_per_nanosecondsq][hartree_per_bohr]                        = to_si*conv[newton][hartree_per_bohr];
  conv[attogram_nanometer_per_nanosecondsq][picogram_micrometer_per_microsecondsq]   = to_si*conv[newton][picogram_micrometer_per_microsecondsq];
  conv[attogram_nanometer_per_nanosecondsq][attogram_nanometer_per_nanosecondsq]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Torque conversion
double get_torque_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[newton_meter][newton_meter]                                                       = 1.0;
  conv[newton_meter][kcal_per_mole]                                                      = 1.0/kcal_per_mole_si;
  conv[newton_meter][ev]                                                                 = 1.0/ev_si;
  conv[newton_meter][dyne_centimeter]                                                    = 1.0/dyne_centimeter_si;
  conv[newton_meter][hartree]                                                            = 1.0/hartree_si;
  conv[newton_meter][picogram_micrometersq_per_microsecondsq]                            = 1.0/picogram_micrometersq_per_microsecondsq_si;
  conv[newton_meter][attogram_nanometersq_per_nanosecondsq]                              = 1.0/attogram_nanometersq_per_nanosecondsq_si;

  to_si = 1.0/conv[newton_meter][kcal_per_mole];
  conv[kcal_per_mole][newton_meter]                                                      = to_si*conv[newton_meter][newton_meter];
  conv[kcal_per_mole][kcal_per_mole]                                                     = 1.0;
  conv[kcal_per_mole][ev]                                                                = to_si*conv[newton_meter][ev];
  conv[kcal_per_mole][dyne_centimeter]                                                   = to_si*conv[newton_meter][dyne_centimeter];
  conv[kcal_per_mole][hartree]                                                           = to_si*conv[newton_meter][hartree];
  conv[kcal_per_mole][picogram_micrometersq_per_microsecondsq]                           = to_si*conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[kcal_per_mole][attogram_nanometersq_per_nanosecondsq]                             = to_si*conv[newton_meter][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[newton_meter][ev];
  conv[ev][newton_meter]                                                                 = to_si*conv[newton_meter][newton_meter];
  conv[ev][kcal_per_mole]                                                                = to_si*conv[newton_meter][kcal_per_mole];
  conv[ev][ev]                                                                           = 1.0;
  conv[ev][dyne_centimeter]                                                              = to_si*conv[newton_meter][dyne_centimeter];
  conv[ev][hartree]                                                                      = to_si*conv[newton_meter][hartree];
  conv[ev][picogram_micrometersq_per_microsecondsq]                                      = to_si*conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[ev][attogram_nanometersq_per_nanosecondsq]                                        = to_si*conv[newton_meter][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[newton_meter][dyne_centimeter];
  conv[dyne_centimeter][newton_meter]                                                    = to_si*conv[newton_meter][newton_meter];
  conv[dyne_centimeter][kcal_per_mole]                                                   = to_si*conv[newton_meter][kcal_per_mole];
  conv[dyne_centimeter][ev]                                                              = to_si*conv[newton_meter][ev];
  conv[dyne_centimeter][dyne_centimeter]                                                 = 1.0;
  conv[dyne_centimeter][hartree]                                                         = to_si*conv[newton_meter][hartree];
  conv[dyne_centimeter][picogram_micrometersq_per_microsecondsq]                         = to_si*conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[dyne_centimeter][attogram_nanometersq_per_nanosecondsq]                           = to_si*conv[newton_meter][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[newton_meter][hartree];
  conv[hartree][newton_meter]                                                            = to_si*conv[newton_meter][newton_meter];
  conv[hartree][kcal_per_mole]                                                           = to_si*conv[newton_meter][kcal_per_mole];
  conv[hartree][ev]                                                                      = to_si*conv[newton_meter][ev];
  conv[hartree][dyne_centimeter]                                                         = to_si*conv[newton_meter][dyne_centimeter];
  conv[hartree][hartree]                                                                 = 1.0;
  conv[hartree][picogram_micrometersq_per_microsecondsq]                                 = to_si*conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[hartree][attogram_nanometersq_per_nanosecondsq]                                   = to_si*conv[newton_meter][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[picogram_micrometersq_per_microsecondsq][newton_meter]                            = to_si*conv[newton_meter][newton_meter];
  conv[picogram_micrometersq_per_microsecondsq][kcal_per_mole]                           = to_si*conv[newton_meter][kcal_per_mole];
  conv[picogram_micrometersq_per_microsecondsq][ev]                                      = to_si*conv[newton_meter][ev];
  conv[picogram_micrometersq_per_microsecondsq][dyne_centimeter]                         = to_si*conv[newton_meter][dyne_centimeter];
  conv[picogram_micrometersq_per_microsecondsq][hartree]                                 = to_si*conv[newton_meter][hartree];
  conv[picogram_micrometersq_per_microsecondsq][picogram_micrometersq_per_microsecondsq] = 1.0;
  conv[picogram_micrometersq_per_microsecondsq][attogram_nanometersq_per_nanosecondsq]   = to_si*conv[newton_meter][attogram_nanometersq_per_nanosecondsq];

  to_si = 1.0/conv[newton_meter][attogram_nanometersq_per_nanosecondsq];
  conv[attogram_nanometersq_per_nanosecondsq][newton_meter]                              = to_si*conv[newton_meter][newton_meter];
  conv[attogram_nanometersq_per_nanosecondsq][kcal_per_mole]                             = to_si*conv[newton_meter][kcal_per_mole];
  conv[attogram_nanometersq_per_nanosecondsq][ev]                                        = to_si*conv[newton_meter][ev];
  conv[attogram_nanometersq_per_nanosecondsq][dyne_centimeter]                           = to_si*conv[newton_meter][dyne_centimeter];
  conv[attogram_nanometersq_per_nanosecondsq][hartree]                                   = to_si*conv[newton_meter][hartree];
  conv[attogram_nanometersq_per_nanosecondsq][picogram_micrometersq_per_microsecondsq]   = to_si*conv[newton_meter][picogram_micrometersq_per_microsecondsq];
  conv[attogram_nanometersq_per_nanosecondsq][attogram_nanometersq_per_nanosecondsq]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Temperature conversion
double get_temperature_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[kelvin][kelvin] = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Pressure conversion
double get_pressure_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[pascal][pascal]                                                               = 1.0;
  conv[pascal][atmosphere]                                                           = 1.0/atmosphere_si;
  conv[pascal][bar]                                                                  = 1.0/bar_si;
  conv[pascal][dyne_per_centimetersq]                                                = 1.0/dyne_per_centimetersq_si;
  conv[pascal][picogram_per_micrometer_microsecondsq]                                = 1.0/picogram_per_micrometer_microsecondsq_si;
  conv[pascal][attogram_per_nanometer_nanosecondsq]                                  = 1.0/attogram_per_nanometer_nanosecondsq_si;

  to_si = 1.0/conv[pascal][atmosphere];
  conv[atmosphere][pascal]                                                           = to_si*conv[pascal][pascal];
  conv[atmosphere][atmosphere]                                                       = 1.0;
  conv[atmosphere][bar]                                                              = to_si*conv[pascal][bar];
  conv[atmosphere][dyne_per_centimetersq]                                            = to_si*conv[pascal][dyne_per_centimetersq];
  conv[atmosphere][picogram_per_micrometer_microsecondsq]                            = to_si*conv[pascal][picogram_per_micrometer_microsecondsq];
  conv[atmosphere][attogram_per_nanometer_nanosecondsq]                              = to_si*conv[pascal][attogram_per_nanometer_nanosecondsq];

  to_si = 1.0/conv[pascal][bar];
  conv[bar][pascal]                                                                  = to_si*conv[pascal][pascal];
  conv[bar][atmosphere]                                                              = to_si*conv[pascal][atmosphere];
  conv[bar][bar]                                                                     = 1.0;
  conv[bar][dyne_per_centimetersq]                                                   = to_si*conv[pascal][dyne_per_centimetersq];
  conv[bar][picogram_per_micrometer_microsecondsq]                                   = to_si*conv[pascal][picogram_per_micrometer_microsecondsq];
  conv[bar][attogram_per_nanometer_nanosecondsq]                                     = to_si*conv[pascal][attogram_per_nanometer_nanosecondsq];

  to_si = 1.0/conv[pascal][dyne_per_centimetersq];
  conv[dyne_per_centimetersq][pascal]                                                = to_si*conv[pascal][pascal];
  conv[dyne_per_centimetersq][atmosphere]                                            = to_si*conv[pascal][atmosphere];
  conv[dyne_per_centimetersq][bar]                                                   = to_si*conv[pascal][bar];
  conv[dyne_per_centimetersq][dyne_per_centimetersq]                                 = 1.0;
  conv[dyne_per_centimetersq][picogram_per_micrometer_microsecondsq]                 = to_si*conv[pascal][picogram_per_micrometer_microsecondsq];
  conv[dyne_per_centimetersq][attogram_per_nanometer_nanosecondsq]                   = to_si*conv[pascal][attogram_per_nanometer_nanosecondsq];

  to_si = 1.0/conv[pascal][picogram_per_micrometer_microsecondsq];
  conv[picogram_per_micrometer_microsecondsq][pascal]                                = to_si*conv[pascal][pascal];
  conv[picogram_per_micrometer_microsecondsq][atmosphere]                            = to_si*conv[pascal][atmosphere];
  conv[picogram_per_micrometer_microsecondsq][bar]                                   = to_si*conv[pascal][bar];
  conv[picogram_per_micrometer_microsecondsq][dyne_per_centimetersq]                 = to_si*conv[pascal][dyne_per_centimetersq];
  conv[picogram_per_micrometer_microsecondsq][picogram_per_micrometer_microsecondsq] = 1.0;
  conv[picogram_per_micrometer_microsecondsq][attogram_per_nanometer_nanosecondsq]   = to_si*conv[pascal][attogram_per_nanometer_nanosecondsq];

  to_si = 1.0/conv[pascal][attogram_per_nanometer_nanosecondsq];
  conv[attogram_per_nanometer_nanosecondsq][pascal]                                  = to_si*conv[pascal][pascal];
  conv[attogram_per_nanometer_nanosecondsq][atmosphere]                              = to_si*conv[pascal][atmosphere];
  conv[attogram_per_nanometer_nanosecondsq][bar]                                     = to_si*conv[pascal][bar];
  conv[attogram_per_nanometer_nanosecondsq][dyne_per_centimetersq]                   = to_si*conv[pascal][dyne_per_centimetersq];
  conv[attogram_per_nanometer_nanosecondsq][picogram_per_micrometer_microsecondsq]   = to_si*conv[pascal][picogram_per_micrometer_microsecondsq];
  conv[attogram_per_nanometer_nanosecondsq][attogram_per_nanometer_nanosecondsq]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Viscosity conversion
double get_viscosity_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[pascal_second][pascal_second]                                             = 1.0;
  conv[pascal_second][poise]                                                     = 1.0/poise_si;
  conv[pascal_second][amu_per_bohr_femtosecond]                                  = 1.0/amu_per_bohr_femtosecond_si;
  conv[pascal_second][picogram_per_micrometer_microsecond]                       = 1.0/picogram_per_micrometer_microsecond_si;
  conv[pascal_second][attogram_per_nanometer_nanosecond]                         = 1.0/attogram_per_nanometer_nanosecond_si;

  to_si = 1.0/conv[pascal_second][poise];
  conv[poise][pascal_second]                                                     = to_si*conv[pascal_second][pascal_second];
  conv[poise][poise]                                                             = 1.0;
  conv[poise][amu_per_bohr_femtosecond]                                          = to_si*conv[pascal_second][amu_per_bohr_femtosecond];
  conv[poise][picogram_per_micrometer_microsecond]                               = to_si*conv[pascal_second][picogram_per_micrometer_microsecond];
  conv[poise][attogram_per_nanometer_nanosecond]                                 = to_si*conv[pascal_second][attogram_per_nanometer_nanosecond];

  to_si = 1.0/conv[pascal_second][amu_per_bohr_femtosecond];
  conv[amu_per_bohr_femtosecond][pascal_second]                                  = to_si*conv[pascal_second][pascal_second];
  conv[amu_per_bohr_femtosecond][poise]                                          = to_si*conv[pascal_second][poise];
  conv[amu_per_bohr_femtosecond][amu_per_bohr_femtosecond]                       = 1.0;
  conv[amu_per_bohr_femtosecond][picogram_per_micrometer_microsecond]            = to_si*conv[pascal_second][picogram_per_micrometer_microsecond];
  conv[amu_per_bohr_femtosecond][attogram_per_nanometer_nanosecond]              = to_si*conv[pascal_second][attogram_per_nanometer_nanosecond];

  to_si = 1.0/conv[pascal_second][picogram_per_micrometer_microsecond];
  conv[picogram_per_micrometer_microsecond][pascal_second]                       = to_si*conv[pascal_second][pascal_second];
  conv[picogram_per_micrometer_microsecond][poise]                               = to_si*conv[pascal_second][poise];
  conv[picogram_per_micrometer_microsecond][amu_per_bohr_femtosecond]            = to_si*conv[pascal_second][amu_per_bohr_femtosecond];
  conv[picogram_per_micrometer_microsecond][picogram_per_micrometer_microsecond] = 1.0;
  conv[picogram_per_micrometer_microsecond][attogram_per_nanometer_nanosecond]   = to_si*conv[pascal_second][attogram_per_nanometer_nanosecond];

  to_si = 1.0/conv[pascal_second][attogram_per_nanometer_nanosecond];
  conv[attogram_per_nanometer_nanosecond][pascal_second]                         = to_si*conv[pascal_second][pascal_second];
  conv[attogram_per_nanometer_nanosecond][poise]                                 = to_si*conv[pascal_second][poise];
  conv[attogram_per_nanometer_nanosecond][amu_per_bohr_femtosecond]              = to_si*conv[pascal_second][amu_per_bohr_femtosecond];
  conv[attogram_per_nanometer_nanosecond][picogram_per_micrometer_microsecond]   = to_si*conv[pascal_second][picogram_per_micrometer_microsecond];
  conv[attogram_per_nanometer_nanosecond][attogram_per_nanometer_nanosecond]     = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Charge conversion
double get_charge_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[coulomb][coulomb]         = 1.0;
  conv[coulomb][echarge]         = 1.0/echarge_si;
  conv[coulomb][statcoulomb]     = 1.0/statcoulomb_si;
  conv[coulomb][picocoulomb]     = 1.0/picocoulomb_si;

  to_si = 1.0/conv[coulomb][echarge];
  conv[echarge][coulomb]         = to_si*conv[coulomb][coulomb];
  conv[echarge][echarge]         = 1.0;
  conv[echarge][statcoulomb]     = to_si*conv[coulomb][statcoulomb];
  conv[echarge][picocoulomb]     = to_si*conv[coulomb][picocoulomb];

  to_si = 1.0/conv[coulomb][statcoulomb];
  conv[statcoulomb][coulomb]     = to_si*conv[coulomb][coulomb];
  conv[statcoulomb][echarge]     = to_si*conv[coulomb][echarge];
  conv[statcoulomb][statcoulomb] = 1.0;
  conv[statcoulomb][picocoulomb] = to_si*conv[coulomb][picocoulomb];

  to_si = 1.0/conv[coulomb][picocoulomb];
  conv[picocoulomb][coulomb]     = to_si*conv[coulomb][coulomb];
  conv[picocoulomb][echarge]     = to_si*conv[coulomb][echarge];
  conv[picocoulomb][statcoulomb] = to_si*conv[coulomb][statcoulomb];
  conv[picocoulomb][picocoulomb] = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Dipole conversion
double get_dipole_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[coulomb_meter][coulomb_meter]                   = 1.0;
  conv[coulomb_meter][electron_angstrom]               = 1.0/electron_angstrom_si;
  conv[coulomb_meter][statcoulomb_centimeter]          = 1.0/statcoulomb_centimeter_si;
  conv[coulomb_meter][debye]                           = 1.0/debye_si;
  conv[coulomb_meter][picocoulomb_micrometer]          = 1.0/picocoulomb_micrometer_si;
  conv[coulomb_meter][electron_nanometer]              = 1.0/electron_nanometer_si;

  to_si = 1.0/conv[coulomb_meter][electron_angstrom];
  conv[electron_angstrom][coulomb_meter]               = to_si*conv[coulomb_meter][coulomb_meter];
  conv[electron_angstrom][electron_angstrom]           = 1.0;
  conv[electron_angstrom][statcoulomb_centimeter]      = to_si*conv[coulomb_meter][statcoulomb_centimeter];
  conv[electron_angstrom][debye]                       = to_si*conv[coulomb_meter][debye];
  conv[electron_angstrom][picocoulomb_micrometer]      = to_si*conv[coulomb_meter][picocoulomb_micrometer];
  conv[electron_angstrom][electron_nanometer]          = to_si*conv[coulomb_meter][electron_nanometer];

  to_si = 1.0/conv[coulomb_meter][statcoulomb_centimeter];
  conv[statcoulomb_centimeter][coulomb_meter]          = to_si*conv[coulomb_meter][coulomb_meter];
  conv[statcoulomb_centimeter][electron_angstrom]      = to_si*conv[coulomb_meter][electron_angstrom];
  conv[statcoulomb_centimeter][statcoulomb_centimeter] = 1.0;
  conv[statcoulomb_centimeter][debye]                  = to_si*conv[coulomb_meter][debye];
  conv[statcoulomb_centimeter][picocoulomb_micrometer] = to_si*conv[coulomb_meter][picocoulomb_micrometer];
  conv[statcoulomb_centimeter][electron_nanometer]     = to_si*conv[coulomb_meter][electron_nanometer];

  to_si = 1.0/conv[coulomb_meter][debye];
  conv[debye][coulomb_meter]                           = to_si*conv[coulomb_meter][coulomb_meter];
  conv[debye][electron_angstrom]                       = to_si*conv[coulomb_meter][electron_angstrom];
  conv[debye][statcoulomb_centimeter]                  = to_si*conv[coulomb_meter][statcoulomb_centimeter];
  conv[debye][debye]                                   = 1.0;
  conv[debye][picocoulomb_micrometer]                  = to_si*conv[coulomb_meter][picocoulomb_micrometer];
  conv[debye][electron_nanometer]                      = to_si*conv[coulomb_meter][electron_nanometer];

  to_si = 1.0/conv[coulomb_meter][picocoulomb_micrometer];
  conv[picocoulomb_micrometer][coulomb_meter]          = to_si*conv[coulomb_meter][coulomb_meter];
  conv[picocoulomb_micrometer][electron_angstrom]      = to_si*conv[coulomb_meter][electron_angstrom];
  conv[picocoulomb_micrometer][statcoulomb_centimeter] = to_si*conv[coulomb_meter][statcoulomb_centimeter];
  conv[picocoulomb_micrometer][debye]                  = to_si*conv[coulomb_meter][debye];
  conv[picocoulomb_micrometer][picocoulomb_micrometer] = 1.0;
  conv[picocoulomb_micrometer][electron_nanometer]     = to_si*conv[coulomb_meter][electron_nanometer];

  to_si = 1.0/conv[coulomb_meter][electron_nanometer];
  conv[electron_nanometer][coulomb_meter]              = to_si*conv[coulomb_meter][coulomb_meter];
  conv[electron_nanometer][electron_angstrom]          = to_si*conv[coulomb_meter][electron_angstrom];
  conv[electron_nanometer][statcoulomb_centimeter]     = to_si*conv[coulomb_meter][statcoulomb_centimeter];
  conv[electron_nanometer][debye]                      = to_si*conv[coulomb_meter][debye];
  conv[electron_nanometer][picocoulomb_micrometer]     = to_si*conv[coulomb_meter][picocoulomb_micrometer];
  conv[electron_nanometer][electron_nanometer]         = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Electric field conversion
double get_efield_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[volt_per_meter][volt_per_meter]                   = 1.0;
  conv[volt_per_meter][volt_per_angstrom]                = 1.0/volt_per_angstrom_si;
  conv[volt_per_meter][statvolt_per_centimeter]          = 1.0/statvolt_per_centimeter_si;
  conv[volt_per_meter][volt_per_centimeter]              = 1.0/volt_per_centimeter_si;
  conv[volt_per_meter][volt_per_micrometer]              = 1.0/volt_per_micrometer_si;
  conv[volt_per_meter][volt_per_nanometer]               = 1.0/volt_per_nanometer_si;

  to_si = 1.0/conv[volt_per_meter][volt_per_angstrom];
  conv[volt_per_angstrom][volt_per_meter]                = to_si*conv[volt_per_meter][volt_per_meter];
  conv[volt_per_angstrom][volt_per_angstrom]             = 1.0;
  conv[volt_per_angstrom][statvolt_per_centimeter]       = to_si*conv[volt_per_meter][statvolt_per_centimeter];
  conv[volt_per_angstrom][volt_per_centimeter]           = to_si*conv[volt_per_meter][volt_per_centimeter];
  conv[volt_per_angstrom][volt_per_micrometer]           = to_si*conv[volt_per_meter][volt_per_micrometer];
  conv[volt_per_angstrom][volt_per_nanometer]            = to_si*conv[volt_per_meter][volt_per_nanometer];

  to_si = 1.0/conv[volt_per_meter][statvolt_per_centimeter];
  conv[statvolt_per_centimeter][volt_per_meter]          = to_si*conv[volt_per_meter][volt_per_meter];
  conv[statvolt_per_centimeter][volt_per_angstrom]       = to_si*conv[volt_per_meter][volt_per_angstrom];
  conv[statvolt_per_centimeter][statvolt_per_centimeter] = 1.0;
  conv[statvolt_per_centimeter][volt_per_centimeter]     = to_si*conv[volt_per_meter][volt_per_centimeter];
  conv[statvolt_per_centimeter][volt_per_micrometer]     = to_si*conv[volt_per_meter][volt_per_micrometer];
  conv[statvolt_per_centimeter][volt_per_nanometer]      = to_si*conv[volt_per_meter][volt_per_nanometer];

  to_si = 1.0/conv[volt_per_meter][volt_per_centimeter];
  conv[volt_per_centimeter][volt_per_meter]              = to_si*conv[volt_per_meter][volt_per_meter];
  conv[volt_per_centimeter][volt_per_angstrom]           = to_si*conv[volt_per_meter][volt_per_angstrom];
  conv[volt_per_centimeter][statvolt_per_centimeter]     = to_si*conv[volt_per_meter][statvolt_per_centimeter];
  conv[volt_per_centimeter][volt_per_centimeter]         = 1.0;
  conv[volt_per_centimeter][volt_per_micrometer]         = to_si*conv[volt_per_meter][volt_per_micrometer];
  conv[volt_per_centimeter][volt_per_nanometer]          = to_si*conv[volt_per_meter][volt_per_nanometer];

  to_si = 1.0/conv[volt_per_meter][volt_per_micrometer];
  conv[volt_per_micrometer][volt_per_meter]              = to_si*conv[volt_per_meter][volt_per_meter];
  conv[volt_per_micrometer][volt_per_angstrom]           = to_si*conv[volt_per_meter][volt_per_angstrom];
  conv[volt_per_micrometer][statvolt_per_centimeter]     = to_si*conv[volt_per_meter][statvolt_per_centimeter];
  conv[volt_per_micrometer][volt_per_centimeter]         = to_si*conv[volt_per_meter][volt_per_centimeter];
  conv[volt_per_micrometer][volt_per_micrometer]         = 1.0;
  conv[volt_per_micrometer][volt_per_nanometer]          = to_si*conv[volt_per_meter][volt_per_nanometer];

  to_si = 1.0/conv[volt_per_meter][volt_per_nanometer];
  conv[volt_per_nanometer][volt_per_meter]               = to_si*conv[volt_per_meter][volt_per_meter];
  conv[volt_per_nanometer][volt_per_angstrom]            = to_si*conv[volt_per_meter][volt_per_angstrom];
  conv[volt_per_nanometer][statvolt_per_centimeter]      = to_si*conv[volt_per_meter][statvolt_per_centimeter];
  conv[volt_per_nanometer][volt_per_centimeter]          = to_si*conv[volt_per_meter][volt_per_centimeter];
  conv[volt_per_nanometer][volt_per_micrometer]          = to_si*conv[volt_per_meter][volt_per_micrometer];
  conv[volt_per_nanometer][volt_per_nanometer]           = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

// Demsity conversion
double get_density_conversion_factor(units from_unit_enum, units to_unit_enum)
{
  map<units, map<units, double> > conv;
  double to_si;

  conv[kilogram_per_metercu][kilogram_per_metercu]           = 1.0;
  conv[kilogram_per_metercu][gram_per_centimetercu]          = 1.0/gram_per_centimetercu_si;
  conv[kilogram_per_metercu][amu_per_bohrcu]                 = 1.0/amu_per_bohrcu_si;
  conv[kilogram_per_metercu][picogram_per_micrometercu]      = 1.0/picogram_per_micrometercu_si;
  conv[kilogram_per_metercu][attogram_per_nanometercu]       = 1.0/attogram_per_nanometercu_si;

  to_si = 1.0/conv[kilogram_per_metercu][gram_per_centimetercu];
  conv[gram_per_centimetercu][kilogram_per_metercu]          = to_si*conv[kilogram_per_metercu][kilogram_per_metercu];
  conv[gram_per_centimetercu][gram_per_centimetercu]         = 1.0;
  conv[gram_per_centimetercu][amu_per_bohrcu]                = to_si*conv[kilogram_per_metercu][amu_per_bohrcu];
  conv[gram_per_centimetercu][picogram_per_micrometercu]     = to_si*conv[kilogram_per_metercu][picogram_per_micrometercu];
  conv[gram_per_centimetercu][attogram_per_nanometercu]      = to_si*conv[kilogram_per_metercu][attogram_per_nanometercu];

  to_si = 1.0/conv[kilogram_per_metercu][amu_per_bohrcu];
  conv[amu_per_bohrcu][kilogram_per_metercu]                 = to_si*conv[kilogram_per_metercu][kilogram_per_metercu];
  conv[amu_per_bohrcu][gram_per_centimetercu]                = to_si*conv[kilogram_per_metercu][gram_per_centimetercu];
  conv[amu_per_bohrcu][amu_per_bohrcu]                       = 1.0;
  conv[amu_per_bohrcu][picogram_per_micrometercu]            = to_si*conv[kilogram_per_metercu][picogram_per_micrometercu];
  conv[amu_per_bohrcu][attogram_per_nanometercu]             = to_si*conv[kilogram_per_metercu][attogram_per_nanometercu];

  to_si = 1.0/conv[kilogram_per_metercu][picogram_per_micrometercu];
  conv[picogram_per_micrometercu][kilogram_per_metercu]      = to_si*conv[kilogram_per_metercu][kilogram_per_metercu];
  conv[picogram_per_micrometercu][gram_per_centimetercu]     = to_si*conv[kilogram_per_metercu][gram_per_centimetercu];
  conv[picogram_per_micrometercu][amu_per_bohrcu]            = to_si*conv[kilogram_per_metercu][amu_per_bohrcu];
  conv[picogram_per_micrometercu][picogram_per_micrometercu] = 1.0;
  conv[picogram_per_micrometercu][attogram_per_nanometercu]  = to_si*conv[kilogram_per_metercu][attogram_per_nanometercu];

  to_si = 1.0/conv[kilogram_per_metercu][attogram_per_nanometercu];
  conv[attogram_per_nanometercu][kilogram_per_metercu]       = to_si*conv[kilogram_per_metercu][kilogram_per_metercu];
  conv[attogram_per_nanometercu][gram_per_centimetercu]      = to_si*conv[kilogram_per_metercu][gram_per_centimetercu];
  conv[attogram_per_nanometercu][amu_per_bohrcu]             = to_si*conv[kilogram_per_metercu][amu_per_bohrcu];
  conv[attogram_per_nanometercu][picogram_per_micrometercu]  = to_si*conv[kilogram_per_metercu][picogram_per_micrometercu];
  conv[attogram_per_nanometercu][attogram_per_nanometercu]   = 1.0;

  return conv[from_unit_enum][to_unit_enum];
}

/* ---------------------------------------------------------------------- */

//  This routine returns the unit conversion factor between the
//  `from_system_enum` to the `to_system_enum` for the `unit_type_enum`.
double get_unit_conversion_factor(unit_type &unit_type_enum,
                                  sys_type from_system_enum,
                                  sys_type to_system_enum)
{
  units from_unit = get_lammps_system_unit(from_system_enum, unit_type_enum);
  units to_unit   = get_lammps_system_unit(to_system_enum, unit_type_enum);
  switch(unit_type_enum) {
    case mass :
      return get_mass_conversion_factor(from_unit, to_unit);
    case distance :
      return get_distance_conversion_factor(from_unit, to_unit);
    case time :
      return get_time_conversion_factor(from_unit, to_unit);
    case energy :
      return get_energy_conversion_factor(from_unit, to_unit);
    case velocity :
      return get_velocity_conversion_factor(from_unit, to_unit);
    case force :
      return get_force_conversion_factor(from_unit, to_unit);
    case torque :
      return get_torque_conversion_factor(from_unit, to_unit);
    case temperature :
      return get_temperature_conversion_factor(from_unit, to_unit);
    case pressure :
      return get_pressure_conversion_factor(from_unit, to_unit);
    case viscosity :
      return get_viscosity_conversion_factor(from_unit, to_unit);
    case charge :
      return get_charge_conversion_factor(from_unit, to_unit);
    case dipole :
      return get_dipole_conversion_factor(from_unit, to_unit);
    case efield :
      return get_efield_conversion_factor(from_unit, to_unit);
    case density :
    default :    // This is here to a prevent a compiler warning
      return get_density_conversion_factor(from_unit, to_unit);
  }
}

} // end of anonymous name space

/* ---------------------------------------------------------------------- */

//  Wrapper to the routine that gets the unit conversion. Translates strings
//  to enumerations and then call get_unit_conversion_factor()
int lammps_unit_conversion(string const &unit_type_str,
                           string const &from_system_str,
                           string const &to_system_str,
                           double &conversion_factor)
{
    // initialize
    conversion_factor = 0.0;
    initialize_dictionaries();

    // convert input to enumeration
    unit_type unit_type_enum;
    {
      map<string, unit_type>::const_iterator itr = unit_dic.find(unit_type_str);
      if (itr != unit_dic.end()) unit_type_enum = itr->second;
      else return 1;  // error
    }
    sys_type from_system_enum;
    {
      map<string, sys_type>::const_iterator
          itr = system_dic.find(from_system_str);
      if (itr != system_dic.end()) from_system_enum = itr->second;
      else return 1;  // error
    }
    sys_type to_system_enum;
    {
      map<string, sys_type>::const_iterator
          itr = system_dic.find(to_system_str);
      if (itr != system_dic.end()) to_system_enum = itr->second;
      else return 1;
    }

    // process unit conversions
    conversion_factor = get_unit_conversion_factor(unit_type_enum,
                                                   from_system_enum,
                                                   to_system_enum);
    return 0;
}


