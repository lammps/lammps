/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio)
------------------------------------------------------------------------- */

#include "elements.h"

#include "utils.h"
#include "error.h"

using namespace LAMMPS_NS;

static int MAX_ATOMIC_NUMBER = 96;

static const std::string Symbol[] = {"X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm"};
static const std::string Name[] = {"X","Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon","Sodium","Magnesium","Aluminum","Silicon","Phosphorus","Sulfur","Chlorine","Argon","Potassium","Calcium","Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt","Nickel","Copper","Zinc","Gallium","Germanium","Arsenic","Selenium","Bromine","Krypton","Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium","Indium","Tin","Antimony","Tellurium","Iodine","Xenon","Cesium","Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium","Europium","Gadolinium","Terbium","Dysprosium","Holmium","Erbium","Thulium","Ytterbium","Lutetium","Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum","Gold","Mercury","Thallium","Lead","Bismuth","Polonium","Astatine","Radon","Francium","Radium","Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium","Americium","Curium"};
static const std::string CPKHexColor[] = {"000000","FFFFFF","D9FFFF","CC80FF","C2FF00","FFB5B5","909090","3050F8","FF0D0D","90E050","B3E3F5","AB5CF2","8AFF00","BFA6A6","F0C8A0","FF8000","FFFF30","1FF01F","80D1E3","8F40D4","3DFF00","E6E6E6","BFC2C7","A6A6AB","8A99C7","9C7AC7","E06633","F090A0","50D050","C88033","7D80B0","C28F8F","668F8F","BD80E3","FFA100","A62929","5CB8D1","702EB0","00FF00","94FFFF","94E0E0","73C2C9","54B5B5","3B9E9E","248F8F","0A7D8C","6985","C0C0C0","FFD98F","A67573","668080","9E63B5","D47A00","940094","429EB0","57178F","00C900","70D4FF","FFFFC7","D9FFC7","C7FFC7","A3FFC7","8FFFC7","61FFC7","45FFC7","30FFC7","1FFFC7","00FF9C","000000","00D452","00BF38","00AB24","4DC2FF","4DA6FF","2194D6","267DAB","266696","175487","D0D0E0","FFD123","B8B8D0","A6544D","575961","9E4FB5","AB5C00","754F45","428296","420066","007D00","70ABFA","00BAFF","00A1FF","008FFF","0080FF","006BFF","545CF2","785CE3"};
static const double AtomicMass[] = {0.0,1.008,4.0026,7.0,9.012183,10.81,12.011,14.007,15.999,18.99840316,20.18,22.9897693,24.305,26.981538,28.085,30.973762,32.07,35.45,39.9,39.0983,40.08,44.95591,47.867,50.9415,51.996,54.93804,55.84,58.93319,58.693,63.55,65.4,69.723,72.63,74.92159,78.97,79.9,83.8,85.468,87.62,88.90584,91.22,92.90637,95.95,96.90636,101.1,102.9055,106.42,107.868,112.41,114.818,118.71,121.76,127.6,126.9045,131.29,132.905452,137.33,138.9055,140.116,140.90766,144.24,144.91276,150.4,151.964,157.2,158.92535,162.5,164.93033,167.26,168.93422,173.05,174.9668,178.49,180.9479,183.84,186.207,190.2,192.22,195.08,196.96657,200.59,204.383,207.0,208.9804,208.98243,209.98715,222.01758,223.01973,226.02541,227.02775,232.038,231.03588,238.0289,237.048172,244.0642,243.06138,247.07035};
static const double VdwRadius[] = {0.0,120.0,140.0,182.0,153.0,192.0,170.0,155.0,152.0,135.0,154.0,227.0,173.0,184.0,210.0,180.0,180.0,175.0,188.0,275.0,231.0,211.0,187.0,179.0,189.0,197.0,194.0,192.0,163.0,140.0,139.0,187.0,211.0,185.0,190.0,183.0,202.0,303.0,249.0,219.0,186.0,207.0,209.0,209.0,207.0,195.0,202.0,172.0,158.0,193.0,217.0,206.0,206.0,198.0,216.0,343.0,268.0,240.0,235.0,239.0,229.0,236.0,229.0,233.0,237.0,221.0,229.0,216.0,235.0,227.0,242.0,221.0,212.0,217.0,210.0,217.0,216.0,202.0,209.0,166.0,209.0,196.0,202.0,207.0,197.0,202.0,220.0,348.0,283.0,260.0,237.0,243.0,240.0,221.0,243.0,244.0,245.0};
static const double CovalentRadius[] =  {0.0,31.0,28.0,128.0,96.0,85.0,76.0,71.0,66.0,57.0,58.0,166.0,141.0,121.0,111.0,107.0,105.0,102.0,106.0,203.0,176.0,170.0,160.0,153.0,139.0,139.0,132.0,126.0,124.0,132.0,122.0,122.0,120.0,119.0,120.0,120.0,116.0,220.0,195.0,190.0,175.0,164.0,154.0,147.0,146.0,142.0,139.0,145.0,144.0,142.0,139.0,139.0,138.0,139.0,140.0,244.0,215.0,207.0,204.0,203.0,201.0,199.0,198.0,198.0,196.0,194.0,192.0,192.0,189.0,190.0,187.0,187.0,175.0,170.0,162.0,151.0,144.0,141.0,136.0,136.0,132.0,145.0,146.0,148.0,140.0,150.0,150.0,260.0,221.0,215.0,206.0,200.0,196.0,190.0,187.0,180.0,169.0};


/* ------------------------------------------------------------------ */

// return "n/a" if atomic_number out of range, otherwise
// elements.cpp:44:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]

std::string elements::symbol(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return "n/a";
  } else
    return Symbol[atomic_number];
}

/* ------------------------------------------------------------------ */

std::string elements::name(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return "n/a";
  } else
    return Name[atomic_number];
}

/* ------------------------------------------------------------------ */

std::string elements::cpkHexColor(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return "n/a";
  } else
    return CPKHexColor[atomic_number];
}

/* ------------------------------------------------------------------ */

// return -1.0 if atomic_number out of range, otherwise
// elements.cpp:44:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]

double elements::atomic_mass(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return -1.0;
  } else
    return AtomicMass[atomic_number];
}

/* ------------------------------------------------------------------ */

double elements::vdw_radius(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return -1.0;
  } else
    return VdwRadius[atomic_number];
}

/* ------------------------------------------------------------------ */

double elements::covalent_radius(int atomic_number, Error *error)
{
  if( atomic_number<0 || atomic_number>MAX_ATOMIC_NUMBER) {
    error->all(FLERR, "atomic_number {} out of range (0-{})", atomic_number, MAX_ATOMIC_NUMBER);
    return -1.0;
  } else
    return CovalentRadius[atomic_number];
}

/* ------------------------------------------------------------------ */

int elements::atomic_number_with_symbol(const std::string &symbol, Error *error)
{
  int i=0;
  
  while( i<=MAX_ATOMIC_NUMBER ) {
    if( LAMMPS_NS::utils::strmatch(Symbol[i], symbol) )
       return i;
    else
       i++;
  }
  
  error->all(FLERR, "symbol {} not found", symbol);
  return -1;
}

/* ------------------------------------------------------------------ */

int elements::atomic_number_with_closest_mass(double mass, Error *error)
{
  int i=0;
  double previous_mass=0.0;
  
  if( mass<0 )
    error->all(FLERR, "atomic mass {} is negative, must be >=0", mass);
      
  while( i<=MAX_ATOMIC_NUMBER ) {
  
    double current_mass=AtomicMass[i];
    
    if( mass<=current_mass ) {
    
      // now we have previous_mass < mass <= current_mass
    
      if( mass-previous_mass < current_mass-mass )
        return i-1;
      else
        return i;
    }
        
    // otherwise keep looking
    previous_mass=current_mass;
    i++;
  }
  
  error->all(FLERR, "atomic mass {} higher than heaviest element {} ({}) available",
                              mass, Name[MAX_ATOMIC_NUMBER], AtomicMass[MAX_ATOMIC_NUMBER]);
                              
  return -1;
}

