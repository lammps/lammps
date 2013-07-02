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

#include "citeme.h"
#include "version.h"

#include <stdio.h>
#include <set>

using namespace LAMMPS_NS;

// the list of publications is below
static const char * const publication[] = {
  /* PLIMPTON_1995 */
  "S.J. Plimpton,\n"
  " Fast Parallel Algorithms for Short-Range Molecular Dynamics,\n"
  "  J. Comp. Phys., 117, 1-19 (1995)\n\n",
  /* PLIMPTON_1997 */
  "S.J. Plimpton, R. Pollock, M. Stevens,\n"
  " Particle-Mesh Ewald and r-RESPA for Parallel Molecular Dynamics Simulations,\n"
  " Proc. of the 8th SIAM Conference on Parallel Processing for Scientific Computing,\n"
  "  Minneapolis, MN (March 1997).\n\n",
  /* AUHL_2003 */
  "R. Auhl, R. Everaers, G.S. Grest, K. Kremer, S.J. Plimpton,\n"
  " Equilibration of long chain polymer melts in computer simulations,\n"
  " J. Chem. Phys., 119, 12718-12728 (2003)\n\n",
  /* JANSSENS_2006 */
  "K.G.F. Janssens, D. Olmsted, E.A. Holm, S.M. Foiles, S.J. Plimpton, P.M. Derlet,\n"
  " Computing the Mobility of Grain Boundaries,\n"
  " Nature Materials, 5, 124-127 (2006).\n\n",
  /* INTVELD_2008 */
  "P.J. in 't Veld, S.J. Plimpton, G.S. Grest,\n"
  " Accurate and Efficient Methods for Modeling Colloidal Mixtures in an\n"
  " Explicit Solvent using Molecular Dynamics,\n"
  " Comp. Phys. Comm., 179, 320-329 (2008).\n\n",
  /* PARKS_2008 */
  "M.L. Parks, R.B. Lehoucq, S.J. Plimpton, S.A. Silling,\n"
  " Implementing peridynamics within a molecular dynamics code,\n"
  " Comp. Phys. Comm., 179, 777-783 (2008)\n\n",
  /* MUKHERJEE_2008 */
  "R.M. Mukherjee, P.S. Crozier, S.J. Plimpton, K.S. Anderson,\n"
  " Substructured molecular dynamics using multibody dynamics algorithms,\n"
  " Intl. J. of Non-Linear Mechanics, 43, 1045-1055 (2008)\n\n",
  /* BROWN_2009 */
  "W.M. Brown, M.K. Petersen, S.J. Plimpton, and G.S. Grest,\n"
  " Liquid crystal nanodroplets in solution,\n"
  " J. Chem. Phys., 130, 044901 (2009)\n\n",
  /* THOMPSON_2009 */
  "A.P. Thompson, S.J. Plimpton, W. Mattson,\n"
  " General formulation of pressure and stress tensor for arbitrary many-body\n"
  " interaction potentials under periodic boundary conditions,\n"
  " J. Chem. Phys., 131, 154107 (2009)\n\n",
  /* PETERSEN_2010 */
  "M.K. Petersen, J.B. Lechman, S.J. Plimpton, G.S. Grest, P.J. in't Veld, P.R. Schunk,\n"
  " Mesoscale Hydrodynamics via Stochastic Rotation Dynamics:\n"
  "  Comparison with Lennard-Jones Fluid,\n"
  " J. Chem. Phys., 132, 174106 (2010)\n\n",
  /* BROWN_2011 */
  "W.M. Brown, P. Wang, S.J. Plimpton, A.N. Tharrington,\n"
  " Implementing Molecular Dynamics on Hybrid High Performance Computers\n"
  " - Short Range Forces,\n"
  " Comp. Phys. Comm., 182, 898-911, (2011)\n\n",
  /* BROWN_2012 */
  "W.M. Brown, A. Kohlmeyer, S.J. Plimpton, A.N. Tharrington,\n"
  " Implementing Molecular Dynamics on Hybrid High Performance Computers\n"
  " - Particle-Particle Particle-Mesh,\n"
  " Comp. Phys. Comm., 183, 449-459 (2012)\n\n",
  /* JARAMILLO_BOTERO_2011 */
  "A. Jaramillo-Botero, J. Su, A. Qi, W.A. Goddard III,\n"
  " Large-Scale, Long-Term Nonadiabatic Electron Molecular Dynamics\n"
  " for Describing Material Properties and Phenomena in Extreme Environments,\n"
  " J. Comp. Chem., 32, 497-512 (2011)\n\n",
  /* KONG_2011 */
  "L.T. Kong,\n"
  " Phonon dispersion measured directly from molecular dynamics simulations,\n"
  " Comp Phys Comm, 182, 2201-2207 (2011)\n\n",
  /* AKTULGA_2012 */
  "H.M. Aktulga, J.C. Fogarty, S.A. Pandit, A.Y. Grama,\n"
  " Parallel reactive molecular dynamics: Numerical methods and algorithmic techniques,\n"
  " Parallel Computing, 38, 245-259 (2012)\n\n",
  /* PLIMPTON_2012 */
  "S.J. Plimpton and A.P. Thompson,\n"
  " Computational Aspects of Many-body Potentials,\n"
  " MRS Bulletin, 37, 513-521 (2012)\n\n",
  /* SIRK_2013 */
  "T.W. Sirk, S. Moore, E.F. Brown,\n"
  " Characteristics of thermal conductivity in classical water models,\n"
  " J. Chem. Phys., 138, 064505 (2013)\n\n",
  NULL 
};

// to simplify the code below
typedef std::set<int> intset;

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp) {

  intset *c = new intset;
  pubs = (void *)c;
}

/* ---------------------------------------------------------------------- */

void CiteMe::add(int ref)
{
  intset *c = (intset *) pubs;
  if ((ref >= 0 ) && (ref < LAST_ENTRY))
    c->insert(ref);
}

/* ---------------------------------------------------------------------- */

static const char cite_header[] = "\n"
  "------------------------------------------------------------------------\n"
  "This simulation made use of algorithms and methodologies described\n"
  "in the following references:\n\n"
  "The LAMMPS Molecular Dynamics Simulator, Version " LAMMPS_VERSION "\n"
  "    http://lammps.sandia.gov\n\n";


CiteMe::~CiteMe(){
  intset *c = (intset *)(pubs);

  if (screen)
    fputs(cite_header,screen);

  if (logfile)
    fputs(cite_header,logfile);

  for (intset::const_iterator i = c->begin(); i != c->end(); ++i) {
    if (screen)
      fputs(publication[*i],screen);
    if (logfile)
      fputs(publication[*i],logfile);
  }

  delete c;
}


