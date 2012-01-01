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

#include "math.h"
#include "compute_gyration_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationMolecule::ComputeGyrationMolecule(LAMMPS *lmp, 
						 int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute gyration/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute gyration/molecule requires molecular atom style");

  vector_flag = 1;
  extvector = 0;

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);
  size_vector = nmolecules;

  memory->create(massproc,nmolecules,"gyration/molecule:massproc");
  memory->create(masstotal,nmolecules,"gyration/molecule:masstotal");
  memory->create(com,nmolecules,3,"gyration/molecule:com");
  memory->create(comall,nmolecules,3,"gyration/molecule:comall");
  memory->create(rg,nmolecules,"gyration/molecule:rg");
  memory->create(rgall,nmolecules,"gyration/molecule:rgall");
  vector = rgall;

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int i,imol;
  double massone;

  for (i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

ComputeGyrationMolecule::~ComputeGyrationMolecule()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(rg);
  memory->destroy(rgall);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute gyration/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationMolecule::compute_vector()
{
  int i,imol;
  double xbox,ybox,zbox,dx,dy,dz;
  double massone;

  invoked_array = update->ntimestep;

  for (i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      com[imol][0] += (x[i][0] + xbox*xprd) * massone;
      com[imol][1] += (x[i][1] + ybox*yprd) * massone;
      com[imol][2] += (x[i][2] + zbox*zprd) * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
		MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmolecules; i++) {
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }

  for (i = 0; i < nmolecules; i++) rg[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      dx = (x[i][0] + xbox*xprd) - comall[imol][0];
      dy = (x[i][1] + ybox*yprd) - comall[imol][1];
      dz = (x[i][2] + zbox*zprd) - comall[imol][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg[imol] += (dx*dx + dy*dy + dz*dz) * massone;
    }

  MPI_Allreduce(rg,rgall,nmolecules,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nmolecules; i++) rgall[i] = sqrt(rgall[i]/masstotal[i]);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGyrationMolecule::memory_usage()
{
  double bytes = 4*nmolecules * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += 2*nmolecules*3 * sizeof(double);
  return bytes;
}
