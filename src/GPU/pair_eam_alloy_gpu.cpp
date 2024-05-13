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

/* ----------------------------------------------------------------------
   Contributing authors: Trung Dac Nguyen (ORNL), W. Michael Brown (ORNL)
------------------------------------------------------------------------- */

#include "pair_eam_alloy_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "suffix.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int eam_alloy_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                       int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                       double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                       double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r,
                       int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                       int &fp_size);
void eam_alloy_gpu_clear();
int **eam_alloy_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int &host_start, int **ilist,
                              int **jnum, const double cpu_time, bool &success, int &inum,
                              void **fp_ptr);
void eam_alloy_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                           double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           int &host_start, const double cpu_time, bool &success, void **fp_ptr);
void eam_alloy_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                                 const bool vatom);
double eam_alloy_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairEAMAlloyGPU::PairEAMAlloyGPU(LAMMPS *lmp) : PairEAM(lmp), gpu_mode(GPU_FORCE)
{
  one_coeff = 1;
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ---------------------------------------------------------------------- */

PairEAMAlloyGPU::~PairEAMAlloyGPU()
{
  eam_alloy_gpu_clear();
}

/* ---------------------------------------------------------------------- */

double PairEAMAlloyGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + eam_alloy_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairEAMAlloyGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // compute density on each atom on GPU

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int inum, host_start, inum_dev;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = eam_alloy_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo,
                                         subhi, atom->tag, atom->nspecial, atom->special, eflag,
                                         vflag, eflag_atom, vflag_atom, host_start, &ilist,
                                         &numneigh, cpu_time, success, inum_dev, &fp_pinned);
  } else {    // gpu_mode == GPU_FORCE
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    eam_alloy_gpu_compute(neighbor->ago, inum, nlocal, nall, atom->x, atom->type, ilist, numneigh,
                          firstneigh, eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                          success, &fp_pinned);
  }

  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  // communicate derivative of embedding function

  comm->forward_comm(this);

  // compute forces on each atom on GPU
  if (gpu_mode != GPU_FORCE)
    eam_alloy_gpu_compute_force(nullptr, eflag, vflag, eflag_atom, vflag_atom);
  else
    eam_alloy_gpu_compute_force(ilist, eflag, vflag, eflag_atom, vflag_atom);
  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMAlloyGPU::init_style()
{

  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i, j);
        cut *= cut;
        if (cut > maxcut) maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int fp_size;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = eam_alloy_gpu_init(
      atom->ntypes + 1, cutforcesq, type2rhor, type2z2r, type2frho, rhor_spline, z2r_spline,
      frho_spline, cutsq, rdr, rdrho, rhomax, nrhor, nrho, nz2r, nfrho, nr, atom->nlocal,
      atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen, fp_size);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
  fp_single = fp_size != sizeof(double);

  embedstep = -1;
}

/* ---------------------------------------------------------------------- */

double PairEAMAlloyGPU::single(int i, int j, int itype, int jtype, double rsq,
                               double /* factor_coul */, double /* factor_lj */, double &fforce)
{
  int m;
  double r, p, rhoip, rhojp, z2, z2p, recip, phi, phip, psip;
  double *coeff;

  r = sqrt(rsq);
  p = r * rdr + 1.0;
  m = static_cast<int>(p);
  m = MIN(m, nr - 1);
  p -= m;
  p = MIN(p, 1.0);

  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0] * p + coeff[1]) * p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];
  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
  z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

  double fp_i, fp_j;
  if (!fp_single) {
    fp_i = ((double *) fp_pinned)[i];
    fp_j = ((double *) fp_pinned)[j];
  } else {
    fp_i = ((float *) fp_pinned)[i];
    fp_j = ((float *) fp_pinned)[j];
  }

  recip = 1.0 / r;
  phi = z2 * recip;
  phip = z2p * recip - phi * recip;
  psip = fp_i * rhojp + fp_j * rhoip + phip;
  fforce = -psip * recip;

  return phi;
}

/* ---------------------------------------------------------------------- */

int PairEAMAlloyGPU::pack_forward_comm(int n, int *list, double *buf, int /* pbc_flag */,
                                       int * /* pbc */)
{
  int i, j, m;

  m = 0;

  if (fp_single) {
    auto fp_ptr = (float *) fp_pinned;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<double>(fp_ptr[j]);
    }
  } else {
    auto fp_ptr = (double *) fp_pinned;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = fp_ptr[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairEAMAlloyGPU::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (fp_single) {
    auto fp_ptr = (float *) fp_pinned;
    for (i = first; i < last; i++) fp_ptr[i] = buf[m++];
  } else {
    auto fp_ptr = (double *) fp_pinned;
    for (i = first; i < last; i++) fp_ptr[i] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMAlloyGPU::coeff(int narg, char **arg)
{
  int i, j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Number of element to type mappings does not match number of atom types");

  // read EAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    memory->destroy(setfl->mass);
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }
  setfl = new Setfl();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i - 2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i], setfl->elements[j]) == 0) break;
    if (j < setfl->nelements)
      map[i - 2] = j;
    else
      error->all(FLERR, "No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++) setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR, i, setfl->mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMAlloyGPU::read_file(char *filename)
{
  Setfl *file = setfl;

  // read potential file
  if (comm->me == 0) {
    PotentialFileReader reader(PairEAM::lmp, filename, "eam/alloy", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY, unit_convert);
    try {
      reader.skip_line();
      reader.skip_line();
      reader.skip_line();

      // extract element names from nelements line
      ValueTokenizer values = reader.next_values(1);
      file->nelements = values.next_int();

      if ((int) values.count() != file->nelements + 1)
        error->one(FLERR, "Incorrect element names in EAM potential file");

      file->elements = new char *[file->nelements];
      for (int i = 0; i < file->nelements; i++)
        file->elements[i] = utils::strdup(values.next_string());

      //

      values = reader.next_values(5);
      file->nrho = values.next_int();
      file->drho = values.next_double();
      file->nr = values.next_int();
      file->dr = values.next_double();
      file->cut = values.next_double();

      if ((file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
        error->one(FLERR, "Invalid EAM potential file");

      memory->create(file->mass, file->nelements, "pair:mass");
      memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
      memory->create(file->rhor, file->nelements, file->nr + 1, "pair:rhor");
      memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");

      for (int i = 0; i < file->nelements; i++) {
        values = reader.next_values(2);
        values.next_int();    // ignore
        file->mass[i] = values.next_double();

        reader.next_dvector(&file->frho[i][1], file->nrho);
        reader.next_dvector(&file->rhor[i][1], file->nr);
        if (unit_convert) {
          for (int j = 1; j < file->nrho; ++j) file->frho[i][j] *= conversion_factor;
        }
      }

      for (int i = 0; i < file->nelements; i++) {
        for (int j = 0; j <= i; j++) {
          reader.next_dvector(&file->z2r[i][j][1], file->nr);
          if (unit_convert) {
            for (int k = 1; k < file->nr; ++k) file->z2r[i][j][k] *= conversion_factor;
          }
        }
      }
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  // broadcast potential information
  MPI_Bcast(&file->nelements, 1, MPI_INT, 0, world);

  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  // allocate memory on other procs
  if (comm->me != 0) {
    file->elements = new char *[file->nelements];
    for (int i = 0; i < file->nelements; i++) file->elements[i] = nullptr;
    memory->create(file->mass, file->nelements, "pair:mass");
    memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
    memory->create(file->rhor, file->nelements, file->nr + 1, "pair:rhor");
    memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");
  }

  // broadcast file->elements string array
  for (int i = 0; i < file->nelements; i++) {
    int n;
    if (comm->me == 0) n = strlen(file->elements[i]) + 1;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (comm->me != 0) file->elements[i] = new char[n];
    MPI_Bcast(file->elements[i], n, MPI_CHAR, 0, world);
  }

  // broadcast file->mass, frho, rhor
  for (int i = 0; i < file->nelements; i++) {
    MPI_Bcast(&file->mass[i], 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&file->frho[i][1], file->nrho, MPI_DOUBLE, 0, world);
    MPI_Bcast(&file->rhor[i][1], file->nr, MPI_DOUBLE, 0, world);
  }

  // broadcast file->z2r
  for (int i = 0; i < file->nelements; i++) {
    for (int j = 0; j <= i; j++) { MPI_Bcast(&file->z2r[i][j][1], file->nr, MPI_DOUBLE, 0, world); }
  }
}

/* ----------------------------------------------------------------------
   copy read-in setfl potential to standard array format
------------------------------------------------------------------------- */

void PairEAMAlloyGPU::file2array()
{
  int i, j, m, n;
  int ntypes = atom->ntypes;

  // set function params directly from setfl file

  nrho = setfl->nrho;
  nr = setfl->nr;
  drho = setfl->drho;
  dr = setfl->dr;
  rhomax = (nrho - 1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho + 1, "pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = setfl->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nr + 1, "pair:rhor");

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements + 1) / 2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nr + 1, "pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}
