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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_cac_eam.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "timer.h"
#include "neigh_list.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"


//#include "math_extra.h"
#define MAXNEIGH1  300
#define MAXNEIGH2  30
#define MAXLINE 1024
#define DELTA 4
#define EXPAND 10
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairCACEAM::PairCACEAM(LAMMPS *lmp) : PairCAC(lmp)
{


  
  restartinfo = 0;
  
  manybody_flag = 1;
  rho = NULL;
  fp = NULL;
  map = NULL;
  type2frho = NULL;

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;
  scale = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;
  outer_neighflag = 1;
 
  inner_neighbor_coords = NULL;
  outer_neighbor_coords = NULL;
  inner_neighbor_types = NULL;
  outer_neighbor_types = NULL;
  rho = NULL;
  fp = NULL;
  nmax = 0;
  

  interior_scales = NULL;
  surface_counts = NULL;

  surface_counts_max[0] = 0;
  surface_counts_max[1] = 0;
  surface_counts_max[2] = 0;
  surface_counts_max_old[0] = 0;
  surface_counts_max_old[1] = 0;
  surface_counts_max_old[2] = 0;
}

/* ---------------------------------------------------------------------- */

PairCACEAM::~PairCACEAM() {

  if (copymode) return;

  memory->destroy(rho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] map;
    delete[] type2frho;
    map = NULL;
    type2frho = NULL;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
    memory->destroy(scale);
    memory->destroy(mass_matrix);

    memory->destroy(surface_counts);
    memory->destroy(interior_scales);
    memory->destroy(inner_neighbor_coords);
    memory->destroy(outer_neighbor_coords);
    memory->destroy(inner_neighbor_types);
    memory->destroy(outer_neighbor_types);
    memory->destroy(rho);
    memory->destroy(fp);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete[] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
    }
    memory->sfree(funcfl);
    funcfl = NULL;
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    delete[] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
    setfl = NULL;
  }

  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete[] fs->elements[i];
    delete[] fs->elements;
    delete[] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
    fs = NULL;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);
  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}



/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  map = new int[n + 1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n + 1];
  memory->create(type2rhor, n + 1, n + 1, "pair:type2rhor");
  memory->create(type2z2r, n + 1, n + 1, "pair:type2z2r");
  memory->create(scale, n + 1, n + 1, "pair:scale");
  memory->create(mass_matrix,max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  memory->create(surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(dof_set, 6, 4, "pairCAC:surf_set");
  memory->create(sort_surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(sort_dof_set, 6, 4, "pairCAC:surf_set");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACEAM::coeff(int narg, char **arg) {
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2], funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *)
      memory->srealloc(funcfl, nfuncfl * sizeof(Funcfl), "pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file, arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i = j

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      if (i == j) {
        setflag[i][i] = 1;
        map[i] = ifuncfl;
        atom->set_mass(FLERR, i, funcfl[ifuncfl].mass);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACEAM::init_one(int i, int j) {

  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file
  double cutmax;
  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax, funcfl[m].cut);
  }
  else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;
  
  cutforcesq = cutmax*cutmax;
  cut_global_s = MAX(cut_global_s, cutmax);
  if (outer_neighflag)
  return 2*cutmax;
  else
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairCACEAM::init_style()
{
  check_existence_flags();
  // convert read-in file(s) to arrays and spline them
  maxneigh_quad_inner = MAXNEIGH2;
  maxneigh_quad_outer = MAXNEIGH1;
  file2array();
  array2spline();

  atom->outer_neigh_flag=1;
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->cac = 1;
  //surface selection array 
  surf_set[0][0] = 1;
  surf_set[0][1] = -1;
  surf_set[1][0] = 1;
  surf_set[1][1] = 1;
  surf_set[2][0] = 2;
  surf_set[2][1] = -1;
  surf_set[3][0] = 2;
  surf_set[3][1] = 1;
  surf_set[4][0] = 3;
  surf_set[4][1] = -1;
  surf_set[5][0] = 3;
  surf_set[5][1] = 1;

  //surface DOF array

  dof_set[0][0] = 0;
  dof_set[0][1] = 3;
  dof_set[0][2] = 4;
  dof_set[0][3] = 7;

  dof_set[1][0] = 1;
  dof_set[1][1] = 2;
  dof_set[1][2] = 5;
  dof_set[1][3] = 6;

  dof_set[2][0] = 0;
  dof_set[2][1] = 1;
  dof_set[2][2] = 4;
  dof_set[2][3] = 5;

  dof_set[3][0] = 2;
  dof_set[3][1] = 3;
  dof_set[3][2] = 6;
  dof_set[3][3] = 7;

  dof_set[4][0] = 0;
  dof_set[4][1] = 1;
  dof_set[4][2] = 2;
  dof_set[4][3] = 3;

  dof_set[5][0] = 4;
  dof_set[5][1] = 5;
  dof_set[5][2] = 6;
  dof_set[5][3] = 7;

  for (int si = 0; si < 6; si++) {
    sort_dof_set[si][0] = dof_set[si][0];
    sort_dof_set[si][1] = dof_set[si][1];
    sort_dof_set[si][2] = dof_set[si][2];
    sort_dof_set[si][3] = dof_set[si][3];
    sort_surf_set[si][0] = surf_set[si][0];
    sort_surf_set[si][1] = surf_set[si][1];
  }
 
}

////////////////////////
void PairCACEAM::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl - 1];

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str, "Cannot open EAM potential file %s", filename);
      error->one(FLERR, str);
    }
  }

  int tmp, nwords;
  if (me == 0) {
    fgets(line, MAXLINE, fptr);
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lg", &tmp, &file->mass);
    fgets(line, MAXLINE, fptr);
    nwords = sscanf(line, "%d %lg %d %lg %lg",
      &file->nrho, &file->drho, &file->nr, &file->dr, &file->cut);
  }

  MPI_Bcast(&nwords, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->mass, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  if ((nwords != 5) || (file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
    error->all(FLERR, "Invalid EAM potential file");

  memory->create(file->frho, (file->nrho + 1), "pair:frho");
  memory->create(file->rhor, (file->nr + 1), "pair:rhor");
  memory->create(file->zr, (file->nr + 1), "pair:zr");

  if (me == 0) grab(fptr, file->nrho, &file->frho[1]);
  MPI_Bcast(&file->frho[1], file->nrho, MPI_DOUBLE, 0, world);

  if (me == 0) grab(fptr, file->nr, &file->zr[1]);
  MPI_Bcast(&file->zr[1], file->nr, MPI_DOUBLE, 0, world);

  if (me == 0) grab(fptr, file->nr, &file->rhor[1]);
  MPI_Bcast(&file->rhor[1], file->nr, MPI_DOUBLE, 0, world);

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
convert read-in funcfl potential(s) to standard array format
interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairCACEAM::file2array()
{
  int i, j, k, m, n;
  int ntypes = atom->ntypes;
  double sixth = 1.0 / 6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr, file->dr);
    drho = MAX(drho, file->drho);
    rmax = MAX(rmax, (file->nr - 1) * file->dr);
    rhomax = MAX(rhomax, (file->nrho - 1) * file->drho);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax / dr + 0.5);
  nrho = static_cast<int> (rhomax / drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho + 1, "pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r, p, cof1, cof2, cof3, cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m - 1)*drho;
      p = r / file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k, file->nrho - 2);
      k = MAX(k, 2);
      p -= k;
      p = MIN(p, 2.0);
      cof1 = -sixth*p*(p - 1.0)*(p - 2.0);
      cof2 = 0.5*(p*p - 1.0)*(p - 2.0);
      cof3 = -0.5*p*(p + 1.0)*(p - 2.0);
      cof4 = sixth*p*(p*p - 1.0);
      frho[n][m] = cof1*file->frho[k - 1] + cof2*file->frho[k] +
        cof3*file->frho[k + 1] + cof4*file->frho[k + 2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho - 1;

    // ------------------------------------------------------------------
    // setup rhor arrays
    // ------------------------------------------------------------------

    // allocate rhor arrays
    // nrhor = # of funcfl files

    nrhor = nfuncfl;
    memory->destroy(rhor);
    memory->create(rhor, nrhor, nr + 1, "pair:rhor");

    // interpolate each file's rhor to a single grid and cutoff

    n = 0;
    for (i = 0; i < nfuncfl; i++) {
      Funcfl *file = &funcfl[i];
      for (m = 1; m <= nr; m++) {
        r = (m - 1)*dr;
        p = r / file->dr + 1.0;
        k = static_cast<int> (p);
        k = MIN(k, file->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p = MIN(p, 2.0);
        cof1 = -sixth*p*(p - 1.0)*(p - 2.0);
        cof2 = 0.5*(p*p - 1.0)*(p - 2.0);
        cof3 = -0.5*p*(p + 1.0)*(p - 2.0);
        cof4 = sixth*p*(p*p - 1.0);
        rhor[n][m] = cof1*file->rhor[k - 1] + cof2*file->rhor[k] +
          cof3*file->rhor[k + 1] + cof4*file->rhor[k + 2];
      }
      n++;
    }

    // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
    // for funcfl files, I,J mapping only depends on I
    // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        type2rhor[i][j] = map[i];

    // ------------------------------------------------------------------
    // setup z2r arrays
    // ------------------------------------------------------------------

    // allocate z2r arrays
    // nz2r = N*(N+1)/2 where N = # of funcfl files

    nz2r = nfuncfl*(nfuncfl + 1) / 2;
    memory->destroy(z2r);
    memory->create(z2r, nz2r, nr + 1, "pair:z2r");

    // create a z2r array for each file against other files, only for I >= J
    // interpolate zri and zrj to a single grid and cutoff

    double zri, zrj;

    n = 0;
    for (i = 0; i < nfuncfl; i++) {
      Funcfl *ifile = &funcfl[i];
      for (j = 0; j <= i; j++) {
        Funcfl *jfile = &funcfl[j];

        for (m = 1; m <= nr; m++) {
          r = (m - 1)*dr;

          p = r / ifile->dr + 1.0;
          k = static_cast<int> (p);
          k = MIN(k, ifile->nr - 2);
          k = MAX(k, 2);
          p -= k;
          p = MIN(p, 2.0);
          cof1 = -sixth*p*(p - 1.0)*(p - 2.0);
          cof2 = 0.5*(p*p - 1.0)*(p - 2.0);
          cof3 = -0.5*p*(p + 1.0)*(p - 2.0);
          cof4 = sixth*p*(p*p - 1.0);
          zri = cof1*ifile->zr[k - 1] + cof2*ifile->zr[k] +
            cof3*ifile->zr[k + 1] + cof4*ifile->zr[k + 2];

          p = r / jfile->dr + 1.0;
          k = static_cast<int> (p);
          k = MIN(k, jfile->nr - 2);
          k = MAX(k, 2);
          p -= k;
          p = MIN(p, 2.0);
          cof1 = -sixth*p*(p - 1.0)*(p - 2.0);
          cof2 = 0.5*(p*p - 1.0)*(p - 2.0);
          cof3 = -0.5*p*(p + 1.0)*(p - 2.0);
          cof4 = sixth*p*(p*p - 1.0);
          zrj = cof1*jfile->zr[k - 1] + cof2*jfile->zr[k] +
            cof3*jfile->zr[k + 1] + cof4*jfile->zr[k + 2];

          z2r[n][m] = 27.2*0.529 * zri*zrj;
        }
        n++;
      }
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

/* ---------------------------------------------------------------------- */

void PairCACEAM::array2spline()
{
  rdr = 1.0 / dr;
  rdrho = 1.0 / drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline, nfrho, nrho + 1, 7, "pair:frho");
  memory->create(rhor_spline, nrhor, nr + 1, 7, "pair:rhor");
  memory->create(z2r_spline, nz2r, nr + 1, 7, "pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho, drho, frho[i], frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr, dr, rhor[i], rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr, dr, z2r[i], z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairCACEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
    spline[m][5] = ((spline[m - 2][6] - spline[m + 2][6]) +
      8.0*(spline[m + 1][6] - spline[m - 1][6])) / 12.0;

  for (int m = 1; m <= n - 1; m++) {
    spline[m][4] = 3.0*(spline[m + 1][6] - spline[m][6]) -
      2.0*spline[m][5] - spline[m + 1][5];
    spline[m][3] = spline[m][5] + spline[m + 1][5] -
      2.0*(spline[m + 1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5] / delta;
    spline[m][1] = 2.0*spline[m][4] / delta;
    spline[m][0] = 3.0*spline[m][3] / delta;
  }
}

/* ----------------------------------------------------------------------
grab n values from file fp and put them in list
values can be several to a line
only called by proc 0
------------------------------------------------------------------------- */

void PairCACEAM::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line, MAXLINE, fptr);
    ptr = strtok(line, " \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL, " \t\n\r\f"))) list[i++] = atof(ptr);
  }
}



//---------------------------------
void PairCACEAM::swap_eam(double *fp_caller, double **fp_caller_hold)
{
  double *tmp = fp;
  fp = fp_caller;
  *fp_caller_hold = tmp;
}

/* ---------------------------------------------------------------------- */

void *PairCACEAM::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *)scale;
  return NULL;
}




//-----------------------------------------------------------------------


void PairCACEAM::force_densities(int iii, double s, double t, double w, double coefficients,
  double &force_densityx, double &force_densityy, double &force_densityz) {

double delx,dely,delz;
double r2inv;
double r6inv;
double shape_func;
double shape_func2;
double unit_cell_mapped[3];
double scanning_unit_cell[3];
double forcelj,factor_lj,fpair;
int *type = atom->type;
double unit_cell[3];
double distancesq;
double current_position[3];
double scan_position[3];
double rcut;

int nodes_per_element;
int *nodes_count_list = atom->nodes_per_element_list;

//equivalent isoparametric cutoff range for a cube of rcut

unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
unit_cell_mapped[2] = 2 / double(current_element_scale[2]);

unit_cell[0] = s;
unit_cell[1] = t;
unit_cell[2] = w;

//scan the surrounding unit cell locations in a cartesian grid
//of isoparametric space until the cutoff is exceeded
//for each grid scan

int distanceflag=0;
    current_position[0]=0;
    current_position[1]=0;
    current_position[2]=0;
  
  if (!atomic_flag) {
    nodes_per_element = nodes_count_list[current_element_type];
    for (int kkk = 0; kkk < nodes_per_element; kkk++) {
      shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kkk + 1);
      //shape_func = (this->*shape_functions[kkk])(unit_cell[0], unit_cell[1], unit_cell[2]);
      current_position[0] += current_nodal_positions[kkk][0] * shape_func;
      current_position[1] += current_nodal_positions[kkk][1] * shape_func;
      current_position[2] += current_nodal_positions[kkk][2] * shape_func;
    }
  }
  else {
    current_position[0] = s;
    current_position[1] = t;
    current_position[2] = w;
  }
  
  rcut = cut_global_s;
  int origin_type = type_array[poly_counter];
  

  int listtype;
  int scan_type, scan_type2;
  int listindex;
  int poly_index;
  double force_contribution[3];
  int element_index;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int neigh_max_inner = inner_quad_lists_counts[iii][neigh_quad_counter];
  int neigh_max_outer = outer_quad_lists_counts[iii][neigh_quad_counter];
  int itype, jtype, ktype;
  double rsq, r, p, rhoip, rhojp, z2, z2p, recip, phip, psip, phi;
  double *coeff;
  int m;
  //if(update->ntimestep==1)
  //timer->stamp(Timer::CAC_INIT);
  if(neigh_max_inner>local_inner_max){
  memory->grow(rho, neigh_max_inner + 1+EXPAND, "Pair_CAC_eam:rho");
  memory->grow(fp, neigh_max_inner + 1+EXPAND, "Pair_CAC_eam:fp");
  memory->grow(inner_neighbor_types, neigh_max_inner+EXPAND, "Pair_CAC_eam:inner_neighbor_types");
  memory->grow(inner_neighbor_coords, neigh_max_inner+EXPAND, 3, "Pair_CAC_eam:inner_neighbor_coords");
  local_inner_max=neigh_max_inner+EXPAND;
  }
  if(neigh_max_outer>local_outer_max){
  memory->grow(outer_neighbor_coords, neigh_max_outer+EXPAND, 3, "Pair_CAC_eam:outer_neighbor_coords");
  memory->grow(outer_neighbor_types, neigh_max_outer+EXPAND, "Pair_CAC_eam:outer_neighbor_types");
  local_outer_max=neigh_max_outer+EXPAND;
  }
  
  for (int l = 0; l < neigh_max_inner+1; l++) {
    rho[l] = 0;
    fp[l] = 0;
  }
  tagint itag, jtag;
  double  rsq1, rsq2;
  double delr1[3], delr2[3], fj[3], fk[3];
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  jlist = firstneigh[iii];
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  origin_type = type_array[poly_counter];
  double inner_scan_position[3];
  int **inner_quad_indices = inner_quad_lists_index[iii][neigh_quad_counter];
  int **outer_quad_indices = outer_quad_lists_index[iii][neigh_quad_counter];
  //precompute virtual neighbor atom locations
    
  
  for (int l = 0; l < neigh_max_inner; l++) {
    //listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
    element_index = inner_quad_indices[l][0];
    poly_index = inner_quad_indices[l][1];
    inner_neighbor_types[l] = node_types[element_index][poly_index];

  }
  for (int l = 0; l < neigh_max_outer; l++) {
    //listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
    element_index = outer_quad_indices[l][0];
    poly_index = outer_quad_indices[l][1];
    outer_neighbor_types[l] = node_types[element_index][poly_index];

  }
  //compute virtual neighbor positions at the current timestep
  interpolation(iii);
  //two body accumulation of electron densities to quadrature site
  for (int l = 0; l < neigh_max_inner; l++) {

    scan_type = inner_neighbor_types[l];
    scan_position[0] = inner_neighbor_coords[l][0];
    scan_position[1] = inner_neighbor_coords[l][1];
    scan_position[2] = inner_neighbor_coords[l][2];


    delx = current_position[0] - scan_position[0];
    dely = current_position[1] - scan_position[1];
    delz = current_position[2] - scan_position[2];
    distancesq = delx*delx + dely*dely + delz*delz;

    
    if (distancesq >= cutforcesq) continue;

    
    p = sqrt(distancesq)*rdr + 1.0;
    m = static_cast<int> (p);
    m = MIN(m, nr - 1);
    p -= m;
    p = MIN(p, 1.0);
    coeff = rhor_spline[type2rhor[scan_type][origin_type]][m];
    rho[0] += ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];
  }
  //two body accumulation of electron densities to inner neighbor list
  //loops over neighbors of inner neighbor list using both the inner and outer neighbor list

  for (int l = 0; l < neigh_max_inner; l++) {

    scan_type = inner_neighbor_types[l];
    inner_scan_position[0] = inner_neighbor_coords[l][0];
    inner_scan_position[1] = inner_neighbor_coords[l][1];
    inner_scan_position[2] = inner_neighbor_coords[l][2];


    delr1[0] = inner_scan_position[0] - current_position[0];
    delr1[1] = inner_scan_position[1] - current_position[1];
    delr1[2] = inner_scan_position[2] - current_position[2];
    

    

    rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
    if (rsq1 >= cutforcesq) continue;
    //add density contribution due to origin atom of the neighborlist
    p = sqrt(rsq1)*rdr + 1.0;
    m = static_cast<int> (p);
    m = MIN(m, nr - 1);
    p -= m;
    p = MIN(p, 1.0);
    coeff = rhor_spline[type2rhor[origin_type][scan_type]][m];
    rho[l+1] += ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];



    for (int k = 0; k < neigh_max_inner; k++) {
      if(l==k) continue;
      scan_type2 = inner_neighbor_types[k];
      scan_position[0] = inner_neighbor_coords[k][0];
      scan_position[1] = inner_neighbor_coords[k][1];
      scan_position[2] = inner_neighbor_coords[k][2];

      delr2[0] = scan_position[0] - inner_scan_position[0];
      delr2[1] = scan_position[1] - inner_scan_position[1];
      delr2[2] = scan_position[2] - inner_scan_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= cutforcesq) continue;
      //add density contribution due to other inner list atoms
      p = sqrt(rsq2)*rdr + 1.0;
      m = static_cast<int> (p);
      m = MIN(m, nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      coeff = rhor_spline[type2rhor[scan_type2][scan_type]][m];
      rho[l + 1] += ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];

    }
    for (int k = 0; k < neigh_max_outer; k++) {
      

      scan_type2 = outer_neighbor_types[k];
      scan_position[0] = outer_neighbor_coords[k][0];
      scan_position[1] = outer_neighbor_coords[k][1];
      scan_position[2] = outer_neighbor_coords[k][2];

      delr2[0] = scan_position[0] - inner_scan_position[0];
      delr2[1] = scan_position[1] - inner_scan_position[1];
      delr2[2] = scan_position[2] - inner_scan_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= cutforcesq) continue;
      //add density contribution due to other outer list atoms
      p = sqrt(rsq2)*rdr + 1.0;
      m = static_cast<int> (p);
      m = MIN(m, nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      coeff = rhor_spline[type2rhor[scan_type2][scan_type]][m];
      rho[l + 1] += ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];
    }

  }
  //compute derivative of the embedding energy for the origin atom
  scan_type = node_types[element_index][poly_index];
  p = rho[0] * rdrho + 1.0;
  m = static_cast<int> (p);
  m = MAX(1, MIN(m, nrho - 1));
  p -= m;
  p = MIN(p, 1.0);
  coeff = frho_spline[type2frho[origin_type]][m];
  fp[0] = (coeff[0] * p + coeff[1])*p + coeff[2];

  if (quad_eflag){
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho[0] > rhomax) phi += fp[0] * (rho[0]-rhomax);
     phi *= scale[origin_type][origin_type];
      quadrature_energy += phi;
  }

  //compute derivative of the embedding energy for all atoms in the inner neighborlist
  for (int l = 0; l < neigh_max_inner; l++) {

    scan_type = inner_neighbor_types[l];;
    p = rho[l+1] * rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1, MIN(m, nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = frho_spline[type2frho[scan_type]][m];
    fp[l+1] = (coeff[0] * p + coeff[1])*p + coeff[2];

  }
  //compute force contribution
  for (int l = 0; l < neigh_max_inner; l++) {


    scan_type = inner_neighbor_types[l];;
    scan_position[0] = inner_neighbor_coords[l][0];
    scan_position[1] = inner_neighbor_coords[l][1];
    scan_position[2] = inner_neighbor_coords[l][2];


    delx = current_position[0] - scan_position[0];
    dely = current_position[1] - scan_position[1];
    delz = current_position[2] - scan_position[2];
    distancesq = delx*delx + dely*dely + delz*delz;


    if (distancesq >= cutforcesq) continue;
    
    r = sqrt(distancesq);
    p = r*rdr + 1.0;
    m = static_cast<int> (p);
    m = MIN(m, nr - 1);
    p -= m;
    p = MIN(p, 1.0);

    // rhoip = derivative of (density at atom j due to atom i)
    // rhojp = derivative of (density at atom i due to atom j)
    // phi = pair potential energy
    // phip = phi'
    // z2 = phi * r
    // z2p = (phi * r)' = (phi' r) + phi
    // psip needs both fp[i] and fp[j] terms since r_ij appears in two
    //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
    //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
    // scale factor can be applied by thermodynamic integration

    coeff = rhor_spline[type2rhor[origin_type][scan_type]][m];
    rhoip = (coeff[0] * p + coeff[1])*p + coeff[2];
    coeff = rhor_spline[type2rhor[scan_type][origin_type]][m];
    rhojp = (coeff[0] * p + coeff[1])*p + coeff[2];
    coeff = z2r_spline[type2z2r[origin_type][scan_type]][m];
    z2p = (coeff[0] * p + coeff[1])*p + coeff[2];
    z2 = ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];

    recip = 1.0 / r;
    phi = z2*recip;
    phip = z2p*recip - phi*recip;
    psip = fp[0] * rhojp + fp[l+1] * rhoip + phip;
    fpair = -scale[origin_type][scan_type] * psip*recip;

    force_densityx += delx*fpair;
    force_densityy += dely*fpair;
    force_densityz += delz*fpair;
    if(atom->CAC_virial){
    virial_density[0] += 0.5*delx*delx*fpair;
    virial_density[1] += 0.5*dely*dely*fpair;
    virial_density[2] += 0.5*delz*delz*fpair;
    virial_density[3] += 0.5*delx*dely*fpair;
    virial_density[4] += 0.5*delx*delz*fpair;
    virial_density[5] += 0.5*dely*delz*fpair;
    }
    if (quad_eflag) 
      quadrature_energy += 0.5*scale[origin_type][scan_type]*phi;

  }
   
//end of scanning loop

}

//--------------------------------------------------------------------------





