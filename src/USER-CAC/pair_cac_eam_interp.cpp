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

#include "pair_cac_eam_interp.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "npair_cac.h"
#include "utils.h"
#include "tokenizer.h"
#include "potential_file_reader.h"


//#include "math_extra.h"
#define PLANE_EPSILON  1e-6 //error tolerance for surface flux calculation
#define MAXNEIGHIN  30
#define MAXLINE 1024
#define DELTA 4
#define EXPAND 10
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairCACEAMInterp::PairCACEAMInterp(LAMMPS *lmp) : PairCAC(lmp)
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
  add_rho = NULL;
  add_fp = NULL;
  nmax = 0;
  
  flux_enable = pre_force_flag=1;
  quad_electron_densities = NULL;
  max_density = 0;
}

/* ---------------------------------------------------------------------- */

PairCACEAMInterp::~PairCACEAMInterp() {

  if (copymode) return;

  memory->destroy(rho);
  memory->destroy(fp);
  if(atom->cac_flux_flag){
    memory->destroy(add_rho);
    memory->destroy(add_fp);
  }

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
    memory->destroy(inner_neighbor_coords);
    memory->destroy(inner_neighbor_types);
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

void PairCACEAMInterp::allocate()
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
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACEAMInterp::coeff(int narg, char **arg) {
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

double PairCACEAMInterp::init_one(int i, int j) {

  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  if (funcfl) {
    cut_global_s = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cut_global_s = MAX(cut_global_s, funcfl[m].cut);
  }
  else if (setfl) cut_global_s = setfl->cut;
  else if (fs) cut_global_s = fs->cut;
    
  cutforcesq = cut_global_s*cut_global_s;

  return cut_global_s;
}

/* ---------------------------------------------------------------------- */


void PairCACEAMInterp::init_style()
{
  PairCAC::init_style();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->ghost_quad_flag = 1;
  atom->sector_flag = sector_flag = 1;
  // convert read-in file(s) to arrays and spline them
  file2array();
  array2spline();
}

/* ---------------------------------------------------------------------- */

void PairCACEAMInterp::read_file(char *filename)
{
    Funcfl *file = &funcfl[nfuncfl-1];

  // read potential file
  if(comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "eam", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);
    try {
      reader.skip_line();

      ValueTokenizer values = reader.next_values(2);
      values.next_int(); // ignore
      file->mass = values.next_double();

      values = reader.next_values(5);
      file->nrho = values.next_int();
      file->drho = values.next_double();
      file->nr   = values.next_int();
      file->dr   = values.next_double();
      file->cut  = values.next_double();

      if ((file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
        error->one(FLERR,"Invalid EAM potential file");

      memory->create(file->frho, (file->nrho+1), "pair:frho");
      memory->create(file->rhor, (file->nr+1), "pair:rhor");
      memory->create(file->zr, (file->nr+1), "pair:zr");

      reader.next_dvector(&file->frho[1], file->nrho);
      reader.next_dvector(&file->zr[1], file->nr);
      reader.next_dvector(&file->rhor[1], file->nr);

      if (unit_convert) {
        const double sqrt_conv = sqrt(conversion_factor);
        for (int i = 1; i <= file->nrho; ++i)
          file->frho[i] *= conversion_factor;
        for (int j = 1; j <= file->nr; ++j)
          file->zr[j] *= sqrt_conv;
      }
    } catch (TokenizerException & e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(&file->mass, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  if(comm->me != 0) {
    memory->create(file->frho, (file->nrho+1), "pair:frho");
    memory->create(file->rhor, (file->nr+1), "pair:rhor");
    memory->create(file->zr, (file->nr+1), "pair:zr");
  }

  MPI_Bcast(&file->frho[1], file->nrho, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->zr[1], file->nr, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->rhor[1], file->nr, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
convert read-in funcfl potential(s) to standard array format
interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairCACEAMInterp::file2array()
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

void PairCACEAMInterp::array2spline()
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

void PairCACEAMInterp::interpolate(int n, double delta, double *f, double **spline)
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

void PairCACEAMInterp::grab(FILE *fptr, int n, double *list)
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
void PairCACEAMInterp::swap_eam(double *fp_caller, double **fp_caller_hold)
{
    double *tmp = fp;
    fp = fp_caller;
    *fp_caller_hold = tmp;
}

/* ---------------------------------------------------------------------- */

void *PairCACEAMInterp::extract(const char *str, int &dim)
{
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *)scale;
    return NULL;
}

//-----------------------------------------------------------------------
//compute nodal and atomic electron densities.
void PairCACEAMInterp::pre_force_densities() {
  double **x = atom->x;
  double ****nodal_positions= atom->nodal_positions;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  //grow array storing electron densities at quadrature points
  atom->npair_cac->allocate_quad_double(quad_electron_densities, 1);

  for (int i = 0; i < atom->nlocal; i++) {
    atomic_flag = 0;
    current_element_type = element_type[i];
    current_element_scale = element_scale[i];
    current_poly_count = poly_count[i];
    type_array = node_types[i];
      
    current_x = x[i];
 
    //determine element type
    if (current_element_type == 0) {
    atomic_flag = 1;
    }
    compute_electron_densities(i);
      
  }
  //test if maxexchange atom is large enough for quadrature count of each element
  comm->increase_max_atom(atom->max_quad_per_element*atom->maxpoly);

  //commmunicate nodal electron densities of atoms/elements that are ghosts of other tasks
  comm->forward_comm_pair(this);
}

//-----------------------------------------------------------------------

void PairCACEAMInterp::compute_electron_densities(int i) {

  int *nodes_count_list = atom->nodes_per_element_list;
  double ****nodal_positions= atom->nodal_positions;
  double coefficients;
  int nodes_per_element;
  double s, t, w;
  double sq, tq, wq;
  nodes_per_element = nodes_count_list[current_element_type];

  //sum over quadrature points to compute force density
  for (int quad_loop=0; quad_loop < quadrature_counts[i] ; quad_loop++){
    if(atom->cac_flux_flag==2&&flux_compute) quad_flux_flag = 1;
    else if(atom->cac_flux_flag==1&&flux_compute){
      //check to make sure this quadrature point is the closest one to a node
      if(quadrature_point_data[qi + quad_loop][7]) quad_flux_flag = 1;
      else quad_flux_flag = 0;
    }
    else quad_flux_flag = 0;

    if(!atomic_flag){
      s = quadrature_point_data[qi][0];
      t = quadrature_point_data[qi][1];
      w = quadrature_point_data[qi][2];
      sq = quadrature_point_data[qi][3];
      tq = quadrature_point_data[qi][4];
      wq = quadrature_point_data[qi][5];
      coefficients = quadrature_point_data[qi][6];
    }
    for (poly_counter = 0; poly_counter < current_poly_count; poly_counter++){
      current_nodal_positions = nodal_positions[i][poly_counter];
      if(!atomic_flag)
        quad_electron_density(i, s, t, w);
      else
        quad_electron_density(i, current_x[0], current_x[1], current_x[2]);
      pqi++;
    }
  qi++;
  }
    
}

//-----------------------------------------------------------------------

void PairCACEAMInterp::quad_electron_density(int i, double s, double t, double w) 
{ 
  double delx,dely,delz;
  int *type = atom->type;
  double distancesq;
  double scan_position[3];
  int nodes_per_element, neigh_max_add, m;
  int *nodes_count_list = atom->nodes_per_element_list;
  int origin_type = type_array[poly_counter];
  int listtype;
  int scan_type, scan_type2;
  int listindex;
  int poly_index;
  int element_index;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int jtype, ktype;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  double rsq, r, p, rhoip, rhojp, z2, z2p, recip, phip, psip, phi;
  double *coeff;
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  origin_type = type_array[poly_counter];

  //allocate arrays specific to this potential here; must be done before allocate_quad_memory.
  if(neigh_max_inner>local_inner_max){
    memory->grow(rho, neigh_max_inner + 1 + EXPAND, "Pair_CAC_eam_interp:rho");
    memory->grow(fp, neigh_max_inner + 1 + EXPAND, "Pair_CAC_eam_interp:fp");
  }
  if(quad_flux_flag){
    neigh_max_add = add_quad_lists_counts[pqi];
    if(neigh_max_add>local_add_max){
    memory->grow(add_rho, neigh_max_add + 1 + EXPAND, "Pair_CAC_eam_interp:add_rho");
    memory->grow(add_fp, neigh_max_add + 1 + EXPAND, "Pair_CAC_eam_interp:add_fp");
  }
  }
  //allocate arrays that store neighbor information around just this quadrature point
  allocate_quad_memory();
  
  rho[0] = 0;
  
  //set virtual neighbor types, etc.
  init_quad_arrays();
  //precompute virtual neighbor atom locations
  //interpolate virtual atom coordinates from shape functions corresponding to unit cells
  interpolation(i,s,t,w);
  
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
   quad_electron_densities[pqi][0] = rho[0];      
  }

//-----------------------------------------------------------------------

void PairCACEAMInterp::force_densities(int iii, double s, double t, double w, double coefficients,
    double &force_densityx, double &force_densityy, double &force_densityz) {

  double delx,dely,delz;
  double fpair;
  int *type = atom->type;
  double distancesq;
  double scan_position[3];
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  int nodes_per_element, neigh_max_add;
  int *nodes_count_list = atom->nodes_per_element_list;
  int origin_type = type_array[poly_counter];
    
  int listtype;
  int scan_type, scan_type2;
  int listindex;
  int poly_index;
  double force_contribution[3];
  int element_index;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int jtype, ktype;
  int quad_index;
  double rsq, r, p, rhoip, rhojp, z2, z2p, recip, phip, psip, phi;
  double *coeff;
  int m;
  
  for (int l = 0; l < neigh_max_inner+1; l++) {
    rho[l] = 0;
    fp[l] = 0;
  }
    
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  origin_type = type_array[poly_counter];

  //assign electron density to quadrature origin
  rho[0] = quad_electron_densities[pqi][0];
  //precompute virtual/real neighbor atom densities
  for (int l = 0; l < neigh_max_inner; l++) {
    quad_index = inner_quad_indices[l][2];
    rho[l+1] = quad_electron_densities[quad_index][0];
  }

  //density of additional neighbors for flux computation
  if(quad_flux_flag){
    neigh_max_add = add_quad_lists_counts[pqi];
    for (int l = 0; l < neigh_max_add; l++) {
      quad_index = inner_quad_lists_index[pqi][l][2];
      add_rho[l] = quad_electron_densities[quad_index][0];
    }
  }

  //set virtual neighbor types, etc.
  init_quad_arrays();
  //interpolate virtual atom coordinates from shape functions corresponding to unit cells
  interpolation(iii,s,t,w);

  //compute derivative of the embedding energy for the origin atom
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

  //compute derivatives for additional neighbors in flux calculation
  if(quad_flux_flag){
    for (int l = 0; l < neigh_max_add; l++) {
      scan_type = add_neighbor_types[l];
      p = add_rho[l] * rdrho + 1.0;
      m = static_cast<int> (p);
      m = MAX(1, MIN(m, nrho - 1));
      p -= m;
      p = MIN(p, 1.0);
      coeff = frho_spline[type2frho[scan_type]][m];
      add_fp[l] = (coeff[0] * p + coeff[1])*p + coeff[2];
    }
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
      
    //cac flux contribution due to current quadrature point and neighbor pair interactions
    if(quad_flux_flag){
      current_quad_flux(l,delx*fpair,dely*fpair,delz*fpair);
    }
  }
  
  //additional cac flux contributions due to neighbors interacting with neighbors
  //  in the vicinity of this quadrature point
  if (quad_flux_flag) {
    //compute_intersections();
    quad_neigh_flux();
  }
}

/* ---------------------------------------------------------------------- 
 Compute the cac flux density due to virtual neighbors around a quadrature point
---------------------------------------------------------------------- */

void PairCACEAMInterp::quad_neigh_flux(){
  int all_neigh, is, isl, normal_flag, sign, scan_type1, scan_type2, index, jindex;
  int m, dim1, dim2, intersection_flag, intersection_count, icontrib;
  double fpair, interaction_forceij[3], interaction_forceji[3], delxa[3];
  double scan_position1[3], scan_position2[3], q1, q2, delx, dely, delz, distancesq, fluxdistsq;  
  double intersection_point[3], proj, lparam, planecoord, plane_limits[2][2];
  double vix, viy, viz, vjx, vjy, vjz, interactionx, interactiony, interactionz;
  double r, p, rho1, rho2, fp1, fp2, *coeff;
  double rhoip, rhojp, z2, z2p, recip, phip, psip, phi;
  double *box_center = atom->box_center;
  double *box_size = atom->box_size;
  int neigh_max = inner_quad_lists_counts[pqi];
  int neigh_add = add_quad_lists_counts[pqi];
  double cut_add = atom->cut_add;
  double fluxcutsq = (cut_global_s+cut_add)*(cut_global_s+cut_add);
  all_neigh = neigh_max + neigh_add;
  //icontrib = 0;

  //determine which of the 6 planes of the atom box are intersected by a given i-j pair
  for(int ineigh=0; ineigh < all_neigh; ineigh++){
    if(ineigh<neigh_max){
      scan_type1 = inner_neighbor_types[ineigh];
      rho1 = rho[ineigh+1];
      fp1 = fp[ineigh+1];
      scan_position1[0] = inner_neighbor_coords[ineigh][0];
      scan_position1[1] = inner_neighbor_coords[ineigh][1];
      scan_position1[2] = inner_neighbor_coords[ineigh][2];
      vix = inner_neighbor_velocities[ineigh][0];
      viy = inner_neighbor_velocities[ineigh][1];
      viz = inner_neighbor_velocities[ineigh][2];
    }
    else{
      index = ineigh-neigh_max;
      scan_type1 = add_neighbor_types[index];
      rho1 = add_rho[index];
      fp1 = add_fp[index];
      scan_position1[0] = add_neighbor_coords[index][0];
      scan_position1[1] = add_neighbor_coords[index][1];
      scan_position1[2] = add_neighbor_coords[index][2];
      vix = add_neighbor_velocities[index][0];
      viy = add_neighbor_velocities[index][1];
      viz = add_neighbor_velocities[index][2];
    }

    delxa[0] = delx = current_position[0] - scan_position1[0];
    delxa[1] = dely = current_position[1] - scan_position1[1];
    delxa[2] = delz = current_position[2] - scan_position1[2];
    fluxdistsq = delx*delx + dely*dely + delz*delz;
    if(fluxdistsq > fluxcutsq) continue;
    
    for(int jneigh=ineigh+1; jneigh < all_neigh; jneigh++){
      intersection_count = 0;
      if(jneigh<neigh_max){
        scan_type2 = inner_neighbor_types[jneigh];
        rho2 = rho[jneigh+1];
        fp2 = fp[jneigh+1];
        scan_position2[0] = inner_neighbor_coords[jneigh][0];
        scan_position2[1] = inner_neighbor_coords[jneigh][1];
        scan_position2[2] = inner_neighbor_coords[jneigh][2];
        vjx = inner_neighbor_velocities[jneigh][0];
        vjy = inner_neighbor_velocities[jneigh][1];
        vjz = inner_neighbor_velocities[jneigh][2];
      }
      else{
        jindex = jneigh-neigh_max;
        scan_type2 = add_neighbor_types[jindex];
        rho2 = add_rho[jindex];
        fp2 = add_fp[jindex];
        scan_position2[0] = add_neighbor_coords[jindex][0];
        scan_position2[1] = add_neighbor_coords[jindex][1];
        scan_position2[2] = add_neighbor_coords[jindex][2];
        vjx = add_neighbor_velocities[jindex][0];
        vjy = add_neighbor_velocities[jindex][1];
        vjz = add_neighbor_velocities[jindex][2];
      }

      delxa[0] = delx = scan_position1[0] - scan_position2[0];
      delxa[1] = dely = scan_position1[1] - scan_position2[1];
      delxa[2] = delz = scan_position1[2] - scan_position2[2];
      distancesq = delx*delx + dely*dely + delz*delz;
      if(distancesq > cutsq[scan_type1][scan_type2]) continue;
      
      //compute pair force
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

      coeff = rhor_spline[type2rhor[scan_type1][scan_type2]][m];
      rhoip = (coeff[0] * p + coeff[1])*p + coeff[2];
      coeff = rhor_spline[type2rhor[scan_type2][scan_type1]][m];
      rhojp = (coeff[0] * p + coeff[1])*p + coeff[2];
      coeff = z2r_spline[type2z2r[scan_type1][scan_type2]][m];
      z2p = (coeff[0] * p + coeff[1])*p + coeff[2];
      z2 = ((coeff[3] * p + coeff[4])*p + coeff[5])*p + coeff[6];

      recip = 1.0 / r;
      phi = z2*recip;
      phip = z2p*recip - phi*recip;
      psip = fp1 * rhojp + fp2 * rhoip + phip;
      fpair = -scale[scan_type1][scan_type2] * psip*recip;

      interaction_forceij[0] = -delx*fpair;
      interaction_forceij[1] = -dely*fpair;
      interaction_forceij[2] = -delz*fpair;
      
      for(int isl=0; isl < 2*domain->dimension; isl++){
        is = isl/2;
        if(is==0){
          dim1 = 1;
          dim2 = 2;
        }
        if(is==1){
          dim1 = 0;
          dim2 = 2;
        }
        if(is==2){
          dim1 = 0;
          dim2 = 1;
        }

        //test negative and positive sides of the box dimension
        if(isl%2==0) planecoord = current_position[is]-box_size[is]/2 + box_center[is];
        else planecoord = current_position[is]+box_size[is]/2 + box_center[is];
        plane_limits[0][0] = current_position[dim1]-box_size[dim1]/2 + box_center[dim1]-PLANE_EPSILON;
        plane_limits[0][1] = current_position[dim1]+box_size[dim1]/2 + box_center[dim1]+PLANE_EPSILON;
        plane_limits[1][0] = current_position[dim2]-box_size[dim2]/2 + box_center[dim2]-PLANE_EPSILON;
        plane_limits[1][1] = current_position[dim2]+box_size[dim2]/2 + box_center[dim2]+PLANE_EPSILON;

        intersection_flag = 1;
        //compute perpendicular projection of the line connecting i and j
        proj = scan_position1[is]-planecoord;

        //test if i-j normal coordinates are on opposing sides of the plane
        if((proj<0&&scan_position2[is]<planecoord)||((proj>0)&&scan_position2[is]>planecoord)) intersection_flag = 0;

        //use the ratio between this projection and the i-j displacement normal to the plane
        //to define the line parameter (0-1) at the point of intersection
        lparam = proj/(delxa[is]);
        if(delxa[is]==0) intersection_flag = 0;

        //use line parameter to extrapolate the possible intersection point between i-j
        intersection_point[dim1] = scan_position2[dim1]+delxa[dim1]*(1-lparam);
        intersection_point[dim2] = scan_position2[dim2]+delxa[dim2]*(1-lparam);

        //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
        if(intersection_point[dim1]<=plane_limits[0][0]||intersection_point[dim1]>plane_limits[0][1]) intersection_flag = 0;
        if(intersection_point[dim2]<=plane_limits[1][0]||intersection_point[dim2]>plane_limits[1][1]) intersection_flag = 0;
        if(intersection_flag){
          intersection_count++;
          if(isl%2==0) normal_flag = 1;
          else normal_flag = -1;

          if(scan_position1[is]<planecoord) sign=-normal_flag;
          else sign=normal_flag;
        
          if(isl==0){
            //flux_contrib[icontrib][0] = -interaction_forceij[0]*sign;
            //flux_contrib[icontrib][1] = ineigh;
            //flux_contrib[icontrib][2] = jneigh;
            //flux_contrib[icontrib][3] = isl;
            //icontrib++;
          }
          flux_density[4*isl] += (interaction_forceij[0]*(vix+vjx) + 
          interaction_forceij[1]*(viy+vjy)+interaction_forceij[2]*(viz+vjz))*sign;
          flux_density[4*isl+1] -= interaction_forceij[0]*sign;
          flux_density[4*isl+2] -= interaction_forceij[1]*sign;
          flux_density[4*isl+3] -= interaction_forceij[2]*sign;
        
        }

        //can intersect with box at most twice
        if(intersection_count==2)
          break;
  
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairCACEAMInterp::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;
  double dx,dy,dz;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;

  m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      for (int quad_index = e2quad_index[j]; quad_index < e2quad_index[j]+quadrature_counts[j]*poly_count[j]; quad_index++)
      {
        buf[m++] = quad_electron_densities[quad_index][0];
      }
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairCACEAMInterp::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (int quad_index = e2quad_index[i]; quad_index < e2quad_index[i]+quadrature_counts[i]*poly_count[i]; quad_index++)
    {
      quad_electron_densities[quad_index][0] = buf[m++];
    }
  }
}
