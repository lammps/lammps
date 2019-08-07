/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxc_lookup.h"
#include <mpi.h>
#include <cstdlib>
#include "reaxc_defs.h"
#include "reaxc_nonbonded.h"
#include "reaxc_tool_box.h"

void Tridiagonal_Solve( const double *a, const double *b,
                        double *c, double *d, double *x, unsigned int n){
  int i;
  double id;

  c[0] /= b[0];        /* Division by zero risk. */
  d[0] /= b[0];        /* Division by zero would imply a singular matrix. */
  for(i = 1; i < n; i++){
    id = (b[i] - c[i-1] * a[i]);  /* Division by zero risk. */
    c[i] /= id;                /* Last value calculated is redundant. */
    d[i] = (d[i] - d[i-1] * a[i])/id;
  }

  x[n - 1] = d[n - 1];
  for(i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i + 1];
}


void Natural_Cubic_Spline( LAMMPS_NS::Error* error_ptr, const double *h, const double *f,
                           cubic_spline_coef *coef, unsigned int n )
{
  int i;
  double *a, *b, *c, *d, *v;

  /* allocate space for the linear system */
  a = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  b = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  c = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  d = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  v = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");

  /* build the linear system */
  a[0] = a[1] = a[n-1] = 0;
  for( i = 2; i < n-1; ++i )
    a[i] = h[i-1];

  b[0] = b[n-1] = 0;
  for( i = 1; i < n-1; ++i )
    b[i] = 2 * (h[i-1] + h[i]);

  c[0] = c[n-2] = c[n-1] = 0;
  for( i = 1; i < n-2; ++i )
    c[i] = h[i];

  d[0] = d[n-1] = 0;
  for( i = 1; i < n-1; ++i )
    d[i] = 6 * ((f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1]);

  v[0] = 0;
  v[n-1] = 0;
  Tridiagonal_Solve( &(a[1]), &(b[1]), &(c[1]), &(d[1]), &(v[1]), n-2 );

  for( i = 1; i < n; ++i ){
    coef[i-1].d = (v[i] - v[i-1]) / (6*h[i-1]);
    coef[i-1].c = v[i]/2;
    coef[i-1].b = (f[i]-f[i-1])/h[i-1] + h[i-1]*(2*v[i] + v[i-1])/6;
    coef[i-1].a = f[i];
  }

  sfree(error_ptr,  a, "cubic_spline:a" );
  sfree(error_ptr,  b, "cubic_spline:b" );
  sfree(error_ptr,  c, "cubic_spline:c" );
  sfree(error_ptr,  d, "cubic_spline:d" );
  sfree(error_ptr,  v, "cubic_spline:v" );
}



void Complete_Cubic_Spline( LAMMPS_NS::Error* error_ptr, const double *h, const double *f, double v0, double vlast,
                            cubic_spline_coef *coef, unsigned int n )
{
  int i;
  double *a, *b, *c, *d, *v;

  /* allocate space for the linear system */
  a = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  b = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  c = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  d = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");
  v = (double*) smalloc(error_ptr,  n * sizeof(double), "cubic_spline:a");

  /* build the linear system */
  a[0] = 0;
  for( i = 1; i < n; ++i )
    a[i] = h[i-1];

  b[0] = 2*h[0];
  for( i = 1; i < n; ++i )
    b[i] = 2 * (h[i-1] + h[i]);

  c[n-1] = 0;
  for( i = 0; i < n-1; ++i )
    c[i] = h[i];

  d[0] = 6 * (f[1]-f[0])/h[0] - 6 * v0;
  d[n-1] = 6 * vlast - 6 * (f[n-1]-f[n-2]/h[n-2]);
  for( i = 1; i < n-1; ++i )
    d[i] = 6 * ((f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1]);

  Tridiagonal_Solve( &(a[0]), &(b[0]), &(c[0]), &(d[0]), &(v[0]), n );

  for( i = 1; i < n; ++i ){
    coef[i-1].d = (v[i] - v[i-1]) / (6*h[i-1]);
    coef[i-1].c = v[i]/2;
    coef[i-1].b = (f[i]-f[i-1])/h[i-1] + h[i-1]*(2*v[i] + v[i-1])/6;
    coef[i-1].a = f[i];
  }

  sfree(error_ptr,  a, "cubic_spline:a" );
  sfree(error_ptr,  b, "cubic_spline:b" );
  sfree(error_ptr,  c, "cubic_spline:c" );
  sfree(error_ptr,  d, "cubic_spline:d" );
  sfree(error_ptr,  v, "cubic_spline:v" );
}


int Init_Lookup_Tables( reax_system *system, control_params *control,
                        storage *workspace, mpi_datatypes *mpi_data, char * /*msg*/ )
{
  int i, j, r;
  int num_atom_types;
  int existing_types[MAX_ATOM_TYPES], aggregated[MAX_ATOM_TYPES];
  double dr;
  double *h, *fh, *fvdw, *fele, *fCEvd, *fCEclmb;
  double v0_vdw, v0_ele, vlast_vdw, vlast_ele;
  LR_lookup_table ** & LR = system->LR;

  /* initializations */
  v0_vdw = 0;
  v0_ele = 0;
  vlast_vdw = 0;
  vlast_ele = 0;

  num_atom_types = system->reax_param.num_atom_types;
  dr = control->nonb_cut / control->tabulate;
  h = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:h");
  fh = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:fh");
  fvdw = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:fvdw");
  fCEvd = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:fCEvd");
  fele = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:fele");
  fCEclmb = (double*)
    smalloc(system->error_ptr,  (control->tabulate+2) * sizeof(double), "lookup:fCEclmb");

  LR = (LR_lookup_table**)
    scalloc(system->error_ptr,  num_atom_types, sizeof(LR_lookup_table*), "lookup:LR");
  for( i = 0; i < num_atom_types; ++i )
    LR[i] = (LR_lookup_table*)
      scalloc(system->error_ptr,  num_atom_types, sizeof(LR_lookup_table), "lookup:LR[i]");

  for( i = 0; i < MAX_ATOM_TYPES; ++i )
    existing_types[i] = 0;
  for( i = 0; i < system->n; ++i )
    existing_types[ system->my_atoms[i].type ] = 1;

  MPI_Allreduce( existing_types, aggregated, MAX_ATOM_TYPES,
                 MPI_INT, MPI_SUM, mpi_data->world );

  for( i = 0; i < num_atom_types; ++i ) {
    if (aggregated[i]) {
      for( j = i; j < num_atom_types; ++j ) {
        if (aggregated[j]) {
          LR[i][j].xmin = 0;
          LR[i][j].xmax = control->nonb_cut;
          LR[i][j].n = control->tabulate + 2;
          LR[i][j].dx = dr;
          LR[i][j].inv_dx = control->tabulate / control->nonb_cut;
          LR[i][j].y = (LR_data*)
            smalloc(system->error_ptr,  LR[i][j].n * sizeof(LR_data), "lookup:LR[i,j].y");
          LR[i][j].H = (cubic_spline_coef*)
            smalloc(system->error_ptr,  LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].H");
          LR[i][j].vdW = (cubic_spline_coef*)
            smalloc(system->error_ptr,  LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].vdW");
          LR[i][j].CEvd = (cubic_spline_coef*)
            smalloc(system->error_ptr,  LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].CEvd");
          LR[i][j].ele = (cubic_spline_coef*)
            smalloc(system->error_ptr,  LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].ele");
          LR[i][j].CEclmb = (cubic_spline_coef*)
            smalloc(system->error_ptr,  LR[i][j].n*sizeof(cubic_spline_coef),
                     "lookup:LR[i,j].CEclmb");

          for( r = 1; r <= control->tabulate; ++r ) {
            LR_vdW_Coulomb( system, workspace, control, i, j, r * dr, &(LR[i][j].y[r]) );
            h[r] = LR[i][j].dx;
            fh[r] = LR[i][j].y[r].H;
            fvdw[r] = LR[i][j].y[r].e_vdW;
            fCEvd[r] = LR[i][j].y[r].CEvd;
            fele[r] = LR[i][j].y[r].e_ele;
            fCEclmb[r] = LR[i][j].y[r].CEclmb;
          }

          // init the start-end points
          h[r] = LR[i][j].dx;
          v0_vdw = LR[i][j].y[1].CEvd;
          v0_ele = LR[i][j].y[1].CEclmb;
          fh[r] = fh[r-1];
          fvdw[r] = fvdw[r-1];
          fCEvd[r] = fCEvd[r-1];
          fele[r] = fele[r-1];
          fCEclmb[r] = fCEclmb[r-1];
          vlast_vdw = fCEvd[r-1];
          vlast_ele = fele[r-1];

          Natural_Cubic_Spline( control->error_ptr, &h[1], &fh[1],
                                &(LR[i][j].H[1]), control->tabulate+1);

          Complete_Cubic_Spline( control->error_ptr, &h[1], &fvdw[1], v0_vdw, vlast_vdw,
                                 &(LR[i][j].vdW[1]), control->tabulate+1);

          Natural_Cubic_Spline( control->error_ptr, &h[1], &fCEvd[1],
                                &(LR[i][j].CEvd[1]), control->tabulate+1);

          Complete_Cubic_Spline( control->error_ptr, &h[1], &fele[1], v0_ele, vlast_ele,
                                 &(LR[i][j].ele[1]), control->tabulate+1);

          Natural_Cubic_Spline( control->error_ptr, &h[1], &fCEclmb[1],
                                &(LR[i][j].CEclmb[1]), control->tabulate+1);
        } else {
          LR[i][j].n = 0;
        }
      }
    }
  }
  free(h);
  free(fh);
  free(fvdw);
  free(fCEvd);
  free(fele);
  free(fCEclmb);

  return 1;
}


void Deallocate_Lookup_Tables( reax_system *system )
{
  int i, j;
  int ntypes;
  LR_lookup_table ** & LR = system->LR;

  ntypes = system->reax_param.num_atom_types;

  for( i = 0; i < ntypes; ++i ) {
    for( j = i; j < ntypes; ++j )
      if (LR[i][j].n) {
        sfree(system->error_ptr,  LR[i][j].y, "LR[i,j].y" );
        sfree(system->error_ptr,  LR[i][j].H, "LR[i,j].H" );
        sfree(system->error_ptr,  LR[i][j].vdW, "LR[i,j].vdW" );
        sfree(system->error_ptr,  LR[i][j].CEvd, "LR[i,j].CEvd" );
        sfree(system->error_ptr,  LR[i][j].ele, "LR[i,j].ele" );
        sfree(system->error_ptr,  LR[i][j].CEclmb, "LR[i,j].CEclmb" );
      }
    sfree(system->error_ptr,  LR[i], "LR[i]" );
  }
  sfree(system->error_ptr,  LR, "LR" );
}
