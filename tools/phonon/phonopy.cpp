#ifdef FFTW3
#include <map>
#include "phonopy.h"
#include "math.h"
#include "kpath.h"
#include "fftw3.h"

/* ----------------------------------------------------------------------------
 * Class Phonopy is designed to interface with phonopy.
 * ---------------------------------------------------------------------------- */
Phonopy::Phonopy(DynMat *dynmat)
{
   dm = dynmat;
   memory = new Memory();
   sysdim = dm->sysdim;
   fftdim = dm->fftdim;
   fftdim2 = fftdim * fftdim;
   nucell = dm->nucell;
   nx = ny = nz = 5;
   write(1);

   char str[MAXLINE];
   if (count_words(fgets(str,MAXLINE,stdin)) >= 3){
      nx = atoi(strtok(str," \t\n\r\f"));
      ny = atoi(strtok(NULL," \t\n\r\f"));
      nz = atoi(strtok(NULL," \t\n\r\f"));
   }
   if (dm->nx == 1) nx = 1;
   if (dm->ny == 1) ny = 1;
   if (dm->nz == 1) nz = 1;
   if (nx < 1) nx = 1;
   if (ny < 1) ny = 1;
   if (nz < 1) nz = 1;
   npt = nx * ny * nz;
   write(2);

   memory->create(mass, nucell, "Phonopy:mass");
   for (int i = 0; i < nucell; ++i){
      double m = 1./dm->M_inv_sqrt[i];
      mass[i] = m * m;
   }

   memory->create(FC_all, npt, fftdim2, "phonopy:FC_all");

   get_my_FC();

   phonopy();

return;

}

/* ----------------------------------------------------------------------------
 * Deconstructor, free the memory used.
 * ---------------------------------------------------------------------------- */
Phonopy::~Phonopy()
{
   memory->destroy(mass);
   memory->destroy(FC_all);
   delete memory;
   dm = NULL;
return;
}

/* ----------------------------------------------------------------------------
 * Subroutine to write the outputs to screen.
 * ---------------------------------------------------------------------------- */
void Phonopy::write(int flag)
{
   if (flag == 1){  // basic information
      for (int ii = 0; ii < 80; ++ii) printf("="); printf("\n");
      printf("Now to prepare the input files for phonopy.\n");
      printf("The dimension of your present supercell is   : %d x %d x %d.\n", dm->nx, dm->ny, dm->nz);
      printf("The size of the force constant matrix will be: %d x %d.\n", dm->npt*3, dm->npt*3);
      printf("Please define your desired dimension [5 5 5] : ");

   } else if (flag == 2){
      printf("\nThe new dimension of the supercell will be   : %d x %d x %d.\n", nx, ny, nz);
      printf("The size of the force constant matrix is then: %d x %d.\n", npt*3, npt*3);

   } else if (flag == 3){
      printf("\nNow to prepare the phonopy FORCE_CONSTANTS ..."); fflush(stdout);

   } else if (flag == 4){
      printf("Done!\nThe force constants information is extracted and written to FORCE_CONSTANTS,\n");
      printf("the primitive cell is written to POSCAR.primitive, and the input file for\n");
      printf("phonopy band evaluation is written to band.conf.\n");
      printf("One should be able to obtain the phonon band structure after correcting\n");
      printf("the element names in POSCAR.primitive and band.conf by running\n");
      printf("`phonopy --readfc -c POSCAR.primitive -p band.conf`.\n");
      for (int ii = 0; ii < 80; ++ii) printf("-");
      printf("\n***          Remember to change the element names.           ***\n");
#ifdef UseSPG
      for (int ii = 0; ii < 80; ++ii) printf("-");
#endif

   } else if (flag == 5){
      for (int ii = 0; ii < 80; ++ii) printf("="); printf("\n");

   }
return;
}

/* ----------------------------------------------------------------------------
 * Driver to obtain the force constant matrix
 * ---------------------------------------------------------------------------- */
void Phonopy::get_my_FC()
{
   double q[3];
   int ipt = 0;
   for (int ix = 0; ix < nx; ++ix)
   for (int iy = 0; iy < ny; ++iy)
   for (int iz = 0; iz < nz; ++iz){
      q[0] = double(ix)/double(nx);
      q[1] = double(iy)/double(ny);
      q[2] = double(iz)/double(nz);
  
      dm->getDMq(q);

      int ndim = 0;
      for (int idim = 0; idim < fftdim; ++idim)
      for (int jdim = 0; jdim < fftdim; ++jdim){
         FC_all[ipt][ndim] = dm->DM_q[idim][jdim];
         double m = sqrt(mass[idim/sysdim] * mass[jdim/sysdim]);
         FC_all[ipt][ndim].r *= m;
         FC_all[ipt][ndim].i *= m;
         ++ndim;
      }
      ++ipt;
   }
return;
}
       
/* ----------------------------------------------------------------------------
 * Method to write out the force constants and related files for
 * postprocessing with phonopy.
 * ---------------------------------------------------------------------------- */
void Phonopy::phonopy()
{
   // output info
   write(3);

   fftw_complex *in, *out;
   double **fc;
   memory->create(in,  npt, "phonopy:in");
   memory->create(out, npt, "phonopy:in");
   memory->create(fc,  npt, fftdim2, "phonopy:in");

   fftw_plan plan = fftw_plan_dft_3d(nx, ny, nz, in, out, -1, FFTW_ESTIMATE);

   double factor = dm->eml2fc / double(npt);
   for (int idim = 0; idim < fftdim2; ++idim){
      for (int i = 0; i < npt; ++i){
         in[i][0] = FC_all[i][idim].r;
         in[i][1] = FC_all[i][idim].i;
      }
      fftw_execute(plan);
      for (int i = 0; i < npt; ++i) fc[i][idim] = out[i][0] * factor;
   }
   fftw_destroy_plan(plan);
   memory->destroy(in);
   memory->destroy(out);

   // in POSCAR, atoms are sorted/aggregated by type, while for LAMMPS there is no such requirment
   int type_id[nucell], num_type[nucell], ntype = 0;
   for (int i = 0; i < nucell; ++i) num_type[i] = 0;
   for (int i = 0; i < nucell; ++i){
      int ip = ntype;
      for (int j = 0; j < ntype; ++j){
         if (dm->attyp[i] == type_id[j]) ip = j;
      }
      if (ip == ntype) type_id[ntype++] = dm->attyp[i];
      num_type[ip]++;
   }
   std::map<int, int> iu_by_type;
   iu_by_type.clear();
   int id_new = 0;
   for (int i = 0; i < ntype; ++i){
      for (int j = 0; j < nucell; ++j){
         if (dm->attyp[j] == type_id[i]) iu_by_type[id_new++] = j;
      }
   }

   // write the FORCE_CONSTANTS file
   FILE *fp = fopen("FORCE_CONSTANTS", "w");
   int natom = npt * nucell;
   fprintf(fp, "%d %d\n", natom, natom);
   for (int i = 0; i < natom; ++i){
      int iu = i / npt;
      int iz = (i % npt) / (nx * ny);
      int iy = (i % (nx *ny) ) / nx;
      int ix = i % nx;
      iu = iu_by_type[iu];

      for (int j = 0; j < natom; ++j){
         int ju = j / npt;
         int jz = (j % npt) / (nx * ny);
         int jy = (j % (nx *ny) ) / nx;
         int jx = j % nx;

         int dx = abs(ix - jx);
         int dy = abs(iy - jy);
         int dz = abs(iz - jz);

         int id = (dx * ny + dy) *  nz + dz;
         ju = iu_by_type[ju];
         fprintf(fp, "%d %d\n", i+1, j+1);
         for (int idim = iu * sysdim; idim < (iu+1)*sysdim; ++idim){
            for (int jdim = ju * sysdim; jdim < (ju+1)*sysdim; ++jdim){
               int dd = idim * fftdim + jdim;
               fprintf(fp, " %lg", fc[id][dd]);
            }
            fprintf(fp, "\n");
         }
      }
   }
   fclose(fp);
   iu_by_type.clear();
   memory->destroy(fc);

   // write the primitive cell in POSCAR format
   fp = fopen("POSCAR.primitive", "w");
   fprintf(fp, "Fix-phonon unit cell");
   for (int ip = 0; ip < ntype; ++ip) fprintf(fp, ", Elem-%d: %lg", type_id[ip], mass[ip]);
   fprintf(fp, "\n1.\n"); 
   int ndim = 0;
   for (int idim = 0; idim < 3; ++idim){
      for (int jdim = 0; jdim < 3; ++jdim) fprintf(fp, "%lg ", dm->basevec[ndim++]);
      fprintf(fp, "\n");
   }
   for (int ip = 0; ip < ntype; ++ip) fprintf(fp, "Elem-%d ", type_id[ip]); fprintf(fp, "\n");
   for (int ip = 0; ip < ntype; ++ip) fprintf(fp, "%d ", num_type[ip]);
   fprintf(fp, "\nDirect\n");
   for (int ip = 0; ip < ntype; ++ip){
      for (int i = 0; i < nucell; ++i){
         if (dm->attyp[i] == type_id[ip]){
            fprintf(fp, "%lg %lg %lg\n", dm->basis[i][0], dm->basis[i][1], dm->basis[i][2]);
         }
      }
   }
   fclose(fp);

#ifdef UseSPG
  // Get high symmetry k-path
   QNodes *q = new QNodes();
   kPath *kp = new kPath(dm, q);
   kp->kpath();
#endif
   
   // band.conf
   fp = fopen("band.conf", "w");
   fprintf(fp, "# From Fix-phonon");
   for (int ip = 0; ip < ntype; ++ip) fprintf(fp, ", Elem-%d: %lg", type_id[ip], mass[ip]);
   fprintf(fp, "\n\nATOM_NAME = ");
   for (int ip = 0; ip < ntype; ++ip) fprintf(fp, "Elem-%d ", type_id[ip]);
   fprintf(fp, "\nDIM = %d %d %d\nBAND = ", nx, ny, nz);
#ifdef UseSPG
   int nsect = q->qs.size();
   for (int i = 0; i < nsect; ++i){
      fprintf(fp, " %lg %lg %lg", q->qs[i][0], q->qs[i][1], q->qs[i][2]);
      if (i+1 < nsect){
         double dq = 0.;
         for (int j = 0; j < 3; ++j) dq += (q->qe[i][j] - q->qs[i+1][j]) * (q->qe[i][j] - q->qs[i+1][j]);
         if (dq > ZERO) {
            fprintf(fp, " %lg %lg %lg,", q->qe[i][0], q->qe[i][1], q->qe[i][2]);
         }
      } else if (i+1 == nsect){
         fprintf(fp, " %lg %lg %lg\n", q->qe[i][0], q->qe[i][1], q->qe[i][2]);
      }
   }
#endif
   fprintf(fp, "\nBAND_POINTS = 21\nBAND_LABELS =");
#ifdef UseSPG
   for (int i = 0; i < q->ndstr.size(); ++i){
      std::size_t found = q->ndstr[i].find("{/Symbol G}");
      if (found != std::string::npos) q->ndstr[i].replace(found, found+11, "$\\Gamma$");
      found = q->ndstr[i].find("/");
      if (found != std::string::npos) q->ndstr[i].replace(found, found, " ");
      fprintf(fp, " %s", q->ndstr[i].c_str());
   }
#endif
   fprintf(fp, "\nFORCE_CONSTANTS = READ\nBAND_CONNECTION = .TRUE.\n");

   // output info
   write(4);
#ifdef UseSPG
   kp->show_path();
   delete kp;
   delete q;
#endif
   write(5);

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int Phonopy::count_words(const char *line)
{
   int n = strlen(line) + 1;
   char *copy;
   memory->create(copy, n, "count_words:copy");
   strcpy(copy,line);
  
   char *ptr;
   if (ptr = strchr(copy,'#')) *ptr = '\0';
  
   if (strtok(copy," \t\n\r\f") == NULL) {
     memory->destroy(copy);
     return 0;
   }
   n = 1;
   while (strtok(NULL," \t\n\r\f")) n++;
  
   memory->destroy(copy);
return n;
}
/*----------------------------------------------------------------------------*/
#endif
