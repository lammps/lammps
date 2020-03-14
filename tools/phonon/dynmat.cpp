#include "dynmat.h"
#include "math.h"
#include "version.h"
#include "global.h"

extern "C" void zheevd_(char *, char *, long int *, doublecomplex *,
                       long int *, double *, doublecomplex *,
                       long int *, double *, long int *, long int *,
                       long int *, long int *);

// to initialize the class
DynMat::DynMat(int narg, char **arg)
{
  attyp = NULL;
  memory = NULL;
  M_inv_sqrt = NULL;
  interpolate = NULL;
  DM_q = DM_all = NULL;
  binfile = funit = dmfile = NULL;

  attyp = NULL;
  basis = NULL;
  flag_reset_gamma = flag_skip = 0;

  // analyze the command line options
  int iarg = 1;
  while (narg > iarg){
    if (strcmp(arg[iarg], "-s") == 0){
      flag_reset_gamma = flag_skip = 1;

    } else if (strcmp(arg[iarg], "-r") == 0){
      flag_reset_gamma = 1;

    } else if (strcmp(arg[iarg], "-h") == 0){
      help();

    } else {
      if (binfile) delete []binfile;
      int n = strlen(arg[iarg]) + 1;
      binfile = new char[n];
      strcpy(binfile, arg[iarg]); 
    }

    iarg++;
  }

  ShowVersion();
  // get the binary file name from user input if not found in command line
  char str[MAXLINE];
  if (binfile == NULL) {
    char *ptr = NULL;
    printf("\n");
    do {
      printf("Please input the binary file name from fix_phonon: ");
      fgets(str,MAXLINE,stdin);
      ptr = strtok(str, " \n\t\r\f");
    } while (ptr == NULL);

    int n = strlen(ptr) + 1;
    binfile = new char[n];
    strcpy(binfile, ptr);
  }

  // open the binary file
  FILE *fp = fopen(binfile, "rb");
  if (fp == NULL) {
    printf("\nFile %s not found! Programe terminated.\n", binfile);
    help();
  }

  // read header info from the binary file
  if ( fread(&sysdim, sizeof(int),    1, fp) != 1) {printf("\nError while reading sysdim from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nx,     sizeof(int),    1, fp) != 1) {printf("\nError while reading nx from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&ny,     sizeof(int),    1, fp) != 1) {printf("\nError while reading ny from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nz,     sizeof(int),    1, fp) != 1) {printf("\nError while reading nz from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nucell, sizeof(int),    1, fp) != 1) {printf("\nError while reading nucell from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&boltz,  sizeof(double), 1, fp) != 1) {printf("\nError while reading boltz from file: %s\n", binfile); fclose(fp); exit(2);}

  fftdim = sysdim*nucell; fftdim2 = fftdim*fftdim;
  npt = nx*ny*nz;

  // display info related to the read file
  printf("\n"); for (int i = 0; i < 80; ++i) printf("="); printf("\n");
  printf("Dynamical matrix is read from file: %s\n", binfile);
  printf("The system size in three dimension: %d x %d x %d\n", nx, ny, nz);
  printf("Number of atoms per unit cell     : %d\n", nucell);
  printf("System dimension                  : %d\n", sysdim);
  printf("Boltzmann constant in used units  : %g\n", boltz);
  for (int i = 0; i < 80; ++i) printf("=");
  printf("\n");
  if (sysdim < 1||sysdim > 3||nx < 1||ny < 1||nz < 1||nucell < 1){
    printf("Wrong values read from header of file: %s, please check the binary file!\n", binfile);
    fclose(fp); exit(3);
  }

  funit = new char[4];
  strcpy(funit, "THz");
  if (boltz == 1.){eml2f = 1.; delete funit; funit = new char[27]; strcpy(funit,"sqrt(epsilon/(m.sigma^2))");}
  else if (boltz == 0.0019872067)  eml2f = 3.256576161;
  else if (boltz == 8.617343e-5)   eml2f = 15.63312493;
  else if (boltz == 1.3806504e-23) eml2f = 1.;
  else if (boltz == 1.3806504e-16) eml2f = 1.591549431e-14;
  else {
    printf("WARNING: Because of float precision, I cannot get the factor to convert sqrt(E/ML^2)\n");
    printf("into THz, instead, I set it to be 1; you should check the unit used by LAMMPS.\n");
    eml2f = 1.;
  }

  // now to allocate memory for DM
  memory = new Memory();
  memory->create(DM_all, npt, fftdim2, "DynMat:DM_all");
  memory->create(DM_q, fftdim,fftdim,"DynMat:DM_q");

  // read all dynamical matrix info into DM_all
  if ( fread(DM_all[0], sizeof(doublecomplex), npt*fftdim2, fp) != size_t(npt*fftdim2)){
    printf("\nError while reading the DM from file: %s\n", binfile);
    fclose(fp);
    exit(1);
  }

  // now try to read unit cell info from the binary file
  memory->create(basis, nucell, sysdim, "DynMat:basis");
  memory->create(attyp, nucell,         "DynMat:attyp");
  memory->create(M_inv_sqrt, nucell,    "DynMat:M_inv_sqrt");
  
  if ( (int) fread(&Tmeasure,      sizeof(double), 1,      fp) != 1     ){printf("\nError while reading temperature from file: %s\n",   binfile); fclose(fp); exit(3);}
  if ( (int) fread(&basevec[0],    sizeof(double), 9,      fp) != 9     ){printf("\nError while reading lattice info from file: %s\n",  binfile); fclose(fp); exit(3);}
  if ( (int) fread(basis[0],       sizeof(double), fftdim, fp) != fftdim){printf("\nError while reading basis info from file: %s\n",    binfile); fclose(fp); exit(3);}
  if ( (int) fread(&attyp[0],      sizeof(int),    nucell, fp) != nucell){printf("\nError while reading atom types from file: %s\n",    binfile); fclose(fp); exit(3);}
  if ( (int) fread(&M_inv_sqrt[0], sizeof(double), nucell, fp) != nucell){printf("\nError while reading atomic masses from file: %s\n", binfile); fclose(fp); exit(3);}
  fclose(fp);

  car2dir();
  real2rec();

  // initialize interpolation
  interpolate = new Interpolate(nx,ny,nz,fftdim2,DM_all);
  if (flag_reset_gamma) interpolate->reset_gamma();

  // Enforcing Austic Sum Rule
  EnforceASR();

  // get the dynamical matrix from force constant matrix: D = 1/M x Phi
  for (int idq = 0; idq < npt; ++idq){
    int ndim =0;
    for (int idim = 0; idim < fftdim; ++idim)
    for (int jdim = 0; jdim < fftdim; ++jdim){
      double inv_mass = M_inv_sqrt[idim/sysdim]*M_inv_sqrt[jdim/sysdim];
      DM_all[idq][ndim].r *= inv_mass;
      DM_all[idq][ndim].i *= inv_mass;
      ndim++;
    }
  }

  // ask for the interpolation method
  interpolate->set_method();

  return;
}

// to destroy the class
DynMat::~DynMat()
{
 // destroy all memory allocated
 if (funit) delete []funit;
 if (dmfile) delete []dmfile;
 if (binfile) delete []binfile;
 if (interpolate) delete interpolate;

 memory->destroy(DM_q);
 memory->destroy(attyp);
 memory->destroy(basis);
 memory->destroy(DM_all);
 memory->destroy(M_inv_sqrt);
 if (memory) delete memory;
}

/* ----------------------------------------------------------------------------
 * method to write DM_q to file, single point
 * ---------------------------------------------------------------------------- */
void DynMat::writeDMq(double *q)
{
  FILE *fp;
  // only ask for file name for the first time
  // other calls will append the result to the file.
  if (dmfile == NULL){
    char str[MAXLINE], *ptr;
    printf("\n");
    while ( 1 ){
      printf("Please input the filename to output the DM at selected q: ");
      fgets(str,MAXLINE,stdin);
      ptr = strtok(str, " \r\t\n\f");
      if (ptr) break;
    }

    int n = strlen(ptr) + 1;
    dmfile = new char[n];
    strcpy(dmfile, ptr);
    fp = fopen(dmfile,"w");

  } else {
    fp = fopen(dmfile,"a");
  }
  fprintf(fp,"# q = [%lg %lg %lg]\n", q[0], q[1], q[2]);

  for (int i = 0; i < fftdim; ++i){
    for (int j = 0; j < fftdim; ++j) fprintf(fp,"%lg %lg\t", DM_q[i][j].r, DM_q[i][j].i);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  fclose(fp);
return;
}

/* ----------------------------------------------------------------------------
 * method to write DM_q to file, dispersion-like
 * ---------------------------------------------------------------------------- */
void DynMat::writeDMq(double *q, const double qr, FILE *fp)
{

  fprintf(fp, "%lg %lg %lg %lg ", q[0], q[1], q[2], qr);

  for (int i = 0; i < fftdim; ++i)
  for (int j = 0; j < fftdim; ++j) fprintf(fp,"%lg %lg\t", DM_q[i][j].r, DM_q[i][j].i);

  fprintf(fp,"\n");
return;
}

/* ----------------------------------------------------------------------------
 * method to evaluate the eigenvalues of current q-point;
 * return the eigenvalues in egv.
 * cLapack subroutine zheevd is employed.
 * ---------------------------------------------------------------------------- */
int DynMat::geteigen(double *egv, int flag)
{
  char jobz, uplo;
  long int n, lda, lwork, lrwork, *iwork, liwork, info;
  doublecomplex *work;
  double *w = &egv[0], *rwork;

  n     = fftdim;
  if (flag) jobz = 'V';
  else jobz = 'N';

  uplo = 'U';
  lwork = (n+2)*n;
  lrwork = 1 + (5+n+n)*n;
  liwork = 3 + 5*n;
  lda    = n;

  memory->create(work,  lwork,  "geteigen:work");
  memory->create(rwork, lrwork, "geteigen:rwork");
  memory->create(iwork, liwork, "geteigen:iwork");

  zheevd_(&jobz, &uplo, &n, DM_q[0], &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
 
  // to get w instead of w^2; and convert w into v (THz hopefully)
  for (int i = 0; i < n; ++i){
    if (w[i]>= 0.) w[i] = sqrt(w[i]);
    else w[i] = -sqrt(-w[i]);

    w[i] *= eml2f;
  }

  memory->destroy(work);
  memory->destroy(rwork);
  memory->destroy(iwork);

return info;
}

/* ----------------------------------------------------------------------------
 * method to get the Dynamical Matrix at q
 * ---------------------------------------------------------------------------- */
void DynMat::getDMq(double *q)
{
  interpolate->execute(q, DM_q[0]);
return;
}

/* ----------------------------------------------------------------------------
 * method to get the Dynamical Matrix at q
 * ---------------------------------------------------------------------------- */
void DynMat::getDMq(double *q, double *wt)
{
  interpolate->execute(q, DM_q[0]);

  if (flag_skip && interpolate->UseGamma ) wt[0] = 0.;
return;
}

/* ----------------------------------------------------------------------------
 * private method to convert the cartisan coordinate of basis into fractional
 * ---------------------------------------------------------------------------- */
void DynMat::car2dir()
{
  double mat[9];
  for (int idim = 0; idim < 9; ++idim) mat[idim] = basevec[idim];
  GaussJordan(3, mat);

  for (int i = 0; i < nucell; ++i){
    double x[3];
    x[0] = x[1] = x[2] = 0.;
    for (int idim = 0; idim < sysdim; idim++) x[idim] = basis[i][idim];
    for (int idim = 0; idim < sysdim; idim++)
      basis[i][idim] = x[0]*mat[idim] + x[1]*mat[3+idim] + x[2]*mat[6+idim];
  }

return;
}

/* ----------------------------------------------------------------------------
 * private method to enforce the acoustic sum rule on force constant matrix at G
 * ---------------------------------------------------------------------------- */
void DynMat::EnforceASR()
{
  char str[MAXLINE];
  int nasr = 20;
  if (nucell <= 1) nasr = 1;
  printf("\n"); for (int i = 0; i < 80; ++i) printf("=");

  // compute and display eigenvalues of Phi at gamma before ASR
  if (nucell > 100){
    printf("\nYour unit cell is rather large, eigenvalue evaluation takes some time...");
    fflush(stdout);
  }

  double egvs[fftdim];
  for (int i = 0; i < fftdim; ++i)
  for (int j = 0; j < fftdim; ++j) DM_q[i][j] = DM_all[0][i*fftdim+j];
  geteigen(egvs, 0);
  printf("\nEigenvalues of Phi at gamma before enforcing ASR:\n");
  for (int i = 0; i < fftdim; ++i){
    printf("%lg ", egvs[i]);
    if (i%10 == 9) printf("\n");
    if (i == 99){ printf("...... (%d more skipped)\n", fftdim-100); break;}
  }
  printf("\n\n");

  // ask for iterations to enforce ASR
  printf("Please input the # of iterations to enforce ASR [%d]: ", nasr);
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str," \t\n\r\f");
  if (ptr) nasr = atoi(ptr);
  if (nasr < 1){
    for (int i=0; i<80; i++) printf("=");
    printf("\n");
    return;
  }

  for (int iit = 0; iit < nasr; ++iit){
    // simple ASR; the resultant matrix might not be symmetric
    for (int a = 0; a < sysdim; ++a)
    for (int b = 0; b < sysdim; ++b){
      for (int k = 0; k < nucell; ++k){
        double sum = 0.;
        for (int kp = 0; kp < nucell; ++kp){
          int idx = (k*sysdim+a)*fftdim+kp*sysdim+b;
          sum += DM_all[0][idx].r;
        }
        sum /= double(nucell);
        for (int kp = 0; kp < nucell; ++kp){
          int idx = (k*sysdim+a)*fftdim+kp*sysdim+b;
          DM_all[0][idx].r -= sum;
        }
      }
    }
   
    // symmetrize
    for (int k  = 0; k  < nucell; ++k)
    for (int kp = k; kp < nucell; ++kp){
      double csum = 0.;
      for (int a = 0; a < sysdim; ++a)
      for (int b = 0; b < sysdim; ++b){
        int idx = (k*sysdim+a)*fftdim+kp*sysdim+b;
        int jdx = (kp*sysdim+b)*fftdim+k*sysdim+a;
        csum = (DM_all[0][idx].r + DM_all[0][jdx].r )*0.5;
        DM_all[0][idx].r = DM_all[0][jdx].r = csum;
      }
    }
  }

  // symmetric ASR
  for (int a = 0; a < sysdim; ++a)
  for (int b = 0; b < sysdim; ++b){
    for (int k = 0; k < nucell; ++k){
      double sum = 0.;
      for (int kp = 0; kp < nucell; ++kp){
        int idx = (k*sysdim+a)*fftdim+kp*sysdim+b;
        sum += DM_all[0][idx].r;
      }
      sum /= double(nucell-k);
      for (int kp = k; kp < nucell; ++kp){
        int idx = (k*sysdim+a)*fftdim+kp*sysdim+b;
        int jdx = (kp*sysdim+b)*fftdim+k*sysdim+a;
        DM_all[0][idx].r -= sum;
        DM_all[0][jdx].r  = DM_all[0][idx].r;
      }
    }
  }

  // compute and display eigenvalues of Phi at gamma after ASR
  for (int i = 0; i < fftdim; ++i)
  for (int j = 0; j < fftdim; ++j) DM_q[i][j] = DM_all[0][i*fftdim+j];
  geteigen(egvs, 0);
  printf("Eigenvalues of Phi at gamma after enforcing ASR:\n");
  for (int i = 0; i < fftdim; ++i){
    printf("%lg ", egvs[i]);
    if (i%10 == 9) printf("\n");
    if (i == 99){ printf("...... (%d more skipped)", fftdim-100); break;}
  }
  printf("\n");
  for (int i = 0; i < 80; ++i) printf("=");
  printf("\n\n");

return;
}

/* ----------------------------------------------------------------------------
 * private method to get the reciprocal lattice vectors from the real space ones
 * ---------------------------------------------------------------------------- */
void DynMat::real2rec()
{
  ibasevec[0] = basevec[4]*basevec[8] - basevec[5]*basevec[7];
  ibasevec[1] = basevec[5]*basevec[6] - basevec[3]*basevec[8];
  ibasevec[2] = basevec[3]*basevec[7] - basevec[4]*basevec[6];

  ibasevec[3] = basevec[7]*basevec[2] - basevec[8]*basevec[1];
  ibasevec[4] = basevec[8]*basevec[0] - basevec[6]*basevec[2];
  ibasevec[5] = basevec[6]*basevec[1] - basevec[7]*basevec[0];

  ibasevec[6] = basevec[1]*basevec[5] - basevec[2]*basevec[4];
  ibasevec[7] = basevec[2]*basevec[3] - basevec[0]*basevec[5];
  ibasevec[8] = basevec[0]*basevec[4] - basevec[1]*basevec[3];

  double vol = 0.;
  for (int i = 0; i < sysdim; ++i) vol += ibasevec[i] * basevec[i];
  vol = 8.*atan(1.)/vol;

  for (int i = 0; i < 9; ++i) ibasevec[i] *= vol;

  printf("\n"); for (int i = 0; i < 80; ++i) printf("=");
  printf("\nBasis vectors of the unit cell in real space:");
  for (int i = 0; i < sysdim; ++i){
    printf("\n     A%d: ", i+1);
    for (int j = 0; j < sysdim; ++j) printf("%8.4f ", basevec[i*3+j]);
  }
  printf("\nBasis vectors of the corresponding reciprocal cell:");
  for (int i = 0; i < sysdim; ++i){
    printf("\n     B%d: ", i+1);
    for (int j = 0; j < sysdim; ++j) printf("%8.4f ", ibasevec[i*3+j]);
  }
  printf("\n"); for (int i = 0; i < 80; ++i) printf("="); printf("\n");

return;
}

/* ----------------------------------------------------------------------
 * private method, to get the inverse of a double matrix by means of
 * Gaussian-Jordan Elimination with full pivoting; square matrix required.
 *
 * Adapted from the Numerical Recipes in Fortran.
 * --------------------------------------------------------------------*/
void DynMat::GaussJordan(int n, double *Mat)
{
  int i,icol=0,irow=0,j,k,l,ll,idr,idc;
  int *indxc,*indxr,*ipiv;
  double big, nmjk;
  double dum, pivinv;

  indxc = new int[n];
  indxr = new int[n];
  ipiv  = new int[n];

  for (i = 0; i < n; ++i) ipiv[i] = 0;
  for (i = 0; i < n; ++i){
    big = 0.;
    for (j = 0; j < n; ++j){
      if (ipiv[j] != 1){
        for (k = 0; k < n; ++k){
          if (ipiv[k] == 0){
            idr = j * n + k;
            nmjk = abs(Mat[idr]);
            if (nmjk >= big){
              big  = nmjk;
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1){
            printf("DynMat: Singular matrix in double GaussJordan!\n"); exit(1);
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l = 0; l < n; ++l){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol * n + icol;
    if (Mat[idr] == 0.){
      printf("DynMat: Singular matrix in double GaussJordan!");
      exit(1);
    }
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.;
    idr = icol*n;
    for (l = 0; l < n; ++l) Mat[idr+l] *= pivinv;
    for (ll = 0; ll < n; ++ll){
      if (ll != icol){
        idc = ll * n + icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l = 0; l < n; ++l) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l = n-1; l >= 0; --l){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k = 0; k < n; ++k){
        idr = k * n + rl;
        idc = k * n + cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }
  delete []indxr;
  delete []indxc;
  delete []ipiv;

return;
}

/* ----------------------------------------------------------------------------
 * Public method to reset the interpolation method
 * ---------------------------------------------------------------------------- */
void DynMat::reset_interp_method()
{
  interpolate->set_method();

return;
}

/* ----------------------------------------------------------------------------
 * Private method to display help info
 * ---------------------------------------------------------------------------- */
void DynMat::help()
{
  ShowVersion();
  printf("\nUsage:\n  phana [options] [file]\n\n");
  printf("Available options:\n");
  printf("  -r          To reset the dynamical matrix at the gamma point by a 4th order\n");
  printf("              polynomial interpolation along the [100] direction; this might be\n");
  printf("              useful for systems with charges. As for charged system, the dynamical\n");
  printf("              matrix at Gamma is far from accurate because of the long range nature\n");
  printf("              of Coulombic interaction. By reset it by interpolation, will partially\n");
  printf("              elliminate the unwanted behavior, but the result is still inaccurate.\n");
  printf("              By default, this is not set; and not expected for uncharged systems.\n\n");
  printf("  -s          This will reset the dynamical matrix at the gamma point, too, but it\n");
  printf("              will also inform the code to skip all q-points that is in the vicinity\n");
  printf("              of the gamma point when evaluating phonon DOS and/or phonon dispersion.\n\n");
  printf("              By default, this is not set; and not expected for uncharged systems.\n\n");
  printf("  -h          To print out this help info.\n\n");
  printf("  file        To define the filename that carries the binary dynamical matrice generated\n");
  printf("              by fix-phonon. If not provided, the code will ask for it.\n");
  printf("\n\n");
  exit(0);
}

/* ----------------------------------------------------------------------------
 * Private method to display the version info
 * ---------------------------------------------------------------------------- */
void DynMat::ShowVersion()
{
  printf("                ____  _   _    __    _  _    __   \n");
  printf("               (  _ \\( )_( )  /__\\  ( \\( )  /__\\  \n");
  printf("                )___/ ) _ (  /(__)\\  )  (  /(__)\\ \n");
  printf("               (__)  (_) (_)(__)(__)(_)\\_)(__)(__)\n");
  printf("\nPHonon ANAlyzer for Fix-Phonon, version 2.%02d, compiled on %s.\n", VERSION, __DATE__);

return;
}
/* --------------------------------------------------------------------*/
