/***************************************************************************                               
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// POD header file 
#include "pod.h"

// LAMMPS header files 
#include "comm.h"
#include "error.h"
#include "memory.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <glob.h>
#include <random>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using namespace LAMMPS_NS;

// constructor 
CPOD::CPOD(LAMMPS* lmp, std::string pod_file, std::string coeff_file) : Pointers(lmp)
{
        // read pod input file to podstruct
    this->read_pod(pod_file);    
        
    // read pod coefficient file to podstruct
    if (coeff_file != "")
        this->read_coeff_file(coeff_file);            
    
    if (pod.snaptwojmax > 0) {
        InitSnap();
    }    
}

// destructor        
CPOD::~CPOD()
{
    pod.freememory(1);
    sna.freememory(1);
}

/****************************************************************************************************************/
void CPOD::print_matrix(const char* desc, int m, int n, double **a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[j][i] );
        printf( "\n" );
    }
}

void CPOD::print_matrix(const char* desc, int m, int n, double* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
        printf( "\n" );
    }
}

void CPOD::print_matrix(const char* desc, int m, int n, int* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
        printf( "\n" );
    }
}

void CPOD::podMatMul(double *c, double *a, double *b, int r1, int c1, int c2) 
{
    int i, j, k;

    for(j = 0; j < c2; j++)        
        for(i = 0; i < r1; i++)
            c[i + r1*j] = 0.0;        
    
    for(j = 0; j < c2; j++)
        for(i = 0; i < r1; i++)        
            for(k = 0; k < c1; k++)            
                c[i + r1*j] += a[i + r1*k] * b[k + c1*j];            
}

void CPOD::podArrayFill(int* output, int start, int length) 
{	
	for (int j = 0; j < length; ++j)	
		output[j] = start + j;
}

double CPOD::podArraySum(double *a, int n)
{
    double e = a[0];
    for (int i=1; i<n; i++)        
        e += a[i];    
    return e;
}

double CPOD::podArrayMin(double *a, int n)
{
    double b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

double CPOD::podArrayMax(double *a, int n)
{
    double b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

int CPOD::podArrayMin(int *a, int n)
{
    int b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

int CPOD::podArrayMax(int *a, int n)
{
    int b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

void CPOD::podKron(double *C, double *A, double *B, double alpha, int M1, int M2)
{            
    int M = M1*M2;        
    for (int idx=0; idx<M; idx++)     
    {
        int ib = idx%M2;
        int ia = (idx-ib)/M2;        
        C[idx] += alpha*A[ia]*B[ib];        
    }
}

void CPOD::podCumsum(int* output, int* input, int length) 
{
	output[0] = 0; 
	for (int j = 1; j < length; ++j)	
		output[j] = input[j - 1] + output[j - 1];	
}

double CPOD::podArrayNorm(double *a, int n)
{
    double e = a[0]*a[0];
    for (int i=1; i<n; i++)        
        e += a[i]*a[i];    
    return sqrt(e);
}

double CPOD::podArrayErrorNorm(double *a, double *b, int n)
{
    double e = (a[0]-b[0])*(a[0]-b[0]);
    for (int i=1; i<n; i++)        
        e += (a[i]-b[i])*(a[i]-b[i]);    
    return sqrt(e);
}

void CPOD::podArraySetValue(double *y, double a, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

void CPOD::podArrayCopy(double *y, double *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = x[i];        
}

void CPOD::rotation_matrix(double *Rmat, double alpha, double beta, double gamma)
{    
    double ca = cos(alpha);
    double cb = cos(beta);
    double cg = cos(gamma);
    double sa = sin(alpha);
    double sb = sin(beta);
    double sg = sin(gamma);

    Rmat[0] = ca*cg*cb - sa*sg;
    Rmat[3] = -ca*cb*sg - sa*cg;
    Rmat[6] = ca*sb;
    
    Rmat[1] = sa*cg*cb + ca*sg;
    Rmat[4] = -sa*cb*sg + ca*cg;
    Rmat[7] = sa*sb;
    
    Rmat[2] = -sb*cg;
    Rmat[5] = sb*sg;
    Rmat[8] = cb;    
}

void CPOD::matrix33_multiplication(double *xrot, double *Rmat, double *x, int natom)
{
    double x1, x2, x3;
    for (int i=0; i < natom; i++) {
        x1 = x[0 + 3*i];
        x2 = x[1 + 3*i];
        x3 = x[2 + 3*i];
        xrot[0 + 3*i] = Rmat[0]*x1 + Rmat[3]*x2 + Rmat[6]*x3;
        xrot[1 + 3*i] = Rmat[1]*x1 + Rmat[4]*x2 + Rmat[7]*x3;
        xrot[2 + 3*i] = Rmat[2]*x1 + Rmat[5]*x2 + Rmat[8]*x3;
    }
}

void CPOD::matrix33_inverse(double *invA, double *A1, double *A2, double *A3)
{                 
    double a11 = A1[0];
    double a21 = A1[1];
    double a31 = A1[2];
    double a12 = A2[0];
    double a22 = A2[1];
    double a32 = A2[2];
    double a13 = A3[0];
    double a23 = A3[1];
    double a33 = A3[2];        
    double detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);

    invA[0] = (a22*a33 - a23*a32)/detA;
    invA[1] = (a23*a31 - a21*a33)/detA;
    invA[2] = (a21*a32 - a22*a31)/detA;
    invA[3] = (a13*a32 - a12*a33)/detA;
    invA[4] = (a11*a33 - a13*a31)/detA;
    invA[5] = (a12*a31 - a11*a32)/detA;
    invA[6] = (a12*a23 - a13*a22)/detA;
    invA[7] = (a13*a21 - a11*a23)/detA;
    invA[8] = (a11*a22 - a12*a21)/detA;            
}

void CPOD::triclinic_lattice_conversion(double *a, double *b, double *c, double *A, double *B, double *C)
{
    double Anorm = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    double Bnorm = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    double Cnorm = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);
    
    double Ahat[3];
    Ahat[0] = A[0]/Anorm; Ahat[1] = A[1]/Anorm; Ahat[2] = A[2]/Anorm;
 
    double ax = Anorm;
    double bx = B[0]*Ahat[0] + B[1]*Ahat[1] + B[2]*Ahat[2]; //dot(B,Ahat);
    double by = sqrt(Bnorm*Bnorm - bx*bx); //sqrt(Bnorm^2 - bx^2);// #norm(cross(Ahat,B));
    double cx = C[0]*Ahat[0] + C[1]*Ahat[1] + C[2]*Ahat[2]; // dot(C,Ahat);
    double cy = (B[0]*C[0] + B[1]*C[1] + B[2]*C[2] - bx*cx)/by; // (dot(B, C) - bx*cx)/by;
    double cz = sqrt(Cnorm*Cnorm - cx*cx - cy*cy); // sqrt(Cnorm^2 - cx^2 - cy^2);    
    
    a[0] = ax; a[1] = 0.0; a[2] = 0.0;
    b[0] = bx; b[1] = by;  b[2] = 0.0;
    c[0] = cx; c[1] = cy;  c[2] = cz;
}

/**************************************************************************************************************/

void podsnapshots(double *rbf, double *xij, double *besselparams, double rin, double rcut, 
        int besseldegree, int inversedegree, int nbesselpars, int N)
{
    double rmax = rcut-rin;
    for (int n=0; n<N; n++) {            
        double dij = xij[n];    

        double r = dij - rin;        
        double y = r/rmax;            
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double fcut = y6/exp(-1.0);
        
        for (int j=0; j<nbesselpars; j++) {            
            double alpha = besselparams[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            
            for (int i=0; i<besseldegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*besseldegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
            }
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            int nij = n + N*p;     
            double a = pow(dij, (double) (i+1.0));
            rbf[nij] = fcut/a;
        }        
    }
}

void podeigenvaluedecomposition(double *Phi, double *Lambda, double *besselparams, double rin, double rcut, 
        int besseldegree, int inversedegree, int nbesselpars, int N)
{
    int ns = besseldegree*nbesselpars + inversedegree;
    
    double *xij = (double *) malloc(N*sizeof(double));
    double *S = (double *) malloc(N*ns*sizeof(double));
    double *Q = (double *) malloc(N*ns*sizeof(double));
    double *A = (double *) malloc(ns*ns*sizeof(double));
    double *b = (double *) malloc(ns*sizeof(double));
    
    for (int i=0; i<N; i++)
        xij[i] = (rin+1e-6) + (rcut-rin-1e-6)*(i*1.0/(N-1));        
    
    podsnapshots(S, xij, besselparams, rin, rcut, besseldegree, inversedegree, nbesselpars, N);
        
    char chn = 'N';
    char cht = 'T';
    char chv = 'V';
    char chu = 'U';
    double alpha = 1.0, beta = 0.0;    
    DGEMM(&cht, &chn, &ns, &ns, &N, &alpha, S, &N, S, &N, &beta, A, &ns);    
        
    for (int i=0; i<ns*ns; i++)
        A[i] = A[i]*(1.0/N);        
    
    // Declaring Function Input for DSYEV
//     char jobz = 'V';    // 'V':  Compute eigenvalues and eigenvectors.
//     char uplo = 'U';    // 'U':  Upper triangle of A is stored.
    int lwork = ns * ns;  // The length of the array work, lwork >= max(1,3*N-1).
    int info = 1;       // = 0:  successful exit
    double work[ns*ns];         
    DSYEV(&chv, &chu, &ns, A, &ns, b, work, &lwork, &info);
    
    // order eigenvalues and eigenvectors from largest to smallest
    for (int j=0; j<ns; j++)
        for (int i=0; i<ns; i++)        
            Phi[i + ns*(ns-j-1)] = A[i + ns*j];        

    for (int i=0; i<ns; i++)        
        Lambda[(ns-i-1)] = b[i];        
            
    DGEMM(&chn, &chn, &N, &ns, &ns, &alpha, S, &N, Phi, &ns, &beta, Q, &N);        
    for (int i=0; i<(N-1); i++)
        xij[i] = xij[i+1] - xij[i];
    double area;
    for (int m=0; m<ns; m++) {
        area = 0.0;
        for (int i=0; i<(N-1); i++)            
            area += 0.5*xij[i]*(Q[i + N*m]*Q[i + N*m] + Q[i+1 + N*m]*Q[i+1 + N*m]);
        for (int i=0; i<ns; i++)            
            Phi[i + ns*m] = Phi[i + ns*m]/sqrt(area);
    }            
        
    free(xij); free(S); free(A); free(b); free(Q);
}


void CPOD::read_pod(std::string pod_file)
{
    pod.nbesselpars = 3;
    pod.besselparams = (double *) malloc(3*sizeof(double));
    pod.pbc = (int *) malloc(3*sizeof(int));
    
    std::ifstream file_in(pod_file);
    if (!file_in) {error->all(FLERR,"Error: POD input file is not found");}
        
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s;
            double d;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
            if (s == "species") {
                std::string element;
                while(ss_line >> element){                                        
                    pod.species.push_back(element);                    
                    pod.nelements += 1;
                }                        
            }
            
            if (s == "pbc") {
                int i, j = 0;                
                while(ss_line >> i) {                                        
                    pod.pbc[j++] = i;                    
                }                        
            }
            
            if ((s != "#") && (s != "species") && (s != "pbc")) {
                ss_line >> d;
                
                if (s == "rin") pod.rin = d;
                if (s == "rcut") pod.rcut = d;
                if (s == "bessel_scaling_parameter1") pod.besselparams[0] = d;
                if (s == "bessel_scaling_parameter2") pod.besselparams[1] = d;
                if (s == "bessel_scaling_parameter3") pod.besselparams[2] = d;                
                if (s == "bessel_polynomial_degree") pod.besseldegree = (int) d;
                if (s == "inverse_polynomial_degree") pod.inversedegree = (int) d;
                if (s == "onebody") pod.onebody = (int) d;
                if (s == "twobody_bessel_polynomial_degree") pod.twobody[0] = (int) d;
                if (s == "twobody_inverse_polynomial_degree") pod.twobody[1] = (int) d;
                if (s == "twobody_number_radial_basis_functions") pod.twobody[2] = (int) d;
                if (s == "threebody_bessel_polynomial_degree") pod.threebody[0] = (int) d;
                if (s == "threebody_inverse_polynomial_degree") pod.threebody[1] = (int) d;
                if (s == "threebody_number_radial_basis_functions") pod.threebody[2] = (int) d;
                if (s == "threebody_number_angular_basis_functions") pod.threebody[3] = (int) (d-1);
                if (s == "fourbody_bessel_polynomial_degree") pod.fourbody[0] = (int) d;
                if (s == "fourbody_inverse_polynomial_degree") pod.fourbody[1] = (int) d;
                if (s == "fourbody_number_radial_basis_functions") pod.fourbody[2] = (int) d;
                if (s == "fourbody_snap_twojmax") pod.snaptwojmax = (int) d;
                if (s == "fourbody_snap_chemflag") pod.snapchemflag = (int) d;
                if (s == "fourbody_snap_rfac0") pod.snaprfac0 = d;
                if (s == "fourbody_snap_neighbor_weight1") pod.snapelementweight[0] = d;
                if (s == "fourbody_snap_neighbor_weight2") pod.snapelementweight[1] = d;
                if (s == "fourbody_snap_neighbor_weight3") pod.snapelementweight[2] = d;
                if (s == "fourbody_snap_neighbor_weight4") pod.snapelementweight[3] = d;
                if (s == "fourbody_snap_neighbor_weight5") pod.snapelementweight[4] = d;
                //if (s == "fourbody_number_spherical_harmonic_basis_functions") pod.fourbody[3] = (int) d;
                if (s == "quadratic22_number_twobody_basis_functions") pod.quadratic22[0] = (int) d;
                if (s == "quadratic22_number_twobody_basis_functions") pod.quadratic22[1] = (int) d;
                if (s == "quadratic23_number_twobody_basis_functions") pod.quadratic23[0] = (int) d;
                if (s == "quadratic23_number_threebody_basis_functions") pod.quadratic23[1] = (int) d;
                if (s == "quadratic24_number_twobody_basis_functions") pod.quadratic24[0] = (int) d;
                if (s == "quadratic24_number_fourbody_basis_functions") pod.quadratic24[1] = (int) d;
                if (s == "quadratic33_number_threebody_basis_functions") pod.quadratic33[0] = (int) d;
                if (s == "quadratic33_number_threebody_basis_functions") pod.quadratic33[1] = (int) d;
                if (s == "quadratic34_number_threebody_basis_functions") pod.quadratic34[0] = (int) d;
                if (s == "quadratic34_number_fourbody_basis_functions") pod.quadratic34[1] = (int) d;
                if (s == "quadratic44_number_fourbody_basis_functions") pod.quadratic44[0] = (int) d;
                if (s == "quadratic44_number_fourbody_basis_functions") pod.quadratic44[1] = (int) d;     
                if (s == "cubic234_number_twobody_basis_functions") pod.cubic234[0] = (int) d;
                if (s == "cubic234_number_threebody_basis_functions") pod.cubic234[1] = (int) d;
                if (s == "cubic234_number_fourbody_basis_functions") pod.cubic234[2] = (int) d;
                if (s == "cubic333_number_threebody_basis_functions") pod.cubic333[0] = (int) d;
                if (s == "cubic444_number_fourbody_basis_functions") pod.cubic444[0] = (int) d;
            }
        }        
    }          
    file_in.close();
        
    pod.twobody[0] = pod.besseldegree;
    pod.twobody[1] = pod.inversedegree;
    pod.threebody[0] = pod.besseldegree;
    pod.threebody[1] = pod.inversedegree;
    
    // number of snapshots
    pod.ns2 = pod.nbesselpars*pod.twobody[0] + pod.twobody[1];
    pod.ns3 = pod.nbesselpars*pod.threebody[0] + pod.threebody[1];
    pod.ns4 = pod.nbesselpars*pod.fourbody[0] + pod.fourbody[1];
    
    for (int i = 0; i < pod.nbesselpars; i++)
        if (fabs(pod.besselparams[i]) < 1e-3) pod.besselparams[i] = 1e-3;
            
    // allocate memory for eigenvectors and eigenvalues
    pod.Phi2 = (double *) malloc(pod.ns2*pod.ns2*sizeof(double));
    pod.Lambda2 = (double *) malloc(pod.ns2*sizeof(double));
    pod.Phi3 = (double *) malloc(pod.ns3*pod.ns3*sizeof(double));
    pod.Lambda3 = (double *) malloc(pod.ns3*sizeof(double));
    pod.Phi4 = (double *) malloc(pod.ns4*pod.ns4*sizeof(double));
    pod.Lambda4 = (double *) malloc(pod.ns4*sizeof(double));    
    
    if (pod.ns2>0) {
        podeigenvaluedecomposition(pod.Phi2, pod.Lambda2, pod.besselparams, pod.rin, pod.rcut, 
            pod.twobody[0], pod.twobody[1], pod.nbesselpars, 2000);                
            
//         /* Print eigenvalues */
//         print_matrix( "Eigenvalues for two-body potential:", 1, pod.ns2, pod.Lambda2, 1 );
// 
//         /* Print eigenvectors */
//         print_matrix( "Eigenvectors for two-body potential:", pod.ns2, pod.ns2, pod.Phi2, pod.ns2);        
    }
    if (pod.ns3>0) {
        podeigenvaluedecomposition(pod.Phi3, pod.Lambda3, pod.besselparams, pod.rin, pod.rcut, 
            pod.threebody[0], pod.threebody[1], pod.nbesselpars, 2000);        
    }
    if (pod.ns4>0) {
        podeigenvaluedecomposition(pod.Phi4, pod.Lambda4, pod.besselparams, pod.rin, pod.rcut, 
            pod.fourbody[0], pod.fourbody[1], pod.nbesselpars, 2000);        
    }
    
    // number of chemical combinations
    pod.nc2 = pod.nelements*(pod.nelements+1)/2;
    pod.nc3 = pod.nelements*pod.nelements*(pod.nelements+1)/2;            
    pod.nc4 = pod.snapchemflag ? pod.nelements*pod.nelements*pod.nelements*pod.nelements : pod.nelements;
            
    // number of basis functions and descriptors for one-body potential
    if (pod.onebody==1) {
        pod.nbf1 = 1;
        pod.nd1 = pod.nelements;
    }
    else {
        pod.nbf1 = 0;
        pod.nd1 = 0;        
    }    
    
    // number of basis functions and descriptors for two-body potential
    pod.nbf2 = pod.twobody[2];
    pod.nd2 = pod.nbf2*pod.nc2;
    
    // number of basis functions and descriptors for three-body potential
    pod.nrbf3 = pod.threebody[2];
    pod.nabf3 = pod.threebody[3];
    pod.nbf3 = pod.nrbf3*(1 + pod.nabf3);
    pod.nd3 = pod.nbf3*pod.nc3;

    // number of basis functions and descriptors for four-body potential 
    int twojmax = pod.snaptwojmax;    
    int idxb_count = 0;    
    if (twojmax > 0) {
        for(int j1 = 0; j1 <= twojmax; j1++)
            for(int j2 = 0; j2 <= j1; j2++)
                for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
                    if (j >= j1) idxb_count++;
    }
    pod.nbf4 = idxb_count;
    pod.nd4 = pod.nbf4*pod.nc4;
        
    pod.quadratic22[0] = PODMIN(pod.quadratic22[0], pod.nbf2);
    pod.quadratic22[1] = PODMIN(pod.quadratic22[1], pod.nbf2);
    pod.quadratic23[0] = PODMIN(pod.quadratic23[0], pod.nbf2);
    pod.quadratic23[1] = PODMIN(pod.quadratic23[1], pod.nbf3);
    pod.quadratic24[0] = PODMIN(pod.quadratic24[0], pod.nbf2);
    pod.quadratic24[1] = PODMIN(pod.quadratic24[1], pod.nbf4);
    pod.quadratic33[0] = PODMIN(pod.quadratic33[0], pod.nbf3);
    pod.quadratic33[1] = PODMIN(pod.quadratic33[1], pod.nbf3);
    pod.quadratic34[0] = PODMIN(pod.quadratic34[0], pod.nbf3);
    pod.quadratic34[1] = PODMIN(pod.quadratic34[1], pod.nbf4);
    pod.quadratic44[0] = PODMIN(pod.quadratic44[0], pod.nbf4);
    pod.quadratic44[1] = PODMIN(pod.quadratic44[1], pod.nbf4);
    
    pod.cubic234[0] = PODMIN(pod.cubic234[0], pod.nbf2);
    pod.cubic234[1] = PODMIN(pod.cubic234[1], pod.nbf3);
    pod.cubic234[2] = PODMIN(pod.cubic234[2], pod.nbf4);    
    pod.cubic333[0] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic333[1] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic333[2] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic444[0] = PODMIN(pod.cubic444[0], pod.nbf4);
    pod.cubic444[1] = PODMIN(pod.cubic444[0], pod.nbf4);
    pod.cubic444[2] = PODMIN(pod.cubic444[0], pod.nbf4);
    
    // number of descriptors for quadratic POD potentials        
    pod.nd22 = pod.quadratic22[0]*pod.quadratic22[1]*pod.nc2*pod.nc2;
    pod.nd23 = pod.quadratic23[0]*pod.quadratic23[1]*pod.nc2*pod.nc3;
    pod.nd24 = pod.quadratic24[0]*pod.quadratic24[1]*pod.nc2*pod.nc4;    
    pod.nd33 = pod.quadratic33[0]*pod.quadratic33[1]*pod.nc3*pod.nc3;
    pod.nd34 = pod.quadratic34[0]*pod.quadratic34[1]*pod.nc3*pod.nc4;
    pod.nd44 = pod.quadratic44[0]*pod.quadratic44[1]*pod.nc4*pod.nc4;
    
    int nq;
    nq = pod.quadratic22[0]*pod.nc2; pod.nd22 = nq*(nq+1)/2;
    nq = pod.quadratic33[0]*pod.nc3; pod.nd33 = nq*(nq+1)/2;
    nq = pod.quadratic44[0]*pod.nc4; pod.nd44 = nq*(nq+1)/2;
    
    // number of descriptors for cubic POD potentials        
    pod.nd234 = pod.cubic234[0]*pod.cubic234[1]*pod.cubic234[2]*pod.nc2*pod.nc3*pod.nc4;
    nq = pod.cubic333[0]*pod.nc3; pod.nd333 = nq*(nq+1)*(nq+2)/6;    
    nq = pod.cubic444[0]*pod.nc4; pod.nd444 = nq*(nq+1)*(nq+2)/6;    
    
    // total number of descriptors for all POD potentials
    pod.nd = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4 + pod.nd22 + pod.nd23 + pod.nd24 + 
             pod.nd33 + pod.nd34 + pod.nd44 + pod.nd234 + pod.nd333 + pod.nd444;             
    pod.nd1234 = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4;
    
    int nelements = pod.nelements;
    pod.elemindex = (int*) malloc (sizeof (int)*(nelements*nelements));     
        
    int k = 1;
    for (int i=0; i < nelements; i++) 
        for (int j=i; j < nelements; j++) {
            pod.elemindex[i + nelements*j] = k;
            pod.elemindex[j + nelements*i] = k;
            k += 1;
        }          
    
    std::cout<<"**************** Begin of POD Potentials ****************"<<std::endl;
    std::cout<<"species: ";
    for (int i=0; i<pod.nelements; i++)
        std::cout<<pod.species[i]<<" ";
    std::cout<<std::endl;      
    std::cout<<"periodic boundary conditions:  "<<pod.pbc[0]<<"  "<<pod.pbc[1]<<"  "<<pod.pbc[2]<<std::endl;
    std::cout<<"inner cut-off radius: "<<pod.rin<<std::endl;
    std::cout<<"outer cut-off radius: "<<pod.rcut<<std::endl;
    std::cout<<"bessel parameters: "<<pod.besselparams[0]<<"  "<<pod.besselparams[1]<<"  "<<pod.besselparams[2]<<std::endl;
    std::cout<<"bessel polynomial degree: "<<pod.besseldegree<<std::endl;
    std::cout<<"inverse polynomial degree: "<<pod.inversedegree<<std::endl;
    std::cout<<"one-body potential: "<<pod.onebody<<std::endl;
    std::cout<<"two-body potential: "<<pod.twobody[0]<<"  "<<pod.twobody[1]<<"  "<<pod.twobody[2]<<std::endl;
    std::cout<<"three-body potential: "<<pod.threebody[0]<<"  "<<pod.threebody[1]<<"  "<<pod.threebody[2]<<"  "<<pod.threebody[3]+1<<std::endl;
    std::cout<<"four-body SNAP potential: "<<pod.snaptwojmax<<"  "<<pod.snapchemflag<<std::endl;
    std::cout<<"three-body quadratic22 potential: "<<pod.quadratic22[0]<<"  "<<pod.quadratic22[1]<<std::endl;
    std::cout<<"four-body quadratic23 potential: "<<pod.quadratic23[0]<<"  "<<pod.quadratic23[1]<<std::endl;
    std::cout<<"five-body quadratic24 potential: "<<pod.quadratic24[0]<<"  "<<pod.quadratic24[1]<<std::endl;
    std::cout<<"five-body quadratic33 potential: "<<pod.quadratic33[0]<<"  "<<pod.quadratic33[1]<<std::endl;
    std::cout<<"six-body quadratic34 potential: "<<pod.quadratic34[0]<<"  "<<pod.quadratic34[1]<<std::endl;
    std::cout<<"seven-body quadratic44 potential: "<<pod.quadratic44[0]<<"  "<<pod.quadratic44[1]<<std::endl;    
    std::cout<<"seven-body cubic234 potential: "<<pod.cubic234[0]<<"  "<<pod.cubic234[1]<<"  "<<pod.cubic234[2]<<std::endl;
    std::cout<<"seven-body cubic333 potential: "<<pod.cubic333[0]<<"  "<<pod.cubic333[1]<<"  "<<pod.cubic333[2]<<std::endl;    
    std::cout<<"ten-body cubic444 potential: "<<pod.cubic444[0]<<"  "<<pod.cubic444[1]<<"  "<<pod.cubic444[2]<<std::endl;    
    std::cout<<"number of snapshots for two-body potential: "<<pod.ns2<<std::endl;
    std::cout<<"number of snapshots for three-body potential: "<<pod.ns3<<std::endl;
    std::cout<<"number of snapshots for four-body potential: "<<pod.ns4<<std::endl;    
    std::cout<<"number of basis functions for one-body potential: "<<pod.nbf1<<std::endl;
    std::cout<<"number of basis functions for two-body potential: "<<pod.nbf2<<std::endl;
    std::cout<<"number of basis functions for three-body potential: "<<pod.nbf3<<std::endl;
    std::cout<<"number of basis functions for four-body potential: "<<pod.nbf4<<std::endl;
    std::cout<<"number of descriptors for one-body potential: "<<pod.nd1<<std::endl;
    std::cout<<"number of descriptors for two-body potential: "<<pod.nd2<<std::endl;
    std::cout<<"number of descriptors for three-body potential: "<<pod.nd3<<std::endl;
    std::cout<<"number of descriptors for four-body potential: "<<pod.nd4<<std::endl;
    std::cout<<"number of descriptors for three-body quadratic22 potential: "<<pod.nd22<<std::endl;
    std::cout<<"number of descriptors for four-body quadratic23 potential: "<<pod.nd23<<std::endl;
    std::cout<<"number of descriptors for five-body quadratic24 potential: "<<pod.nd24<<std::endl;
    std::cout<<"number of descriptors for five-body quadratic33 potential: "<<pod.nd33<<std::endl;
    std::cout<<"number of descriptors for six-body quadratic34 potential: "<<pod.nd34<<std::endl;
    std::cout<<"number of descriptors for seven-body quadratic44 potential: "<<pod.nd44<<std::endl;
    std::cout<<"number of descriptors for seven-body cubic234 potential: "<<pod.nd333<<std::endl;
    std::cout<<"number of descriptors for seven-body cubic333 potential: "<<pod.nd234<<std::endl;
    std::cout<<"number of descriptors for ten-body cubic444 potential: "<<pod.nd444<<std::endl;
    std::cout<<"total number of descriptors for all POD potentials: "<<pod.nd<<std::endl;    
    std::cout<<"**************** End of POD Potentials ****************"<<std::endl<<std::endl;
}

void CPOD::read_coeff_file(std::string coeff_file)
{
    std::ifstream file_in(coeff_file);
    if (!file_in) {error->all(FLERR,"Error: Coefficient input file is not found");}
    
    int ncoeff=0;
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s;
            int n;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
                        
            if (s == "POD_coefficients:") {
                ss_line >> n;           
                ncoeff = n;   
                break;
            }
        }
    }        
    
    pod.coeff = (double *) malloc(ncoeff*sizeof(double));
    
    int k = 0;
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            double d;            
            std::istringstream ss_line(line);                                    
            ss_line >> d;                                    
            pod.coeff[k] = d;
            k += 1;
        }
    }        
    
    file_in.close();
    
//     std::cout<<"**************** Begin of Coefficient File ****************"<<std::endl;    
//     std::cout<<"number of POD coefficients: "<<ncoeff<<std::endl;
//     print_matrix( "POD coefficient vector:", 1, ncoeff, pod.coeff, 1); 
//     std::cout<<"**************** End of Coefficient File ****************"<<std::endl<<std::endl;
}

/*********************************************************************************************************/

void CPOD::linear_descriptors(double *gd, double *efatom, double *y, double *tmpmem, int *atomtype, 
            int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij)
{
    int dim = 3;    
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nrbf2 = pod.nbf2;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;    
    int *pdegree2 = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi2 = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    double *fatom1 = &efatom[0];
    double *fatom2 = &efatom[dim*natom*(nd1)];
    double *fatom3 = &efatom[dim*natom*(nd1+nd2)];    
    double *fatom4 = &efatom[dim*natom*(nd1+nd2+nd3)];    
    double *eatom1 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)];
    double *eatom2 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*nd1];
    double *eatom3 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2)];
    double *eatom4 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2+nd3)];    
        
    podArraySetValue(fatom1, 0.0, (1+dim)*natom*(nd1+nd2+nd3+nd4));    
    
//     // peratom descriptors for one-body, two-body, and three-body linear potentials
//     this->poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, y, Phi2, besselparams, 
//             tmpmem, rin, rcut, atomtype, alist, pairlist, pairnum, pairnumsum, 
//             elemindex, pdegree2, tmpint, nbesselpars, nrbf2, nrbf3, nabf3, 
//             nelements, Nij, natom);                    
//     
//     if (pod.snaptwojmax>0) 
//         this->snapdesc(eatom4, fatom4, y, tmpmem, atomtype, alist, 
//                 pairlist, pairnumsum, tmpint, natom, Nij);            

    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    this->podNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    this->poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, rij, Phi2, besselparams, 
            &tmpmem[3*Nij], rin, rcut, pairnumsum, atomtype, ai, aj, ti, tj, elemindex, pdegree2, 
            nbesselpars, nrbf2, nrbf3, nabf3, nelements, Nij, natom);                    
    
    if (pod.snaptwojmax>0) 
        this->snapdesc(eatom4, fatom4, rij, &tmpmem[3*Nij], atomtype, ai, aj, ti, tj, natom, Nij);            
    
    // global descriptors for one-body, two-body, three-body, and four-bodt linear potentials    
    podArraySetValue(tmpmem, 1.0, natom);
    
    char cht = 'T';
    double one = 1.0, zero = 0.0;    
    int inc1 = 1;
    DGEMV(&cht, &natom, &nd1234, &one, eatom1, &natom, tmpmem, &inc1, &zero, gd, &inc1);            
}

void CPOD::quadratic_descriptors(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
        int M2, int M3, int N)
{
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = 0; m2<M2; m2++)
        {
            int m = m2 + M2*m3;
            d23[m] = d2[m2]*d3[m3];                
            for (int n=0; n<N; n++)
                dd23[n + N*m] = d2[m2]*dd3[n + N*m3] + dd2[n + N*m2]*d3[m3];
        }
}

void CPOD::quadratic_descriptors(double* d33, double *dd33, double *d3, double *dd3, int M3, int N)
{
    int m = 0;
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = m3; m2<M3; m2++)
        {            
            d33[m] = d3[m2]*d3[m3];                
            for (int n=0; n<N; n++)
                dd33[n + N*m] = d3[m2]*dd3[n + N*m3] + dd3[n + N*m2]*d3[m3];
            m += 1;
        }
}

void CPOD::cubic_descriptors(double* d234, double *dd234, double* d2, double *d3, double *d4, 
        double* dd2, double *dd3, double *dd4, int M2, int M3, int M4, int N)
{
    for (int m4 = 0; m4<M4; m4++)
        for (int m3 = 0; m3<M3; m3++)
            for (int m2 = 0; m2<M2; m2++)
            {
                int m = m2 + M2*m3 + M2*M3*m4;
                d234[m] = d2[m2]*d3[m3]*d4[m4];                
                for (int n=0; n<N; n++)
                    dd234[n + N*m] = d2[m2]*d3[m3]*dd4[n + N*m4] + 
                                     d2[m2]*dd3[n + N*m3]*d4[m4] + 
                                     dd2[n + N*m2]*d3[m3]*d4[m4];
            }
}

void CPOD::cubic_descriptors(double* d333, double *Dd333, double *d3, double *Dd3, int M3, int N)
{
    int m = 0;
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = m3; m2<M3; m2++)
            for (int m1 = m2; m1<M3; m1++)
            {            
                d333[m] = d3[m1]*d3[m2]*d3[m3];                
                for (int n=0; n<N; n++)
                    Dd333[n + N*m] = d3[m1]*d3[m2]*Dd3[n + N*m3] + d3[m1]*Dd3[n + N*m2]*d3[m3] + Dd3[n + N*m1]*d3[m2]*d3[m3];
                m += 1;
            }
}

double CPOD::quadratic_coefficients(double *c2, double *c3, double *d2, double *d3, 
        double *coeff23, int *quadratic, int nc2, int nc3)
{        
    int nd2 = quadratic[0]*nc2;
    int nd3 = quadratic[1]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=0; k< nd2; k++) {
            energy += coeff23[m]*d3[j]*d2[k];  
            c2[k] += coeff23[m]*d3[j];
            c3[j] += coeff23[m]*d2[k];
            m += 1;        
        }
            
    return energy;
}

double CPOD::quadratic_coefficients(double *c3, double *d3, double *coeff33, 
        int *quadratic, int nc3)
{        
    int nd3 = quadratic[0]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=j; k< nd3; k++) {
            energy += coeff33[m]*d3[j]*d3[k];  
            c3[k] += coeff33[m]*d3[j];
            c3[j] += coeff33[m]*d3[k];
            m += 1;        
        }
            
    return energy;
}

double CPOD::cubic_coefficients(double *c2, double *c3, double *c4, double *d2, double *d3, double *d4, 
        double *coeff234, int *cubic, int nc2, int nc3, int nc4)
{        
    int nd2 = cubic[0]*nc2;
    int nd3 = cubic[1]*nc3;
    int nd4 = cubic[2]*nc4;
    
    double energy = 0.0;    
    int m = 0;
    for (int i=0; i< nd4; i++)
        for (int j=0; j< nd3; j++)
            for (int k=0; k< nd2; k++) {
                energy += coeff234[m]*d4[i]*d3[j]*d2[k];  
                c2[k] += coeff234[m]*d4[i]*d3[j];
                c3[j] += coeff234[m]*d4[i]*d2[k];
                c4[i] += coeff234[m]*d3[j]*d2[k];                
                m += 1;        
            }
            
    return energy;
}

double CPOD::cubic_coefficients(double *c3, double *d3, double *coeff333, int *cubic, int nc3)
{        
    int nd3 = cubic[0]*nc3;
        
    double energy = 0.0;
    
    int m = 0;
    for (int i=0; i< nd3; i++)
        for (int j=i; j< nd3; j++)
            for (int k=j; k< nd3; k++) {
                energy += coeff333[m]*d3[i]*d3[j]*d3[k];  
                c3[k] += coeff333[m]*d3[i]*d3[j];
                c3[j] += coeff333[m]*d3[i]*d3[k];
                c3[i] += coeff333[m]*d3[j]*d3[k];                
                m += 1;        
            }
                                
    return energy;
}

double CPOD::quadratic_coefficients(double *ce2, double *ce3, double *c2, double *c3, double *d2, double *d3, 
        double *coeff23, int *quadratic, int nc2, int nc3)
{        
    int nd2 = quadratic[0]*nc2;
    int nd3 = quadratic[1]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=0; k< nd2; k++) {
            energy += coeff23[m]*d3[j]*d2[k];  
            c2[k] += coeff23[m]*d3[j];
            c3[j] += coeff23[m]*d2[k];
            ce2[k] += coeff23[m]*d3[j]/2.0;
            ce3[j] += coeff23[m]*d2[k]/2.0;
            m += 1;        
        }
            
    return energy;
}

double CPOD::quadratic_coefficients(double *ce3, double *c3, double *d3, double *coeff33, 
        int *quadratic, int nc3)
{        
    int nd3 = quadratic[0]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=j; k< nd3; k++) {
            energy += coeff33[m]*d3[j]*d3[k];  
            c3[k] += coeff33[m]*d3[j];
            c3[j] += coeff33[m]*d3[k];
            ce3[k] += coeff33[m]*d3[j];
            ce3[j] += coeff33[m]*d3[k];
            m += 1;        
        }
            
    return energy;
}

double CPOD::cubic_coefficients(double *ce2, double *ce3, double *ce4, double *c2, double *c3, double *c4, 
        double *d2, double *d3, double *d4, double *coeff234, int *cubic, int nc2, int nc3, int nc4)
{        
    int nd2 = cubic[0]*nc2;
    int nd3 = cubic[1]*nc3;
    int nd4 = cubic[2]*nc4;
    
    double energy = 0.0;    
    int m = 0;
    for (int i=0; i< nd4; i++)
        for (int j=0; j< nd3; j++)
            for (int k=0; k< nd2; k++) {
                energy += coeff234[m]*d4[i]*d3[j]*d2[k];  
                c2[k] += coeff234[m]*d4[i]*d3[j];
                c3[j] += coeff234[m]*d4[i]*d2[k];
                c4[i] += coeff234[m]*d3[j]*d2[k];                
                ce2[k] += coeff234[m]*d4[i]*d3[j]/3.0;
                ce3[j] += coeff234[m]*d4[i]*d2[k]/3.0;
                ce4[i] += coeff234[m]*d3[j]*d2[k]/3.0;                
                m += 1;        
            }
            
    return energy;
}

double CPOD::cubic_coefficients(double *ce3, double *c3, double *d3, double *coeff333, int *cubic, int nc3)
{        
    int nd3 = cubic[0]*nc3;
        
    double energy = 0.0;
    
    int m = 0;
    for (int i=0; i< nd3; i++)
        for (int j=i; j< nd3; j++)
            for (int k=j; k< nd3; k++) {
                energy += coeff333[m]*d3[i]*d3[j]*d3[k];  
                c3[k] += coeff333[m]*d3[i]*d3[j];
                c3[j] += coeff333[m]*d3[i]*d3[k];
                c3[i] += coeff333[m]*d3[j]*d3[k];                
                ce3[k] += coeff333[m]*d3[i]*d3[j];
                ce3[j] += coeff333[m]*d3[i]*d3[k];
                ce3[i] += coeff333[m]*d3[j]*d3[k];                
                m += 1;        
            }
                                
    return energy;
}

double CPOD::calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, double *tmp, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;        
    int nd22 = pod.nd22;
    int nd23 = pod.nd23;
    int nd24 = pod.nd24;
    int nd33 = pod.nd33;
    int nd34 = pod.nd34;
    int nd44 = pod.nd44;
    int nd234 = pod.nd234;
    int nd333 = pod.nd333;
    int nd444 = pod.nd444;    
    int nc2 = pod.nc2;
    int nc3 = pod.nc3;
    int nc4 = pod.nc4;
        
    // two-body, three-body, and four-body descriptors
    double *d2 = &gd[nd1];
    double *d3 = &gd[nd1+nd2];
    double *d4 = &gd[nd1+nd2+nd3];
        
    // quadratic and cubic POD coefficients
    double *coeff22 = &coeff[nd1234];
    double *coeff23 = &coeff[nd1234+nd22];
    double *coeff24 = &coeff[nd1234+nd22+nd23];
    double *coeff33 = &coeff[nd1234+nd22+nd23+nd24];
    double *coeff34 = &coeff[nd1234+nd22+nd23+nd24+nd33];
    double *coeff44 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34];  
    double *coeff234 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44];  
    double *coeff333 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234];  
    double *coeff444 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234+nd333];  
            
    // effective POD coefficients for calculating force 
    double *c1 = &tmp[0];
    double *c2 = &tmp[nd1];  
    double *c3 = &tmp[nd1+nd2];
    double *c4 = &tmp[nd1+nd2+nd3];
        
    // calculate energy for linear potentials
    double energy = 0.0;
    for (int i=0; i< nd1234; i++) {
        c1[i] = 0.0;    
        energy += coeff[i]*gd[i];    
    }
    
    // calculate energy for quadratic22 potential
    if (nd22>0) energy += this->quadratic_coefficients(c2, d2, coeff22, pod.quadratic22, nc2);
    
    // calculate energy for quadratic23 potential
    if (nd23>0) energy += this->quadratic_coefficients(c2, c3, d2, d3, coeff23, pod.quadratic23, nc2, nc3);
    
    // calculate energy for quadratic24 potential
    if (nd24>0) energy += this->quadratic_coefficients(c2, c4, d2, d4, coeff24, pod.quadratic24, nc2, nc4);
    
    // calculate energy for quadratic33 potential
    if (nd33>0) energy += this->quadratic_coefficients(c3, d3, coeff33, pod.quadratic33, nc3);
    
    // calculate energy for quadratic34 potential
    if (nd34>0) energy += this->quadratic_coefficients(c3, c4, d3, d4, coeff34, pod.quadratic34, nc3, nc4);
    
    // calculate energy for quadratic44 potential
    if (nd44>0) energy += this->quadratic_coefficients(c4, d4, coeff44, pod.quadratic44, nc4);
    
    // calculate energy for cubic234 potential
    if (nd234>0) energy += this->cubic_coefficients(c2, c3, c4, d2, d3, d4, coeff234, pod.cubic234, nc2, nc3, nc4);
    
    // calculate energy for cubic333 potential
    if (nd333>0) energy += this->cubic_coefficients(c3, d3, coeff333, pod.cubic333, nc3);

    // calculate energy for cubic444 potential
    if (nd444>0) energy += this->cubic_coefficients(c4, d4, coeff444, pod.cubic444, nc4);
    
    // calculate effective POD coefficients
    for (int i=0; i< nd1234; i++) c1[i] += coeff[i];    
    
    // calculate force = gdd * c1
    char chn = 'N';
    double one = 1.0, zero = 0.0;    
    int inc1 = 1;
    DGEMV(&chn, &nforce, &nd1234, &one, gdd, &nforce, c1, &inc1, &zero, force, &inc1);        
        
    return energy;
}

double CPOD::energyforce_calculation(double *force, double *gd, double *gdd, double *coeff, double *y, 
    int *atomtype, int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij)         
{
    int dim = 3;
    int nd1234 = pod.nd1+pod.nd2+pod.nd3+pod.nd4;        
    double *tmpmem = &gdd[dim*natom*nd1234+natom*nd1234];
    
    // calculate POD and SNAP descriptors and their derivatives
    this->linear_descriptors(gd, gdd, y, tmpmem, atomtype, alist, 
            pairlist, pairnum, pairnumsum, tmpint, natom, Nij);                
    
    // calculate energy and force
    double energy = 0.0;
    energy = this->calculate_energyforce(force, gd, gdd, coeff, &gdd[dim*natom*nd1234], natom);
            
    return energy;
}

void CPOD::podNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *atomtype, int *alist, int inum, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ii;       // atom i
        int itype = atomtype[i];        
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start; // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = pairlist[l + start];  // atom j              
            int k = start + l;                                     
            ai[k]        = i;
            aj[k]        = alist[j];          
            ti[k]        = itype;       
            tj[k]        = atomtype[alist[j]];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
        }
    }    
};

void CPOD::podradialbasis(double *rbf, double *drbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
    for (int n=0; n<N; n++) {    
        double xij1 = xij[0+3*n];
        double xij2 = xij[1+3*n];
        double xij3 = xij[2+3*n];

        double dij = pow(xij1*xij1 + xij2*xij2 + xij3*xij3, 0.5);    
        double dr1 = xij1/dij;    
        double dr2 = xij2/dij;    
        double dr3 = xij3/dij;    

        double r = dij - rin;        
        double y = r/rmax;    
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double y7 = pow(y4, 1.5);
        double fcut = y6/exp(-1.0);
        double dfcut = ((3.0/(rmax*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;

        for (int j=0; j<nbesselpars; j++) {            
            double alpha = besselparams[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            double dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(1.0 - exp(-alpha));        

            for (int i=0; i<besseldegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*besseldegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
                double drbfdr = b*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r*r) + a*cos(a*x)*fcut*dx/r);
                drbf[0 + 3*nij] = drbfdr*dr1;
                drbf[1 + 3*nij] = drbfdr*dr2;
                drbf[2 + 3*nij] = drbfdr*dr3;
            }
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            int nij = n + N*p;     
            double a = pow(dij, (double) (i+1.0));
            rbf[nij] = fcut/a;
            double drbfdr = dfcut/a - (i+1.0)*fcut/(a*dij);  
            drbf[0 + 3*nij] = drbfdr*dr1;
            drbf[1 + 3*nij] = drbfdr*dr2;
            drbf[2 + 3*nij] = drbfdr*dr3;
        }
    }       
}

void CPOD::podtally2b(double *eatom, double *fatom, double *eij, double *fij, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        //int mij = (elemindex[typei + typej*nelements]-1)*nbf;
        for (int m=0; m<nbf; m++) {               
            //int im = i1 + natom*(m + mij);
            //int jm = j1 + natom*(m + mij);
            int im =  i1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int jm =  j1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int nm = n + N*m;
            eatom[im] += eij[nm];
            fatom[0 + 3*im] += fij[0 + 3*nm];
            fatom[1 + 3*im] += fij[1 + 3*nm];
            fatom[2 + 3*im] += fij[2 + 3*nm];
            fatom[0 + 3*jm] -= fij[0 + 3*nm];
            fatom[1 + 3*jm] -= fij[1 + 3*nm];
            fatom[2 + 3*jm] -= fij[2 + 3*nm];          
        }
    }
}

void CPOD::pod1body(double *eatom, double *fatom, int *atomtype, int nelements, int natom)
{
    for (int m=1; m<=nelements; m++)       
        for (int i=0; i<natom; i++)         
            eatom[i + natom*(m-1)] = (atomtype[i] == m) ? 1.0 : 0.0;
        
    for (int i=0; i<3*natom*nelements; i++)  
        fatom[i] = 0.0;
}

void CPOD::pod3body(double *eatom, double *fatom, double *yij, double *e2ij, double *f2ij, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, nijk, nijk3, typei, typej, typek, ij, ik, i, j, k;    
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    double uj, uk, rbf, drbf1, drbf2, drbf3, drbf4, drbf5, drbf6;
    double eijk, fj1, fj2, fj3, fk1, fk2, fk3;
            
    double *abf = &tmpmem[0];
    double *dabf1 = &tmpmem[nabf1];
    double *dabf2 = &tmpmem[2*nabf1];
    double *dabf3 = &tmpmem[3*nabf1];
    double *dabf4 = &tmpmem[4*nabf1];
    double *dabf5 = &tmpmem[5*nabf1];
    double *dabf6 = &tmpmem[6*nabf1];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            ij = lj + s;
            i = ai[ij];  // atom i                        
            j = aj[ij];  // atom j    
            typei = ti[ij] - 1;           
            typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);                         
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                ik = lk + s;
                k = aj[ik];  // atom k                       
                typek = tj[ik] - 1;                         
                xik1 = yij[0+dim*ik];  // xk - xi           
                xik2 = yij[1+dim*ik];  // xk - xi           
                xik3 = yij[2+dim*ik];  // xk - xi           
                riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
                rik = pow(riksq, 0.5); 

                xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
                costhe = xdot/(rij*rik);    
                costhe = costhe > 1.0 ? 1.0 : costhe;
                costhe = costhe < -1.0 ? -1.0 : costhe;
                xdot = costhe*(rij*rik);

                sinthe = pow(1.0 - costhe*costhe, 0.5);
                sinthe = sinthe > 1e-12 ? sinthe : 1e-12;    
                theta = acos(costhe);            
                dtheta = -1.0/sinthe; 

                tm1 = pow(rijsq,1.5)*rik;
                tm2 = rij*pow(riksq,1.5);
                tm1 = 1.0/tm1;
                tm2 = 1.0/tm2;
                dct1 = (xik1*rijsq - xij1*xdot)*tm1; 
                dct2 = (xik2*rijsq - xij2*xdot)*tm1;
                dct3 = (xik3*rijsq - xij3*xdot)*tm1;
                dct4 = (xij1*riksq - xik1*xdot)*tm2;
                dct5 = (xij2*riksq - xik2*xdot)*tm2;
                dct6 = (xij3*riksq - xik3*xdot)*tm2;

                for (int p=0; p <nabf1; p++) {
                    abf[p] = cos(p*theta);                
                    tm = -p*sin(p*theta)*dtheta;                    
                    dabf1[p] = tm*dct1;
                    dabf2[p] = tm*dct2;
                    dabf3[p] = tm*dct3;
                    dabf4[p] = tm*dct4;
                    dabf5[p] = tm*dct5;
                    dabf6[p] = tm*dct6;        
                }

                for (int m=0; m<nrbf; m++) {
                    uj = e2ij[lj + s + Nij*m];
                    uk = e2ij[lk + s + Nij*m];
                    rbf = uj*uk;
                    drbf1 = f2ij[0 + dim*(lj + s) + dim*Nij*m]*uk;
                    drbf2 = f2ij[1 + dim*(lj + s) + dim*Nij*m]*uk;
                    drbf3 = f2ij[2 + dim*(lj + s) + dim*Nij*m]*uk;                                                        
                    drbf4 = f2ij[0 + dim*(lk + s) + dim*Nij*m]*uj;
                    drbf5 = f2ij[1 + dim*(lk + s) + dim*Nij*m]*uj;
                    drbf6 = f2ij[2 + dim*(lk + s) + dim*Nij*m]*uj;     

                    for (int p=0; p <nabf1; p++) {
                        eijk = rbf*abf[p];
                        fj1 = drbf1*abf[p] + rbf*dabf1[p];
                        fj2 = drbf2*abf[p] + rbf*dabf2[p];
                        fj3 = drbf3*abf[p] + rbf*dabf3[p];
                        fk1 = drbf4*abf[p] + rbf*dabf4[p];
                        fk2 = drbf5*abf[p] + rbf*dabf5[p];
                        fk3 = drbf6*abf[p] + rbf*dabf6[p];

                        n = p + (nabf1)*m;
                        nijk = natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n);
                        eatom[i + nijk] += eijk;

                        nijk3 = 3*i + 3*nijk;                        
                        fatom[0 + nijk3] += fj1 + fk1;
                        fatom[1 + nijk3] += fj2 + fk2;
                        fatom[2 + nijk3] += fj3 + fk3;

                        nijk3 = 3*j + 3*nijk;    
                        fatom[0 + nijk3] -= fj1;
                        fatom[1 + nijk3] -= fj2;
                        fatom[2 + nijk3] -= fj3;

                        nijk3 = 3*k + 3*nijk;    
                        fatom[0 + nijk3] -= fk1;   
                        fatom[1 + nijk3] -= fk2;   
                        fatom[2 + nijk3] -= fk3;                           
                    }                    
                }
            }
        }
    }
}

void CPOD::poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
            double *fatom3, double *rij, double *Phi, double *besselparams, double *tmpmem, double rin, 
            double rcut, int *pairnumsum, int *atomtype, int *ai, int *aj, int *ti, int *tj, int *elemindex, 
            int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom)
{       
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *f2ij = &tmpmem[Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    

    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    this->podMatMul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);

    // one-body descriptors
    this->pod1body(eatom1, fatom1, atomtype, nelements, natom);

    // two-body descriptors
    this->podtally2b(eatom2, fatom2, e2ij, f2ij, ai, aj, ti, tj, elemindex, nelements, nrbf2, natom, Nij);   
    
    // three-body descriptors
    this->pod3body(eatom3, fatom3, rij, e2ij, f2ij, &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, 
             ai, aj, ti, tj, nrbf3, nabf, nelements, natom, Nij);            
}

/****************************************************************************************************************/

void snapBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax)
{
  // index list for cglist

  int jdim = twojmax + 1;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j + j2*jdim + j1*jdim*jdim] = idxcg_count;  
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idx_max[0] = idxcg_count;
          
  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  //idxu_max = idxu_count;
  idx_max[1] = idxu_count;
  
  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  int idxb_max = idxb_count;
  idx_max[2] = idxb_max;

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count*3 + 0] = j1;
          idxb[idxb_count*3 + 1] = j2;
          idxb[idxb_count*3 + 2] = j;  
          idxb_count++;
        }
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j + j2*jdim + j1*jdim*jdim] = idxb_count;    
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  int idxz_max = idxz_count;
  idx_max[3] = idxz_max;

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j + j2*jdim + j1*jdim*jdim] = idxz_count;    

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {

            idxz[idxz_count*10 + 0] = j1;
            idxz[idxz_count*10 + 1] = j2;
            idxz[idxz_count*10 + 2] = j;
            idxz[idxz_count*10 + 3] = PODMAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 4] = (2 * ma - j - (2 * idxz[idxz_count*10 + 3] - j1) + j2) / 2;
            idxz[idxz_count*10 + 5] = PODMIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count*10 + 3] + 1;
            idxz[idxz_count*10 + 6] = PODMAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 7] = (2 * mb - j - (2 * idxz[idxz_count*10 + 6] - j1) + j2) / 2;
            idxz[idxz_count*10 + 8] = PODMIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count*10 + 6] + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count*10 + 9] = jju;
              
            idxz_count++;
          }
      }
};

void snapInitRootpqArray(double *rootpqarray, int twojmax)
{
  int jdim = twojmax + 1;
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p*jdim + q] = sqrt(((double) p)/q);
};

double snapDeltacg(double *factorial, int j1, int j2, int j)
{
  double sfaccg = factorial[(j1 + j2 + j) / 2 + 1];
  return sqrt(factorial[(j1 + j2 - j) / 2] *
              factorial[(j1 - j2 + j) / 2] *
              factorial[(-j1 + j2 + j) / 2] / sfaccg);
};

void snapInitClebschGordan(double *cglist, double *factorial, int twojmax)
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = PODMAX(0, PODMAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= PODMIN((j1 + j2 - j) / 2,
                          PODMIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial[z] *
                 factorial[(j1 + j2 - j) / 2 - z] *
                 factorial[(j1 - aa2) / 2 - z] *
                 factorial[(j2 + bb2) / 2 - z] *
                 factorial[(j - j2 + aa2) / 2 + z] *
                 factorial[(j - j1 - bb2) / 2 + z]);
            }

            cc2 = 2 * m - j;
            dcg = snapDeltacg(factorial, j1, j2, j);
            sfaccg = sqrt(factorial[(j1 + aa2) / 2] *
                          factorial[(j1 - aa2) / 2] *
                          factorial[(j2 + bb2) / 2] *
                          factorial[(j2 - bb2) / 2] *
                          factorial[(j  + cc2) / 2] *
                          factorial[(j  - cc2) / 2] *
                          (j + 1));

            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}

void snapInitSna(double *rootpqarray, double *cglist, double *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax)
{
    snapBuildIndexList(idx_max, idxz, idxz_block, idxb, 
            idxb_block, idxu_block, idxcg_block, twojmax);
    
    snapInitRootpqArray(rootpqarray, twojmax);
    snapInitClebschGordan(cglist, factorial, twojmax);        
}

void CPOD::snapSetup(int twojmax, int ntypes)
{
    int backend = 1;
    sna.twojmax = twojmax;
    sna.ntypes = ntypes;
    
    int jdim = twojmax + 1;    
    int jdimpq = twojmax + 2;  
    
    TemplateMalloc(&sna.map, ntypes+1, backend);
    TemplateMalloc(&sna.idxcg_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxz_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxb_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxu_block, jdim, backend);   
    TemplateMalloc(&sna.idx_max, 5, backend);     
    
    int idxb_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) idxb_count++;
    
    int idxz_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
          for (int mb = 0; 2*mb <= j; mb++)
            for (int ma = 0; ma <= j; ma++)
              idxz_count++;
    
    int idxcg_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
                for (int m1 = 0; m1 <= j1; m1++)
                  for (int m2 = 0; m2 <= j2; m2++)
                    idxcg_count++;
            }
   
    TemplateMalloc(&sna.idxz, idxz_count*10, backend);        
    TemplateMalloc(&sna.idxb, idxb_count*3, backend);        
    
    TemplateMalloc(&sna.rcutsq, (ntypes+1)*(ntypes+1), backend);
    TemplateMalloc(&sna.radelem, ntypes+1, backend);
    TemplateMalloc(&sna.wjelem, ntypes+1, backend);
        
    TemplateMalloc(&sna.rootpqarray, jdimpq*jdimpq, backend);      
    TemplateMalloc(&sna.cglist, idxcg_count, backend);        
    TemplateMalloc(&sna.bzero, jdim, backend); 
    TemplateMalloc(&sna.fac, 168, backend); 
        
    for (int i=0; i<jdimpq*jdimpq; i++)
        sna.rootpqarray[i] = 0;
    
    double facn = 1.0;
    for (int i=0; i<168; i++) {
        sna.fac[i] = facn;
        facn = facn*(i+1);
    }
    
    snapInitSna(sna.rootpqarray, sna.cglist, sna.fac, sna.idx_max, sna.idxz, 
        sna.idxz_block, sna.idxb, sna.idxb_block, sna.idxu_block, sna.idxcg_block, sna.twojmax);    

    sna.idxcg_max = sna.idx_max[0];
    sna.idxu_max = sna.idx_max[1];
    sna.idxb_max = sna.idx_max[2];
    sna.idxz_max = sna.idx_max[3];        
}

void CPOD::InitSnap()
{
  double *elemradius = pod.snapelementradius;
  double *elemweight = pod.snapelementweight;
  double rcutfac = pod.rcut;
  double rmin0 = 0.0;  
  double rfac0 = pod.snaprfac0;
  int twojmax = pod.snaptwojmax;
  int ntypes = pod.nelements;
  int chemflag = pod.snapchemflag;
        
  int backend=1;  
  int bzeroflag = 0;  
  int switchflag = 1;
  int wselfallflag = 0; 
  int bnormflag = chemflag; 
  double wself=1.0;  
              
  // Calculate maximum cutoff for all elements
  double rcutmax = 0.0;
  for (int ielem = 0; ielem < ntypes; ielem++)
    rcutmax = PODMAX(2.0*elemradius[ielem]*rcutfac,rcutmax);
      
  this->snapSetup(twojmax, ntypes);  
  TemplateCopytoDevice(&sna.radelem[1], elemradius, ntypes, backend);
  TemplateCopytoDevice(&sna.wjelem[1], elemweight, ntypes, backend);  
  podArrayFill(&sna.map[1], (int) 0, ntypes);
  
  double cutsq[100];
  for (int i=0; i<ntypes; i++)
      for (int j=0; j<ntypes; j++) {
          double cut = (elemradius[i] + elemradius[j])*rcutfac;
          cutsq[j+1 + (i+1)*(ntypes+1)] = cut*cut;
      }
  TemplateCopytoDevice(sna.rcutsq, cutsq, (ntypes+1)*(ntypes+1), backend);  
  
  if (bzeroflag) {
    double www = wself*wself*wself;
    double bzero[100];
    for (int j = 0; j <= twojmax; j++)
      if (bnormflag)
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
    TemplateCopytoDevice(sna.bzero, bzero, twojmax+1, backend);
  }
  
  int nelements = ntypes;
  if (!chemflag)    
    nelements = 1;
  
  sna.nelements = nelements;    
  sna.ndoubles = nelements*nelements;   // number of multi-element pairs
  sna.ntriples = nelements*nelements*nelements;   // number of multi-element triplets      
  sna.bnormflag = bnormflag;
  sna.chemflag = chemflag;    
  sna.switchflag = switchflag;
  sna.bzeroflag = bzeroflag;
  sna.wselfallflag = wselfallflag;
  sna.wself = wself;
  sna.rmin0 = rmin0;
  sna.rfac0 = rfac0;
  sna.rcutfac = rcutfac;
  sna.rcutmax = rcutmax;      
  sna.ncoeff = sna.idxb_max*sna.ntriples;
}

void CPOD::snapComputeUlist(double *Sr, double *Si, double *dSr, double *dSi, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
{    
  double *Srx = &dSr[0];  
  double *Sry = &dSr[idxu_max*ijnum];  
  double *Srz = &dSr[2*idxu_max*ijnum];  
  double *Six = &dSi[0];  
  double *Siy = &dSi[idxu_max*ijnum];  
  double *Siz = &dSi[2*idxu_max*ijnum];  

  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);
    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan(theta0);                
    double dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
            
    double sfac = 0.0, dsfac = 0.0;        
    if (switch_flag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switch_flag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else if(r > rcutij) {
            sfac = 0.0;
            dsfac = 0.0;
        }
        else {
            double rcutfac0 = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac0) + 1.0);   
            dsfac = -0.5 * sin((r - rmin0) * rcutfac0) * rcutfac0;                    
//             double y = (r - rmin0)/(rcutij - rmin0);    
//             double y2 = y*y;
//             double y3 = 1.0 - y2*y;
//             double y4 = y3*y3 + 1e-6;
//             double y5 = pow(y4, 0.5);
//             double y6 = exp(-1.0/y5);
//             double y7 = pow(y4, 1.5);
//             sfac = y6/exp(-1.0);
//             dsfac = ((3.0/((rcutij - rmin0)*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;            
        }
    } 
    sfac *= wjelem[tj[ij]];
    dsfac *= wjelem[tj[ij]];

    //sfac = 1.0; 
    //dsfac = 0.0;
    
    double r0inv, dr0invdr;
    double a_r, a_i, b_r, b_i;
    double da_r[3], da_i[3], db_r[3], db_i[3];
    double dz0[3], dr0inv[3];
    double rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    for (int k = 0; k < 3; k++) {
        da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
        da_i[k] = -z * dr0inv[k];
    }
    da_i[2] += -r0inv;

    for (int k = 0; k < 3; k++) {
        db_r[k] = y * dr0inv[k];
        db_i[k] = -x * dr0inv[k];
    }
    db_i[0] += -r0inv;
    db_r[1] += r0inv;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    Srx[ij+0*ijnum] = 0.0;
    Six[ij+0*ijnum] = 0.0;
    Sry[ij+0*ijnum] = 0.0;
    Siy[ij+0*ijnum] = 0.0;
    Srz[ij+0*ijnum] = 0.0;
    Siz[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            Srx[ij+jju*ijnum] = 0.0;
            Six[ij+jju*ijnum] = 0.0;
            Sry[ij+jju*ijnum] = 0.0;
            Siy[ij+jju*ijnum] = 0.0;
            Srz[ij+jju*ijnum] = 0.0;
            Siz[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                double u_r = Sr[njjup];
                double u_i = Si[njjup];
                double ux_r = Srx[njjup];
                double ux_i = Six[njjup];
                double uy_r = Sry[njjup];
                double uy_i = Siy[njjup];
                double uz_r = Srz[njjup];
                double uz_i = Siz[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
                Srx[njju] += rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
                Six[njju] += rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
                Sry[njju] += rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
                Siy[njju] += rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
                Srz[njju] += rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
                Siz[njju] += rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                Srx[njju1] = -rootpq * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
                Six[njju1] = -rootpq * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
                Sry[njju1] = -rootpq * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
                Siy[njju1] = -rootpq * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
                Srz[njju1] = -rootpq * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
                Siz[njju1] = -rootpq * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
                jju++;
                jjup++;
            }
            jju++;
        }
                   
        jju = idxu_block[j];
        jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb <= j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                if (mapar == 1) {
                    Sr[njjup] = Sr[njju];
                    Si[njjup] = -Si[njju];
                    if (j%2==1 && mb==(j/2)) {                    
                    Srx[njjup] =  Srx[njju];
                    Six[njjup] = -Six[njju];
                    Sry[njjup] =  Sry[njju];
                    Siy[njjup] = -Siy[njju];
                    Srz[njjup] =  Srz[njju];
                    Siz[njjup] = -Siz[njju];
                    }
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                    if (j%2==1 && mb==(j/2)) {
                    Srx[njjup] = -Srx[njju];
                    Six[njjup] =  Six[njju];
                    Sry[njjup] = -Sry[njju];
                    Siy[njjup] =  Siy[njju];
                    Srz[njjup] = -Srz[njju];
                    Siz[njjup] =  Siz[njju];                    
                    }
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        

    for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            int ijk = ij+jju*ijnum;               
            Srx[ijk] = dsfac * Sr[ijk] * ux + sfac * Srx[ijk]; 
            Six[ijk] = dsfac * Si[ijk] * ux + sfac * Six[ijk]; 
            Sry[ijk] = dsfac * Sr[ijk] * uy + sfac * Sry[ijk]; 
            Siy[ijk] = dsfac * Si[ijk] * uy + sfac * Siy[ijk]; 
            Srz[ijk] = dsfac * Sr[ijk] * uz + sfac * Srz[ijk]; 
            Siz[ijk] = dsfac * Si[ijk] * uz + sfac * Siz[ijk];                  
            jju++;
          }
    }
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};

void CPOD::snapZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum)
{
    int N1 = inum;
    int N2 = N1*(twojmax+1);
    int N3 = N2*nelements;                                
    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;  // inum*(twojmax+1)
        int ii = l%N1;    // inum
        int j = (l-ii)/N1; // (twojmax+1)
        int jelem = (idx-l)/N2; // nelements   
        //int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        int ielem = (chemflag) ? map[type[ii]]: 0;                
        int jju = idxu_block[j];        
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                int n = ii + inum*jju + inum*idxu_max*jelem;        
                Stotr[n] = 0.0;
                Stoti[n] = 0.0;
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)
                        Stotr[n] = wself; ///// double check this
                jju++;
            }
        }
        
    }                    
};

 void snapKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *map, int *ai, int *tj, int inum, int ijnum, int N1, int N2, int chemflag)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max   
        int jelem = (chemflag) ? map[tj[ij]] : 0;     
        int i = ai[ij] + inum*jju + N1*jelem;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
 void snapKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *ai, int inum, int ijnum, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max        
        int i = ai[ij] + inum*jju;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
void CPOD::snapAddUarraytot(double *Stotr, double *Stoti, double *Sr, 
        double *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag)
{   
    int N1 = inum*idxu_max;    
    int N2 = ijnum*idxu_max;    
    if (chemflag==0) {
        snapKernelAddUarraytot(Stotr, Stoti, Sr, Si, ai, inum, ijnum, N2);          
    } else
        snapKernelAddUarraytot(Stotr, Stoti, Sr, Si, map, ai, tj, inum, ijnum, N1, N2, chemflag);  
};

void CPOD::snapComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum)
{
    int jdim = twojmax + 1;    
    int N1 = inum;    
    int N2 = N1*idxz_max;
    int N3 = N2*nelements*nelements;                                
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;   //  inum*idxz_max
        int ii = l%inum;    // inum
        int jjz = (l-ii)/inum; // idxz_max
        int ielem = (idx-l)/N2;  // nelements*nelements  
        int elem2 = ielem%nelements; // nelements
        int elem1 = (ielem-elem2)/nelements; // nelements
              
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];
        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;

        const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
        double qr = 0.0;
        double qi = 0.0;          
        for (int ib = 0; ib < nb; ib++) {
            double suma1_r = 0.0;
            double suma1_i = 0.0;

            // Stotr: inum*idxu_max*nelements  
            const double *u1_r = &Stotr[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u1_i = &Stoti[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u2_r = &Stotr[ii + inum*jju2 + inum*idxu_max*elem2];
            const double *u2_i = &Stoti[ii + inum*jju2 + inum*idxu_max*elem2];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
                suma1_r += cgblock[icga] * (u1_r[inum*ma1] * u2_r[inum*ma2] - u1_i[inum*ma1] * u2_i[inum*ma2]);
                suma1_i += cgblock[icga] * (u1_r[inum*ma1] * u2_i[inum*ma2] + u1_i[inum*ma1] * u2_r[inum*ma2]);
                ma1++;
                ma2--;
                icga += j2;
            } // end loop over ia

            qr += cgblock[icgb] * suma1_r;
            qi += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
        } // end loop over ib
        
        if (bnorm_flag) {
            qr /= (j+1);
            qi /= (j+1);
        }        
        
        zlist_r[idx] = qr;
        zlist_i[idx] = qi;          
    }
};

void snapKernelComputeBi1(double *blist, double *zlist_r, double *zlist_i, 
        double *Stotr, double *Stoti, int *idxb, int *idxu_block, int *idxz_block, int jdim,         
        int nelements, int nelemsq, int nz_max, int nu_max, int inum, int N2, int N3)
{    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int idu;
        int idz;
        double sumzu = 0.0;
        for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);
        } // end if jeven

        blist[idx] = 2.0 * sumzu;                
    }
}
void snapKernelComputeBi2(double *blist, double *bzero,int *ilist, int *type,
       int *map, int *idxb, int nelements, int nb_max, int inum, int N2, int chemflag)
{        
    for (int idx=0; idx < N2; idx++) {        
        int ii = idx%inum;        
        int jjb = (idx-ii)/inum;    
        
        int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
        int itriple = (ielem*nelements+ielem)*nelements+ielem;

        const int j = idxb[jjb*3 + 2];  
        blist[ii + inum*jjb + nb_max*itriple] -= bzero[j];                
    }
}
void snapKernelComputeBi4(double *blist, double *bzero,
       int *idxb, int inum, int N2, int N3)
{        
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;
        int ii = l%inum;        
        int jjb = (l-ii)/inum;    
        int j = idxb[jjb*3 + 2];  
        blist[idx] -= bzero[j];        
    }
}
void CPOD::snapComputeBi1(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        int *idxb, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    //int nb_max = idxb_max*inum;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;
    int jdim = twojmax+1;

    snapKernelComputeBi1(blist, zlist_r, zlist_i, Stotr, Stoti, idxb, idxu_block, idxz_block, 
            jdim, nelements, nelemsq, nz_max, nu_max, inum, N2, N3);

};
void snapComputeBi2(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;
    int jdim = twojmax+1;

    snapKernelComputeBi1(blist, zlist_r, zlist_i, Stotr, Stoti, 
            idxb, idxu_block, idxz_block, jdim, nelements, nelemsq, nz_max, nu_max, inum, N2, N3);

    if (bzero_flag) {
        if (!wselfall_flag) {
            snapKernelComputeBi2(blist, bzero, ilist, type, map, 
                    idxb, nelements, nb_max, inum, N2, chemflag);
        }
        else {
            snapKernelComputeBi4(blist, bzero, idxb, inum, N2, N3);            
        }
    }
};

void CPOD::snapComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, 
        double *dulist_r, double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum)
{                    
    int nz_max = idxz_max*inum;
    int nb_max = idxb_max*ijnum;    
    int nu_max = idxu_max*ijnum;    
    int N2 = ijnum*idxb_max;
    int jdim = twojmax+1;

    for (int i=0; i<nb_max*3*nelements*nelements*nelements; i++)
        dblist[i] = 0.0;
    //snapArraySetValue(dblist, (double) 0.0, nb_max*3*nelements*nelements*nelements);
        
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;              
        int jjb = (idx-ij)/ijnum;                              
        int elem3 = (chemflag) ? map[tj[ij]] : 0;//(chemflag) ? map[type[alist[aj[ij]]]] : 0;
        int i = ai[ij]; // atom i
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

       // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
        for(int elem1 = 0; elem1 < nelements; elem1++)
            for(int elem2 = 0; elem2 < nelements; elem2++) {

            int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
            int jju = idxu_block[j];
            int idouble = elem1*nelements+elem2;
            int itriple = (elem1*nelements+elem2)*nelements+elem3;
            int nimax = nz_max*idouble;                      

            double *dbdr = &dblist[nb_max*3*itriple];
            double sumzdu_r[3];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j; mb++)
              for (int ma = 0; ma <= j; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j even, handle middle column

            if (j % 2 == 0) {
              int mb = j / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              }

              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if jeven
                                
            for (int k = 0; k < 3; k++)
              dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
            
            // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
            double j1fac = (j + 1) / (j1 + 1.0);
            idouble = elem1*nelements+elem2;
            itriple = (elem3*nelements+elem2)*nelements+elem1;            
            //jjz = idxz_block[j][j2][j1];
            jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
            jju = idxu_block[j1];

            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j1; mb++)
              for (int ma = 0; ma <= j1; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j1 even, handle middle column

            if (j1 % 2 == 0) {
              int mb = j1 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              }
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j1even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j1fac;

            // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)
            double j2fac = (j + 1) / (j2 + 1.0);
            idouble = elem2*nelements+elem1;
            itriple = (elem1*nelements+elem3)*nelements+elem2;
            //jjz = idxz_block[j][j1][j2];
            jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
            jju = idxu_block[j2];
            
            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j2; mb++)
              for (int ma = 0; ma <= j2; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                            
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j2 even, handle middle column

            if (j2 % 2 == 0) {
              int mb = j2 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                                                 
                jjz++;
                jju++;
              }
              
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j2even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j2fac;
          }        
    }
}

void snapTallyBispectrumDeriv(double *db, double *dbdr, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int ntype)
{   
    //snapArraySetValue(db, (double) 0.0, inum*3*nperdim*ntype);
    for (int i=0; i<inum*3*ncoeff*ntype; i++)
        db[i] = 0.0;
        
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = ncoeff*(itype-1);        
        int nii = inum*3*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        //printf("%i %i %i %i %i %i %i %i %i %i %i \n", idx, ij, icoeff, ii, i, j, itype, n, nii, nij, quadraticflag);

        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];        
        db[0 + 3*i + nii] += bix; 
        db[1 + 3*i + nii] += biy;
        db[2 + 3*i + nii] += biz;
        db[0 + 3*j + nii] -= bix;
        db[1 + 3*j + nii] -= biy;
        db[2 + 3*j + nii] -= biz;        
    }
}

void CPOD::snapdesc(double *blist, double *bd, double *rij, double *tmpmem, int *atomtype, int *ai, 
        int *aj, int *ti, int *tj, int natom, int Nij)            
{    
    int dim = 3;    
    //int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    //int ncoeffall = sna.ncoeffall;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    //int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    //int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    //int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    //double rcutmax = sna.rcutmax;        
    //double *bzero = sna.bzero;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    //double *rcutsq = sna.rcutsq;    
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
        
//     int *ai = &tmpint[0];     // Nij
//     int *aj = &tmpint[Nij];   // Nij 
//     int *ti = &tmpint[2*Nij]; // Nij
//     int *tj = &tmpint[3*Nij]; // Nij
    
    int ne = 0;
    //double *rij = &tmpmem[ne]; // Nij*dim    
    //ne += Nij*dim; 
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *dUr = &tmpmem[ne];
    ne += idxu_max*dim*Nij;
    double *dUi = &tmpmem[ne];
    ne += idxu_max*dim*Nij;    
    //double *blist = &tmpmem[ne]; // idxb_max*ntriples*natom          
    //ne += idxb_max*ntriples*natom;    
    double *dblist = &tmpmem[ne]; // idxb_max*ntriples*dim*Nij          
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    snapComputeUlist(Ur, Ui, dUr, dUi, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
    
    snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, ai, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, ai, tj, idxu_max, natom, Nij, chemflag);    
        
    snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);

    snapComputeBi1(blist, Zr, Zi, Utotr, Utoti, idxb, idxu_block, idxz_block, twojmax, idxb_max, 
            idxu_max, idxz_max, nelem, natom);           
            
    snapComputeDbidrj(dblist, Zr, Zi, dUr, dUi, idxb, idxu_block, idxz_block, map, ai, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, natom, Nij);
    
    snapTallyBispectrumDeriv(bd, dblist, ai, aj, ti, natom, Nij, ncoeff, ntypes);                        
}

/******************************************************************************/

void CPOD::podNeighPairs(double *rij, double *x, int *idxi, int *ai, int *aj,  int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *jlist, int *alist, int inum)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        //int gi = ilist[ii];       // atom i
        int gi = ii;       // atom i
        int itype = atomtype[gi];
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;
            int gj = jlist[k];  // atom j                                                               
            idxi[k]      = ii;
            ai[k]        = alist[gi];
            aj[k]        = alist[gj];          
            ti[k]        = itype;       
            tj[k]        = atomtype[aj[k]];              
            rij[k*3+0]   = x[gj*3+0] -  x[gi*3+0];  // xj - xi            
            rij[k*3+1]   = x[gj*3+1] -  x[gi*3+1];  // xj - xi            
            rij[k*3+2]   = x[gj*3+2] -  x[gi*3+2];  // xj - xi            
        }
    }            
};

int CPOD::lammpsNeighPairs(double *rij, double **x, double rcutsq, int *idxi, int *ai, int *aj,  int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *numneigh, int *ilist, int **jlist, int inum)
{  
    
    int ninside = 0;
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int gi = ilist[ii];       // atom i
        int itype = atomtype[gi];
        int m = numneigh[gi];
        pairnumsum[ii+1] = 0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int gj = jlist[gi][l];  // atom j     
            double delx   = x[gj][0] -  x[gi][0];  // xj - xi            
            double dely   = x[gj][1] -  x[gi][1];  // xj - xi            
            double delz   = x[gj][2] -  x[gi][2];  // xj - xi            
            double rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < rcutsq && rsq > 1e-20) {
                rij[ninside*3 + 0] = delx;
                rij[ninside*3 + 1] = dely;
                rij[ninside*3 + 2] = delz;
                idxi[ninside]      = ii;
                ai[ninside]        = gi;
                aj[ninside]        = gj;          
                ti[ninside]        = itype;       
                tj[ninside]        = atomtype[gj];       
                ninside++;
                pairnumsum[ii+1] += 1;
            }                        
        }
    }    
    
    pairnumsum[0] = 0;
    for (int ii=0; ii<inum; ii++)
        pairnumsum[ii+1] = pairnumsum[ii+1] + pairnumsum[ii];
    
    
    return ninside;
};

void CPOD::podradialbasis(double *rbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
    for (int n=0; n<N; n++) {    
        double xij1 = xij[0+3*n];
        double xij2 = xij[1+3*n];
        double xij3 = xij[2+3*n];

        double dij = pow(xij1*xij1 + xij2*xij2 + xij3*xij3, 0.5);    
        double r = dij - rin;        
        double y = r/rmax;    
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double fcut = y6/exp(-1.0);

        for (int j=0; j<nbesselpars; j++) {            
            double x =  (1.0 - exp(-besselparams[j]*r/rmax))/(1.0-exp(-besselparams[j]));
            for (int i=0; i<besseldegree; i++)                 
                rbf[n + N*i + N*besseldegree*j] = ((sqrt(2.0/(rmax))/(i+1)))*fcut*sin((i+1)*M_PI*x)/r;            
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            double a = pow(dij, (double) (i+1.0));
            rbf[n + N*p] = fcut/a;
        }
    }                    
}

void CPOD::podtally2b(double *eatom, double *eij, int *idxi, int *ti, int *tj, int *elemindex, 
        int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = idxi[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        for (int m=0; m<nbf; m++) {               
            int im =  i1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int nm = n + N*m;
            eatom[im] += eij[nm];
        }
    }
}

void CPOD::pod1body(double *eatom, int *atomtype, int nelements, int natom)
{
    for (int m=1; m<=nelements; m++)       
        for (int i=0; i<natom; i++)         
            eatom[i + natom*(m-1)] = (atomtype[i] == m) ? 1.0 : 0.0;        
}

void CPOD::pod3body(double *eatom, double *yij, double *e2ij, double *tmpmem, int *elemindex, int *pairnumsum, 
        int *idxi, int *ti, int *tj, int nrbf, int nabf, int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, nijk, typei, typej, typek, ij, ik, i, j, k;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, theta; 
    double uj, uk, rbf;
            
    double *abf = &tmpmem[0];
    
//     for (int m=0; m<nrbf; m++)
//         for (int p=0; p <nabf1; p++) 
//             for (int ii=0; ii<natom; ii++) {
//                 int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
//                 int s = pairnumsum[ii];        
//                 for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
//                     ij = lj + s;
//                     i = idxi[ij];  // atom i                        
//                     typei = ti[ij] - 1;           
//                     typej = tj[ij] - 1;                   
//                     xij1 = yij[0+dim*ij];  // xj - xi           
//                     xij2 = yij[1+dim*ij];  // xj - xi           
//                     xij3 = yij[2+dim*ij];  // xj - xi           
//                     rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
//                     rij = pow(rijsq, 0.5);                         
//                     for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
//                         ik = lk + s;
//                         typek = tj[ik] - 1;                         
//                         xik1 = yij[0+dim*ik];  // xk - xi           
//                         xik2 = yij[1+dim*ik];  // xk - xi           
//                         xik3 = yij[2+dim*ik];  // xk - xi           s
//                         riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
//                         rik = pow(riksq, 0.5); 
// 
//                         xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
//                         costhe = xdot/(rij*rik);    
//                         costhe = costhe > 1.0 ? 1.0 : costhe;
//                         costhe = costhe < -1.0 ? -1.0 : costhe;
//                         xdot = costhe*(rij*rik);
//                         theta = acos(costhe);            
// 
//                         uj = e2ij[lj + s + Nij*m];
//                         uk = e2ij[lk + s + Nij*m];
//                         rbf = uj*uk;
// 
//                         n = p + (nabf1)*m;
//                         nijk = natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n);
//                         eatom[i + nijk] += rbf*cos(p*theta);
//                     }
//                 }
//             }
    
    double *etm = &tmpmem[nabf1];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        
        for (int m=0; m<nrbf*nabf1*nelements2*nelements; m++)
            etm[m] = 0.0;
        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            ij = lj + s;
            i = idxi[ij];  // atom i                        
            typei = ti[ij] - 1;           
            typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);                         
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                ik = lk + s;
                typek = tj[ik] - 1;                         
                xik1 = yij[0+dim*ik];  // xk - xi           
                xik2 = yij[1+dim*ik];  // xk - xi           
                xik3 = yij[2+dim*ik];  // xk - xi           s
                riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
                rik = pow(riksq, 0.5); 

                xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
                costhe = xdot/(rij*rik);    
                costhe = costhe > 1.0 ? 1.0 : costhe;
                costhe = costhe < -1.0 ? -1.0 : costhe;
                xdot = costhe*(rij*rik);
                theta = acos(costhe);            

                for (int p=0; p <nabf1; p++) 
                    abf[p] = cos(p*theta);                

                for (int m=0; m<nrbf; m++) {
                    uj = e2ij[lj + s + Nij*m];
                    uk = e2ij[lk + s + Nij*m];
                    rbf = uj*uk;
                    for (int p=0; p <nabf1; p++) {                        
                        n = p + (nabf1)*m;
                        nijk = (elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n;                        
                        etm[nijk] += rbf*abf[p];
                        //nijk = natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n);
                        //eatom[i + nijk] += rbf*abf[p];
                    }                    
                }
            }
        }
        for (int m=0; m<nrbf*nabf1*nelements2*nelements; m++)
            eatom[ii + natom*m] += etm[m];
    }
}


void CPOD::poddesc_ij(double *eatom1, double *eatom2, double *eatom3, double *rij, double *Phi, double *besselparams, 
            double *tmpmem, double rin, double rcut, int *pairnumsum, int *atomtype, int *idxi, int *ti, int *tj, 
            int *elemindex, int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom)
{       
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *e2ijt = &tmpmem[Nij*nrbf]; // Nij*ns

//     struct timeval tv1, tv2;    
//     gettimeofday(&tv1, NULL); 
    
    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
      
    // one-body descriptors
    this->pod1body(eatom1, atomtype, nelements, natom);

    this->podtally2b(eatom2, e2ij, idxi, ti, tj, elemindex, nelements, nrbf2, natom, Nij);   

//     gettimeofday(&tv2, NULL);            
//     printf("\nExecution time (in millisec) for pod2body %g\n", 
//         (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
//         (double)(tv2.tv_sec -tv1.tv_sec )*1000);    
//     
//     gettimeofday(&tv1, NULL); 
    
    // three-body descriptors
    this->pod3body(eatom3, rij, e2ij, &tmpmem[Nij*nrbf], elemindex, pairnumsum, 
             idxi, ti, tj, nrbf3, nabf, nelements, natom, Nij);     
        
//     gettimeofday(&tv2, NULL);            
//     printf("\nExecution time (in millisec) for pod3body %g\n", 
//         (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
//         (double)(tv2.tv_sec -tv1.tv_sec )*1000);    
    
}

void CPOD::snapComputeUij(double *Sr, double *Si, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
{    
  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan(theta0);                
            
    double sfac = 0.0; 
    if (switch_flag == 0) {
        sfac = 1.0;
    }
    else if (switch_flag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
        }
        else if(r > rcutij) {
            sfac = 0.0;
        }
        else {
            double rcutfac0 = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac0) + 1.0);   
        }
    } 
    sfac *= wjelem[tj[ij]];
    
    double r0inv;
    double a_r, a_i, b_r, b_i;
    double rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                double u_r = Sr[njjup];
                double u_i = Si[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                jju++;
                jjup++;
            }
            jju++;
        }
                   
        jju = idxu_block[j];
        jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb <= j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                if (mapar == 1) {
                    Sr[njjup] = Sr[njju];
                    Si[njjup] = -Si[njju];
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};

void CPOD::snapdesc_ij(double *blist, double *rij, double *tmpmem, int *atomtype, int *idxi, 
        int *ti, int *tj, int natom, int Nij)            
{    
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    //int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    //int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
            
    int ne = 0;
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    this->snapComputeUij(Ur, Ui, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
        
    this->snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, idxi, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    this->snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, idxi, tj, idxu_max, natom, Nij, chemflag);    
        
    this->snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);

    this->snapComputeBi1(blist, Zr, Zi, Utotr, Utoti, idxb, idxu_block, idxz_block, twojmax, idxb_max, 
            idxu_max, idxz_max, nelem, natom);                        
}

void CPOD::linear_descriptors_ij(double *gd, double *eatom, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *idxi, int *ti, int *tj, int natom, int Nij)
{
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nrbf2 = pod.nbf2;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;     
    int *pdegree2 = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi2 = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    double *eatom1 = &eatom[0];
    double *eatom2 = &eatom[0+natom*nd1];
    double *eatom3 = &eatom[0+natom*(nd1+nd2)];
    double *eatom4 = &eatom[0+natom*(nd1+nd2+nd3)];    
        
    podArraySetValue(eatom1, 0.0, natom*nd1234);    
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    this->poddesc_ij(eatom1, eatom2, eatom3, rij, Phi2, besselparams, 
            tmpmem, rin, rcut, pairnumsum, atomtype, idxi, ti, tj, elemindex, pdegree2, 
            nbesselpars, nrbf2, nrbf3, nabf3, nelements, Nij, natom);                    
    
    // peratom snap descriptors
    if (pod.snaptwojmax>0) 
        this->snapdesc_ij(eatom4, rij, tmpmem, atomtype, idxi, ti, tj, natom, Nij);                
    
    // global descriptors for one-body, two-body, three-body, and four-bodt linear potentials    
    podArraySetValue(tmpmem, 1.0, natom);
    
    char cht = 'T';
    double one = 1.0;    
    int inc1 = 1;
    DGEMV(&cht, &natom, &nd1234, &one, eatom1, &natom, tmpmem, &inc1, &one, gd, &inc1);                        
}

double CPOD::calculate_energy(double *effectivecoeff, double *gd, double *coeff)
{        
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;        
    int nd22 = pod.nd22;
    int nd23 = pod.nd23;
    int nd24 = pod.nd24;
    int nd33 = pod.nd33;
    int nd34 = pod.nd34;
    int nd44 = pod.nd44;
    int nd234 = pod.nd234;
    int nd333 = pod.nd333;
    int nd444 = pod.nd444;    
    int nc2 = pod.nc2;
    int nc3 = pod.nc3;
    int nc4 = pod.nc4;
        
    // two-body, three-body, and four-body descriptors
    double *d2 = &gd[nd1];
    double *d3 = &gd[nd1+nd2];
    double *d4 = &gd[nd1+nd2+nd3];
        
    // quadratic and cubic POD coefficients
    double *coeff22 = &coeff[nd1234];
    double *coeff23 = &coeff[nd1234+nd22];
    double *coeff24 = &coeff[nd1234+nd22+nd23];
    double *coeff33 = &coeff[nd1234+nd22+nd23+nd24];
    double *coeff34 = &coeff[nd1234+nd22+nd23+nd24+nd33];
    double *coeff44 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34];  
    double *coeff234 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44];  
    double *coeff333 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234];  
    double *coeff444 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234+nd333];  
                    
    // calculate energy for linear potentials
    double energy = 0.0;
    for (int i=0; i< nd1234; i++) {
        effectivecoeff[i] = 0.0;    
        energy += coeff[i]*gd[i];    
    }
    
    // effective POD coefficients for calculating force 
    //double *c1 = &effectivecoeff[0];
    double *c2 = &effectivecoeff[nd1];  
    double *c3 = &effectivecoeff[nd1+nd2];
    double *c4 = &effectivecoeff[nd1+nd2+nd3];
    
    // calculate energy for quadratic22 potential
    if (nd22>0) energy += this->quadratic_coefficients(c2, d2, coeff22, pod.quadratic22, nc2);
    
    // calculate energy for quadratic23 potential
    if (nd23>0) energy += this->quadratic_coefficients(c2, c3, d2, d3, coeff23, pod.quadratic23, nc2, nc3);
    
    // calculate energy for quadratic24 potential
    if (nd24>0) energy += this->quadratic_coefficients(c2, c4, d2, d4, coeff24, pod.quadratic24, nc2, nc4);
    
    // calculate energy for quadratic33 potential
    if (nd33>0) energy += this->quadratic_coefficients(c3, d3, coeff33, pod.quadratic33, nc3);
    
    // calculate energy for quadratic34 potential
    if (nd34>0) energy += this->quadratic_coefficients(c3, c4, d3, d4, coeff34, pod.quadratic34, nc3, nc4);
    
    // calculate energy for quadratic44 potential
    if (nd44>0) energy += this->quadratic_coefficients(c4, d4, coeff44, pod.quadratic44, nc4);
    
    // calculate energy for cubic234 potential
    if (nd234>0) energy += this->cubic_coefficients(c2, c3, c4, d2, d3, d4, coeff234, pod.cubic234, nc2, nc3, nc4);
    
    // calculate energy for cubic333 potential
    if (nd333>0) energy += this->cubic_coefficients(c3, d3, coeff333, pod.cubic333, nc3);

    // calculate energy for cubic444 potential
    if (nd444>0) energy += this->cubic_coefficients(c4, d4, coeff444, pod.cubic444, nc4);
    
    // calculate effective POD coefficients
    for (int i=0; i< nd1234; i++) effectivecoeff[i] += coeff[i];    
                
    return energy;
}

double CPOD::calculate_energy(double *energycoeff, double *forcecoeff, double *gd, double *coeff)
{        
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;        
    int nd22 = pod.nd22;
    int nd23 = pod.nd23;
    int nd24 = pod.nd24;
    int nd33 = pod.nd33;
    int nd34 = pod.nd34;
    int nd44 = pod.nd44;
    int nd234 = pod.nd234;
    int nd333 = pod.nd333;
    int nd444 = pod.nd444;    
    int nc2 = pod.nc2;
    int nc3 = pod.nc3;
    int nc4 = pod.nc4;
        
    // two-body, three-body, and four-body descriptors
    double *d2 = &gd[nd1];
    double *d3 = &gd[nd1+nd2];
    double *d4 = &gd[nd1+nd2+nd3];
        
    // quadratic and cubic POD coefficients
    double *coeff22 = &coeff[nd1234];
    double *coeff23 = &coeff[nd1234+nd22];
    double *coeff24 = &coeff[nd1234+nd22+nd23];
    double *coeff33 = &coeff[nd1234+nd22+nd23+nd24];
    double *coeff34 = &coeff[nd1234+nd22+nd23+nd24+nd33];
    double *coeff44 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34];  
    double *coeff234 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44];  
    double *coeff333 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234];  
    double *coeff444 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234+nd333];  
                    
    // calculate energy for linear potentials
    double energy = 0.0;
    for (int i=0; i< nd1234; i++) {
        energycoeff[i] = 0.0;    
        forcecoeff[i] = 0.0;    
        energy += coeff[i]*gd[i];    
    }
    
    // effective POD coefficients for calculating force 
    //double *c1 = &forcecoeff[0];
    double *c2 = &forcecoeff[nd1];  
    double *c3 = &forcecoeff[nd1+nd2];
    double *c4 = &forcecoeff[nd1+nd2+nd3];

    double *ce2 = &energycoeff[nd1];  
    double *ce3 = &energycoeff[nd1+nd2];
    double *ce4 = &energycoeff[nd1+nd2+nd3];

    // calculate energy for quadratic22 potential
    if (nd22>0) energy += this->quadratic_coefficients(ce2, c2, d2, coeff22, pod.quadratic22, nc2);

    
    // calculate energy for quadratic23 potential
    if (nd23>0) energy += this->quadratic_coefficients(ce2, ce3, c2, c3, d2, d3, coeff23, pod.quadratic23, nc2, nc3);
    
    // calculate energy for quadratic24 potential
    if (nd24>0) energy += this->quadratic_coefficients(ce2, ce4, c2, c4, d2, d4, coeff24, pod.quadratic24, nc2, nc4);
    
    
    // calculate energy for quadratic33 potential
    if (nd33>0) energy += this->quadratic_coefficients(ce3, c3, d3, coeff33, pod.quadratic33, nc3);
    
    // calculate energy for quadratic34 potential
    if (nd34>0) energy += this->quadratic_coefficients(ce3, ce4, c3, c4, d3, d4, coeff34, pod.quadratic34, nc3, nc4);
    
    // calculate energy for quadratic44 potential
    if (nd44>0) energy += this->quadratic_coefficients(ce4, c4, d4, coeff44, pod.quadratic44, nc4);
    
    // calculate energy for cubic234 potential
    if (nd234>0) energy += this->cubic_coefficients(ce2, ce3, ce4, c2, c3, c4, d2, d3, d4, coeff234, pod.cubic234, nc2, nc3, nc4);
    
    // calculate energy for cubic333 potential
    if (nd333>0) energy += this->cubic_coefficients(ce3, c3, d3, coeff333, pod.cubic333, nc3);

    // calculate energy for cubic444 potential
    if (nd444>0) energy += this->cubic_coefficients(ce4, c4, d4, coeff444, pod.cubic444, nc4);
    
    // calculate effective POD coefficients
    for (int i=0; i< nd1234; i++) {
        energycoeff[i] += coeff[i];    
        forcecoeff[i] += coeff[i];    
    }
                
    return energy;
}

void CPOD::pod2body_force(double *force, double *fij, double *coeff2, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        for (int m=0; m<nbf; m++) {               
            int im =  3*i1;
            int jm =  3*j1;
            int nm = n + N*m;
            int km = (elemindex[typei + typej*nelements] - 1) + nelements2*m;
            double ce = coeff2[km];
            force[0 + im] += fij[0 + 3*nm]*ce;
            force[1 + im] += fij[1 + 3*nm]*ce;
            force[2 + im] += fij[2 + 3*nm]*ce;
            force[0 + jm] -= fij[0 + 3*nm]*ce;
            force[1 + jm] -= fij[1 + 3*nm]*ce;
            force[2 + jm] -= fij[2 + 3*nm]*ce;          
        }
    }
}

void CPOD::pod3body_force(double *force, double *yij, double *e2ij, double *f2ij, double *coeff3, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, c, nijk3, typei, typej, typek, ij, ik, i, j, k;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    //double uj, uk, rbf, drbf1, drbf2, drbf3, drbf4, drbf5, drbf6;
    //double fj1, fj2, fj3, fk1, fk2, fk3;
            
    double *abf = &tmpmem[0];
    double *dabf1 = &tmpmem[nabf1];
    double *dabf2 = &tmpmem[2*nabf1];
    double *dabf3 = &tmpmem[3*nabf1];
    double *dabf4 = &tmpmem[4*nabf1];
    double *dabf5 = &tmpmem[5*nabf1];
    double *dabf6 = &tmpmem[6*nabf1];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            ij = lj + s;
            i = ai[ij];  // atom i                        
            j = aj[ij];  // atom j    
            typei = ti[ij] - 1;           
            typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);            
            
            double fixtmp,fiytmp,fiztmp;
            fixtmp = fiytmp = fiztmp = 0.0;      
            double fjxtmp,fjytmp,fjztmp;
            fjxtmp = fjytmp = fjztmp = 0.0;                  
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                ik = lk + s;
                k = aj[ik];  // atom k                       
                typek = tj[ik] - 1;                         
                xik1 = yij[0+dim*ik];  // xk - xi           
                xik2 = yij[1+dim*ik];  // xk - xi           
                xik3 = yij[2+dim*ik];  // xk - xi           s
                riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
                rik = pow(riksq, 0.5); 

                xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
                costhe = xdot/(rij*rik);    
                costhe = costhe > 1.0 ? 1.0 : costhe;
                costhe = costhe < -1.0 ? -1.0 : costhe;
                xdot = costhe*(rij*rik);

                sinthe = pow(1.0 - costhe*costhe, 0.5);
                sinthe = sinthe > 1e-12 ? sinthe : 1e-12;    
                theta = acos(costhe);            
                dtheta = -1.0/sinthe; 

                tm1 = pow(rijsq,1.5)*rik;
                tm2 = rij*pow(riksq,1.5);
                tm1 = 1.0/tm1;
                tm2 = 1.0/tm2;
                dct1 = (xik1*rijsq - xij1*xdot)*tm1; 
                dct2 = (xik2*rijsq - xij2*xdot)*tm1;
                dct3 = (xik3*rijsq - xij3*xdot)*tm1;
                dct4 = (xij1*riksq - xik1*xdot)*tm2;
                dct5 = (xij2*riksq - xik2*xdot)*tm2;
                dct6 = (xij3*riksq - xik3*xdot)*tm2;

                for (int p=0; p <nabf1; p++) {
                    abf[p] = cos(p*theta);                
                    tm = -p*sin(p*theta)*dtheta;                    
                    dabf1[p] = tm*dct1;
                    dabf2[p] = tm*dct2;
                    dabf3[p] = tm*dct3;
                    dabf4[p] = tm*dct4;
                    dabf5[p] = tm*dct5;
                    dabf6[p] = tm*dct6;        
                }

                double fjx = 0.0, fjy = 0.0, fjz = 0.0;
                double fkx = 0.0, fky = 0.0, fkz = 0.0;
                
                for (int m=0; m<nrbf; m++) {
                    double uj = e2ij[lj + s + Nij*m];
                    double uk = e2ij[lk + s + Nij*m];
                    double rbf = uj*uk;
                    double drbf1 = f2ij[0 + dim*(lj + s) + dim*Nij*m]*uk;
                    double drbf2 = f2ij[1 + dim*(lj + s) + dim*Nij*m]*uk;
                    double drbf3 = f2ij[2 + dim*(lj + s) + dim*Nij*m]*uk;                                                        
                    double drbf4 = f2ij[0 + dim*(lk + s) + dim*Nij*m]*uj;
                    double drbf5 = f2ij[1 + dim*(lk + s) + dim*Nij*m]*uj;
                    double drbf6 = f2ij[2 + dim*(lk + s) + dim*Nij*m]*uj;     

                    for (int p=0; p <nabf1; p++) {
                        tm = abf[p];
                        double fj1 = drbf1*tm + rbf*dabf1[p];
                        double fj2 = drbf2*tm + rbf*dabf2[p];
                        double fj3 = drbf3*tm + rbf*dabf3[p];
                        double fk1 = drbf4*tm + rbf*dabf4[p];
                        double fk2 = drbf5*tm + rbf*dabf5[p];
                        double fk3 = drbf6*tm + rbf*dabf6[p];

                        n = p + (nabf1)*m;
                        c = (elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n;
                        tm = coeff3[c];
                        
                        fjx += fj1*tm;
                        fjy += fj2*tm;
                        fjz += fj3*tm;
                        fkx += fk1*tm;   
                        fky += fk2*tm;   
                        fkz += fk3*tm;                           
                        
//                         nijk3 = 3*i;                        
//                         force[0 + nijk3] += (fj1 + fk1)*tm;
//                         force[1 + nijk3] += (fj2 + fk2)*tm;
//                         force[2 + nijk3] += (fj3 + fk3)*tm;
// 
//                         nijk3 = 3*j;    
//                         force[0 + nijk3] -= fj1*tm;
//                         force[1 + nijk3] -= fj2*tm;
//                         force[2 + nijk3] -= fj3*tm;
// 
//                         nijk3 = 3*k;    
//                         force[0 + nijk3] -= fk1*tm;   
//                         force[1 + nijk3] -= fk2*tm;   
//                         force[2 + nijk3] -= fk3*tm;                           
                    }                    
                }
                nijk3 = 3*k;    
                force[0 + nijk3] -= fkx;   
                force[1 + nijk3] -= fky;   
                force[2 + nijk3] -= fkz;                                           
                fjxtmp += fjx;
                fjytmp += fjy;
                fjztmp += fjz;                
                fixtmp += fjx+fkx;
                fiytmp += fjy+fky;
                fiztmp += fjz+fkz;                                
            }
            nijk3 = 3*j;    
            force[0 + nijk3] -= fjxtmp;
            force[1 + nijk3] -= fjytmp;
            force[2 + nijk3] -= fjztmp;            
            nijk3 = 3*i;                        
            force[0 + nijk3] += fixtmp;
            force[1 + nijk3] += fiytmp;
            force[2 + nijk3] += fiztmp;                                
        }
    }
}

void CPOD::snapTallyForce(double *force, double *dbdr, double *coeff4,
        int *ai, int *aj, int *ti, int ijnum, int ncoeff, int ntype)
{           
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = ncoeff*(itype-1);        
        int nij = ijnum*3*icoeff;
        
        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];      
        double ce = coeff4[icoeff + n];
        
        force[0 + 3*i] += bix*ce; 
        force[1 + 3*i] += biy*ce;
        force[2 + 3*i] += biz*ce;
        force[0 + 3*j] -= bix*ce;
        force[1 + 3*j] -= biy*ce;
        force[2 + 3*j] -= biz*ce;        
    }
}

void CPOD::pod4body_force(double *force, double *rij, double *coeff4, double *tmpmem, int *atomtype, 
        int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)            
{    
    int dim = 3;    
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
            
    int ne = 0;
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *dUr = &tmpmem[ne];
    ne += idxu_max*dim*Nij;
    double *dUi = &tmpmem[ne];
    ne += idxu_max*dim*Nij;    
    double *dblist = &tmpmem[ne]; // idxb_max*ntriples*dim*Nij          
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    this->snapComputeUlist(Ur, Ui, dUr, dUi, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
    
    this->snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, idxi, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    this->snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, idxi, tj, idxu_max, natom, Nij, chemflag);    
        
    this->snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);
            
    this->snapComputeDbidrj(dblist, Zr, Zi, dUr, dUi, idxb, idxu_block, idxz_block, map, idxi, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, natom, Nij);
    
    this->snapTallyForce(force, dblist, coeff4, ai, aj, ti, Nij, ncoeff, ntypes);           
}

void CPOD::calculate_force(double *force, double *effectivecoeff, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)
{        
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nrbf2 = pod.nbf2;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int *pdegree = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    // effective POD coefficients for calculating force 
    double *coeff2 = &effectivecoeff[nd1];  
    double *coeff3 = &effectivecoeff[nd1+nd2];
    double *coeff4 = &effectivecoeff[nd1+nd2+nd3];    
    
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *f2ij = &tmpmem[Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    

    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    this->podMatMul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);
       
    this->pod2body_force(force, f2ij, coeff2, ai, aj, ti, tj, elemindex, nelements, nrbf2, natom, Nij);
    
    this->pod3body_force(force, rij, e2ij, f2ij, coeff3, &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, ai, aj, 
            ti, tj, nrbf3, nabf3, nelements, natom, Nij);
        
    if (pod.snaptwojmax>0) 
        this->pod4body_force(force, rij, coeff4, tmpmem, atomtype, idxi, ai, aj, ti, tj, natom, Nij);                   
}

double CPOD::energyforce_calculation(double *force, double *podcoeff, double *effectivecoeff, double *gd, double *rij, 
        double *tmpmem, int *pairnumsum, int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)
{       
    int nd1234 = pod.nd1+pod.nd2+pod.nd3+pod.nd4;
    double *eatom = &tmpmem[0];
                
    podArraySetValue(gd, 0.0, nd1234);    
    this->linear_descriptors_ij(gd, eatom, rij, &tmpmem[natom*nd1234], pairnumsum, atomtype, idxi, ti, tj, natom, Nij);
    
    // Need to do MPI_Allreduce on gd for paralell
    
    double energy = this->calculate_energy(effectivecoeff, gd, podcoeff);    
        
    podArraySetValue(force, 0.0, 3*natom);    

    this->calculate_force(force, effectivecoeff, rij, tmpmem, pairnumsum, atomtype, idxi, ai, aj, ti, tj, natom, Nij);
    
    return energy;
}


void CPOD::pod2body_force(double **force, double *fij, double *coeff2, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        for (int m=0; m<nbf; m++) {               
            //int im =  3*i1;
            //int jm =  3*j1;
            int nm = n + N*m;
            int km = (elemindex[typei + typej*nelements] - 1) + nelements2*m;
            double ce = coeff2[km];
            force[i1][0] += fij[0 + 3*nm]*ce;
            force[i1][1] += fij[1 + 3*nm]*ce;
            force[i1][2] += fij[2 + 3*nm]*ce;
            force[j1][0] -= fij[0 + 3*nm]*ce;
            force[j1][1] -= fij[1 + 3*nm]*ce;
            force[j1][2] -= fij[2 + 3*nm]*ce;          
        }
    }
}

void CPOD::pod3body_force(double **force, double *yij, double *e2ij, double *f2ij, double *coeff3, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, c, nijk3, typei, typej, typek, ij, ik, i, j, k;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    //double uj, uk, rbf, drbf1, drbf2, drbf3, drbf4, drbf5, drbf6;
    //double fj1, fj2, fj3, fk1, fk2, fk3;
            
    double *abf = &tmpmem[0];
    double *dabf1 = &tmpmem[nabf1];
    double *dabf2 = &tmpmem[2*nabf1];
    double *dabf3 = &tmpmem[3*nabf1];
    double *dabf4 = &tmpmem[4*nabf1];
    double *dabf5 = &tmpmem[5*nabf1];
    double *dabf6 = &tmpmem[6*nabf1];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            ij = lj + s;
            i = ai[ij];  // atom i                        
            j = aj[ij];  // atom j    
            typei = ti[ij] - 1;           
            typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);            
            
            double fixtmp,fiytmp,fiztmp;
            fixtmp = fiytmp = fiztmp = 0.0;      
            double fjxtmp,fjytmp,fjztmp;
            fjxtmp = fjytmp = fjztmp = 0.0;                  
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                ik = lk + s;
                k = aj[ik];  // atom k                       
                typek = tj[ik] - 1;                         
                xik1 = yij[0+dim*ik];  // xk - xi           
                xik2 = yij[1+dim*ik];  // xk - xi           
                xik3 = yij[2+dim*ik];  // xk - xi           s
                riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
                rik = pow(riksq, 0.5); 

                xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
                costhe = xdot/(rij*rik);    
                costhe = costhe > 1.0 ? 1.0 : costhe;
                costhe = costhe < -1.0 ? -1.0 : costhe;
                xdot = costhe*(rij*rik);

                sinthe = pow(1.0 - costhe*costhe, 0.5);
                sinthe = sinthe > 1e-12 ? sinthe : 1e-12;    
                theta = acos(costhe);            
                dtheta = -1.0/sinthe; 

                tm1 = pow(rijsq,1.5)*rik;
                tm2 = rij*pow(riksq,1.5);
                tm1 = 1.0/tm1;
                tm2 = 1.0/tm2;
                dct1 = (xik1*rijsq - xij1*xdot)*tm1; 
                dct2 = (xik2*rijsq - xij2*xdot)*tm1;
                dct3 = (xik3*rijsq - xij3*xdot)*tm1;
                dct4 = (xij1*riksq - xik1*xdot)*tm2;
                dct5 = (xij2*riksq - xik2*xdot)*tm2;
                dct6 = (xij3*riksq - xik3*xdot)*tm2;

                for (int p=0; p <nabf1; p++) {
                    abf[p] = cos(p*theta);                
                    tm = -p*sin(p*theta)*dtheta;                    
                    dabf1[p] = tm*dct1;
                    dabf2[p] = tm*dct2;
                    dabf3[p] = tm*dct3;
                    dabf4[p] = tm*dct4;
                    dabf5[p] = tm*dct5;
                    dabf6[p] = tm*dct6;        
                }

                double fjx = 0.0, fjy = 0.0, fjz = 0.0;
                double fkx = 0.0, fky = 0.0, fkz = 0.0;
                
                for (int m=0; m<nrbf; m++) {
                    double uj = e2ij[lj + s + Nij*m];
                    double uk = e2ij[lk + s + Nij*m];
                    double rbf = uj*uk;
                    double drbf1 = f2ij[0 + dim*(lj + s) + dim*Nij*m]*uk;
                    double drbf2 = f2ij[1 + dim*(lj + s) + dim*Nij*m]*uk;
                    double drbf3 = f2ij[2 + dim*(lj + s) + dim*Nij*m]*uk;                                                        
                    double drbf4 = f2ij[0 + dim*(lk + s) + dim*Nij*m]*uj;
                    double drbf5 = f2ij[1 + dim*(lk + s) + dim*Nij*m]*uj;
                    double drbf6 = f2ij[2 + dim*(lk + s) + dim*Nij*m]*uj;     

                    for (int p=0; p <nabf1; p++) {
                        tm = abf[p];
                        double fj1 = drbf1*tm + rbf*dabf1[p];
                        double fj2 = drbf2*tm + rbf*dabf2[p];
                        double fj3 = drbf3*tm + rbf*dabf3[p];
                        double fk1 = drbf4*tm + rbf*dabf4[p];
                        double fk2 = drbf5*tm + rbf*dabf5[p];
                        double fk3 = drbf6*tm + rbf*dabf6[p];

                        n = p + (nabf1)*m;
                        c = (elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n;
                        tm = coeff3[c];
                        
                        fjx += fj1*tm;
                        fjy += fj2*tm;
                        fjz += fj3*tm;
                        fkx += fk1*tm;   
                        fky += fk2*tm;   
                        fkz += fk3*tm;                           
                        
//                         nijk3 = 3*i;                        
//                         force[0 + nijk3] += (fj1 + fk1)*tm;
//                         force[1 + nijk3] += (fj2 + fk2)*tm;
//                         force[2 + nijk3] += (fj3 + fk3)*tm;
// 
//                         nijk3 = 3*j;    
//                         force[0 + nijk3] -= fj1*tm;
//                         force[1 + nijk3] -= fj2*tm;
//                         force[2 + nijk3] -= fj3*tm;
// 
//                         nijk3 = 3*k;    
//                         force[0 + nijk3] -= fk1*tm;   
//                         force[1 + nijk3] -= fk2*tm;   
//                         force[2 + nijk3] -= fk3*tm;                           
                    }                    
                }
                nijk3 = k;    
                force[nijk3][0] -= fkx;   
                force[nijk3][1] -= fky;   
                force[nijk3][2] -= fkz;                                           
                fjxtmp += fjx;
                fjytmp += fjy;
                fjztmp += fjz;                
                fixtmp += fjx+fkx;
                fiytmp += fjy+fky;
                fiztmp += fjz+fkz;                                
            }
            nijk3 = j;    
            force[nijk3][0] -= fjxtmp;
            force[nijk3][1] -= fjytmp;
            force[nijk3][2] -= fjztmp;            
            nijk3 = i;                        
            force[nijk3][0] += fixtmp;
            force[nijk3][1] += fiytmp;
            force[nijk3][2] += fiztmp;                                
        }
    }
}

void CPOD::snapTallyForce(double **force, double *dbdr, double *coeff4,
        int *ai, int *aj, int *ti, int ijnum, int ncoeff, int ntype)
{           
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = ncoeff*(itype-1);        
        int nij = ijnum*3*icoeff;
        
        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];      
        double ce = coeff4[icoeff + n];
        
        force[i][0] += bix*ce; 
        force[i][1] += biy*ce;
        force[i][2] += biz*ce;
        force[j][0] -= bix*ce;
        force[j][1] -= biy*ce;
        force[j][2] -= biz*ce;        
    }
}

void CPOD::pod4body_force(double **force, double *rij, double *coeff4, double *tmpmem, int *atomtype, 
        int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)            
{    
    int dim = 3;    
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
            
    int ne = 0;
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *dUr = &tmpmem[ne];
    ne += idxu_max*dim*Nij;
    double *dUi = &tmpmem[ne];
    ne += idxu_max*dim*Nij;    
    double *dblist = &tmpmem[ne]; // idxb_max*ntriples*dim*Nij          
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    this->snapComputeUlist(Ur, Ui, dUr, dUi, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
    
    this->snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, idxi, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    this->snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, idxi, tj, idxu_max, natom, Nij, chemflag);    
        
    this->snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);
            
    this->snapComputeDbidrj(dblist, Zr, Zi, dUr, dUi, idxb, idxu_block, idxz_block, map, idxi, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, natom, Nij);
    
    this->snapTallyForce(force, dblist, coeff4, ai, aj, ti, Nij, ncoeff, ntypes);           
}

void CPOD::calculate_force(double **force, double *effectivecoeff, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)
{        
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nrbf2 = pod.nbf2;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int *pdegree = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    // effective POD coefficients for calculating force 
    double *coeff2 = &effectivecoeff[nd1];  
    double *coeff3 = &effectivecoeff[nd1+nd2];
    double *coeff4 = &effectivecoeff[nd1+nd2+nd3];    
    
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *f2ij = &tmpmem[Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    

    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    this->podMatMul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);
       
    this->pod2body_force(force, f2ij, coeff2, ai, aj, ti, tj, elemindex, nelements, nrbf2, natom, Nij);
    
    this->pod3body_force(force, rij, e2ij, f2ij, coeff3, &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, ai, aj, 
            ti, tj, nrbf3, nabf3, nelements, natom, Nij);
        
    if (pod.snaptwojmax>0) 
        this->pod4body_force(force, rij, coeff4, tmpmem, atomtype, idxi, ai, aj, ti, tj, natom, Nij);                   
}


