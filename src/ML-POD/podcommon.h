#ifndef PODCOMMON_H
#define PODCOMMON_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGETRF dgetrf_
#define DGETRI dgetri_
#define DSYEV dsyev_
#define DPOSV dposv_

#define PODMIN(a,b) ((a) < (b) ? (a) : (b))
#define PODMAX(a,b) ((a) > (b) ? (a) : (b))

using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != NULL) {                                                      \
        free(x);                                                          \
        x = NULL;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(int*,double*,int*);
    double DDOT(int*,double*,int*,double*,int*);
    void DAXPY(int*,double*,double*,int*,double*,int*);
    void DGEMV(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);  
    void DGEMM(char*,char*,int*,int*,int*,double*,double*,int*,
             double*,int*,double*,double*,int*);        
    void DGETRF(int*,int*,double*,int*,int*,int*);
    void DGETRI(int*,double*,int*,int*,double*,int*,int*);
    void DTRSM(char *, char*, char*, char *, int *, int *, double*, double*, int*,
             double*, int*);
    
    void DSYEV( char* jobz, char* uplo, int* n, double* a, int* lda,
        double* w, double* work, int* lwork, int* info );   
    
    void DPOSV( char* uplo, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, int* info );
}
template <typename T> static void TemplateMalloc(T **data, int n, int backend)
{
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));     
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));    
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CUDA_CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

template <typename T> static void TemplateCopytoDevice(T *d_data, T *h_data, int n, int backend)
{
    if (backend == 0)       
        //podArrayCopy(d_data, h_data, n);
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoDevice(d_data, h_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoDevice(d_data, h_data, n);
#endif                      
}

template <typename T> static void TemplateCopytoHost(T *h_data, T *d_data, int n, int backend)
{
    if (backend == 0)       
        //podArrayCopy(h_data, d_data, n);
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoHost(h_data, d_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoHost(h_data, d_data, n);
#endif                      
}

template <typename T> static void TemplateFree(T *data, int backend)
{
    if (backend == 0)       
        CPUFREE(data);
    if (backend == 1)         
        CPUFREE(data);
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        GPUFREE(data);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP          
        HIPFREE(data);
#endif                      
}

// template <typename T> void readarrayfromfile(const char* filename, T **a, int N)
// {
//     if (N>0) {
//         // Open file to read
//         ifstream in(filename, ios::in | ios::binary);
// 
//         if (!in) {
//             pod_error("Unable to open file ");
//         }
// 
//         if (in) {
//             *a = (T*) malloc (sizeof (T)*N);
//             in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
//         }
// 
//         in.close();
//     }
// }
// 
// template <typename T> void writearray2file(const char* filename, T *a, int N, int backend)
// {
//     if (N>0) {        
//         // Open file to read
//         ofstream out(filename, ios::out | ios::binary);
// 
//         if (!out) {
//             pod_error("Unable to open file ");
//         }
// 
//         if (backend==2) { //GPU
// #ifdef  USE_CUDA                        
//             T *a_host;            
//             a_host = (T*) malloc (sizeof (T)*N);            
//             
//             // transfer data from GPU to CPU to save in a file
//             cudaMemcpy(&a_host[0], &a[0], N*sizeof(T), cudaMemcpyDeviceToHost);    
//             
//             out.write( reinterpret_cast<char*>( &a_host[0] ), sizeof(T) * N );
//             
//             free(a_host);
// #endif            
//         }
//         else 
//             out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );                    
//         
//         out.close();
//     }
// }


// void print_matrix(const char* desc, int m, int n, double* a, int lda ) 
// {
//     int i, j;
//     printf( "\n %s\n", desc );
// 
//     for( i = 0; i < m; i++ ) 
//     {
//         for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
//         printf( "\n" );
//     }
// }
// 
// void print_matrix(const char* desc, int m, int n, int* a, int lda ) 
// {
//     int i, j;
//     printf( "\n %s\n", desc );
// 
//     for( i = 0; i < m; i++ ) 
//     {
//         for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
//         printf( "\n" );
//     }
// }

// void podKron(double *C, double *A, double *B, double alpha, int M1, int M2)
// {            
//     int M = M1*M2;        
//     for (int idx=0; idx<M; idx++)     
//     {
//         int ib = idx%M2;
//         int ia = (idx-ib)/M2;        
//         C[idx] += alpha*A[ia]*B[ib];        
//     }
// }
// 
// void podCumsum(int* output, int* input, int length) 
// {
// 	output[0] = 0; 
// 	for (int j = 1; j < length; ++j)	
// 		output[j] = input[j - 1] + output[j - 1];	
// }
// 
// void podArrayFill(int* output, int start, int length) 
// {	
// 	for (int j = 0; j < length; ++j)	
// 		output[j] = start + j;
// }
// 
// double podArrayErrorNorm(double *a, double *b, int n)
// {
//     double e = (a[0]-b[0])*(a[0]-b[0]);
//     for (int i=1; i<n; i++)        
//         e += (a[i]-b[i])*(a[i]-b[i]);    
//     return sqrt(e);
// }
// 
// double podArrayMaxError(double *a, double *b, int n)
// {
//     double e = fabs(a[0]-b[0]);
//     for (int i=1; i<n; i++)
//         if (fabs(a[i]-b[i])>e)
//             e = fabs(a[i]-b[i]);    
//     return e;
// }
// 
// double podArrayNorm(double *a, int n)
// {
//     double e = a[0]*a[0];
//     for (int i=1; i<n; i++)        
//         e += a[i]*a[i];    
//     return sqrt(e);
// }
// 
// double podArrayMin(double *a, int n)
// {
//     double b = a[0];
//     for (int i=1; i<n; i++)
//         if (a[i]<b)
//             b = a[i];    
//     return b;
// }
// 
// double podArrayMax(double *a, int n)
// {
//     double b = a[0];
//     for (int i=1; i<n; i++)
//         if (a[i]>b)
//             b = a[i];    
//     return b;
// }
// 
// int podArrayMin(int *a, int n)
// {
//     int b = a[0];
//     for (int i=1; i<n; i++)
//         if (a[i]<b)
//             b = a[i];    
//     return b;
// }
// 
// int podArrayMax(int *a, int n)
// {
//     int b = a[0];
//     for (int i=1; i<n; i++)
//         if (a[i]>b)
//             b = a[i];    
//     return b;
// }
// 
// void podArraySetValue(double *y, double a, int n)
// {    
//     for (int i=0; i<n; i++) 
//         y[i] = a;        
// }
// 
// void podArrayCopy(double *y, double *x, int n)
// {    
//     for (int i=0; i<n; i++) 
//         y[i] = x[i];        
// }

// struct podstruct {     
//     std::vector<std::string> species;    
//     int *pbc=NULL; //[3] = {1,1,1};
//     int *elemindex=NULL;
//     
//     int nelements = 0;
//     int onebody = 1;
//     int besseldegree = 3;
//     int inversedegree = 6;
//     int twobody[3] = {5,10,10};
//     int threebody[4] = {4,8,8,5}; 
//     int fourbody[4] = {0,0,0,0};    
//     
//     int quadratic22[2] = {0,0};
//     int quadratic23[2] = {0,0};
//     int quadratic24[2] = {0,0};
//     int quadratic33[2] = {0,0};
//     int quadratic34[2] = {0,0};
//     int quadratic44[2] = {0,0};        
//     int cubic234[3] = {0,0,0};
//     int cubic333[3] = {0,0,0};
//     int cubic444[3] = {0,0,0};
//     
//     double rin = 0.5;
//     double rcut = 4.6;
//     double *besselparams=NULL; //[3] = {0.0, 2.0, 4.0};        
//     double *Phi2=NULL, *Phi3=NULL, *Phi4=NULL, *Lambda2=NULL, *Lambda3=NULL, *Lambda4=NULL;    
//     double *coeff=NULL;
//         
//     int nbesselpars = 3;    
//     int ns2, ns3, ns4;       // number of snapshots for radial basis functions for linear POD potentials      
//     int nc2, nc3, nc4;       // number of chemical  combinations for linear POD potentials      
//     int nbf1, nbf2, nbf3, nbf4; // number of basis functions for linear POD potentials      
//     int nd1, nd2, nd3, nd4;     // number of descriptors for linear POD potentials 
//     int nd22, nd23, nd24, nd33, nd34, nd44; // number of descriptors for quadratic POD potentials    
//     int nd234, nd333, nd444; // number of descriptors for cubic POD potentials    
//     int nrbf3, nabf3, nrbf4, nabf4;    
//     int nd;
//     
//     int snaptwojmax = 0;
//     int snapchemflag = 0;
//     double snaprfac0 = 0.99363;
//     double snapelementradius[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
//     double snapelementweight[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//     
//     void allocatememory(int backend)
//     {
//         TemplateMalloc(&pbc, 3, backend);
//         TemplateMalloc(&besselparams, 3, backend);
//     }    
//     
//     void freememory(int backend)
//     {
//         TemplateFree(pbc, backend);    
//         TemplateFree(elemindex, backend);    
//         TemplateFree(besselparams, backend);        
//         TemplateFree(Phi2, backend);        
//         TemplateFree(Phi3, backend);        
//         TemplateFree(Phi4, backend);        
//         TemplateFree(Lambda2, backend);        
//         TemplateFree(Lambda3, backend);        
//         TemplateFree(Lambda4, backend);    
//         TemplateFree(coeff, backend);    
//     }        
// };
// 
// struct snastruct {        
//     int twojmax;
//     int ncoeff;
//     int idxb_max;
//     int idxu_max;
//     int idxz_max;
//     int idxcg_max;
//     int ntypes;
//     int nelements;    
//     int ndoubles;   // number of multi-element pairs
//     int ntriples;   // number of multi-element triplets      
//     int bnormflag;
//     int chemflag;    
//     int switchflag;
//     int bzeroflag;
//     int wselfallflag;
//     
//     double wself;
//     double rmin0;
//     double rfac0;
//     double rcutfac;
//     double rcutmax;    
//         
//     int *map=NULL;  // map types to [0,nelements)    
//     int *idx_max=NULL; 
//     int *idxz=NULL;
//     int *idxz_block=NULL;
//     int *idxb=NULL;
//     int *idxb_block=NULL;
//     int *idxu_block=NULL;
//     int *idxcg_block=NULL;
//     
//     double *rcutsq=NULL;    
//     double *radelem=NULL;
//     double *wjelem=NULL; 
//     double *bzero=NULL;
//     double *fac=NULL;
//     double *rootpqarray=NULL; 
//     double *cglist=NULL;
//     
//     void printout()
//     {
//         printf("twojmax %d \n", twojmax); 
//         printf("ncoeff %d \n", ncoeff);         
//         printf("idxb_max %d \n", idxb_max);         
//         printf("idxu_max %d \n", idxu_max);         
//         printf("idxz_max %d \n", idxz_max); 
//         printf("idxcg_max %d \n", idxcg_max);
//         printf("ntypes %d \n", ntypes);
//         printf("nelements %d \n", nelements);
//         printf("ndoubles %d \n", ndoubles);
//         printf("ntriples %d \n", ntriples);
//         printf("bnormflag %d \n", bnormflag);
//         printf("chemflag %d \n", chemflag);
//         printf("switchflag %d \n", switchflag);
//         printf("bzeroflag %d \n", bzeroflag);
//         printf("wselfallflag %d \n", wselfallflag);        
//         printf("rfac0 %g \n", rfac0);
//         printf("rmin0 %g \n", rmin0);
//         printf("rcutfac %g \n", rcutfac);
//         printf("rcutmax %g \n", rcutmax);    
//         print_matrix( "map:", 1, ntypes+1, map, 1); 
//         print_matrix( "radelem:", 1, ntypes+1, radelem, 1); 
//         print_matrix( "wjelem:", 1, ntypes+1, wjelem, 1); 
//         print_matrix( "rcutsq:", ntypes+1, ntypes+1, rcutsq, ntypes+1); 
//         print_matrix( "bzero:", 1, twojmax+1, bzero, 1);
//         print_matrix( "fac:", 1, 20, fac, 1);
//         print_matrix( "rootpqarray:", twojmax+1, twojmax+1, rootpqarray, (twojmax+1));
//         print_matrix( "cglist:", 1, idxcg_max, cglist, 1);            
//     }
//     
//     void freememory(int backend)
//     {   
//         TemplateFree(map, backend);
//         TemplateFree(idx_max, backend);
//         TemplateFree(idxz, backend);
//         TemplateFree(idxb, backend);
//         TemplateFree(idxb_block, backend);
//         TemplateFree(idxu_block, backend);
//         TemplateFree(idxz_block, backend);
//         TemplateFree(idxcg_block, backend);
//         
//         TemplateFree(rootpqarray, backend);
//         TemplateFree(cglist, backend);
//         TemplateFree(fac, backend);
//         TemplateFree(bzero, backend);
//         TemplateFree(wjelem, backend);
//         TemplateFree(radelem, backend);
//         TemplateFree(rcutsq, backend);
//     }                         
// };
// 
// struct datastruct {     
//     std::string file_format;
//     std::string file_extension;    
//     std::string data_path;    
//     std::vector<std::string> data_files;     
//     std::vector<std::string> filenames;     
//     
//     std::vector<int> num_atom;
//     std::vector<int> num_atom_cumsum;
//     std::vector<int> num_atom_each_file;
//     std::vector<int> num_config;
//     std::vector<int> num_config_cumsum;
//     int num_atom_sum; 
//     int num_atom_min; 
//     int num_atom_max; 
//     int num_config_sum;    
//     
//     double *lattice=NULL;
//     double *energy=NULL; 
//     double *stress=NULL;
//     double *position=NULL;
//     double *force=NULL;
//     int *atomtype=NULL;
//     
//     int training = 1;
//     int normalizeenergy = 1;
//     int training_analysis = 1;
//     int test_analysis = 1;
//     int training_calculation = 0;
//     int test_calculation = 0;
//     int randomize = 1;    
//     double percentage = 1.0;    
//         
//     double fitting_weights[12] = {0.0, 0.0, 0.0, 1, 1, 0, 0, 1, 1, 1, 1, 0};
//     
//     void copydatainfo(datastruct &data) {
//         data.data_path = data_path;        
//         data.file_format = file_format;
//         data.file_extension = file_extension;              
//         data.data_files = data_files;
//         data.filenames = filenames;
//         data.training_analysis = training_analysis;
//         data.test_analysis = test_analysis;
//         data.training_calculation = training_calculation;
//         data.test_calculation = test_calculation;
//         data.percentage = percentage;
//         data.randomize = randomize;
//         data.training = training;        
//         data.normalizeenergy = normalizeenergy;
//         for (int i = 0; i < 12; i++)
//             data.fitting_weights[i] = fitting_weights[i];
//     }                
//     
//     void freememory(int backend)
//     {
//         TemplateFree(lattice, backend);        
//         TemplateFree(energy, backend);        
//         TemplateFree(stress, backend);        
//         TemplateFree(position, backend);        
//         TemplateFree(force, backend);        
//         TemplateFree(atomtype, backend);        
//     }            
// };
// 
// // struct tempstruct {     
// //     double *tmpmem;
// //     int *tmpint;        
// //     int szd;
// //     int szi;
// // };
// 
// struct neighborstruct {
//     int *alist=NULL;
//     int *pairnum=NULL;
//     int *pairnum_cumsum=NULL;
//     int *pairlist=NULL;
//     double *y=NULL;        
//     
//     int natom;    
//     int nalist;    
//     int natom_max;
//     int sze;
//     int sza;
//     int szy;    
//     int szp;    
//     
//     void freememory(int backend)
//     {
//         TemplateFree(alist, backend);        
//         TemplateFree(pairnum, backend);        
//         TemplateFree(pairnum_cumsum, backend);        
//         TemplateFree(pairlist, backend);        
//         TemplateFree(y, backend);        
//     }                
// };
//         
// struct descriptorstruct {
//     double *gd=NULL;  // global descriptors    
//     double *gdd=NULL; // derivatives of global descriptors and peratom descriptors
//     double *A=NULL;  // least-square matrix for all descriptors    
//     double *b=NULL;  // least-square vector for all descriptors    
//     double *c=NULL;  // coefficents of descriptors
//     int *tmpint=NULL;        
//     int szd;
//     int szi;
//     
//     void freememory(int backend)
//     {
//         TemplateFree(gd, backend);     
//         TemplateFree(gdd, backend);     
//         TemplateFree(A, backend);        
//         TemplateFree(b, backend);       
//         TemplateFree(c, backend);       
//         TemplateFree(tmpint, backend);        
//     }                
// };

#endif  
