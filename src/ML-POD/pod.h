/***************************************************************************                                                   
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __POD_H__
#define __POD_H__

#include "podcommon.h"
#include "pointers.h"

namespace LAMMPS_NS {
    
class CPOD : protected Pointers {    

private:
    // ***********************  implemented in podinputfiles.cpp **************************/
    void read_pod(std::string pod_file);
    
    void read_coeff_file(std::string coeff_file);
    // ******************************************************************************/
    
    // ***********************  implemented in poddescriptors.cpp **************************/
    void podradialbasis(double *rbf, double *drbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N);
               
    void pod1body(double *eatom, double *fatom, int *atomtype, int nelements, int natom);
            
    void podtally2b(double *eatom, double *fatom, double *eij, double *fij, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N);
                    
    void pod3body(double *eatom, double *fatom, double *rij, double *e2ij, double *f2ij, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij);
                
    void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
            double *fatom3, double *rij, double *Phi, double *besselparams, double *tmpmem, double rin, 
            double rcut, int *pairnumsum, int *atomtype, int *ai, int *aj, int *ti, int *tj, int *elemindex, 
            int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom);
    
    double quadratic_coefficients(double *c2, double *c3, double *d2, double *d3, 
        double *coeff23, int *quadratic, int nc2, int nc3);
    
    double quadratic_coefficients(double *c3, double *d3, double *coeff33, 
        int *quadratic, int nc3);

    double cubic_coefficients(double *c2, double *c3, double *c4, double *d2, double *d3, double *d4, 
        double *coeff234, int *cubic, int nc2, int nc3, int nc4);

    double cubic_coefficients(double *c3, double *d3, double *coeff333, int *cubic, int nc3);   
    
    double quadratic_coefficients(double *ce2, double *ce3, double *c2, double *c3, double *d2, double *d3, 
        double *coeff23, int *quadratic, int nc2, int nc3);
    
    double quadratic_coefficients(double *ce3, double *c3, double *d3, double *coeff33, 
        int *quadratic, int nc3);

    double cubic_coefficients(double *ce2, double *ce3, double *ce4, double *c2, double *c3, double *c4, 
            double *d2, double *d3, double *d4, double *coeff234, int *cubic, int nc2, int nc3, int nc4);

    double cubic_coefficients(double *ce3, double *c3, double *d3, double *coeff333, int *cubic, int nc3);           
    // ******************************************************************************/
    
    // ***********************  implemented in podsnap.cpp **************************/
    void snapSetup(int twojmax, int ntypes);
    
    void InitSnap();
    
    void snapComputeUlist(double *Sr, double *Si, double *dSr, double *dSi, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);
    
    void snapZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum);

    void snapAddUarraytot(double *Stotr, double *Stoti, double *Sr, 
        double *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag);
    
    void snapComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum);
    
    void snapComputeBi1(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        int *idxb, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int inum);
 
    void snapComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, 
        double *dulist_r, double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum);
    
    void snapdesc(double *blist, double *bd, double *rij, double *tmpmem, int *atomtype, int *ai, 
        int *aj, int *ti, int *tj, int natom, int Nij);           
    // ******************************************************************************/
    
    // ***********************  implemented in podenergyforce.cpp **************************/
    void podradialbasis(double *rbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N);
    
    void pod1body(double *eatom, int *atomtype, int nelements, int natom);
                        
    void podtally2b(double *eatom, double *eij, int *ai, int *ti, int *tj, int *elemindex, 
        int nelements, int nbf, int natom, int N);

    void pod3body(double *eatom, double *yij, double *e2ij, double *tmpmem, int *elemindex, int *pairnumsum, 
        int *ai, int *ti, int *tj, int nrbf, int nabf, int nelements, int natom, int Nij);
        
    void poddesc_ij(double *eatom1, double *eatom2, double *eatom3, double *rij, double *Phi, double *besselparams, 
                double *tmpmem, double rin, double rcut, int *pairnumsum, int *atomtype, int *ai, int *ti, int *tj, 
                int *elemindex, int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom);        
    
    void snapComputeUij(double *Sr, double *Si, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);
    
    void snapdesc_ij(double *blist, double *rij, double *tmpmem, int *atomtype, int *ai, 
        int *ti, int *tj, int natom, int Nij);        
    
    void pod2body_force(double *force, double *fij, double *coeff2, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int Nij);
    
    void pod3body_force(double *force, double *yij, double *e2ij, double *f2ij, double *coeff3, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij);
        
    void snapTallyForce(double *force, double *dbdr, double *coeff4, int *ai, int *aj, int *ti, int ijnum, 
            int ncoeff, int ntype);
    
    void pod4body_force(double *force, double *rij, double *coeff4, double *tmpmem, int *atomtype, 
        int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);           
    
    void pod2body_force(double **force, double *fij, double *coeff2, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int Nij);
    
    void pod3body_force(double **force, double *yij, double *e2ij, double *f2ij, double *coeff3, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij);
        
    void snapTallyForce(double **force, double *dbdr, double *coeff4, int *ai, int *aj, int *ti, int ijnum, 
            int ncoeff, int ntype);
    
    void pod4body_force(double **force, double *rij, double *coeff4, double *tmpmem, int *atomtype, 
        int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);                        
    // ******************************************************************************/    
    
public:          
    struct podstruct {     
        std::vector<std::string> species;    
        int *pbc=NULL; //[3] = {1,1,1};
        int *elemindex=NULL;

        int nelements = 0;
        int onebody = 1;
        int besseldegree = 3;
        int inversedegree = 6;
        int twobody[3] = {5,10,10};
        int threebody[4] = {4,8,8,5}; 
        int fourbody[4] = {0,0,0,0};    

        int quadraticpod = 0;
        int quadratic22[2] = {0,0};
        int quadratic23[2] = {0,0};
        int quadratic24[2] = {0,0};
        int quadratic33[2] = {0,0};
        int quadratic34[2] = {0,0};
        int quadratic44[2] = {0,0};        
        int cubic234[3] = {0,0,0};
        int cubic333[3] = {0,0,0};
        int cubic444[3] = {0,0,0};

        double rin = 0.5;
        double rcut = 4.6;
        double *besselparams=NULL; //[3] = {0.0, 2.0, 4.0};        
        double *Phi2=NULL, *Phi3=NULL, *Phi4=NULL, *Lambda2=NULL, *Lambda3=NULL, *Lambda4=NULL;    
        double *coeff=NULL;

        int nbesselpars = 3;    
        int ns2, ns3, ns4;       // number of snapshots for radial basis functions for linear POD potentials      
        int nc2, nc3, nc4;       // number of chemical  combinations for linear POD potentials      
        int nbf1, nbf2, nbf3, nbf4; // number of basis functions for linear POD potentials      
        int nd1, nd2, nd3, nd4;     // number of descriptors for linear POD potentials 
        int nd22, nd23, nd24, nd33, nd34, nd44; // number of descriptors for quadratic POD potentials    
        int nd234, nd333, nd444; // number of descriptors for cubic POD potentials    
        int nrbf3, nabf3, nrbf4, nabf4;    
        int nd, nd1234;

        int snaptwojmax = 0;
        int snapchemflag = 0;
        double snaprfac0 = 0.99363;
        double snapelementradius[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        double snapelementweight[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        void allocatememory(int backend)
        {
            TemplateMalloc(&pbc, 3, backend);
            TemplateMalloc(&besselparams, 3, backend);
        }    

        void freememory(int backend)
        {
            TemplateFree(pbc, backend);    
            TemplateFree(elemindex, backend);    
            TemplateFree(besselparams, backend);        
            TemplateFree(Phi2, backend);        
            TemplateFree(Phi3, backend);        
            TemplateFree(Phi4, backend);        
            TemplateFree(Lambda2, backend);        
            TemplateFree(Lambda3, backend);        
            TemplateFree(Lambda4, backend);    
            TemplateFree(coeff, backend);    
        }        
    };

    struct snastruct {        
        int twojmax;
        int ncoeff;
        int idxb_max;
        int idxu_max;
        int idxz_max;
        int idxcg_max;
        int ntypes;
        int nelements;    
        int ndoubles;   // number of multi-element pairs
        int ntriples;   // number of multi-element triplets      
        int bnormflag;
        int chemflag;    
        int switchflag;
        int bzeroflag;
        int wselfallflag;

        double wself;
        double rmin0;
        double rfac0;
        double rcutfac;
        double rcutmax;    

        int *map=NULL;  // map types to [0,nelements)    
        int *idx_max=NULL; 
        int *idxz=NULL;
        int *idxz_block=NULL;
        int *idxb=NULL;
        int *idxb_block=NULL;
        int *idxu_block=NULL;
        int *idxcg_block=NULL;

        double *rcutsq=NULL;    
        double *radelem=NULL;
        double *wjelem=NULL; 
        double *bzero=NULL;
        double *fac=NULL;
        double *rootpqarray=NULL; 
        double *cglist=NULL;

        void printout()
        {
            printf("twojmax %d \n", twojmax); 
            printf("ncoeff %d \n", ncoeff);         
            printf("idxb_max %d \n", idxb_max);         
            printf("idxu_max %d \n", idxu_max);         
            printf("idxz_max %d \n", idxz_max); 
            printf("idxcg_max %d \n", idxcg_max);
            printf("ntypes %d \n", ntypes);
            printf("nelements %d \n", nelements);
            printf("ndoubles %d \n", ndoubles);
            printf("ntriples %d \n", ntriples);
            printf("bnormflag %d \n", bnormflag);
            printf("chemflag %d \n", chemflag);
            printf("switchflag %d \n", switchflag);
            printf("bzeroflag %d \n", bzeroflag);
            printf("wselfallflag %d \n", wselfallflag);        
            printf("rfac0 %g \n", rfac0);
            printf("rmin0 %g \n", rmin0);
            printf("rcutfac %g \n", rcutfac);
            printf("rcutmax %g \n", rcutmax);    
        }

        void freememory(int backend)
        {   
            TemplateFree(map, backend);
            TemplateFree(idx_max, backend);
            TemplateFree(idxz, backend);
            TemplateFree(idxb, backend);
            TemplateFree(idxb_block, backend);
            TemplateFree(idxu_block, backend);
            TemplateFree(idxz_block, backend);
            TemplateFree(idxcg_block, backend);

            TemplateFree(rootpqarray, backend);
            TemplateFree(cglist, backend);
            TemplateFree(fac, backend);
            TemplateFree(bzero, backend);
            TemplateFree(wjelem, backend);
            TemplateFree(radelem, backend);
            TemplateFree(rcutsq, backend);
        }                         
    };

    podstruct pod;
    snastruct sna;
    
    // constructor 
    CPOD(LAMMPS *, std::string pod_file, std::string coeff_file); 
    
    CPOD(LAMMPS *lmp) : Pointers(lmp){};
            
    // destructor        
    ~CPOD() override;             
            
    // ***********************  implemented in podarrayfunctions.cpp **************************/        
    void print_matrix(const char* desc, int m, int n, int* a, int lda ); 

    void print_matrix(const char* desc, int m, int n, double* a, int lda ); 
    
    void print_matrix(const char* desc, int m, int n, double **a, int lda ); 
    
    void podMatMul(double *c, double *a, double *b, int r1, int c1, int c2);    

    void podCumsum(int* output, int* input, int length);
    
    double podArrayNorm(double *a, int n);
    
    double podArrayErrorNorm(double *a, double *b, int n);
    
    void podArraySetValue(double *y, double a, int n);

    void podArrayCopy(double *y, double *x, int n);    
    
    void podArrayFill(int* output, int start, int length); 

    double podArrayMin(double *a, int n);

    double podArrayMax(double *a, int n);
    
    double podArraySum(double *a, int n);
    
    int podArrayMin(int *a, int n);

    int podArrayMax(int *a, int n);

    void podKron(double *C, double *A, double *B, double alpha, int M1, int M2);        
    
    void rotation_matrix(double *Rmat, double alpha, double beta, double gamma);    
    void triclinic_lattice_conversion(double *a, double *b, double *c, double *A, double *B, double *C);
    void matrix33_multiplication(double *xrot, double *Rmat, double *x, int natom);
    void matrix33_inverse(double *invA, double *A1, double *A2, double *A3);    
    // ******************************************************************************/
    
    // ***********************  implemented in poddescriptors.cpp **************************/
    void podNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *atomtype, int *alist, int inum, int dim);
        
    void linear_descriptors(double *gd, double *efatom, double *y, double *tmpmem, int *atomtype, 
            int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij);            

    void quadratic_descriptors(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
        int M2, int M3, int N);
    
    void quadratic_descriptors(double* d33, double *dd33, double *d3, double *dd3, int M3, int N);

    void cubic_descriptors(double* d234, double *dd234, double* d2, double *d3, double *d4, 
        double* dd2, double *dd3, double *dd4, int M2, int M3, int M4, int N);
    
    void cubic_descriptors(double* d333, double *Dd333, double *d3, double *Dd3, int M3, int N);

    double calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, double *tmp, int natom);
    
    double energyforce_calculation(double *f, double *gd, double *gdd, double *coeff, double *y, int *atomtype, 
            int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij);         
    // ******************************************************************************/
        
    // ***********************  implemented in podenergyforce.cpp **************************/
    void podNeighPairs(double *rij, double *x, int *idxi, int *ai, int *aj, int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *jlist, int *alist, int inum);
        
    //podptr->podNeighPairs(rij, nb.y, idxi, ai, aj, ti, tj, nb.pairnum_cumsum, atomtype, nb.pairlist, nb.alist, natom);
    
    int lammpsNeighPairs(double *rij, double **x, double rcutsq, int *idxi, int *ai, int *aj,  int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *numneigh, int *ilist, int **jlist, int inum);
        
    void linear_descriptors_ij(double *gd, double *eatom, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *ai, int *ti, int *tj, int natom, int Nij);
    
    double calculate_energy(double *effectivecoeff, double *gd, double *coeff);
        
    double calculate_energy(double *energycoeff, double *forcecoeff, double *gd, double *coeff);
    
    void calculate_force(double *force, double *effectivecoeff, double *rij, double *tmpmem, int *pairnumsum,
            int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);    

    void calculate_force(double **force, double *effectivecoeff, double *rij, double *tmpmem, int *pairnumsum,
            int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);    
    
    double energyforce_calculation(double *force, double *podcoeff, double *effectivecoeff, double *gd, double *rij, 
        double *tmpmem, int *pairnumsum, int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);            
    // ******************************************************************************/

};

}    // namespace LAMMPS_NS

#endif        

