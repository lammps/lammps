#ifndef CPUPOD
#define CPUPOD

#include "cpuPOD.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int neighborlist(int *ai, int *aj, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
{
    int k = 0;
    for (int i = 0; i<nx; i++) {
        double *ri = &r[i*dim];
        int inc = 0;
        for (int j=0; j<N; j++) {
            double *rj = &r[dim*j];                        
            double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
            if  ((rijsq > 1e-12) && (rijsq <= rcutsq))  { 
                inc += 1;                                
                ai[k] = i+1;
                aj[k] = j+1;          
                k += 1;                                                  
            }
        }
        numneigh[i] = inc; 
    }
    return k; 
}

void cpuCumsum(int* output, int* input, int length) 
{
	output[0] = 0; 
	for (int j = 1; j < length; ++j)	
		output[j] = input[j - 1] + output[j - 1];	
}

void podhybrid23(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
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

int cpuNeighPairList(int *pairnum, int *pairnumsum, int *pairlist, double *x, double rcutsq, 
        int *neighlist, int *neighnumsum, int inum,  int dim)
{    
    int sumcount = 0;
    for (int ii=0; ii<inum; ii++) {
        int i = ii;       // atom i
        int n1 = neighnumsum[i];    
        int m = neighnumsum[i+1] - n1;
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[n1 + l];         
            // distance between atom i and atom j                                    
            double xij0 = x[j*dim] - x[i*dim];  // xj - xi
            double xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            double xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            double dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij < rcutsq && dij>1e-20) {
                pairlist[count + sumcount] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
        sumcount += count; 
    }    

    cpuCumsum(pairnumsum, pairnum, inum+1); 

    return sumcount; 
};

void cpuNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
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

void cpuTripletnum(int* tripletnum, int* pairnum, int length) 
{	
    // pairnum = Int32.(pairnumsum[2:end]-pairnumsum[1:end-1]);
    // tripletnum  = Int32.((pairnum.-1).*pairnum/2);
	for (int ii = 0; ii < length; ++ii)	
		tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       
}

void cpuNeighTripletList(int *tripletlist, int *tripletnum, int *tripletnumsum, int *pairlist, 
     int *pairnum, int *pairnumsum, int inum)
{        
	for (int ii = 0; ii < inum; ++ii)	
		tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       

    cpuCumsum(tripletnumsum, tripletnum, inum+1); 

    int Nijk = tripletnumsum[inum];

    for (int ii=0; ii<inum; ii++) {
        //int i = ii;       // atom i
        int p = pairnum[ii];      // number of pairs (i,j) around i         
        int q = tripletnumsum[ii];
        int s = pairnumsum[ii];
        int count = 0;
        for (int lj=0; lj<p ; lj++) {   // loop over each atom j around atom i
            int gj = pairlist[lj + s];  // atom j           
            for (int lk=lj+1; lk<p; lk++) { // loop over each atom k around atom i (k > j)
                int gk = pairlist[lk + s];  // atom k
                //tripletlist[0 + 2*(count + q)] = gj;
                //tripletlist[1 + 2*(count + q)] = gk;
                tripletlist[(count + q)] = gj;
                tripletlist[(count + q) + Nijk] = gk;                
                count += 1;
            }
        }                        
    }
}

void cpuNeighTriplets(double *xij, double *xik, double *x, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletlist, int *tripletnumsum, 
      int *alist,  int *atomtype, int inum, int dim)
{        
    int Nijk = tripletnumsum[inum];
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ii;       // atom i
        int itype = atomtype[i];        
        int start = tripletnumsum[ii];   
        //std::cout<<" "<<ii<<" "<<i<<std::endl;                              
        int m = tripletnumsum[ii+1]-tripletnumsum[ii];        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom pair (j,k) around atom i
            // int gj = tripletlist[0 + 2*(l + start)];  // ghost index of atom j  
            // int gk = tripletlist[1 + 2*(l + start)];  // ghost index of atom k  
            int gj = tripletlist[(l + start)];  // ghost index of atom j  
            int gk = tripletlist[(l + start) + Nijk];  // ghost index of atom k  
            int j = alist[gj];  // atom j
            int k = alist[gk];  // atom k
            int n = start + l;       
            //std::cout<<" "<<n<<std::endl;                              
            ai[n]        = i;
            aj[n]        = j;    
            ak[n]        = k;    
            ti[n]        = itype;       
            tj[n]        = atomtype[j];     
            tk[n]        = atomtype[k];     
            for (int d=0; d<dim; d++) {
                xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
                xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
            }
        }
    }    
}

void bispectrumderiv(double *dbsa, double *Ar, double *Ai, double *dAr, double *dAi, double *cg, 
    int *pairnum, int *indl, int *indm, int *rowm, int K, int J, int Q, int M, int N, int Nij)
{
    for (int k=0; k<K; k++) {
        double *pAr = &Ar[N*M*k];
        double *pAi = &Ai[N*M*k];
        for (int j=0; j<J; j++) {
            int l2 = indl[j + 0*J];
            int l1 = indl[j + 1*J];
            int l  = indl[j + 2*J];                                     
            int nm = rowm[j+1]-rowm[j];       
            double *pdb = &dbsa[3*Nij*(j + J*k)];
            for (int i=0; i<nm; i++) {
                int q = rowm[j]+i;
                double c = cg[q];
                int m2 = indm[q + 0*Q];
                int m1 = indm[q + 1*Q];
                int m  = indm[q + 2*Q];
                int ml = m + l + (l*l);          
                int ml1 = m1 + l1 + (l1*l1);
                int ml2 = m2 + l2 + (l2*l2);
                for (int n=0; n<N; n++)  {
                    int n0 = pairnum[n];
                    int nj = pairnum[n+1]-n0;

                    double a1 = pAr[n + N*ml];
                    double b1 = pAi[n + N*ml];
                    double a2 = pAr[n + N*ml1];
                    double b2 = pAi[n + N*ml1];
                    double a3 = pAr[n + N*ml2];
                    double b3 = pAi[n + N*ml2];

                    double a2a3 = a2*a3;
                    double a1a3 = a1*a3;
                    double a1a2 = a1*a2;
                    double b1b2 = b1*b2;
                    double b1b3 = b1*b3;
                    double b2b3 = b2*b3;
                    double a1b2 = a1*b2;
                    double a1b3 = a1*b3;
                    double a2b3 = a2*b3;
                    double a2b1 = a2*b1;
                    double a3b1 = a3*b1;
                    double a3b2 = a3*b2;

                    double c1 = c*(a2a3-b2b3);
                    double c2 = c*(a1a3+b1b3);
                    double c3 = c*(a1a2+b1b2);
                    double c4 = c*(a2b3+a3b2);
                    double c5 = c*(a3b1-a1b3);
                    double c6 = c*(a2b1-a1b2);
                    
                    double *pdA1 = &dAr[3*Nij*(ml + M*k)];
                    double *pdA2 = &dAr[3*Nij*(ml1 + M*k)];
                    double *pdA3 = &dAr[3*Nij*(ml2 + M*k)];
                    double *pdA4 = &dAi[3*Nij*(ml + M*k)];                     
                    double *pdA5 = &dAi[3*Nij*(ml1 + M*k)];                     
                    double *pdA6 = &dAi[3*Nij*(ml2 + M*k)]; 

                    for (int ij=0; ij<nj; ij++) {
                        int n1 = n0 + ij;     

                        int d = 0 + 3*n1; 
                        pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       

                        d = 1 + 3*n1; 
                        pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       

                        d = 2 + 3*n1;                  
                        pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       
                    }
                }                                
            }
        }
    }
}

void bispectrumderiv2(double *dbsa, double *Ar, double *Ai, double *dAr, double *dAi, double *cg, 
    int *pairnum, int *indl, int *indm, int *rowm, int *tj, int K, int J, int Q, int S, int M, int N, int Nij)
{
    for (int k=0; k<K; k++) {
        double *pAr = &Ar[N*M*S*k];
        double *pAi = &Ai[N*M*S*k];
        for (int j=0; j<J; j++) {
            int l2 = indl[j + 0*J];
            int l1 = indl[j + 1*J];
            int l  = indl[j + 2*J];                                     
            int nm = rowm[j+1]-rowm[j];       
            double *pdb = &dbsa[3*Nij*(j + J*k)];
            for (int i=0; i<nm; i++) {
                int q = rowm[j]+i;
                double c = cg[q];
                int m2 = indm[q + 0*Q];
                int m1 = indm[q + 1*Q];
                int m  = indm[q + 2*Q];
                int ml = m + l + (l*l);          
                int ml1 = m1 + l1 + (l1*l1);
                int ml2 = m2 + l2 + (l2*l2);
                double *pdA1 = &dAr[3*Nij*(ml + M*k)];
                double *pdA2 = &dAr[3*Nij*(ml1 + M*k)];
                double *pdA3 = &dAr[3*Nij*(ml2 + M*k)];
                double *pdA4 = &dAi[3*Nij*(ml + M*k)];                     
                double *pdA5 = &dAi[3*Nij*(ml1 + M*k)];                     
                double *pdA6 = &dAi[3*Nij*(ml2 + M*k)]; 
                for (int n=0; n<N; n++)  {
                    int n0 = pairnum[n];
                    int nj = pairnum[n+1]-n0;

                    for (int sj=0; sj<S; sj++) {
                        double a1 = pAr[n + N*ml + N*M*sj];
                        double b1 = pAi[n + N*ml + N*M*sj];
                        double a2 = pAr[n + N*ml1 + N*M*sj];
                        double b2 = pAi[n + N*ml1 + N*M*sj];
                        double a3 = pAr[n + N*ml2 + N*M*sj];
                        double b3 = pAi[n + N*ml2 + N*M*sj];

                        double a2a3 = a2*a3;
                        double a1a3 = a1*a3;
                        double a1a2 = a1*a2;
                        double b1b2 = b1*b2;
                        double b1b3 = b1*b3;
                        double b2b3 = b2*b3;
                        double a1b2 = a1*b2;
                        double a1b3 = a1*b3;
                        double a2b3 = a2*b3;
                        double a2b1 = a2*b1;
                        double a3b1 = a3*b1;
                        double a3b2 = a3*b2;

                        double c1 = c*(a2a3-b2b3);
                        double c2 = c*(a1a3+b1b3);
                        double c3 = c*(a1a2+b1b2);
                        double c4 = c*(a2b3+a3b2);
                        double c5 = c*(a3b1-a1b3);
                        double c6 = c*(a2b1-a1b2);
                        
                        for (int ij=0; ij<nj; ij++) {                        
                            int n1 = n0 + ij;     
                            if (tj[n1] == sj) {
                                int d = 0 + 3*n1; 
                                pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       

                                d = 1 + 3*n1; 
                                pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       

                                d = 2 + 3*n1;                  
                                pdb[d] += pdA1[d]*c1 + pdA2[d]*c2 + pdA3[d]*c3 + pdA4[d]*c4 + pdA5[d]*c5 + pdA6[d]*c6;       
                            }
                        }
                    }
                }                                
            }
        }
    }
}

void makeindjk(int *indj, int *indk, int n)
{    
    int k1 = 1;
    for (int i=1; i<=(n-1); i++) {
        for (int j=k1; j<=(k1+n-i-1); j++) {        
            indj[j-1] = i-1;
            indk[j-1] = i+1+j-k1-1;
        }
        k1 = k1 + (n-i);
    }         
}

void makejk0(double *uij, double *uik, double *e2ij, double *ei, int *pairnum, int *tripletnum, 
        int *indj, int *indk, int Nj, int M, int inum, int Nij, int Nijk)
{
    //int count3 = 0;
    //int count2 = 0;
    for (int ii =0; ii<inum; ii++) {    
        int i = ii;     
        int n2 = pairnum[i+1] - pairnum[i]; 
        int n3 = tripletnum[i+1] - tripletnum[i];         

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n2; i2++) {        
                int ind1 = i2 + Nj*i1;    
                int ind2 = pairnum[i] + i2 + Nij*i1;
                ei[ind1] = e2ij[ind2];
            }
            
        makeindjk(indj, indk, n2);

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n3; i2++) {        
                int jj = indj[i2] + Nj*i1;             
                int kk = indk[i2] + Nj*i1;                    
                int ind3 = tripletnum[i] + i2 + Nijk*i1;                
                uij[ind3] = ei[jj];
                uik[ind3] = ei[kk];
            }
    }
}

void makejk(double *uij, double *uik, double *wij, double *wik, double *e2ij, double *f2ij, 
        double *ei, double *f1, double *f2, double *f3, int *pairnum, int *tripletnum, 
        int *indj, int *indk, int Nj, int M, int inum, int Nij, int Nijk)
{
    //int count3 = 0;
    //int count2 = 0;
    for (int ii =0; ii<inum; ii++) {    
        int i = ii;     
        int n2 = pairnum[i+1] - pairnum[i]; 
        int n3 = tripletnum[i+1] - tripletnum[i];         

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n2; i2++) {        
                int ind1 = i2 + Nj*i1;    
                int ind2 = pairnum[i] + i2 + Nij*i1;
                ei[ind1] = e2ij[ind2];
                f1[ind1] = f2ij[0 + 3*ind2];
                f2[ind1] = f2ij[1 + 3*ind2];
                f3[ind1] = f2ij[2 + 3*ind2];
            }
            
        //std::cout<<pairnum[i]<<" "<<tripletnum[i]<<" "<<n2<<" "<<n3<<std::endl;    
        makeindjk(indj, indk, n2);

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n3; i2++) {        
                int jj = indj[i2] + Nj*i1;             
                int kk = indk[i2] + Nj*i1;                    
                int ind3 = tripletnum[i] + i2 + Nijk*i1;                
                uij[ind3] = ei[jj];
                uik[ind3] = ei[kk];
                wij[0 + 3*ind3] = f1[jj];
                wij[1 + 3*ind3] = f2[jj];
                wij[2 + 3*ind3] = f3[jj];
                wik[0 + 3*ind3] = f1[kk];                       
                wik[1 + 3*ind3] = f2[kk];                       
                wik[2 + 3*ind3] = f3[kk];         
            }
        //count3 = count3 + n3;                   
        //count2 = count2 + n2;                        
    }
}

void radialbasis0(double *rbf, double *xij, double *scalefac,
                double rin, double rmax, int pdegree, int K, int M, int N)
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

        for (int j=0; j<M; j++) {            
            double alpha = scalefac[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));

            for (int i=0; i<pdegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*pdegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
            }
        }

        for (int i=0; i<K; i++) {
            int p = pdegree*M + i;
            int nij = n + N*p;     
            double a = pow(dij, (double) (i+1.0));
            rbf[nij] = fcut/a;
        }
    }
}

void radialbasis(double *rbf, double *drbf, double *xij, double *scalefac,
                double rin, double rmax, int pdegree, int K, int M, int N)
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
        //fcut = exp(-1.0/sqrt((1.0 - y^3)^2 + epsil))/exp(-1.0);
        //dfcut = ((3.0/(rmax*exp(-1.0)))*(y^2)*exp(-1.0/(epsil + (y^3 - 1.0)^2)^(1/2))*(y^3 - 1.0))/(epsil + (y^3 - 1.0)^2)^(3/2);
        //if (n==0) std::cout<<rin<<" "<<rmax<<" "<<r<<" "<<y<<" "<<y5<<" "<<y6<<" "<<fcut<<" "<<dfcut<<std::endl;

        for (int j=0; j<M; j++) {            
            double alpha = scalefac[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            double dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(1.0 - exp(-alpha));        

            for (int i=0; i<pdegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*pdegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
                double drbfdr = b*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r*r) + a*cos(a*x)*fcut*dx/r);
                drbf[0 + 3*nij] = drbfdr*dr1;
                drbf[1 + 3*nij] = drbfdr*dr2;
                drbf[2 + 3*nij] = drbfdr*dr3;
            }
        }

        for (int i=0; i<K; i++) {
            int p = pdegree*M + i;
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

void cosinbasis0(double *abf, double *xij, double *xik, int nabf, int N)
{
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, theta; 

    for (int n=0; n<N; n++) {
        xij1 = xij[0+3*n];
        xij2 = xij[1+3*n];
        xij3 = xij[2+3*n];
        xik1 = xik[0+3*n];
        xik2 = xik[1+3*n];
        xik3 = xik[2+3*n];

        xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;
        rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
        riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;    
        rij = pow(rijsq, 0.5); 
        rik = pow(riksq, 0.5); 

        costhe = xdot/(rij*rik);    
        costhe = costhe > 1.0 ? 1.0 : costhe;
        costhe = costhe < -1.0 ? -1.0 : costhe;
        theta = acos(costhe);       

        for (int p=0; p <nabf; p++) 
            abf[n + N*p] = cos((p+1)*theta);                        
    }
}

void cosinbasis(double *abf, double *dabf, double *xij, double *xik, int nabf, int N)
{
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    int np6;

    for (int n=0; n<N; n++) {
        xij1 = xij[0+3*n];
        xij2 = xij[1+3*n];
        xij3 = xij[2+3*n];
        xik1 = xik[0+3*n];
        xik2 = xik[1+3*n];
        xik3 = xik[2+3*n];

        xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;
        rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
        riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;    
        rij = pow(rijsq, 0.5); 
        rik = pow(riksq, 0.5); 

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

        for (int p=0; p <nabf; p++) {
            abf[n + N*p] = cos((p+1)*theta);                
            tm = -(p+1)*sin((p+1)*theta)*dtheta;
            np6 = 6*n + 6*N*p;
            dabf[0 + np6] = tm*dct1;
            dabf[1 + np6] = tm*dct2;
            dabf[2 + np6] = tm*dct3;
            dabf[3 + np6] = tm*dct4;
            dabf[4 + np6] = tm*dct5;
            dabf[5 + np6] = tm*dct6;            
        }
    }
}

void podtally2(double *eatom, double *fatom, double *vatom, double *rij, double *eij, double *fij, 
               int *ai, int *aj, int *ti, int *tj, int *ind, int S, int natom, int M, int N)
{
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        int mij = (ind[typei + typej*S]-1)*M;
        for (int m=0; m<M; m++) {   
            int im = i1 + natom*(m + mij);
            int jm = j1 + natom*(m + mij);
            int nm = n + N*m;
            eatom[im] += eij[nm];
            fatom[0 + 3*im] += fij[0 + 3*nm];
            fatom[1 + 3*im] += fij[1 + 3*nm];
            fatom[2 + 3*im] += fij[2 + 3*nm];
            fatom[0 + 3*jm] -= fij[0 + 3*nm];
            fatom[1 + 3*jm] -= fij[1 + 3*nm];
            fatom[2 + 3*jm] -= fij[2 + 3*nm];

            double v0 = rij[0+3*n]*fij[0 + 3*nm];        
            double v1 = rij[1+3*n]*fij[1 + 3*nm];        
            double v2 = rij[2+3*n]*fij[2 + 3*nm];        
            double v3 = rij[1+3*n]*fij[2 + 3*nm];        
            double v4 = rij[0+3*n]*fij[2 + 3*nm];        
            double v5 = rij[0+3*n]*fij[1 + 3*nm];      
            vatom[0 + 6*im] += v0; 
            vatom[1 + 6*im] += v1;
            vatom[2 + 6*im] += v2; 
            vatom[3 + 6*im] += v3;
            vatom[4 + 6*im] += v4; 
            vatom[5 + 6*im] += v5;        
            vatom[0 + 6*jm] += v0; 
            vatom[1 + 6*jm] += v1;
            vatom[2 + 6*jm] += v2; 
            vatom[3 + 6*jm] += v3;
            vatom[4 + 6*jm] += v4; 
            vatom[5 + 6*jm] += v5;        
        }
    }
}

void podtally2b(double *eatom, double *fatom, double *rij, double *eij, double *fij, 
               int *ai, int *aj, int *ti, int *tj, int *ind, int S, int natom, int M, int N)
{
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        int mij = (ind[typei + typej*S]-1)*M;
        for (int m=0; m<M; m++) {   
            int im = i1 + natom*(m + mij);
            int jm = j1 + natom*(m + mij);
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

void podtally3b0(double *eatom, double *uij, double *uik, double *uijk, int *ai, 
            int *ti, int *tj, int *tk, int *ind, int nrbf, int nabf, int natom, int N, int S)
{
    double eijk, uj, uk, ujk;            
    int S2 = S*(S+1)/2;
    int M = nrbf*(1+nabf);

    int K = -1;
    for (int m =0; m<nrbf; m++) {                    
        K += 1;
        for (int n=0; n<N; n++) {        
            int nm = n + N*m;
            uj = uij[nm];
            uk = uik[nm];
            ujk = uj*uk;
            eijk = uj*uk;

            int typei = ti[n]-1;
            int typej = tj[n]-1;
            int typek = tk[n]-1;
            int nijk =  natom*(K + (ind[typej + typek*S] - 1)*M + (typei)*M*S2);

            int i1 = ai[n];
            int k = i1 + nijk;            
            eatom[k] += eijk;            
        }

        for (int p=0; p<nabf; p++) {               
            K = K + 1;                
            for (int n=0; n<N; n++) {       
                int nm = n + N*m;
                int np = n + N*p;

                uj = uij[nm];
                uk = uik[nm];
                ujk = uj*uk;                
                eijk = ujk*uijk[np];                
           
                int typei = ti[n]-1;
                int typej = tj[n]-1;
                int typek = tk[n]-1;
                int nijk = natom*(K + (ind[typej + typek*S] - 1)*M + (typei)*M*S2);

                int i1 = ai[n];
                int k = i1 + nijk;
                eatom[k] += eijk;
           }
        }
    }
}

void podtally3b(double *eatom, double *fatom, double *xij, double *xik, double *uij, double *uik, 
                double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj, int *ak, 
                int *ti, int *tj, int *tk, int *ind, int nrbf, int nabf, int natom, int N, int S)
{
    double fij0, fij1, fij2; // xij0, xij1, xij2;
    double fik0, fik1, fik2; // xik0, xik1, xik2;
    //double v0, v1, v2, v3, v4, v5;
    double eijk, uj, uk, ujk, wij0, wij1, wij2, wik0, wik1, wik2;
            
    int S2 = S*(S+1)/2;
    int M = nrbf*(1+nabf);

    int K = -1;
    for (int m =0; m<nrbf; m++) {                    
        K += 1;
        for (int n=0; n<N; n++) {        
            int nm = n + N*m;
            uj = uij[nm];
            uk = uik[nm];
            ujk = uj*uk;
            // double *xij3 = &xij[3*n];
            // double *xik3 = &xik[3*n];
            double *wij3 = &wij[3*nm];
            double *wik3 = &wik[3*nm];
            // xij0 = xij3[0];
            // xij1 = xij3[1];
            // xij2 = xij3[2];
            // xik0 = xik3[0];
            // xik1 = xik3[1];
            // xik2 = xik3[2];
            wij0 = wij3[0];
            wij1 = wij3[1];
            wij2 = wij3[2];
            wik0 = wik3[0];
            wik1 = wik3[1];
            wik2 = wik3[2];

            eijk = uj*uk;
            fij0 = wij0*uk;
            fij1 = wij1*uk;
            fij2 = wij2*uk;
            fik0 = uj*wik0;
            fik1 = uj*wik1;
            fik2 = uj*wik2;                        

            int typei = ti[n]-1;
            int typej = tj[n]-1;
            int typek = tk[n]-1;
            int nijk =  natom*(K + (ind[typej + typek*S] - 1)*M + (typei)*M*S2);
            int nK3 = 3*nijk;

            int i1 = ai[n];
            int k = i1 + nijk;            
            int k0 = 0 + 3*i1 + nK3;
            int k1 = 1 + 3*i1 + nK3;
            int k2 = 2 + 3*i1 + nK3;            

            eatom[k] += eijk;
            fatom[k0] += fij0 + fik0;
            fatom[k1] += fij1 + fik1;
            fatom[k2] += fij2 + fik2;              
            
            i1 = aj[n];
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;            
            fatom[k0] -= fij0;
            fatom[k1] -= fij1;
            fatom[k2] -= fij2;

            i1 = ak[n];
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;
            fatom[k0] -= fik0;   
            fatom[k1] -= fik1;   
            fatom[k2] -= fik2;   
        }

        for (int p=0; p<nabf; p++) {               
            K = K + 1;                
            for (int n=0; n<N; n++) {       
                int nm = n + N*m;
                int np = n + N*p;
                int np0 = 0 + 6*n + 6*N*p;
                int np1 = 1 + 6*n + 6*N*p;
                int np2 = 2 + 6*n + 6*N*p;
                int np3 = 3 + 6*n + 6*N*p;
                int np4 = 4 + 6*n + 6*N*p;
                int np5 = 5 + 6*n + 6*N*p;      

                uj = uij[nm];
                uk = uik[nm];
                ujk = uj*uk;
                // double *xij3 = &xij[3*n];
                // double *xik3 = &xik[3*n];
                double *wij3 = &wij[3*nm];
                double *wik3 = &wik[3*nm];
                // xij0 = xij3[0];
                // xij1 = xij3[1];
                // xij2 = xij3[2];
                // xik0 = xik3[0];
                // xik1 = xik3[1];
                // xik2 = xik3[2];
                wij0 = wij3[0];
                wij1 = wij3[1];
                wij2 = wij3[2];
                wik0 = wik3[0];
                wik1 = wik3[1];
                wik2 = wik3[2];

                double u = uijk[np];   
                eijk = ujk*u;                
                fij0 = wij0*uk*u + ujk*wijk[np0];
                fij1 = wij1*uk*u + ujk*wijk[np1];
                fij2 = wij2*uk*u + ujk*wijk[np2];
                fik0 = uj*wik0*u + ujk*wijk[np3];
                fik1 = uj*wik1*u + ujk*wijk[np4];
                fik2 = uj*wik2*u + ujk*wijk[np5];
           
                int typei = ti[n]-1;
                int typej = tj[n]-1;
                int typek = tk[n]-1;
                int nijk = natom*(K + (ind[typej + typek*S] - 1)*M + (typei)*M*S2);
                int nK3 = 3*nijk;

                int i1 = ai[n];
                int k = i1 + nijk;
                int k0 = 0 + 3*i1 + nK3;
                int k1 = 1 + 3*i1 + nK3;
                int k2 = 2 + 3*i1 + nK3;

                eatom[k] += eijk;
                fatom[k0] += fij0 + fik0;
                fatom[k1] += fij1 + fik1;
                fatom[k2] += fij2 + fik2;
                
                i1 = aj[n];
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                fatom[k0] -= fij0;
                fatom[k1] -= fij1;
                fatom[k2] -= fij2;

                i1 = ak[n];
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                fatom[k0] -= fik0;   
                fatom[k1] -= fik1;   
                fatom[k2] -= fik2;   
           }
        }
    }
}

void matmul(double *c, double *a, double *b, int r1, int c1, int c2) 
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

void pod3body(double *eatom, double *fatom, double *y, double *phi, double *gamma0, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *ind, int *pairlist, int *pairnum, 
             int *pairnumsum, int *atomtype, int *alist, int *pdegree, int ngm, int nbf, 
             int nrbf, int nabf, int nelements, int natom, int Nj, int Nij, int Nijk)
{
    int dim = 3;
    
    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    cpuNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);

    int *tripletnum = &tmpint[4*Nij];  // natom
    int *tripletnumsum = &tmpint[4*Nij+natom]; // natom+1
    int *tripletlist = &tmpint[4*Nij+2*natom+1]; // 2*Nijk
    cpuNeighTripletList(tripletlist, tripletnum, tripletnumsum, pairlist, 
            pairnum, pairnumsum, natom);        

    double *e2ij = &tmpmem[3*Nij]; // Nij*nrbf
    double *f2ij = &tmpmem[3*Nij+Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[3*Nij+4*Nij*nrbf]; // Nij*nbf
    double *f2ijt = &tmpmem[3*Nij+4*Nij*nrbf+Nij*nbf]; // dim*Nij*nbf    
    // for (int i = 0; i<3; i++) {                
    //     for (int j = 0; j<Nij; j++) {
    //         std::cout<<rij[i+3*j]<<"  ";
    //     }
    //     std::cout<<std::endl;
    // }
    radialbasis(e2ijt, f2ijt, rij, gamma0, rin, rcut-rin, pdegree[0], pdegree[1], ngm, Nij);
    // for (int i = 0; i<nbf; i++) {                
    //     std::cout<<std::endl;
    //     for (int j = 0; j<nrbf; j++) {
    //         std::cout<<phi[i+nbf*j]<<"  ";
    //     }
    //     std::cout<<std::endl;
    // }
    matmul(e2ij, e2ijt, phi, Nij, nbf, nrbf);
    matmul(f2ij, f2ijt, phi, 3*Nij, nbf, nrbf);

    int n = 3*Nij+ (1+dim)*Nij*nrbf;
    double *uij = &tmpmem[n]; // Nijk*nrbf
    double *uik = &tmpmem[n+Nijk*nrbf]; // Nijk*nrbf
    double *wij = &tmpmem[n+2*Nijk*nrbf]; // dim*Nijk*nrbf
    double *wik = &tmpmem[n+(2+dim)*Nijk*nrbf]; // dim*Nijk*nrbf

    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf;
    double *ei = &tmpmem[n]; // Nj*nrbf 
    double *f1 = &tmpmem[n + Nj*nrbf]; // Nj*nrbf
    double *f2 = &tmpmem[n + 2*Nj*nrbf]; // Nj*nrbf
    double *f3 = &tmpmem[n + 3*Nj*nrbf]; // Nj*nrbf
    
    int m = 4*Nij+2*natom+1+2*Nijk;
    int *indj = &tmpint[m]; //(Nj-1)*Nj/2 
    int *indk = &tmpint[m+(Nj-1)*Nj/2]; //(Nj-1)*Nj/2
    makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnumsum, tripletnumsum, 
            indj, indk, Nj, nrbf, natom, Nij, Nijk);

    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf;
    double *xij = &tmpmem[n]; // dim*Nijk
    double *xik = &tmpmem[n+dim*Nijk]; // dim*Nijk

    m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj;
    ai = &tmpint[m];     // Nijk
    aj = &tmpint[m+Nijk];   // Nijk 
    ti = &tmpint[m+2*Nijk]; // Nijk
    tj = &tmpint[m+3*Nijk]; // Nijk
    int *ak = &tmpint[m+4*Nijk]; // Nijk
    int *tk = &tmpint[m+5*Nijk]; // Nijk
    cpuNeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
        alist, atomtype,  natom, dim);                    
    
    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf + 2*dim*Nijk;
    double *uijk = &tmpmem[n]; // Nijk*nabf
    double *wijk = &tmpmem[n+Nijk*nabf]; // 6*Nijk*nabf
    cosinbasis(uijk, wijk, xij, xik, nabf, Nijk);

    podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
        ti, tj, tk, ind, nrbf, nabf, natom, Nijk, nelements);   

    //n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf + 2*dim*Nijk + 7*Nijk*nabf;
    //m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;
}

void pod3body0(double *eatom, double *y, double *phi, double *gamma0, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *ind, int *pairlist, int *pairnum, 
             int *pairnumsum, int *atomtype, int *alist, int *pdegree, int ngm, int nbf, 
             int nrbf, int nabf, int nelements, int natom, int Nj, int Nij, int Nijk)
{
    int dim = 3;
    
    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    cpuNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);

    int *tripletnum = &tmpint[4*Nij];  // natom
    int *tripletnumsum = &tmpint[4*Nij+natom]; // natom+1
    int *tripletlist = &tmpint[4*Nij+2*natom+1]; // 2*Nijk
    cpuNeighTripletList(tripletlist, tripletnum, tripletnumsum, pairlist, 
            pairnum, pairnumsum, natom);        

    double *e2ij = &tmpmem[3*Nij]; // Nij*nrbf
    double *e2ijt = &tmpmem[3*Nij+Nij*nrbf]; // Nij*nbf
    radialbasis0(e2ijt, rij, gamma0, rin, rcut-rin, pdegree[0], pdegree[1], ngm, Nij);
    matmul(e2ij, e2ijt, phi, Nij, nbf, nrbf);

    int n = 3*Nij+ Nij*nrbf;
    double *uij = &tmpmem[n]; // Nijk*nrbf
    double *uik = &tmpmem[n+Nijk*nrbf]; // Nijk*nrbf

    n = 3*Nij+ Nij*nrbf + 2*Nijk*nrbf;
    double *ei = &tmpmem[n]; // Nj*nrbf 
    
    int m = 4*Nij+2*natom+1+2*Nijk;
    int *indj = &tmpint[m]; //(Nj-1)*Nj/2 
    int *indk = &tmpint[m+(Nj-1)*Nj/2]; //(Nj-1)*Nj/2
    makejk0(uij, uik, e2ij, ei, pairnumsum, tripletnumsum, 
            indj, indk, Nj, nrbf, natom, Nij, Nijk);

    n = 3*Nij + Nij*nrbf + 2*Nijk*nrbf + Nj*nrbf;
    double *xij = &tmpmem[n]; // dim*Nijk
    double *xik = &tmpmem[n+dim*Nijk]; // dim*Nijk

    m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj;
    ai = &tmpint[m];     // Nijk
    aj = &tmpint[m+Nijk];   // Nijk 
    ti = &tmpint[m+2*Nijk]; // Nijk
    tj = &tmpint[m+3*Nijk]; // Nijk
    int *ak = &tmpint[m+4*Nijk]; // Nijk
    int *tk = &tmpint[m+5*Nijk]; // Nijk
    cpuNeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
        alist, atomtype,  natom, dim);                    
    
    n = 3*Nij + Nij*nrbf + 2*Nijk*nrbf + Nj*nrbf + 2*dim*Nijk;
    double *uijk = &tmpmem[n]; // Nijk*nabf
    cosinbasis0(uijk, xij, xik, nabf, Nijk);

    podtally3b0(eatom,  uij, uik, uijk, ai, ti, tj, tk, ind, nrbf, nabf, natom, Nijk, nelements);   
}

void podtally3(double *eatom, double *fatom, double *vatom, double *xij, double *xik, double *uij, 
             double *uik, double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj,
             int *ak, int nrbf, int nabf, int natom, int N)
{
    // int K = -1;
    // for (int m =0; m<nrbf; m++) {         
    //     K += 1;
    //     double *uj = &uij[N*m];
    //     double *uk = &uik[N*m];
    //     double *wj = &wij[3*N*m];
    //     double *wk = &wik[3*N*m];
    //     int nK = natom*K;
    //     int Nm3 = 3*N*m;        
    //     int nK3 = 3*natom*K;
    //     int nK6 = 6*natom*K;        
    //     for (int n=0; n<N; n++) {                      
    //         double *wj3 = &wj[3*n];
    //         double *wk3 = &wk[3*n];
    //         double fij0 = wj3[0]*uk[n];
    //         double fij1 = wj3[1]*uk[n];
    //         double fij2 = wj3[2]*uk[n];
    //         double fik0 = uj[n]*wk3[0];
    //         double fik1 = uj[n]*wk3[1];
    //         double fik2 = uj[n]*wk3[2];         
            
    //         int i1 = ai[n];                      
    //         int k = i1 + nK;            
    //         eatom[k] += uj[n]*uk[n];

    //         k = 3*i1 + nK3;
    //         fatom[k+0] += fij0 + fik0;
    //         fatom[k+1] += fij1 + fik1;
    //         fatom[k+2] += fij2 + fik2;

    //         i1 = aj[n];
    //         k = 3*i1 + nK3;
    //         fatom[k+0] -= fij0;
    //         fatom[k+1] -= fij1;
    //         fatom[k+2] -= fij2;

    //         i1 = ak[n];
    //         k = 3*i1 + nK3;
    //         fatom[k+0] -= fik0;   
    //         fatom[k+1] -= fik1;   
    //         fatom[k+2] -= fik2;   
    //     }
    //     for (int p=0; p<nabf; p++) {               
    //         K = K + 1;                
    //         double *u = &uijk[N*p];
    //         double *w = &wijk[6*N*p];
    //         for (int n=0; n<N; n++) {      
    //             //int k = ai[n] + natom*K;                   
    //             double ujk = uj[n]*uk[n];
    //             double ujn = uj[n]*u[n];
    //             double ukn = uk[n]*u[n];
    //             //eatom[k] += ujk*u[n];             
    //             double *wj3 = &wj[3*n];
    //             double *wk3 = &wk[3*n];
    //             double *w6 = &w[6*n];
    //             double fij0 = wj3[0]*ukn + ujk*w6[0];
    //             double fij1 = wj3[1]*ukn + ujk*w6[1];
    //             double fij2 = wj3[2]*ukn + ujk*w6[2];
    //             double fik0 = wk3[0]*ujn + ujk*w6[3];
    //             double fik1 = wk3[1]*ujn + ujk*w6[4];
    //             double fik2 = wk3[2]*ujn + ujk*w6[5];      

    //             int i1 = ai[n];                      
    //             int k = i1 + nK;            
    //             eatom[k] += ujk*u[n];
                
    //             k = 3*i1 + nK3;
    //             fatom[k+0] += fij0 + fik0;
    //             fatom[k+1] += fij1 + fik1;
    //             fatom[k+2] += fij2 + fik2;

    //             i1 = aj[n];
    //             k = 3*i1 + nK3;
    //             fatom[k+0] -= fij0;
    //             fatom[k+1] -= fij1;
    //             fatom[k+2] -= fij2;

    //             i1 = ak[n];
    //             k = 3*i1 + nK3;
    //             fatom[k+0] -= fik0;   
    //             fatom[k+1] -= fik1;   
    //             fatom[k+2] -= fik2;                         
    //         }
    //     }
    // }

    double fij0, fij1, fij2, xij0, xij1, xij2;
    double fik0, fik1, fik2, xik0, xik1, xik2;
    //double v0, v1, v2, v3, v4, v5;
    double eijk, uj, uk, ujk, wij0, wij1, wij2, wik0, wik1, wik2;
            
    int K = -1;
    for (int m =0; m<nrbf; m++) {                    
        K += 1;
        //int Nm3 = 3*N*m;        
        int nK3 = 3*natom*K;
        int nK6 = 6*natom*K;        
        for (int n=0; n<N; n++) {        
            int nm = n + N*m;
            // int nm0 = 0 + 3*n + Nm3;
            // int nm1 = 1 + 3*n + Nm3;
            // int nm2 = 2 + 3*n + Nm3;
            // int n0 = 0 + 3*n;
            // int n1 = 1 + 3*n;
            // int n2 = 2 + 3*n;

            uj = uij[nm];
            uk = uik[nm];
            ujk = uj*uk;
            double *xij3 = &xij[3*n];
            double *xik3 = &xik[3*n];
            double *wij3 = &wij[3*nm];
            double *wik3 = &wik[3*nm];
            xij0 = xij3[0];
            xij1 = xij3[1];
            xij2 = xij3[2];
            xik0 = xik3[0];
            xik1 = xik3[1];
            xik2 = xik3[2];
            wij0 = wij3[0];
            wij1 = wij3[1];
            wij2 = wij3[2];
            wik0 = wik3[0];
            wik1 = wik3[1];
            wik2 = wik3[2];

            eijk = uj*uk;
            fij0 = wij0*uk;
            fij1 = wij1*uk;
            fij2 = wij2*uk;
            fik0 = uj*wik0;
            fik1 = uj*wik1;
            fik2 = uj*wik2;                        

            // v0 = xij0*fij0 + xik0*fik0;
            // v1 = xij1*fij1 + xik1*fik1;        
            // v2 = xij2*fij2 + xik2*fik2;
            // v3 = xij1*fij2 + xik1*fik2;
            // v4 = xij0*fij2 + xik0*fik2;            
            // v5 = xij0*fij1 + xik0*fik1;                                

            int i1 = ai[n];
            int k = i1 + natom*K;
            int k0 = 0 + 3*i1 + nK3;
            int k1 = 1 + 3*i1 + nK3;
            int k2 = 2 + 3*i1 + nK3;
            // int l0 = 0 + 6*i1 + nK6;
            // int l1 = 1 + 6*i1 + nK6;
            // int l2 = 2 + 6*i1 + nK6;
            // int l3 = 3 + 6*i1 + nK6;
            // int l4 = 4 + 6*i1 + nK6;
            // int l5 = 5 + 6*i1 + nK6;

            eatom[k] += eijk;
            fatom[k0] += fij0 + fik0;
            fatom[k1] += fij1 + fik1;
            fatom[k2] += fij2 + fik2;
            // vatom[l0] += v0; 
            // vatom[l1] += v1;
            // vatom[l2] += v2; 
            // vatom[l3] += v3;
            // vatom[l4] += v4; 
            // vatom[l5] += v5;        
            
            i1 = aj[n];
            k = i1 + natom*K;
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;
            // l0 = 0 + 6*i1 + nK6;
            // l1 = 1 + 6*i1 + nK6;
            // l2 = 2 + 6*i1 + nK6;
            // l3 = 3 + 6*i1 + nK6;
            // l4 = 4 + 6*i1 + nK6;
            // l5 = 5 + 6*i1 + nK6;
            fatom[k0] -= fij0;
            fatom[k1] -= fij1;
            fatom[k2] -= fij2;
            // vatom[l0] += v0; 
            // vatom[l1] += v1;
            // vatom[l2] += v2; 
            // vatom[l3] += v3;
            // vatom[l4] += v4; 
            // vatom[l5] += v5;        

            i1 = ak[n];
            k = i1 + natom*K;
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;
            // l0 = 0 + 6*i1 + nK6;
            // l1 = 1 + 6*i1 + nK6;
            // l2 = 2 + 6*i1 + nK6;
            // l3 = 3 + 6*i1 + nK6;
            // l4 = 4 + 6*i1 + nK6;
            // l5 = 5 + 6*i1 + nK6;
            fatom[k0] -= fik0;   
            fatom[k1] -= fik1;   
            fatom[k2] -= fik2;   
            // vatom[l0] += v0; 
            // vatom[l1] += v1;
            // vatom[l2] += v2; 
            // vatom[l3] += v3;
            // vatom[l4] += v4; 
            // vatom[l5] += v5;        
        }

        for (int p=0; p<nabf; p++) {               
            K = K + 1;                
            nK3 = 3*natom*K;
            nK6 = 6*natom*K;
            for (int n=0; n<N; n++) {       
                int nm = n + N*m;
                // int nm0 = 0 + 3*n + Nm3;
                // int nm1 = 1 + 3*n + Nm3;
                // int nm2 = 2 + 3*n + Nm3;
                // int n0 = 0 + 3*n;
                // int n1 = 1 + 3*n;
                // int n2 = 2 + 3*n;
                int np = n + N*p;
                int np0 = 0 + 6*n + 6*N*p;
                int np1 = 1 + 6*n + 6*N*p;
                int np2 = 2 + 6*n + 6*N*p;
                int np3 = 3 + 6*n + 6*N*p;
                int np4 = 4 + 6*n + 6*N*p;
                int np5 = 5 + 6*n + 6*N*p;      

                uj = uij[nm];
                uk = uik[nm];
                ujk = uj*uk;
                double *xij3 = &xij[3*n];
                double *xik3 = &xik[3*n];
                double *wij3 = &wij[3*nm];
                double *wik3 = &wik[3*nm];
                xij0 = xij3[0];
                xij1 = xij3[1];
                xij2 = xij3[2];
                xik0 = xik3[0];
                xik1 = xik3[1];
                xik2 = xik3[2];
                wij0 = wij3[0];
                wij1 = wij3[1];
                wij2 = wij3[2];
                wik0 = wik3[0];
                wik1 = wik3[1];
                wik2 = wik3[2];

                double u = uijk[np];   
                eijk = ujk*u;                
                fij0 = wij0*uk*u + ujk*wijk[np0];
                fij1 = wij1*uk*u + ujk*wijk[np1];
                fij2 = wij2*uk*u + ujk*wijk[np2];
                fik0 = uj*wik0*u + ujk*wijk[np3];
                fik1 = uj*wik1*u + ujk*wijk[np4];
                fik2 = uj*wik2*u + ujk*wijk[np5];
           
                // v0 = xij0*fij0 + xik0*fik0;
                // v1 = xij1*fij1 + xik1*fik1;        
                // v2 = xij2*fij2 + xik2*fik2;
                // v3 = xij1*fij2 + xik1*fik2;
                // v4 = xij0*fij2 + xik0*fik2;            
                // v5 = xij0*fij1 + xik0*fik1;                                

                int i1 = ai[n];
                int k = i1 + natom*K;
                int k0 = 0 + 3*i1 + nK3;
                int k1 = 1 + 3*i1 + nK3;
                int k2 = 2 + 3*i1 + nK3;
                // int l0 = 0 + 6*i1 + nK6;
                // int l1 = 1 + 6*i1 + nK6;
                // int l2 = 2 + 6*i1 + nK6;
                // int l3 = 3 + 6*i1 + nK6;
                // int l4 = 4 + 6*i1 + nK6;
                // int l5 = 5 + 6*i1 + nK6;
                eatom[k] += eijk;
                fatom[k0] += fij0 + fik0;
                fatom[k1] += fij1 + fik1;
                fatom[k2] += fij2 + fik2;
                // vatom[l0] += v0; 
                // vatom[l1] += v1;
                // vatom[l2] += v2; 
                // vatom[l3] += v3;
                // vatom[l4] += v4; 
                // vatom[l5] += v5;        
                
                i1 = aj[n];
                k = i1 + natom*K;
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                // l0 = 0 + 6*i1 + nK6;
                // l1 = 1 + 6*i1 + nK6;
                // l2 = 2 + 6*i1 + nK6;
                // l3 = 3 + 6*i1 + nK6;
                // l4 = 4 + 6*i1 + nK6;
                // l5 = 5 + 6*i1 + nK6;
                fatom[k0] -= fij0;
                fatom[k1] -= fij1;
                fatom[k2] -= fij2;
                // vatom[l0] += v0; 
                // vatom[l1] += v1;
                // vatom[l2] += v2; 
                // vatom[l3] += v3;
                // vatom[l4] += v4; 
                // vatom[l5] += v5;      

                i1 = ak[n];
                k = i1 + natom*K;
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                // l0 = 0 + 6*i1 + nK6;
                // l1 = 1 + 6*i1 + nK6;
                // l2 = 2 + 6*i1 + nK6;
                // l3 = 3 + 6*i1 + nK6;
                // l4 = 4 + 6*i1 + nK6;
                // l5 = 5 + 6*i1 + nK6;
                fatom[k0] -= fik0;   
                fatom[k1] -= fik1;   
                fatom[k2] -= fik2;   
                // vatom[l0] += v0; 
                // vatom[l1] += v1;
                // vatom[l2] += v2; 
                // vatom[l3] += v3;
                // vatom[l4] += v4; 
                // vatom[l5] += v5;            
           }
        }
    }
    // for (int n=0; n<N; n++) {
    //     int K = -1;
    //     for (int m =0; m<nrbf; m++) {
    //         K += 1;
    //         int nm = n + N*m;
    //         int nm0 = 0 + 3*n + 3*N*m;
    //         int nm1 = 1 + 3*n + 3*N*m;
    //         int nm2 = 2 + 3*n + 3*N*m;
    //         int n0 = 0 + 3*n;
    //         int n1 = 1 + 3*n;
    //         int n2 = 2 + 3*n;

    //         double uj = uij[nm];
    //         double uk = uik[nm];
    //         double ujk = uj*uk;
    //         double wij0 = wij[nm0];
    //         double wij1 = wij[nm1];
    //         double wij2 = wij[nm2];
    //         double wik0 = wik[nm0];
    //         double wik1 = wik[nm1];
    //         double wik2 = wik[nm2];

    //         double eijk = uj*uk;
    //         fij[0] = wij0*uk;
    //         fij[1] = wij1*uk;
    //         fij[2] = wij2*uk;
    //         fik[0] = uj*wik0;
    //         fik[1] = uj*wik1;
    //         fik[2] = uj*wik2;            

    //         double v0 = xij[n0]*fij[0] + xik[n0]*fik[0];
    //         double v1 = xij[n1]*fij[1] + xik[n1]*fik[1];        
    //         double v2 = xij[n2]*fij[2] + xik[n2]*fik[2];
    //         double v3 = xij[n1]*fij[2] + xik[n1]*fik[2];
    //         double v4 = xij[n0]*fij[2] + xik[n0]*fik[2];            
    //         double v5 = xij[n0]*fij[1] + xik[n0]*fik[1];                                

    //         int i1 = ai[n];
    //         int k = i1 + natom*K;
    //         int k0 = 0 + 3*i1 + 3*natom*K;
    //         int k1 = 1 + 3*i1 + 3*natom*K;
    //         int k2 = 2 + 3*i1 + 3*natom*K;
    //         int l0 = 0 + 6*i1 + 6*natom*K;
    //         int l1 = 1 + 6*i1 + 6*natom*K;
    //         int l2 = 2 + 6*i1 + 6*natom*K;
    //         int l3 = 3 + 6*i1 + 6*natom*K;
    //         int l4 = 4 + 6*i1 + 6*natom*K;
    //         int l5 = 5 + 6*i1 + 6*natom*K;

    //         eatom[k] += eijk;
    //         fatom[k0] += fij[0] + fik[0];
    //         fatom[k1] += fij[1] + fik[1];
    //         fatom[k2] += fij[2] + fik[2];
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        
            
    //         i1 = aj[n];
    //         k = i1 + natom*K;
    //         k0 = 0 + 3*i1 + 3*natom*K;
    //         k1 = 1 + 3*i1 + 3*natom*K;
    //         k2 = 2 + 3*i1 + 3*natom*K;
    //         l0 = 0 + 6*i1 + 6*natom*K;
    //         l1 = 1 + 6*i1 + 6*natom*K;
    //         l2 = 2 + 6*i1 + 6*natom*K;
    //         l3 = 3 + 6*i1 + 6*natom*K;
    //         l4 = 4 + 6*i1 + 6*natom*K;
    //         l5 = 5 + 6*i1 + 6*natom*K;
    //         fatom[k0] -= fij[0];
    //         fatom[k1] -= fij[1];
    //         fatom[k2] -= fij[2];
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        

    //         i1 = ak[n];
    //         k = i1 + natom*K;
    //         k0 = 0 + 3*i1 + 3*natom*K;
    //         k1 = 1 + 3*i1 + 3*natom*K;
    //         k2 = 2 + 3*i1 + 3*natom*K;
    //         l0 = 0 + 6*i1 + 6*natom*K;
    //         l1 = 1 + 6*i1 + 6*natom*K;
    //         l2 = 2 + 6*i1 + 6*natom*K;
    //         l3 = 3 + 6*i1 + 6*natom*K;
    //         l4 = 4 + 6*i1 + 6*natom*K;
    //         l5 = 5 + 6*i1 + 6*natom*K;
    //         fatom[k0] -= fik[0];   
    //         fatom[k1] -= fik[1];   
    //         fatom[k2] -= fik[2];   
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        

    //         for (int p=0; p<nabf; p++) {               
    //             K = K + 1;                
    //             int np = n + N*p;
    //             int np0 = 0 + 6*n + 6*N*p;
    //             int np1 = 1 + 6*n + 6*N*p;
    //             int np2 = 2 + 6*n + 6*N*p;
    //             int np3 = 3 + 6*n + 6*N*p;
    //             int np4 = 4 + 6*n + 6*N*p;
    //             int np5 = 5 + 6*n + 6*N*p;            
    //             double u = uijk[np];   
    //             eijk = ujk*u;                
    //             fij[0] = wij0*uk*u + ujk*wijk[np0];
    //             fij[1] = wij1*uk*u + ujk*wijk[np1];
    //             fij[2] = wij2*uk*u + ujk*wijk[np2];
    //             fik[0] = uj*wik0*u + ujk*wijk[np3];
    //             fik[1] = uj*wik1*u + ujk*wijk[np4];
    //             fik[2] = uj*wik2*u + ujk*wijk[np5];
           
    //             v0 = xij[n0]*fij[0] + xik[n0]*fik[0];
    //             v1 = xij[n1]*fij[1] + xik[n1]*fik[1];        
    //             v2 = xij[n2]*fij[2] + xik[n2]*fik[2];
    //             v3 = xij[n1]*fij[2] + xik[n1]*fik[2];
    //             v4 = xij[n0]*fij[2] + xik[n0]*fik[2];            
    //             v5 = xij[n0]*fij[1] + xik[n0]*fik[1];                                

    //             i1 = ai[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             eatom[k] += eijk;
    //             fatom[k0] += fij[0] + fik[0];
    //             fatom[k1] += fij[1] + fik[1];
    //             fatom[k2] += fij[2] + fik[2];
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;        
                
    //             i1 = aj[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             fatom[k0] -= fij[0];
    //             fatom[k1] -= fij[1];
    //             fatom[k2] -= fij[2];
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;      

    //             i1 = ak[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             fatom[k0] -= fik[0];   
    //             fatom[k1] -= fik[1];   
    //             fatom[k2] -= fik[2];   
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;            
    //        }
    //     }
    // }
}

#endif

