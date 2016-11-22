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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_psiN_gift.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "complex"
#include "vector"
#include "algorithm"
#include "iterator"




using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */



ComputePsiNGift::ComputePsiNGift(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  //if (domain->dimension != 2) error->all(FLERR,"PsiN is only for 2d systems for xy plane");  
  if (narg != 6) error->all(FLERR,"Illegal compute coord/PsiN command");

  
  Npsi = atof(arg[3]);
  
  sphere_step=atof(arg[4]);
  
  double cutoff = atof(arg[5]);
  cutsq = cutoff*cutoff;
  
  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;
  PsiN = NULL;
  
}

/* ---------------------------------------------------------------------- */

ComputePsiNGift::~ComputePsiNGift()
{

 memory->destroy(PsiN);

}

/* ---------------------------------------------------------------------- */

void ComputePsiNGift::init()
{
  if (force->pair == NULL) 
    error->all(FLERR,"Compute coord/PsiN requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce) 
    error->all(FLERR,"Compute coord/PsiN cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"coord/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute coord/atom");
}

/* ---------------------------------------------------------------------- */

void ComputePsiNGift::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePsiNGift::compute_peratom()
{
    
  int i,j,ii,jj,inum,jnum,n,nn;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  
  double cosA,sinA,phi;
  
  std::complex<double> psi = std::complex<double>(0.0);
  std::complex<double> j1 = std::complex<double>(0.0,1.0);
  
  std::vector<int> hull;  // id of atoms of closest hull
  std::vector<Point> cvx_hull;
  std::vector<int>::iterator it;

  int *image = atom->image;
  nn=0;
  int *ilist,*jlist,*numneigh,**firstneigh;
 

  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nlocal > nmax) {
   memory->destroy(PsiN);
   nmax = atom->nmax;
   memory->create(PsiN,nmax,4,"PsiNGift/atom:PsiN");
   array_atom = PsiN;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute coordination number for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;
  double PI= 4*atan(1); //Pi
  //Reset counters
  for (i=0;i<inum;i++) PsiN[i]==0;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    std::vector<int> neib_keys;
    std::vector<double> distance;
    std::vector<double> angle;
    std::vector<ComputePsiNGift::Point> neib_points;
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      n = 0;
      for (jj = 0; jj < jnum; jj++) {
	    j = jlist[jj];
	    if (j >= nall) j %= nall;
        
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely;
        if (rsq < cutsq)
        {
           neib_keys.push_back(j);
           distance.push_back(sqrt(rsq));
           angle.push_back(calculatePhi(delx,dely));
           Point p;
           p.key=j;
           if (delx !=0) {
           p.x=delx/rsq;
           } else{
           p.x=0;
           };
           if (dely !=0){
           p.y=dely/rsq;
           } else {
           p.y=0;
           }
           neib_points.push_back(p);
        }
      }
      
      //Pishem suda
       //hull = createHull(i,neib_keys,distance,angle);
       cvx_hull=convex_hull(neib_points);

       std::complex<double> psi = std::complex<double>(0.0);

       /*for (int kk =0; kk < cvx_hull.size();kk++)
        {
          printf("%i:%i\n",j,cvx_hull[kk].key);
        }*/

   
        for(int bb=0;bb<cvx_hull.size();bb++)
        {  
         
          delx = x[cvx_hull[bb].key][0] - xtmp;
          dely = x[cvx_hull[bb].key][1] - ytmp; 
                   
          phi=calculatePhi(delx,dely);
          psi=psi+exp(j1*phi*Npsi);
          nn=nn+1;
        }
      
        
        PsiN[i][0] = abs(psi);  // module
        PsiN[i][1] = psi.real(); // real part
        PsiN[i][2] = psi.imag(); // image part
        PsiN[i][3] = cvx_hull.size(); //near atoms
     
     }
    }  
     
  }
  
//typedef double coord_t;   // coordinate type
//typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2


// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
ComputePsiNGift::coord2_t ComputePsiNGift::cross(const Point &O, const Point &A, const Point &B)
{
	return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
std::vector<ComputePsiNGift::Point> ComputePsiNGift::convex_hull(std::vector<ComputePsiNGift::Point> P)
{
	int n = P.size(), k = 0;
	std::vector<Point> H(2*n);

	// Sort points lexicographically
	std::sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k-1);
	return H;
}


std::vector<int> ComputePsiNGift::createHull(int i,std::vector<int> keys,std::vector<double> distance,double** &x)
{
    
    


}



std::vector<int> ComputePsiNGift::createHull(int i,std::vector<int> keys,std::vector<double> distance, std::vector<double> angle)
{
   
 std::vector<int> hull;
 sort2arraysbyvalue(keys,distance);

 

 if(distance.size()>0){
     
 double prev_distance=distance[0];
 double max_dist=2*distance[0];
      
    for (int ii=0;ii<distance.size();ii++)
    {   
        //if((distance[ii]-prev_distance)/prev_distance<sphere_step & hull.size()<=Npsi+1)
        if(distance[ii]<max_dist)
       // if(hull.size()<=Npsi-1)
        {
            hull.push_back(keys[ii]);   
            prev_distance=distance[ii];                
        } else
        {
            break;
        }
        
        
    }
}
 
       return hull;
       
}


std::vector<double> ComputePsiNGift::sort(std::vector<double> array)
{
   for (int ii=0;ii<array.size()-1;ii++)
   {
       if(array[ii]>array[ii+1])
       {
           
           double tmp_array;
           tmp_array=array[ii];
           array[ii]=array[ii+1];
           array[ii+1]=tmp_array;
           printf("afterswitch=%f,%f\n",array[ii],array[ii+1]);
           ii=0;
       }
   }
   return array;
}




void ComputePsiNGift::sort2arraysbyvalue(std::vector<int> &keys,std::vector<double> &distance)
{
    if (keys.size()>-1){
    
    for (int ii=0;ii<keys.size()-1;ii++)
   {
          
       if(distance[ii]>distance[ii+1])
       {    
           
           int tmp_key=0;
      
           tmp_key= keys[ii];
           keys[ii]=keys[ii+1];
           keys[ii+1]=tmp_key;
           
           double tmp_distance=0;
           tmp_distance=distance[ii];
           distance[ii]=distance[ii+1];
           distance[ii+1]=tmp_distance;
           ii=-1;
       }
       
   }
   }
   /* 
    printf("Sorted :");
    
    for (int ii=0;ii<keys.size();ii++)
    {
        printf("%i-%f\n,",keys[ii],distance[ii]);
    }
    
     printf("\n:");
    */
}


double ComputePsiNGift::calculatePhi(double delx,double dely)
{
 double phi=0;
 double PI= 4*atan(1); //Pi
 // printf("delx=%f,dely=%f\n",delx,dely);  
 double cosA=delx/sqrt(delx*delx+dely*dely);
 double sinA=dely/sqrt(delx*delx+dely*dely);
  
    if (sinA>=0 & cosA>=0  ) phi=asin(sinA);
    if (sinA>=0 & cosA < 0 ) phi=acos(cosA);
    if (sinA < 0 & cosA < 0) phi = -asin(sinA) + PI;
    if (sinA < 0 & cosA >= 0) phi = asin(sinA) + 2*PI; 
 
 return phi;   
    
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */



double ComputePsiNGift::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  return bytes;
}
