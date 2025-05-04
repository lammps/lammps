/* -*- c++ -*- ---------------------------------------------------------- */
#include "compute_esp_grid.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "update.h"
#include "utils.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeESPGrid::ComputeESPGrid(LAMMPS *lmp,int narg,char **arg) :
  Compute(lmp,narg,arg),
  spacing(1.0),nx(0),ny(0),nz(0),
  bcut_acks2(nullptr),esp(nullptr),reference(nullptr),reaxflag(0)
{
  if(narg<4) error->all(FLERR,"Illegal compute esp/grid command");

  int iarg=3;
  while(iarg<narg){
    if(strcmp(arg[iarg],"spacing")==0){
      if(iarg+1>=narg) utils::missing_cmd_args(FLERR,"compute esp/grid",error);
      spacing=utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg+=2;
    } else error->all(FLERR,"Unknown keyword in compute esp/grid");
  }

  scalar_flag=1;
  extscalar=1;
  vector_flag=1;
  extvector=0;
  pergrid_flag=1;
  invoked_pergrid=-1; // make sure invoked with run 0

  allocate_grid();

}

/* ---------------------------------------------------------------------- */

ComputeESPGrid::~ComputeESPGrid()
{
  deallocate_grid();
  if(!reaxflag) memory->destroy(bcut_acks2);
}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::allocate_grid()
{
  // Get current box dimensions
  xlo=domain->boxlo[0];
  ylo=domain->boxlo[1];
  zlo=domain->boxlo[2];

  // Calculate grid dimensions
  nx=std::max(1,int(std::ceil(domain->prd[0]/spacing)));
  ny=std::max(1,int(std::ceil(domain->prd[1]/spacing)));
  nz=std::max(1,int(std::ceil(domain->prd[2]/spacing)));

  utils::logmesg(lmp, "*** xlo ylo zlo {} {} {} prd {} {} {} spacing {} nx ny nz {} {} {}\n",
    xlo, ylo, zlo, domain->prd[0], domain->prd[1], domain->prd[2], spacing, nx, ny, nz);

  // Allocate memory
  memory->create3d_offset(esp,0,nz-1,0,ny-1,0,nx-1,"esp/grid:esp");
  memory->create3d_offset(reference,0,nz-1,0,ny-1,0,nx-1,"esp/grid:reference");
  memory->create3d_offset(weight,0,nz-1,0,ny-1,0,nx-1,"esp/grid:weight");
  
  // Set up vector interface
  size_vector = nx * ny * nz;
  vector = &(reference[0][0][0]);

}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::deallocate_grid()
{
  memory->destroy(esp);
  memory->destroy(reference);
  memory->destroy(weight);
}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::reset_grid()
{
  deallocate_grid();
  allocate_grid();
}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::init()
{
  int flag=0;
  Pair *pair=force->pair_match("^reaxff",0);
  if(pair){
    reaxflag=1;
    bcut_acks2=static_cast<double *>(pair->extract("bcut_acks2",flag));
    if(!flag) error->all(FLERR,"compute esp/grid could not extract bcut_acks2");
  } else {
    reaxflag=0;
    memory->create(bcut_acks2,1,"esp/bcut_dummy");
    bcut_acks2[0]=5.0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::compute_pergrid()
{
  invoked_pergrid = update->ntimestep;

  double **x=atom->x;
  double *q =atom->q;
  int   *type=atom->type;
  int    nlocal=atom->nlocal;

  for(int iz=0;iz<nz;++iz){
    double gz=zlo+(iz+0.5)*spacing;
    for(int iy=0;iy<ny;++iy){
      double gy=ylo+(iy+0.5)*spacing;
      for(int ix=0;ix<nx;++ix){
        double gx=xlo+(ix+0.5)*spacing;

        double rmin2=1e60;
        int tnear=-1;
        for(int i=0;i<nlocal;++i){
          double dx=gx-x[i][0],dy=gy-x[i][1],dz=gz-x[i][2];
          double r2=dx*dx+dy*dy+dz*dz;
          if(r2<rmin2){rmin2=r2;tnear=type[i];}
        }
        if(tnear<0) continue;

        double rcut=reaxflag?bcut_acks2[tnear]:bcut_acks2[0];
        double r=std::sqrt(rmin2);
        esp[iz][iy][ix]=0.0;
        if(r<1.4||r>rcut) {
          weight[iz][iy][ix]=0.0;
          continue;
        }
        weight[iz][iy][ix]=compute_weight(r,rcut);

        for(int i=0;i<nlocal;++i){
          double dx=gx-x[i][0],dy=gy-x[i][1],dz=gz-x[i][2];
          double r2=dx*dx+dy*dy+dz*dz+1e-12;
          esp[iz][iy][ix]+=q[i]/std::sqrt(r2);
        }

      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeESPGrid::compute_vector()
{
  invoked_vector = update->ntimestep;
  // vector pointer is already set to reference in allocate_grid()
}

/* ---------------------------------------------------------------------- */

double ComputeESPGrid::compute_scalar()
{
  if (invoked_pergrid != update->ntimestep) compute_pergrid();
  invoked_scalar=update->ntimestep;

  double loss_sum=0.0,weight_sum=0.0;

  for(int iz=0;iz<nz;++iz){
    for(int iy=0;iy<ny;++iy){
      for(int ix=0;ix<nx;++ix){
        double w=weight[iz][iy][ix];
        if(w==0.0) continue;
        double diff=esp[iz][iy][ix]-reference[iz][iy][ix];
        loss_sum+=w*diff*diff;
        weight_sum+=w;
      }
    }
  }

  scalar = (weight_sum>0.0)?loss_sum/weight_sum:0.0;
  utils::logmesg(lmp, "*** scalar {}\n", scalar);
  return scalar;

}

inline double ComputeESPGrid::compute_weight(double r,double rcut) const
{
  double w=1.0/(r*r);
  w*=1.0-std::exp(-std::pow((r-1.4)/0.3,6));
  w*=std::exp(-std::pow(r/rcut,6));
  return w;
}

int ComputeESPGrid::get_grid_by_name(const std::string &name, int &dim)
{
  if (name == "grid")
    return 0;
  else if (name == "reference")
    return 1;

  return -1;
}

void *ComputeESPGrid::get_grid_by_index(int index)
{
  if (index == 0)
    return esp_grid;
  else if (index == 1)
    return reference_grid;

  return nullptr;
}

int ComputeESPGrid::get_griddata_by_name(int igrid, const std::string &name, int &ncol)
{
  if (((igrid == 0) || (igrid == 1)) && (name == "data")) {
    ncol = 0;
    return 0;
  }

  return -1;
}

void *ComputeESPGrid::get_griddata_by_index(int index)
{
  if (index == 0)
    return esp;
  else if (index == 1)
    return reference;

  return nullptr;
}
