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

ComputeESPGrid::ComputeESPGrid(LAMMPS *lmp,int narg,char **arg) :
  Compute(lmp,narg,arg),
  spacing(0.3),nx(0),ny(0),nz(0),
  esp_ref(nullptr),bcut_acks2(nullptr),reaxflag(0)
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
}

ComputeESPGrid::~ComputeESPGrid()
{
  memory->destroy(esp_ref);
  if(!reaxflag) memory->destroy(bcut_acks2);
}

void ComputeESPGrid::init()
{
  xlo=domain->boxlo[0];
  ylo=domain->boxlo[1];
  zlo=domain->boxlo[2];

  nx=std::max(1,int(std::ceil(domain->prd[0]/spacing)));
  ny=std::max(1,int(std::ceil(domain->prd[1]/spacing)));
  nz=std::max(1,int(std::ceil(domain->prd[2]/spacing)));

  memory->create3d_offset(esp_ref,0,nz-1,0,ny-1,0,nx-1,"esp/ref");

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

double ComputeESPGrid::compute_scalar()
{
  invoked_scalar=update->ntimestep;

  double **x=atom->x;
  double *q =atom->q;
  int   *type=atom->type;
  int    nlocal=atom->nlocal;

  double loss_sum=0.0,weight_sum=0.0;

  for(int iz=0;iz<nz;++iz){
    double gz=zlo+(iz+0.5)*spacing;
    for(int iy=0;iy<ny;++iy){
      double gy=ylo+(iy+0.5)*spacing;
      for(int ix=0;ix<nx;++ix){
        double gx=xlo+(ix+0.5)*spacing;

        double rmin2=1e60;
        int tnear=-1;
        for(int a=0;a<nlocal;++a){
          double dx=gx-x[a][0],dy=gy-x[a][1],dz=gz-x[a][2];
          double r2=dx*dx+dy*dy+dz*dz;
          if(r2<rmin2){rmin2=r2;tnear=type[a];}
        }
        if(tnear<0) continue;

        double rcut=reaxflag?bcut_acks2[tnear]:bcut_acks2[0];
        double r=std::sqrt(rmin2);
        if(r<1.4||r>rcut) continue;

        double Vmodel=0.0;
        for(int a=0;a<nlocal;++a){
          double dx=gx-x[a][0],dy=gy-x[a][1],dz=gz-x[a][2];
          double r2=dx*dx+dy*dy+dz*dz+1e-12;
          Vmodel+=q[a]/std::sqrt(r2);
        }

        double diff=Vmodel-esp_ref[iz][iy][ix];
        double w=weight(r,rcut);

        loss_sum+=w*diff*diff;
        weight_sum+=w;
      }
    }
  }
  return (weight_sum>0.0)?loss_sum/weight_sum:0.0;
}

inline double ComputeESPGrid::weight(double r,double rcut) const
{
  double w=1.0/(r*r);
  w*=1.0-std::exp(-std::pow((r-1.4)/0.3,6));
  w*=std::exp(-std::pow(r/rcut,6));
  return w;
}

void *ComputeESPGrid::extract_reference()
{
  return static_cast<void *>(esp_ref);
}
