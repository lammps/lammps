/* -*- c++ -*- ---------------------------------------------------------- */
#include "compute_esp_grid.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid3d.h"
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
  spacing(1.0),nx(0),ny(0),nz(0),bcut_acks2(nullptr),reaxflag(0),
  esp_grid(nullptr),esp(nullptr),reference(nullptr),weight(nullptr)
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

  // Allocate grid
  esp_grid = new Grid3d(lmp, world, nx, ny, nz);
  esp_grid->setup_grid(
    ixlo, ixhi, iylo, iyhi, izlo, izhi,
    oxlo, oxhi, oylo, oyhi, ozlo, ozhi
  );

  // Allocate memory
  memory->create3d_offset(esp, izlo, izhi, iylo, iyhi, ixlo, ixhi, "esp/grid:esp");
  memory->create3d_offset(reference, izlo, izhi, iylo, iyhi, ixlo, ixhi, "esp/grid:reference");
  memory->create3d_offset(weight, izlo, izhi, iylo, iyhi, ixlo, ixhi, "esp/grid:weight");

  // Set up vector interface
  int local_size = (izhi-izlo+1) * (iyhi-iylo+1) * (ixhi-ixlo+1);
  size_vector = local_size;
  vector = &(reference[izlo][iylo][ixlo]);

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

void ComputeESPGrid::compute_pergrid()
{
  invoked_pergrid = update->ntimestep;

  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;
  
  // FIXED: Explicitly calculate grid spacing - do not compute inside the loop
  double dx = (oxhi - oxlo) / std::max(1, ixhi - ixlo);
  double dy = (oyhi - oylo) / std::max(1, iyhi - iylo);
  double dz = (ozhi - ozlo) / std::max(1, izhi - izlo);
  
  // Debug grid spacing
  // utils::logmesg(lmp, "*** Grid spacing: dx={:.6f}, dy={:.6f}, dz={:.6f}\n", dx, dy, dz);

  for (int iz = izlo; iz <= izhi; ++iz) {
    double gz = ozlo + (iz - izlo + 0.5) * dz;
    for (int iy = iylo; iy <= iyhi; ++iy) {
      double gy = oylo + (iy - iylo + 0.5) * dy;
      for (int ix = ixlo; ix <= ixhi; ++ix) {
        double gx = oxlo + (ix - ixlo + 0.5) * dx;

        // Initialize ESP and weight values for this cell
        esp[iz][iy][ix] = 0.0;
        weight[iz][iy][ix] = 0.0;
        

        // Find nearest atom type including ghost atoms
        double rmin2 = 1e60;
        int tnear = -1;
        for (int i = 0; i < ntotal; ++i) {
          double dx = gx - x[i][0];
          double dy = gy - x[i][1];
          double dz = gz - x[i][2];
          
          // Apply minimum image convention for periodic boundaries
          domain->minimum_image(dx, dy, dz);
          
          double r2 = dx*dx + dy*dy + dz*dz;
          if (r2 < rmin2) {rmin2 = r2; tnear = type[i];}
        }
        
        if (tnear < 0 || tnear > atom->ntypes) continue;

        double rcut = reaxflag ? bcut_acks2[tnear] : bcut_acks2[0];
        double r = std::sqrt(rmin2);

        if (r < 1.4 || r > rcut) {
          weight[iz][iy][ix] = 0.0;
          continue;
        }
        weight[iz][iy][ix] = compute_weight(r, rcut);

        // Calculate ESP including contributions from all atoms
        for (int i = 0; i < ntotal; ++i) {
          double dx = gx - x[i][0];
          double dy = gy - x[i][1];
          double dz = gz - x[i][2];
          
          // Apply minimum image convention for periodic boundaries
          domain->minimum_image(dx, dy, dz);
          
          double r2 = dx*dx + dy*dy + dz*dz + 1e-12;
          esp[iz][iy][ix] += q[i] / std::sqrt(r2);
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
  invoked_scalar = update->ntimestep;

  double loss_sum = 0.0, weight_sum = 0.0;

  for (int iz = izlo; iz <= izhi; ++iz) {
    for (int iy = iylo; iy <= iyhi; ++iy) {
      for (int ix = ixlo; ix <= ixhi; ++ix) {
        double w = weight[iz][iy][ix];
        if (w == 0.0) continue;
        double diff = esp[iz][iy][ix] - reference[iz][iy][ix];
        loss_sum += w * diff * diff;
        weight_sum += w;
        //utils::logmesg(lmp, "*** iz iy ix {} {} {} weight {} esp {} reference {}\n", iz, iy, ix, w, esp[iz][iy][ix], reference[iz][iy][ix] );
      }
    }
  }

  double local_result[2], global_result[2];
  local_result[0] = loss_sum;
  local_result[1] = weight_sum;
  
  MPI_Allreduce(local_result, global_result, 2, MPI_DOUBLE, MPI_SUM, world);
  
  scalar = (global_result[1] > 0.0) ? global_result[0] / global_result[1] : 0.0;
  if (comm->me == 0) utils::logmesg(lmp, "*** scalar {}\n", scalar);
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
  if (name == "esp")
    return 0;
  else if (name == "reference")
    return 1;

  return -1;
}

void *ComputeESPGrid::get_grid_by_index(int index)
{
  if( (index == 0) || (index == 1) )
    return esp_grid;

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
