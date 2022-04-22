#include "compute_deeptensor_atom.h"
#include <cstring>
#include <algorithm>
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"

using namespace LAMMPS_NS;

#ifdef HIGH_PREC
#define VALUETYPE double
#else
#define VALUETYPE float
#endif

/* ---------------------------------------------------------------------- */

ComputeDeeptensorAtom::ComputeDeeptensorAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  dp(lmp),
  tensor(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute deeptensor/atom command");

  // parse args
  std::string model_file = std::string(arg[3]);

  // initialize deeptensor
  int gpu_rank = dp.get_node_rank();
  std::string model_file_content = dp.get_file_content(model_file);
  dt.init(model_file, gpu_rank);
  sel_types = dt.sel_types();
  std::sort(sel_types.begin(), sel_types.end());

  peratom_flag = 1;
  size_peratom_cols = dt.output_dim();
  pressatomflag = 0;
  timeflag = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeDeeptensorAtom::~ComputeDeeptensorAtom()
{
  memory->destroy(tensor);
}

/* ---------------------------------------------------------------------- */

void ComputeDeeptensorAtom::init()
{
  // need an occasional full neighbor list

#if LAMMPS_VERSION_NUMBER>=20220324
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
#else
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
#endif
}

void ComputeDeeptensorAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeDeeptensorAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local tensor array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(tensor);
    nmax = atom->nmax;
    memory->create(tensor, nmax, size_peratom_cols, "deeptensor/atom:tensor");
    array_atom = tensor;
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_pair = force->newton_pair;

  std::vector<VALUETYPE > dcoord (nall * 3, 0.);
  std::vector<VALUETYPE > dbox (9, 0) ;
  std::vector<int > dtype (nall);
  // get type
  for (int ii = 0; ii < nall; ++ii){
    dtype[ii] = type[ii] - 1;
  }
  // get box
  dbox[0] = domain->h[0];	// xx
  dbox[4] = domain->h[1];	// yy
  dbox[8] = domain->h[2];	// zz
  dbox[7] = domain->h[3];	// zy
  dbox[6] = domain->h[4];	// zx
  dbox[3] = domain->h[5];	// yx
  // get coord
  for (int ii = 0; ii < nall; ++ii){
    for (int dd = 0; dd < 3; ++dd){
      dcoord[ii*3+dd] = x[ii][dd] - domain->boxlo[dd];
    }
  }

  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);
  deepmd::InputNlist lmp_list (list->inum, list->ilist, list->numneigh, list->firstneigh);
  
  // declare outputs
  std::vector<VALUETYPE > gtensor, force, virial, atensor, avirial;

  // compute tensors
  dt.compute (gtensor, force, virial, atensor, avirial,
	      dcoord, dtype, dbox, nghost, lmp_list);
  
  // store the result in tensor
  int iter_tensor = 0;
  for(int ii = 0; ii < nlocal; ++ii){
    std::vector<int>::iterator _it =
	std::find(sel_types.begin(), sel_types.end(), dtype[ii]);
    bool selected = (_it != sel_types.end());    
    bool ingroup = (mask[ii] & groupbit);    
    // record when selected and in group
    if (selected && ingroup){
      for(int jj = 0; jj < size_peratom_cols; ++jj){
	tensor[ii][jj] = atensor[iter_tensor+jj];
      }
    }
    // if not selected or not in group set to 0.
    else{
      for(int jj = 0; jj < size_peratom_cols; ++jj){
	tensor[ii][jj] = 0.0;
      }
    }    
    if (selected) {
      iter_tensor += size_peratom_cols;
    }
  }
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDeeptensorAtom::memory_usage()
{
  double bytes = nmax*size_peratom_cols * sizeof(double);
  return bytes;
}
