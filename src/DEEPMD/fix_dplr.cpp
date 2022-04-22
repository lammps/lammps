#include <iostream>
#include <iomanip>
#include <limits>
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "fix.h"
#include "fix_dplr.h"
#include "pppm_dplr.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

static bool 
is_key (const string& input) 
{
  vector<string> keys ;
  keys.push_back("model");
  keys.push_back("type_associate");
  keys.push_back("bond_type");
  keys.push_back("efield");
  for (int ii = 0; ii < keys.size(); ++ii){
    if (input == keys[ii]) {
      return true;
    }
  }
  return false;
}


FixDPLR::FixDPLR(LAMMPS *lmp, int narg, char **arg) 
    :Fix(lmp, narg, arg), 
     efield(3, 0.0), 
     efield_fsum(4, 0.0), 
     efield_fsum_all(4, 0.0), 
     efield_force_flag(0)
{
#if LAMMPS_VERSION_NUMBER>=20210210
  // lammps/lammps#2560
  energy_global_flag = 1;
  virial_global_flag = 1;
#else
  virial_flag = 1;
#endif

  if (strcmp(update->unit_style,"metal") != 0) {
    error->all(FLERR,"Pair deepmd requires metal unit, please set it by \"units metal\"");
  }
  
  int iarg = 3;
  vector<int> map_vec;
  bond_type.clear();
  while (iarg < narg) {
    if (! is_key(arg[iarg])) {
      error->all(FLERR,"Illegal pair_style command\nwrong number of parameters\n");
    }
    if (string(arg[iarg]) == string("model")) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix adapt command");
      model = string(arg[iarg+1]);
      iarg += 2;
    }
    else if (string(arg[iarg]) == string("efield")) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command, efield should be provided 3 float numbers");
      efield[0] = atof(arg[iarg+1]);
      efield[1] = atof(arg[iarg+2]);
      efield[2] = atof(arg[iarg+3]);
      iarg += 4;
    }
    else if (string(arg[iarg]) == string("type_associate")) {
      int iend = iarg+1;
      while (iend < narg && (! is_key(arg[iend]) )) {
	map_vec.push_back(atoi(arg[iend])-1);
	iend ++;
      }
      iarg = iend;
    }
    else if (string(arg[iarg]) == string("bond_type")) {
      int iend = iarg+1;
      while (iend < narg && (! is_key(arg[iend]) )) {
	bond_type.push_back(atoi(arg[iend])-1);
	iend ++;
      }
      sort(bond_type.begin(), bond_type.end());
      iarg = iend;
    }
    else {
      break;
    }
  }
  assert(map_vec.size() % 2 == 0), "number of ints provided by type_associate should be even";
  for (int ii = 0; ii < map_vec.size()/2; ++ii){
    type_asso[map_vec[ii*2+0]] = map_vec[ii*2+1];
    bk_type_asso[map_vec[ii*2+1]] = map_vec[ii*2+0];
  }

  // dpt.init(model);
  // dtm.init("frozen_model.pb");
  dpt.init(model, 0, "dipole_charge");
  dtm.init(model, 0, "dipole_charge");

  sel_type = dpt.sel_types();
  sort(sel_type.begin(), sel_type.end());
  dpl_type.clear();
  for (int ii = 0; ii < sel_type.size(); ++ii){
    dpl_type.push_back(type_asso[sel_type[ii]]);
  }

  pair_deepmd = (PairDeepMD *) force->pair_match("deepmd",1);
  if (!pair_deepmd) {
    error->all(FLERR,"pair_style deepmd should be set before this fix\n");
  }

  // set comm size needed by this fix
  comm_reverse = 3;
}

int FixDPLR::setmask()
{
  int mask = 0;
#if LAMMPS_VERSION_NUMBER<20210210
  // THERMO_ENERGY removed in lammps/lammps#2560
  mask |= THERMO_ENERGY;
#endif
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

void FixDPLR::init()
{
  // double **xx = atom->x;
  // double **vv = atom->v;
  // int nlocal = atom->nlocal;
  // for (int ii = 0; ii < nlocal; ++ii){
  //   cout << xx[ii][0] << " " 
  // 	 << xx[ii][1] << " " 
  // 	 << xx[ii][2] << "   " 
  // 	 << vv[ii][0] << " " 
  // 	 << vv[ii][1] << " " 
  // 	 << vv[ii][2] << " " 
  // 	 << endl;
  // }
}

void FixDPLR::setup(int vflag)
{
  // if (strstr(update->integrate_style,"verlet"))
  //   post_force(vflag);
  // else {
  //   error->all(FLERR, "respa is not supported by this fix");
  // }
  if (vflag) {
    v_setup(vflag);
  }
  else {
    evflag = 0;
  }
}


void
FixDPLR::get_valid_pairs(vector<pair<int,int> >& pairs)
{
  pairs.clear();
  
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  vector<int > dtype (nall);
  // get type
  {
    int *type = atom->type;
    for (int ii = 0; ii < nall; ++ii){
      dtype[ii] = type[ii] - 1;
    }
  }

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  for (int ii = 0; ii < nbondlist; ++ii){
    int idx0=-1, idx1=-1;
    int bd_type = bondlist[ii][2] - 1;
    if ( ! binary_search(bond_type.begin(), bond_type.end(), bd_type) ){
      continue;
    }
    if (binary_search(sel_type.begin(), sel_type.end(), dtype[bondlist[ii][0]]) 
	&& 
	binary_search(dpl_type.begin(), dpl_type.end(), dtype[bondlist[ii][1]])
	){
      idx0 = bondlist[ii][0];
      idx1 = bondlist[ii][1];
    }
    else if (binary_search(sel_type.begin(), sel_type.end(), dtype[bondlist[ii][1]])
	     &&
	     binary_search(dpl_type.begin(), dpl_type.end(), dtype[bondlist[ii][0]])
	){
      idx0 = bondlist[ii][1];
      idx1 = bondlist[ii][0];
    }
    else {
      error->all(FLERR, "find a bonded pair the types of which are not associated");
    }
    if ( ! (idx0 < nlocal && idx1 < nlocal) ){
      error->all(FLERR, "find a bonded pair that is not on the same processor, something should not happen");
    }
    pairs.push_back(pair<int,int>(idx0, idx1));
  }
}

void FixDPLR::post_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  vector<pair<int,int> > valid_pairs;
  get_valid_pairs(valid_pairs);  
  
  for (int ii = 0; ii < valid_pairs.size(); ++ii){
    int idx0 = valid_pairs[ii].first;
    int idx1 = valid_pairs[ii].second;
    for (int dd = 0; dd < 3; ++dd){
      x[idx1][dd] = x[idx0][dd] ;
      v[idx1][dd] = v[idx0][dd] ;
      // v[idx1][dd] = 0.0;
    }
  }
}

void FixDPLR::pre_force(int vflag)
{
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  // if (eflag_atom) {
  //   error->all(FLERR,"atomic energy calculation is not supported by this fix\n");
  // }
  
  // declear inputs
  vector<int > dtype (nall);
  vector<FLOAT_PREC > dbox (9, 0) ;
  vector<FLOAT_PREC > dcoord (nall * 3, 0.);
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
  // get lammps nlist
  NeighList * list = pair_deepmd->list;
  deepmd::InputNlist lmp_list (list->inum, list->ilist, list->numneigh, list->firstneigh);
  // declear output
  vector<FLOAT_PREC> tensor;
  // compute
  dpt.compute(tensor, dcoord, dtype, dbox, nghost, lmp_list);
  // cout << "tensor of size " << tensor.size() << endl;
  // cout << "nghost " << nghost << endl;
  // cout << "nall " << dtype.size() << endl;
  // cout << "nloc " << nlocal << endl;
  // for (int ii = 0; ii < tensor.size(); ++ii){
  //   if (ii%3 == 0){
  //     cout << endl;
  //   }
  //   cout << tensor[ii] << "\t";
  // }
  // cout << endl;
  // for (int ii = 0; ii < nlocal * 3; ++ii){
  //   if (ii%3 == 0){
  //     cout << endl;
  //   }
  //   cout << dcoord[ii] << "\t";
  // }
  // int max_type = 0;
  // for (int ii = 0; ii < dtype.size(); ++ii){
  //   if (dtype[ii] > max_type) {
  //     max_type = dtype[ii];
  //   }
  // }

  // selected type
  vector<int> dpl_type;
  for (int ii = 0; ii < sel_type.size(); ++ii){
    dpl_type.push_back(type_asso[sel_type[ii]]);
  }
  vector<int> sel_fwd, sel_bwd;
  int sel_nghost;
  deepmd::select_by_type(sel_fwd, sel_bwd, sel_nghost, dcoord, dtype, nghost, sel_type);
  int sel_nall = sel_bwd.size();
  int sel_nloc = sel_nall - sel_nghost;
  vector<int> sel_type(sel_bwd.size());
  deepmd::select_map<int>(sel_type, dtype, sel_fwd, 1);
  
  // Yixiao: because the deeptensor already return the correct order, the following map is no longer needed
  // deepmd::AtomMap<FLOAT_PREC> atom_map(sel_type.begin(), sel_type.begin() + sel_nloc);
  // const vector<int> & sort_fwd_map(atom_map.get_fwd_map());

  vector<pair<int,int> > valid_pairs;
  get_valid_pairs(valid_pairs);  
  
  int odim = dpt.output_dim();
  assert(odim == 3);
  dipole_recd.resize(nall * 3);
  fill(dipole_recd.begin(), dipole_recd.end(), 0.0);
  for (int ii = 0; ii < valid_pairs.size(); ++ii){
    int idx0 = valid_pairs[ii].first;
    int idx1 = valid_pairs[ii].second;
    assert(idx0 < sel_fwd.size()); // && sel_fwd[idx0] < sort_fwd_map.size());
    // Yixiao: the sort map is no longer needed
    // int res_idx = sort_fwd_map[sel_fwd[idx0]];
    int res_idx = sel_fwd[idx0];
    // int ret_idx = dpl_bwd[res_idx];
    for (int dd = 0; dd < 3; ++dd){
      x[idx1][dd] = x[idx0][dd] + tensor[res_idx * 3 + dd];
      // res_buff[idx1 * odim + dd] = tensor[res_idx * odim + dd];
      dipole_recd[idx0*3+dd] = tensor[res_idx * 3 + dd];
    }
  }
  // cout << "-------------------- fix/dplr: pre force " << endl;
  // for (int ii = 0; ii < nlocal; ++ii){
  //   cout << ii << "    ";
  //   for (int dd = 0; dd < 3; ++dd){
  //     cout << x[ii][dd] << " " ;
  //   }
  //   cout << endl;
  // }
}


void FixDPLR::post_force(int vflag)
{
  if (vflag) {
    v_setup(vflag);
  }
  else {
    evflag = 0;
  }
  if (vflag_atom) {
    error->all(FLERR,"atomic virial calculation is not supported by this fix\n");
  }

  PPPMDPLR * pppm_dplr = (PPPMDPLR*) force->kspace_match("pppm/dplr", 1);
  if (!pppm_dplr) {
    error->all(FLERR,"kspace_style pppm/dplr should be set before this fix\n");
  }
  const vector<double > & dfele_(pppm_dplr->get_fele());
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  vector<FLOAT_PREC> dcoord(nall*3, 0.0), dbox(9, 0.0), dfele(nlocal*3, 0.0);
  vector<int> dtype(nall, 0);
  // set values for dcoord, dbox, dfele
  {
    int *type = atom->type;
    for (int ii = 0; ii < nall; ++ii){
      dtype[ii] = type[ii] - 1;
    }
    dbox[0] = domain->h[0];	// xx
    dbox[4] = domain->h[1];	// yy
    dbox[8] = domain->h[2];	// zz
    dbox[7] = domain->h[3];	// zy
    dbox[6] = domain->h[4];	// zx
    dbox[3] = domain->h[5];	// yx
    // get coord
    double ** x = atom->x;
    for (int ii = 0; ii < nall; ++ii){
      for (int dd = 0; dd < 3; ++dd){
	dcoord[ii*3+dd] = x[ii][dd] - domain->boxlo[dd];
      }
    }
    assert(dfele_.size() == nlocal * 3);
    // revise force according to efield
    for (int ii = 0; ii < nlocal*3; ++ii){
      dfele[ii] = dfele_[ii];
    }
    // revise force and virial according to efield
    double * q = atom->q;
    imageint *image = atom->image;
    double unwrap[3];
    double v[6];
    efield_fsum[0] = efield_fsum[1] = efield_fsum[2] = efield_fsum[3] = 0.0;
    efield_force_flag = 0;
    for (int ii = 0; ii < nlocal; ++ii){
      double tmpf[3];
      for (int dd = 0; dd < 3; ++dd){
	tmpf[dd] = q[ii] * efield[dd];
      }
      for (int dd = 0; dd < 3; ++dd){
	dfele[ii*3+dd] += tmpf[dd];
      }
      domain->unmap(x[ii],image[ii],unwrap);
      efield_fsum[0] -= tmpf[0]*unwrap[0]+tmpf[1]*unwrap[1]+tmpf[2]*unwrap[2];
      efield_fsum[1] += tmpf[0];
      efield_fsum[2] += tmpf[1];
      efield_fsum[3] += tmpf[2];
      if (evflag) {
	v[0] = tmpf[0] *unwrap[0];
	v[1] = tmpf[1] *unwrap[1];
	v[2] = tmpf[2] *unwrap[2];
	v[3] = tmpf[0] *unwrap[1];
	v[4] = tmpf[0] *unwrap[2];
	v[5] = tmpf[1] *unwrap[2];
	v_tally(ii, v);
      }
    }
  }
  // lmp nlist
  NeighList * list = pair_deepmd->list;
  deepmd::InputNlist lmp_list (list->inum, list->ilist, list->numneigh, list->firstneigh);
  // bonded pairs
  vector<pair<int,int> > valid_pairs;
  get_valid_pairs(valid_pairs);  
  // output vects
  vector<FLOAT_PREC> dfcorr, dvcorr;
  // compute
  dtm.compute(dfcorr, dvcorr, dcoord, dtype, dbox, valid_pairs, dfele, nghost, lmp_list);
  assert(dfcorr.size() == dcoord.size());
  assert(dfcorr.size() == nall * 3);
  // backward communication of fcorr
  dfcorr_buff.resize(dfcorr.size());
  copy(dfcorr.begin(), dfcorr.end(), dfcorr_buff.begin());
#if LAMMPS_VERSION_NUMBER>=20220324
  comm->reverse_comm(this,3);
#else
  comm->reverse_comm_fix(this,3);
#endif
  copy(dfcorr_buff.begin(), dfcorr_buff.end(), dfcorr.begin());
  // // check and print
  // cout << "-------------------- fix/dplr: post force " << endl;
  // cout << "dfcorr.size() " << dfcorr.size() << endl;
  // cout << "dcoord.size() " << dcoord.size() << endl;
  // for (int ii = 0; ii < nlocal; ++ii){
  //   cout << ii << "\t x: ";
  //   for (int dd = 0; dd < 3; ++dd){
  //     cout << dcoord[ii*3+dd] << " \t " ;
  //   }
  //   cout << ii << "\t f: ";
  //   for (int dd = 0; dd < 3; ++dd){
  //     cout << dfcorr[ii*3+dd] << " \t " ;
  //   }
  //   cout << endl;
  // }
  // apply the force correction
  double ** f = atom->f;
  for (int ii = 0; ii < nlocal; ++ii){
    for(int dd = 0; dd < 3; ++dd){
      f[ii][dd] += dfcorr[ii*3+dd];
    }
  }
  // cout << "virial corr1 ";
  // for (int ii = 0; ii < 9; ++ii){
  //   cout << dvcorr[ii] << " " ;
  // }
  // cout << endl;
  for (int ii = 0; ii < valid_pairs.size(); ++ii){
    int idx0 = valid_pairs[ii].first;
    int idx1 = valid_pairs[ii].second;
    for (int dd0 = 0; dd0 < 3; ++dd0){
      for (int dd1 = 0; dd1 < 3; ++dd1){
	dvcorr[dd0*3+dd1] -= dfele[idx1*3+dd0] * dipole_recd[idx0*3+dd1];
      }
    }    
  }
  // cout << "virial corr2 ";
  // for (int ii = 0; ii < 9; ++ii){
  //   cout << dvcorr[ii] << " " ;
  // }
  // cout << endl;
  if (evflag){
    double vv[6] = {0.0};
    vv[0] += dvcorr[0];
    vv[1] += dvcorr[4];
    vv[2] += dvcorr[8];
    vv[3] += dvcorr[3];
    vv[4] += dvcorr[6];
    vv[5] += dvcorr[7];
    v_tally(0, vv);
  }
}


int FixDPLR::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = dfcorr_buff[3*i+0];
    buf[m++] = dfcorr_buff[3*i+1];
    buf[m++] = dfcorr_buff[3*i+2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixDPLR::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    dfcorr_buff[3*j+0] += buf[m++];
    dfcorr_buff[3*j+1] += buf[m++];
    dfcorr_buff[3*j+2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixDPLR::compute_scalar(void)
{
  if (efield_force_flag == 0) {
    MPI_Allreduce(&efield_fsum[0],&efield_fsum_all[0],4,MPI_DOUBLE,MPI_SUM,world);
    efield_force_flag = 1;
  }
  return efield_fsum_all[0];
}

/* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

double FixDPLR::compute_vector(int n)
{
  if (efield_force_flag == 0) {
    MPI_Allreduce(&efield_fsum[0],&efield_fsum_all[0],4,MPI_DOUBLE,MPI_SUM,world);
    efield_force_flag = 1;
  }
  return efield_fsum_all[n+1];
}
