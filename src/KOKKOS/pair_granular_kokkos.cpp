// clang-format off
/* ----------------------------------------------------------------------
   Kokkos implementation of the modular pair_style granular.

   This first implementation supports the non-cohesive Hooke/Hertz normal
   families, standard damping laws, and linear/Mindlin tangential models.
------------------------------------------------------------------------- */

#include "pair_granular_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_neigh_history_kokkos.h"
#include "force.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_rolling.h"
#include "gran_sub_mod_tangential.h"
#include "gran_sub_mod_twisting.h"
#include "gran_sub_mod_heat.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <string>

using namespace LAMMPS_NS;
using namespace LAMMPS_NS::Granular_NS;

template<class DeviceType>
PairGranularKokkos<DeviceType>::PairGranularKokkos(LAMMPS *lmp) : PairGranular(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK |
    TYPE_MASK | MASK_MASK | ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  fix_historyKK = nullptr;
  use_history_kk = size_history_kk = 0;
}

template<class DeviceType>
PairGranularKokkos<DeviceType>::~PairGranularKokkos()
{
  if (copymode) return;
  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::create_kokkos_history()
{
  int size_max[NSUBMODELS] = {0};
  use_history = 0;
  size_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;

  for (int n = 0; n < nmodels; n++) {
    auto *model = models_list[n];
    if (model->beyond_contact) {
      beyond_contact = 1;
      use_history = 1;
    }
    if (model->size_history) use_history = 1;
    if (model->nondefault_history_transfer) nondefault_history_transfer = 1;
    for (int m = 0; m < NSUBMODELS; m++)
      size_max[m] = MAX(size_max[m],model->sub_models[m]->size_history);
  }

  if (!use_history) return;
  for (int m = 0; m < NSUBMODELS; m++) size_history += size_max[m];
  size_history = MAX(size_history,1);
  for (int n = 0; n < nmodels; n++) {
    int next = 0;
    for (int m = 0; m < NSUBMODELS; m++) {
      models_list[n]->sub_models[m]->history_index = next;
      next += size_max[m];
    }
  }

  if (id_history == nullptr)
    id_history = utils::strdup(std::string("NEIGH_HISTORY_GRANULAR") +
                               std::to_string(instance_index()));

  if (fix_history == nullptr) {
    std::string cmd = std::string(id_history) + " all ";
    if (execution_space == Device) cmd += "NEIGH_HISTORY/KK/DEVICE ";
    else cmd += "NEIGH_HISTORY/KK/HOST ";
    cmd += std::to_string(size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix(id_dummy,cmd,1));
    fix_history->pair = this;
  } else {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id(id_history));
    if (!fix_history) error->all(FLERR,"Could not find Kokkos granular history fix");
  }
  fix_historyKK = dynamic_cast<FixNeighHistoryKokkos<DeviceType> *>(fix_history);
  if (!fix_historyKK) error->all(FLERR,"Could not create Kokkos granular history fix");
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::init_style()
{
  create_kokkos_history();
  PairGranular::init_style();

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Pair granular/kk requires a half or half/thread neighbor list");
  if (use_history && force->newton_pair)
    error->all(FLERR,"Pair granular/kk with contact history requires newton off");
  if (fix_rigid)
    error->all(FLERR,"Pair granular/kk does not yet support fix rigid masses");

  use_history_kk = use_history;
  size_history_kk = size_history;
}

template<class DeviceType>
double PairGranularKokkos<DeviceType>::init_one(int i, int j)
{
  const double cut = PairGranular::init_one(i,j);
  pack_models();
  return cut;
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::pack_models()
{
  d_model_i = decltype(d_model_i)("pair_granular:model_i",nmodels,MI_COUNT);
  d_model_f = decltype(d_model_f)("pair_granular:model_f",nmodels,MF_COUNT);
  d_type_model = decltype(d_type_model)("pair_granular:type_model",atom->ntypes+1,
                                        atom->ntypes+1);
  auto hmi = Kokkos::create_mirror_view(d_model_i);
  auto hmf = Kokkos::create_mirror_view(d_model_f);
  auto htm = Kokkos::create_mirror_view(d_type_model);

  for (int n = 0; n < nmodels; n++) {
    auto *gm = models_list[n];
    const std::string &normal = gm->normal_model->name;
    const std::string &damping = gm->damping_model->name;
    const std::string &tangential = gm->tangential_model->name;

    if (normal == "hooke") hmi(n,MI_NORMAL) = N_HOOKE;
    else if (normal == "hertz" || normal == "hertz/material")
      hmi(n,MI_NORMAL) = N_HERTZ;
    else
      error->all(FLERR,"Pair granular/kk does not yet support normal model {}",normal);

    if (damping == "velocity") hmi(n,MI_DAMPING) = D_VELOCITY;
    else if (damping == "mass_velocity") hmi(n,MI_DAMPING) = D_MASS_VELOCITY;
    else if (damping == "viscoelastic") hmi(n,MI_DAMPING) = D_VISCOELASTIC;
    else if (damping == "tsuji" || damping == "coeff_restitution")
      hmi(n,MI_DAMPING) = D_TSUJI;
    else
      error->all(FLERR,"Pair granular/kk does not yet support damping model {}",damping);

    if (tangential == "linear_nohistory")
      hmi(n,MI_TANGENTIAL) = T_LINEAR_NOHISTORY;
    else if (tangential == "linear_history")
      hmi(n,MI_TANGENTIAL) = T_LINEAR_HISTORY;
    else if (tangential == "linear_history_classic")
      hmi(n,MI_TANGENTIAL) = T_LINEAR_CLASSIC;
    else if (tangential == "mindlin_classic")
      hmi(n,MI_TANGENTIAL) = T_MINDLIN_CLASSIC;
    else if (tangential == "mindlin" || tangential == "mindlin/force")
      hmi(n,MI_TANGENTIAL) = T_MINDLIN;
    else
      error->all(FLERR,"Pair granular/kk does not yet support tangential model {}",tangential);

    if (gm->rolling_model->name != "none" || gm->twisting_model->name != "none" ||
        gm->heat_model->name != "none")
      error->all(FLERR,"Pair granular/kk does not yet support rolling, twisting, or heat models");
    if (gm->synchronized_verlet)
      error->all(FLERR,"Pair granular/kk does not yet support synchronized_verlet");
    if (gm->tangential_model->get_mindlin_rescale())
      error->all(FLERR,"Pair granular/kk does not yet support Mindlin rescale history");

    hmi(n,MI_HISTORY_INDEX) = gm->tangential_model->history_index;
    hmi(n,MI_LIMIT_DAMPING) = gm->limit_damping;
    hmi(n,MI_MINDLIN_FORCE) = gm->tangential_model->get_mindlin_force();
    hmi(n,MI_MINDLIN_RESCALE) = gm->tangential_model->get_mindlin_rescale();
    hmf(n,MF_NORMAL_K) = gm->normal_model->get_k();
    hmf(n,MF_DAMP) = gm->damping_model->get_damp();
    hmf(n,MF_TANGENTIAL_K) = gm->tangential_model->get_k();
    hmf(n,MF_XT) = gm->tangential_model->get_xt();
    hmf(n,MF_MU) = gm->tangential_model->get_mu();
  }
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++) htm(i,j) = types_indices[i][j];

  Kokkos::deep_copy(d_model_i,hmi);
  Kokkos::deep_copy(d_model_f,hmf);
  Kokkos::deep_copy(d_type_model,htm);
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;
  eflag = eflag_in;
  vflag = vflag_in;
  ev_init(eflag,vflag,0);
  const int history_update = update->setupflag == 0;

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  for (int k = 0; k < 4; k++) special_lj[k] = force->special_lj[k];
  dt_kk = update->dt;

  const int inum = list->inum;
  auto *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  if (use_history_kk) {
    fix_historyKK->k_firstflag.template sync<DeviceType>();
    fix_historyKK->k_firstvalue.template sync<DeviceType>();
    d_firsttouch = fix_historyKK->k_firstflag.template view<DeviceType>();
    d_firsthistory = fix_historyKK->k_firstvalue.template view<DeviceType>();
    Kokkos::deep_copy(d_firsttouch,0);
  }

  EV_FLOAT ev;
#define RUN_GRAN(NF,NP,VF,HU) \
  if (VF) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairGranularCompute<NF,NP,1,HU>>(0,inum),*this,ev); \
  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairGranularCompute<NF,NP,0,HU>>(0,inum),*this)

  if (neighflag == HALF) {
    if (newton_pair) {
      if (history_update) { if (vflag_either) RUN_GRAN(HALF,1,1,1); else RUN_GRAN(HALF,1,0,1); }
      else { if (vflag_either) RUN_GRAN(HALF,1,1,0); else RUN_GRAN(HALF,1,0,0); }
    } else {
      if (history_update) { if (vflag_either) RUN_GRAN(HALF,0,1,1); else RUN_GRAN(HALF,0,0,1); }
      else { if (vflag_either) RUN_GRAN(HALF,0,1,0); else RUN_GRAN(HALF,0,0,0); }
    }
  } else {
    if (newton_pair) {
      if (history_update) { if (vflag_either) RUN_GRAN(HALFTHREAD,1,1,1); else RUN_GRAN(HALFTHREAD,1,0,1); }
      else { if (vflag_either) RUN_GRAN(HALFTHREAD,1,1,0); else RUN_GRAN(HALFTHREAD,1,0,0); }
    } else {
      if (history_update) { if (vflag_either) RUN_GRAN(HALFTHREAD,0,1,1); else RUN_GRAN(HALFTHREAD,0,0,1); }
      else { if (vflag_either) RUN_GRAN(HALFTHREAD,0,1,0); else RUN_GRAN(HALFTHREAD,0,0,0); }
    }
  }
#undef RUN_GRAN

  if (use_history_kk) {
    fix_historyKK->k_firstflag.template modify<DeviceType>();
    fix_historyKK->k_firstvalue.template modify<DeviceType>();
  }
  if (eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
  if (vflag_global) for (int k = 0; k < 6; k++) virial[k] += (double) ev.v[k];
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  copymode = 0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::operator()(
    TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,HISTORYUPDATE>,
    const int ii, EV_FLOAT &ev) const
{
  using AtomicView = Kokkos::View<KK_ACC_FLOAT*[3],typename DAT::t_kkacc_1d_3::array_layout,
    typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>>;
  AtomicView a_f = f;
  AtomicView a_torque = torque;
  const int i = d_ilist[ii];
  const KK_FLOAT xi=x(i,0), yi=x(i,1), zi=x(i,2);
  const KK_FLOAT vxi=v(i,0), vyi=v(i,1), vzi=v(i,2);
  const KK_FLOAT wxi=omega(i,0), wyi=omega(i,1), wzi=omega(i,2);
  const KK_FLOAT ri=radius[i], mi=rmass[i];
  const int itype=type[i], maski=mask[i], jnum=d_numneigh[i];
  KK_ACC_FLOAT fxi=0, fyi=0, fzi=0, txi=0, tyi=0, tzi=0;

  for (int jj=0; jj<jnum; jj++) {
    int j=d_neighbors(i,jj);
    KK_FLOAT factor=(KK_FLOAT)special_lj[sbmask(j)];
    j &= NEIGHMASK;
    if (factor == 0) continue;
    const KK_FLOAT dx=xi-x(j,0), dy=yi-x(j,1), dz=zi-x(j,2);
    const KK_FLOAT rsq=dx*dx+dy*dy+dz*dz;
    const KK_FLOAT rj=radius[j], radsum=ri+rj;
    const int model=d_type_model(itype,type[j]);
    if (rsq >= radsum*radsum) {
      if (use_history_kk) {
        for (int k=0; k<size_history_kk; k++)
          d_firsthistory(i,size_history_kk*jj+k)=0;
      }
      continue;
    }
    if (use_history_kk) d_firsttouch(i,jj)=1;

    const KK_FLOAT r=sqrt(rsq), rinv=1/r;
    const KK_FLOAT nx=dx*rinv, ny=dy*rinv, nz=dz*rinv;
    const KK_FLOAT delta=radsum-r;
    const KK_FLOAT reff=ri*rj/radsum;
    const KK_FLOAT contact_radius=sqrt(delta*reff);
    const KK_FLOAT vrx=vxi-v(j,0), vry=vyi-v(j,1), vrz=vzi-v(j,2);
    const KK_FLOAT vnnr=vrx*nx+vry*ny+vrz*nz;
    const KK_FLOAT vtx=vrx-vnnr*nx, vty=vry-vnnr*ny, vtz=vrz-vnnr*nz;
    const KK_FLOAT wrx=ri*wxi+rj*omega(j,0);
    const KK_FLOAT wry=ri*wyi+rj*omega(j,1);
    const KK_FLOAT wrz=ri*wzi+rj*omega(j,2);
    const KK_FLOAT vtrx=vtx-(wry*nz-wrz*ny);
    const KK_FLOAT vtry=vty-(wrz*nx-wrx*nz);
    const KK_FLOAT vtrz=vtz-(wrx*ny-wry*nx);
    const KK_FLOAT vrel=sqrt(vtrx*vtrx+vtry*vtry+vtrz*vtrz);

    KK_FLOAT meff=mi*rmass[j]/(mi+rmass[j]);
    if (maski & freeze_group_bit) meff=rmass[j];
    if (mask[j] & freeze_group_bit) meff=mi;

    KK_FLOAT fnormal;
    if (d_model_i(model,MI_NORMAL)==N_HOOKE)
      fnormal=d_model_f(model,MF_NORMAL_K)*delta;
    else
      fnormal=d_model_f(model,MF_NORMAL_K)*contact_radius*delta;

    KK_FLOAT damp_pref;
    const int damping=d_model_i(model,MI_DAMPING);
    if (damping==D_VELOCITY) damp_pref=d_model_f(model,MF_DAMP);
    else if (damping==D_MASS_VELOCITY) damp_pref=d_model_f(model,MF_DAMP)*meff;
    else if (damping==D_VISCOELASTIC)
      damp_pref=d_model_f(model,MF_DAMP)*meff*contact_radius;
    else {
      const KK_FLOAT arg=delta>0 ? meff*fnormal/delta : 0;
      damp_pref=d_model_f(model,MF_DAMP)*sqrt(arg>0?arg:0);
    }
    KK_FLOAT fntot=fnormal-damp_pref*vnnr;
    if (d_model_i(model,MI_LIMIT_DAMPING) && fntot<0) fntot=0;
    const KK_FLOAT fscrit=d_model_f(model,MF_MU)*fabs(fntot);
    const int tangential=d_model_i(model,MI_TANGENTIAL);
    const KK_FLOAT tdamp=d_model_f(model,MF_XT)*damp_pref;
    KK_FLOAT fsx=0,fsy=0,fsz=0;

    if (tangential==T_LINEAR_NOHISTORY) {
      const KK_FLOAT ft=(vrel>0) ? (fscrit<tdamp*vrel?fscrit:tdamp*vrel)/vrel : 0;
      fsx=-ft*vtrx; fsy=-ft*vtry; fsz=-ft*vtrz;
    } else {
      const int hi=d_model_i(model,MI_HISTORY_INDEX)+size_history_kk*jj;
      KK_FLOAT hx=d_firsthistory(i,hi), hy=d_firsthistory(i,hi+1), hz=d_firsthistory(i,hi+2);
      const KK_FLOAT k0=d_model_f(model,MF_TANGENTIAL_K);
      KK_FLOAT kscaled=k0;
      if (tangential==T_MINDLIN || tangential==T_MINDLIN_CLASSIC)
        kscaled*=contact_radius;
      if (HISTORYUPDATE) {
        if (tangential==T_LINEAR_CLASSIC || tangential==T_MINDLIN_CLASSIC) {
          hx+=dt_kk*vtrx; hy+=dt_kk*vtry; hz+=dt_kk*vtrz;
          const KK_FLOAT proj=hx*nx+hy*ny+hz*nz;
          hx-=proj*nx; hy-=proj*ny; hz-=proj*nz;
        } else {
          const KK_FLOAT proj=hx*nx+hy*ny+hz*nz;
          const KK_FLOAT hmag2=hx*hx+hy*hy+hz*hz;
          const KK_FLOAT tang2=hmag2-proj*proj;
          if (tang2>0 && fabs(proj)*kscaled>1.0e-10*fscrit) {
            const KK_FLOAT scale=sqrt(hmag2/tang2);
            hx=(hx-proj*nx)*scale; hy=(hy-proj*ny)*scale; hz=(hz-proj*nz)*scale;
          }
          if (d_model_i(model,MI_MINDLIN_FORCE)) {
            hx-=kscaled*dt_kk*vtrx; hy-=kscaled*dt_kk*vtry; hz-=kscaled*dt_kk*vtrz;
          } else {
            hx+=dt_kk*vtrx; hy+=dt_kk*vtry; hz+=dt_kk*vtrz;
          }
        }
      }
      if (d_model_i(model,MI_MINDLIN_FORCE)) { fsx=hx; fsy=hy; fsz=hz; }
      else { fsx=-kscaled*hx; fsy=-kscaled*hy; fsz=-kscaled*hz; }
      const KK_FLOAT fdx=-tdamp*vtrx, fdy=-tdamp*vtry, fdz=-tdamp*vtrz;
      fsx+=fdx; fsy+=fdy; fsz+=fdz;
      const KK_FLOAT fsmag=sqrt(fsx*fsx+fsy*fsy+fsz*fsz);
      const KK_FLOAT hmag=sqrt(hx*hx+hy*hy+hz*hz);
      if (fsmag>fscrit) {
        if (hmag>0) {
          const KK_FLOAT scale=fscrit/fsmag;
          fsx*=scale; fsy*=scale; fsz*=scale;
          hx=fsx-fdx; hy=fsy-fdy; hz=fsz-fdz;
          if (!d_model_i(model,MI_MINDLIN_FORCE)) {
            hx/=-kscaled; hy/=-kscaled; hz/=-kscaled;
          }
        } else fsx=fsy=fsz=0;
      }
      if (HISTORYUPDATE) {
        d_firsthistory(i,hi)=hx; d_firsthistory(i,hi+1)=hy; d_firsthistory(i,hi+2)=hz;
      }
    }

    KK_FLOAT fx=(fntot*nx+fsx)*factor;
    KK_FLOAT fy=(fntot*ny+fsy)*factor;
    KK_FLOAT fz=(fntot*nz+fsz)*factor;
    fxi+=fx; fyi+=fy; fzi+=fz;
    const KK_FLOAT cx=ny*fsz-nz*fsy;
    const KK_FLOAT cy=nz*fsx-nx*fsz;
    const KK_FLOAT cz=nx*fsy-ny*fsx;
    const KK_FLOAT di=(ri-0.5*delta)*factor;
    const KK_FLOAT dj=(rj-0.5*delta)*factor;
    txi-=di*cx; tyi-=di*cy; tzi-=di*cz;
    if (NEWTON_PAIR || j<nlocal) {
      a_f(j,0)-=fx; a_f(j,1)-=fy; a_f(j,2)-=fz;
      a_torque(j,0)-=dj*cx; a_torque(j,1)-=dj*cy; a_torque(j,2)-=dj*cz;
    }
    if (VFLAG) ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,fx,fy,fz,dx,dy,dz);
  }
  a_f(i,0)+=fxi; a_f(i,1)+=fyi; a_f(i,2)+=fzi;
  a_torque(i,0)+=txi; a_torque(i,1)+=tyi; a_torque(i,2)+=tzi;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::operator()(
    TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,HISTORYUPDATE> tag,
    const int ii) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,VFLAG,HISTORYUPDATE>(tag,ii,ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz, KK_FLOAT dx, KK_FLOAT dy, KK_FLOAT dz) const
{
  using VView = Kokkos::View<KK_ACC_FLOAT*[6],typename DAT::t_kkacc_1d_6::array_layout,
    typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>>;
  VView va=d_vatom;
  if (!(vflag_global || vflag_atom)) return;
  const KK_ACC_FLOAT vv[6]={(KK_ACC_FLOAT)(dx*fx),(KK_ACC_FLOAT)(dy*fy),
    (KK_ACC_FLOAT)(dz*fz),(KK_ACC_FLOAT)(dx*fy),(KK_ACC_FLOAT)(dx*fz),
    (KK_ACC_FLOAT)(dy*fz)};
  if (vflag_global) {
    if (NEWTON_PAIR) for (int k=0;k<6;k++) ev.v[k]+=vv[k];
    else {
      if (i<nlocal) for (int k=0;k<6;k++) ev.v[k]+=0.5*vv[k];
      if (j<nlocal) for (int k=0;k<6;k++) ev.v[k]+=0.5*vv[k];
    }
  }
  if (vflag_atom) {
    if (NEWTON_PAIR || i<nlocal) for (int k=0;k<6;k++) va(i,k)+=0.5*vv[k];
    if (NEWTON_PAIR || j<nlocal) for (int k=0;k<6;k++) va(j,k)+=0.5*vv[k];
  }
}

namespace LAMMPS_NS {
template class PairGranularKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGranularKokkos<LMPHostType>;
#endif
}
