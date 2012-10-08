/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef PPPM_CUDA_CU_H_
#define PPPM_CUDA_CU_H_

extern "C" void pppm_device_init(void* cu_density_brick, void* cu_vdx_brick, void* cu_vdy_brick, void* cu_vdz_brick, void* cu_density_fft, void* cu_energy, void* cu_virial
                                 , void* cu_work1, void* cu_work2, void* cu_work3, void* cu_greensfn, void* cu_fkx, void* cu_fky, void* cu_fkz, void* cu_vg
                                 , int nxlo_in, int nxhi_in, int nylo_in, int nyhi_in, int nzlo_in, int nzhi_in, int nxlo_out, int nxhi_out, int nylo_out, int nyhi_out, int nzlo_out, int nzhi_out, int nx_pppm, int ny_pppm, int nz_pppm
                                 , int cu_nxlo_fft, int cu_nxhi_fft, int cu_nylo_fft, int cu_nyhi_fft, int cu_nzlo_fft, int cu_nzhi_fft, void* cu_gf_b
                                 , double cu_qqrd2e, int cu_order, void* cu_rho_coeff, void* cu_debugdata, void* cu_density_brick_lock, int slabflag
                                );
extern "C" void pppm_device_init_setup(cuda_shared_data* sdata, PPPM_FLOAT shiftone, PPPM_FLOAT delxinv, PPPM_FLOAT delyinv, PPPM_FLOAT delzinv, int nlower, int nupper);
extern "C" void Cuda_PPPM_Setup_fkxyz_vg(int nx_pppma, int ny_pppma, int nz_pppma, PPPM_FLOAT unitkx, PPPM_FLOAT unitky, PPPM_FLOAT unitkz, PPPM_FLOAT g_ewald);
extern "C" void Cuda_PPPM_setup_greensfn(int nx_pppma, int ny_pppma, int nz_pppma, PPPM_FLOAT unitkx, PPPM_FLOAT unitky, PPPM_FLOAT unitkz, PPPM_FLOAT g_ewald,
    int nbx, int nby, int nbz, PPPM_FLOAT xprd, PPPM_FLOAT yprd, PPPM_FLOAT zprd_slab);

extern "C" void pppm_device_update(cuda_shared_data* sdata, void* cu_part2grid, int nlocala, int nmaxa);
extern "C" void pppm_update_nlocal(int nlocala);
extern "C" void poisson_scale(int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_xgrad(int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_ygrad(int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_zgrad(int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_vdx_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_vdy_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_vdz_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppm, int ny_pppm, int nz_pppm);
extern "C" void poisson_energy(int nxlo_fft, int nxhi_fft, int nylo_fft, int nyhi_fft, int nzlo_fft, int nzhi_fft, int vflag);
extern "C" ENERGY_FLOAT sum_energy(void* cu_virial, void* cu_energy, int nx_pppma, int ny_pppma, int nz_pppma, int vflag, ENERGY_FLOAT* cpu_virial);
extern "C" int cuda_particle_map(cuda_shared_data* sdata, void* flag);
extern "C" void cuda_make_rho(cuda_shared_data* sdata, void* flag, PPPM_FLOAT* cu_density_intScale, int ihi, int ilo, int jhi, int jlo, int khi, int klo, void* cu_density_brick, void* cu_density_brick_int);
extern "C" void cuda_fieldforce(cuda_shared_data* sdata, void* flag);
extern "C" double cuda_slabcorr_energy(cuda_shared_data* sdata, ENERGY_FLOAT* buf, ENERGY_FLOAT* dev_buf);
extern "C" void cuda_slabcorr_force(cuda_shared_data* sdata, F_FLOAT ffact);
extern "C" void pppm_initfftdata(cuda_shared_data* sdata, PPPM_FLOAT* in, FFT_FLOAT* out);
#endif /*PPPM_CUDA_CU_H_*/
