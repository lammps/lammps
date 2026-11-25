/* Stand-in for the header of the MBX class bblock::System (src/bblock/system.h of MBX), only
   for testing that the MBX plugin of LAMMPS configures, compiles, links and loads.  It declares
   only the members that the MBX package uses, with the signatures of MBX 1.4.0, and none of
   them does any work.  A member that the package starts to use must be added here. */

#ifndef MOCK_MBX_SYSTEM_H
#define MOCK_MBX_SYSTEM_H

// the standard headers that the real header includes, which LAMMPS may rely on
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(HAVE_MPI) && HAVE_MPI > 0
#include <mpi.h>
#elif !defined(MPI_VERSION) && !defined(MPI_Comm)
typedef int MPI_Comm;
#endif

namespace bblock {
class System {
 public:
  System();
  ~System();
  static constexpr const char *get_mbx_version() { return "0.0.0 (stand-in)"; }

  void AddMonomer(std::vector<double> xyz, std::vector<std::string> atoms, std::string id,
                  size_t islocal = 1, int tag = 0);
  void AddMolecule(std::vector<size_t> molec);
  void SetUpFromJson(char *json_file = 0);
  void SetUpFromJson(std::string json_text);
  void SetMPI(MPI_Comm comm, int nx, int ny, int nz);
  void SetPBC(std::vector<double> box = {});
  void SetPeriodicity(bool periodic);
  void SetBoxPMElocal(std::vector<double> box);
  void SetRealXyz(std::vector<double> xyz);
  void SetRealspaceCutoff(double cutoff_realspace);
  void SetExternalChargesAndPositions(std::vector<double> chg, std::vector<double> xyz);
  void SetExternalChargesAndPositions(std::vector<double> chg, std::vector<double> xyz,
                                      std::vector<size_t> islocal, std::vector<int> tag);
  void SetNumDipoleHistory(size_t num_hist);
  void SetDipoleHistory(size_t indx, std::vector<double> mu_hist);
  void ResetDipoleHistory();

  double Get2bCutoff();
  double Get3bCutoff();
  double Get4bCutoff();
  double GetRealspaceCutoff();
  std::string GetCurrentSystemConfig();
  std::string GetDipoleMethod();
  double GetNumDipoleHistory();
  std::vector<double> GetDipoleHistory(size_t indx);
  void GetMolecularDipoles(std::vector<double> &mu_perm, std::vector<double> &mu_ind);
  void GetEwaldParamsElectrostatics(double &alpha, double &grid_density, size_t &spline_order);
  void GetEwaldParamsDispersion(double &alpha, double &grid_density, size_t &spline_order);
  std::vector<int> GetFFTDimensionElectrostatics(int box_id = 0);
  std::vector<int> GetFFTDimensionDispersion(int box_id = 0);
  std::vector<double> GetRealGrads();
  std::vector<double> GetExternalChargesGradients();
  std::vector<double> GetVirial();
  double GetPermanentElectrostaticEnergy();
  double GetInducedElectrostaticEnergy();
  std::vector<size_t> GetInfoElectrostaticsCounts();
  std::vector<double> GetInfoElectrostaticsTimings();
  std::vector<size_t> GetInfoDispersionCounts();
  std::vector<double> GetInfoDispersionTimings();

  double OneBodyEnergy(bool do_grads);
  double TwoBodyEnergy(bool do_grads, bool use_ghost = 0);
  double ThreeBodyEnergy(bool do_grads, bool use_ghost = 0);
  double FourBodyEnergy(bool do_grads, bool use_ghost = 0);
  double ElectrostaticsMPIlocal(bool do_grads, bool use_ghost = 0);
  double Dispersion(bool do_grads, bool use_ghost = 0);
  double DispersionPMElocal(bool do_grads, bool use_ghost = 0);
  double Buckingham(bool do_grads, bool use_ghost = 0);
  double LennardJones(bool do_grads, bool use_ghost = 0);
};
}    // namespace bblock

#endif
