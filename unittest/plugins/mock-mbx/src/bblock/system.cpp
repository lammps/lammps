// Stand-in for MBX: every member does nothing and returns zero or empty values

#include "bblock/system.h"

namespace bblock {
System::System() = default;
System::~System() = default;

void System::AddMonomer(std::vector<double>, std::vector<std::string>, std::string, size_t, int) {}
void System::AddMolecule(std::vector<size_t>) {}
void System::SetUpFromJson(char *) {}
void System::SetUpFromJson(std::string) {}
void System::SetMPI(MPI_Comm, int, int, int) {}
void System::SetPBC(std::vector<double>) {}
void System::SetPeriodicity(bool) {}
void System::SetBoxPMElocal(std::vector<double>) {}
void System::SetRealXyz(std::vector<double>) {}
void System::SetRealspaceCutoff(double) {}
void System::SetExternalChargesAndPositions(std::vector<double>, std::vector<double>) {}
void System::SetExternalChargesAndPositions(std::vector<double>, std::vector<double>,
                                            std::vector<size_t>, std::vector<int>)
{
}
void System::SetNumDipoleHistory(size_t) {}
void System::SetDipoleHistory(size_t, std::vector<double>) {}
void System::ResetDipoleHistory() {}

double System::Get2bCutoff() { return 0.0; }
double System::Get3bCutoff() { return 0.0; }
double System::Get4bCutoff() { return 0.0; }
double System::GetRealspaceCutoff() { return 0.0; }
std::string System::GetCurrentSystemConfig() { return {}; }
std::string System::GetDipoleMethod() { return {}; }
double System::GetNumDipoleHistory() { return 0.0; }
std::vector<double> System::GetDipoleHistory(size_t) { return {}; }
void System::GetMolecularDipoles(std::vector<double> &, std::vector<double> &) {}
void System::GetEwaldParamsElectrostatics(double &, double &, size_t &) {}
void System::GetEwaldParamsDispersion(double &, double &, size_t &) {}
std::vector<int> System::GetFFTDimensionElectrostatics(int) { return {}; }
std::vector<int> System::GetFFTDimensionDispersion(int) { return {}; }
std::vector<double> System::GetRealGrads() { return {}; }
std::vector<double> System::GetExternalChargesGradients() { return {}; }
std::vector<double> System::GetVirial() { return {}; }
double System::GetPermanentElectrostaticEnergy() { return 0.0; }
double System::GetInducedElectrostaticEnergy() { return 0.0; }
std::vector<size_t> System::GetInfoElectrostaticsCounts() { return {}; }
std::vector<double> System::GetInfoElectrostaticsTimings() { return {}; }
std::vector<size_t> System::GetInfoDispersionCounts() { return {}; }
std::vector<double> System::GetInfoDispersionTimings() { return {}; }

double System::OneBodyEnergy(bool) { return 0.0; }
double System::TwoBodyEnergy(bool, bool) { return 0.0; }
double System::ThreeBodyEnergy(bool, bool) { return 0.0; }
double System::FourBodyEnergy(bool, bool) { return 0.0; }
double System::ElectrostaticsMPIlocal(bool, bool) { return 0.0; }
double System::Dispersion(bool, bool) { return 0.0; }
double System::DispersionPMElocal(bool, bool) { return 0.0; }
double System::Buckingham(bool, bool) { return 0.0; }
double System::LennardJones(bool, bool) { return 0.0; }
}    // namespace bblock
