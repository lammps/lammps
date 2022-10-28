// unit tests for getting LAMMPS configuration through the Fortran wrapper

#include "lammps.h"
#include "library.h"
#include "info.h"

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
int f_lammps_version();
int f_lammps_os_info(const char*);
int f_lammps_mpi_support();
int f_lammps_gzip_support();
int f_lammps_png_support();
int f_lammps_jpeg_support();
int f_lammps_ffmpeg_support();
int f_lammps_has_exceptions();
int f_lammps_has_package(const char*);
int f_lammps_package_count();
char* f_lammps_package_name(int);
char* f_lammps_installed_packages(int);
int f_lammps_config_accelerator(const char*, const char*, const char*);
int f_lammps_has_gpu();
char* f_lammps_get_gpu_info(size_t);
int f_lammps_has_style(const char*, const char*);
int f_lammps_style_count(const char*);
int f_lammps_style_name();
}
namespace LAMMPS_NS {

using ::testing::ContainsRegex;

class LAMMPS_configuration : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        ::testing::internal::CaptureStdout();
        lmp = (LAMMPS *)f_lammps_with_args();

        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }

    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_configuration, version)
{
    EXPECT_LT(20200917, f_lammps_version());
    EXPECT_EQ(lmp->num_ver, f_lammps_version());
};

TEST_F(LAMMPS_configuration, os_info)
{
    std::string str = platform::os_info();
    EXPECT_EQ(f_lammps_os_info(str.c_str()), 1);
};

TEST_F(LAMMPS_configuration, MPI_support)
{
#ifdef MPI_STUBS
    EXPECT_EQ(f_lammps_mpi_support(), 0);
#else
    EXPECT_EQ(f_lammps_mpi_support(), 1);
#endif
};

TEST_F(LAMMPS_configuration, gzip_support)
{
    EXPECT_EQ(f_lammps_gzip_support(), Info::has_gzip_support());
}

TEST_F(LAMMPS_configuration, png_support)
{
    EXPECT_EQ(f_lammps_png_support(), Info::has_png_support());
}

TEST_F(LAMMPS_configuration, jpeg_support)
{
    EXPECT_EQ(f_lammps_jpeg_support(), Info::has_jpeg_support());
}

TEST_F(LAMMPS_configuration, ffmpeg_support)
{
    EXPECT_EQ(f_lammps_ffmpeg_support(), Info::has_ffmpeg_support());
}

TEST_F(LAMMPS_configuration, has_exceptions)
{
    EXPECT_EQ(f_lammps_has_exceptions(), Info::has_exceptions());
}

TEST_F(LAMMPS_configuration, has_package)
{
    std::vector<std::string> pkg_name = {"ADIOS","ASPHERE","ATC","AWPMD","BOCS","BODY","BPM",
      "BROWNIAN","CG-DNA","CLASS2","COLLOID","COLVARS","COMPRESS","CORESHELL","DEPEND",
      "DIELECTRIC","DIFFRACTION","DIPOLE","DPD-BASIC","DPD-MESO","DPD-REACT","DPD-SMOOTH","DRUDE",
      "EFF","ELECTRODE","EXTRA-COMPUTE","EXTRA-DUMP","EXTRA-FIX","EXTRA-MOLECULE","EXTRA-PAIR",
      "FEP","GPU","GRANULAR","H5MD","INTEL","INTEL/TEST","INTERLAYER","KIM","KOKKOS","KSPACE",
      "LATBOLTZ","LATTE","MACHDYN","MAKE","MAKE/MACHINES","MAKE/MINE","MAKE/OPTIONS","MANIFOLD",
      "MANYBODY","MC","MDI","MEAM","MESONT","MGPT","MISC","ML-HDNNP","ML-IAP","ML-PACE","ML-QUIP",
      "ML-RANN","ML-SNAP","MOFFF","MOLECULE","MOLFILE","MPIIO","MSCG","NETCDF","OPENMP","OPT",
      "ORIENT","PERI","PHONON","PLUGIN","PLUMED","POEMS","PTM","PYTHON","QEQ","QMMM","QTB",
      "REACTION","REAXFF","REPLICA","RIGID","SCAFACOS","SHOCK","SMTBQ","SPH","SPIN","SRD","STUBS",
      "TALLY","UEF","VORONOI","VTK","YAFF","CG-SPICA","AMOEBA"};

   for (int i = 0; i < pkg_name.size(); i++)
      EXPECT_EQ(f_lammps_has_package(pkg_name[i].c_str()),
         Info::has_package(pkg_name[i]));
}

TEST_F(LAMMPS_configuration, package_count)
{
   int package_count = 0;
   while (LAMMPS::installed_packages[package_count] != nullptr)
      package_count++;
   EXPECT_EQ(f_lammps_package_count(), package_count);
};

TEST_F(LAMMPS_configuration, package_name)
{
   int i = 0;
   while (LAMMPS::installed_packages[i] != nullptr)
   {
      char* name = f_lammps_package_name(i+1); // 1 in Fortran is 0 in C
      if (name)
         EXPECT_STREQ(LAMMPS::installed_packages[i], name);
      else
         EXPECT_STREQ(LAMMPS::installed_packages[i], "NOT FOUND");
      std::free(name);
      i++;
   }
};

TEST_F(LAMMPS_configuration, installed_packages)
{
    const char *package_name;
    int npackages = lammps_config_package_count();
    char *pkg;
    for (int i=0; i < npackages; i++) {
      package_name = LAMMPS::installed_packages[i];
      pkg = f_lammps_installed_packages(i+1);
      EXPECT_STREQ(package_name, pkg);
      if (pkg) std::free(pkg);
    }
};

TEST_F(LAMMPS_configuration, config_accelerator)
{
    const int npackage = 4;
    const int ncategory = 2;
    const int nsetting_api = 7;
    const int nsetting_precision = 3;
    const std::string package[] = {"GPU","KOKKOS","INTEL","OPENMP"};
    const std::string category[] = {"api","precision"};
    const std::string setting_api[] = {"cuda","hip","phi","pthreads","opencl",
      "openmp","serial"};
    const std::string setting_precision[] = {"double","mixed","single"};

    for (int p=0; p < npackage; p++)
    {
      for (int c=0; c < ncategory; c++)
      {
        if (category[c] == "api")
        {
          for (int s=0; s < nsetting_api; s++)
            EXPECT_EQ(f_lammps_config_accelerator(package[p].c_str(), category[c].c_str(),
                setting_api[s].c_str()),
              Info::has_accelerator_feature(package[p], category[c], setting_api[s]));
        }
        else if (category[c] == "precision")
        {
          for (int s=0; s < nsetting_precision; s++)
            EXPECT_EQ(f_lammps_config_accelerator(package[p].c_str(), category[c].c_str(),
                setting_precision[s].c_str()),
              Info::has_accelerator_feature(package[p],category[c],setting_precision[s]));
        }
      }
   }
};

TEST_F(LAMMPS_configuration, has_gpu)
{
    EXPECT_EQ(Info::has_gpu_device(), f_lammps_has_gpu());
};

TEST_F(LAMMPS_configuration, get_gpu_info)
{
    if (!Info::has_gpu_device()) GTEST_SKIP();
    size_t n;
    std::string cpp_info = Info::get_gpu_device_info();
    n = cpp_info.size();
    char* f_string;
    f_string = f_lammps_get_gpu_info(n);
    EXPECT_STREQ(f_string, cpp_info.c_str());
    std::free(f_string);

    if (n > 80)
    {
      f_string = f_lammps_get_gpu_info(80);
      cpp_info.resize(80);
      EXPECT_STREQ(f_string, cpp_info.c_str());
      std::free(f_string);
    }
};

TEST_F(LAMMPS_configuration, has_style)
{
    std::vector<std::string> category = {"atom","integrate","minimize","pair",
      "bond","angle","dihedral","improper","kspace","fix","compute","region",
      "dump","command"};
    Info info(lmp);
    for (int c = 0; c < category.size(); c++)
    {
        std::vector<std::string> name = info.get_available_styles(category[c]);
        for (int s = 0; s < name.size(); s++)
        {
            EXPECT_EQ(f_lammps_has_style(category[c].c_str(),
              name[s].c_str()), info.has_style(category[c], name[s]));
        }
    }
    EXPECT_EQ(f_lammps_has_style("atom","none"), 0);
};

TEST_F(LAMMPS_configuration, style_count)
{
    Info info(lmp);
    std::vector<std::string> category = {"atom","integrate","minimize","pair",
      "bond","angle","dihedral","improper","kspace","fix","compute","region",
      "dump","command"};
    for (int i = 0; i < category.size(); i++)
{
        EXPECT_EQ(f_lammps_style_count(category[i].c_str()),
          info.get_available_styles(category[i].c_str()).size());
}
};

} // namespace LAMMPS_NS
