// unit tests for getting LAMMPS configuration through the Fortran wrapper

#include "info.h"
#include "lammps.h"
#include "library.h"

#ifdef LMP_PLUGIN
#include "plugin.h"
#endif

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
int f_lammps_version();
int f_lammps_os_info(const char *);
int f_lammps_mpi_support();
int f_lammps_gzip_support();
int f_lammps_png_support();
int f_lammps_jpeg_support();
int f_lammps_ffmpeg_support();
int f_lammps_has_exceptions();
int f_lammps_has_package(const char *);
int f_lammps_package_count();
char *f_lammps_package_name(int);
char *f_lammps_installed_packages(int);
int f_lammps_config_accelerator(const char *, const char *, const char *);
int f_lammps_has_gpu();
char *f_lammps_get_gpu_info(size_t);
int f_lammps_has_style(const char *, const char *);
int f_lammps_style_count(const char *);
char *f_lammps_style_name(const char *, int);
void f_setup_has_id();
int f_lammps_has_id(const char *, const char *);
int f_lammps_id_count(const char *);
char *f_lammps_id_name(const char *, int);
int f_lammps_plugin_count();
int f_lammps_plugin_name();
}
namespace LAMMPS_NS {

using ::testing::ContainsRegex;

class LAMMPS_configuration : public ::testing::Test {
protected:
    LAMMPS *lmp;
    const std::vector<std::string> style_category = {
        "atom",     "integrate", "minimize", "pair",    "bond",   "angle", "dihedral",
        "improper", "kspace",    "fix",      "compute", "region", "dump",  "command"};

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
    // clang-format off
    std::vector<std::string> pkg_name = {
        "ADIOS", "AMOEBA", "ASPHERE", "ATC", "AWPMD", "BOCS", "BODY", "BPM", "BROWNIAN", "CG-DNA",
        "CG-SPICA", "CLASS2", "COLLOID", "COLVARS", "COMPRESS", "CORESHELL", "DEPEND", "DIELECTRIC",
        "DIFFRACTION", "DIPOLE", "DPD-BASIC", "DPD-MESO", "DPD-REACT", "DPD-SMOOTH", "DRUDE",
        "EFF", "ELECTRODE", "EXTRA-COMPUTE", "EXTRA-DUMP", "EXTRA-FIX", "EXTRA-MOLECULE",
        "EXTRA-PAIR", "FEP", "GPU", "GRANULAR", "H5MD", "INTEL", "INTEL/TEST", "INTERLAYER", "KIM",
        "KOKKOS", "KSPACE", "LATBOLTZ", "LATTE", "LEPTON", "MACHDYN", "MAKE", "MAKE/MACHINES",
        "MAKE/MINE", "MAKE/OPTIONS", "MANIFOLD", "MANYBODY", "MC", "MDI", "MEAM", "MESONT", "MGPT",
        "MISC", "ML-HDNNP", "ML-IAP", "ML-PACE", "ML-QUIP", "ML-RANN", "ML-SNAP", "MOFFF",
        "MOLECULE", "MOLFILE", "MPIIO", "MSCG", "NETCDF", "OPENMP", "OPT", "ORIENT", "PERI",
        "PHONON", "PLUGIN", "PLUMED", "POEMS", "PTM", "PYTHON", "QEQ", "QMMM", "QTB", "REACTION",
        "REAXFF", "REPLICA", "RIGID", "SCAFACOS", "SHOCK", "SMTBQ", "SPH", "SPIN", "SRD", "STUBS",
        "TALLY", "UEF", "VORONOI", "VTK", "YAFF"
    };
    // clang-format on

    for (std::size_t i = 0; i < pkg_name.size(); i++)
        EXPECT_EQ(f_lammps_has_package(pkg_name[i].c_str()), Info::has_package(pkg_name[i]));
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
    while (LAMMPS::installed_packages[i] != nullptr) {
        char *name = f_lammps_package_name(i + 1); // 1 in Fortran is 0 in C
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
    for (int i = 0; i < npackages; i++) {
        package_name = LAMMPS::installed_packages[i];
        pkg          = f_lammps_installed_packages(i + 1);
        EXPECT_STREQ(package_name, pkg);
        if (pkg) std::free(pkg);
    }
};

TEST_F(LAMMPS_configuration, config_accelerator)
{
    const std::vector<std::string> package           = {"GPU", "KOKKOS", "INTEL", "OPENMP"};
    const std::vector<std::string> category          = {"api", "precision"};
    const std::vector<std::string> setting_api       = {"cuda",   "hip",    "phi",   "pthreads",
                                                        "opencl", "openmp", "serial"};
    const std::vector<std::string> setting_precision = {"double", "mixed", "single"};

    for (const auto &p : package) {
        for (const auto &c : category) {
            if (c == "api") {
                for (const auto &s : setting_api)
                    EXPECT_EQ(f_lammps_config_accelerator(p.c_str(), c.c_str(), s.c_str()),
                              Info::has_accelerator_feature(p, c, s));
            } else if (c == "precision") {
                for (const auto &s : setting_precision)
                    EXPECT_EQ(f_lammps_config_accelerator(p.c_str(), c.c_str(), s.c_str()),
                              Info::has_accelerator_feature(p, c, s));
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
    auto cpp_info  = Info::get_gpu_device_info();
    char *f_string = f_lammps_get_gpu_info(cpp_info.size());

    EXPECT_STREQ(utils::trim(f_string).c_str(), utils::trim(cpp_info).c_str());
    std::free(f_string);

    if (cpp_info.size() > 80) {
        f_string = f_lammps_get_gpu_info(80);
        cpp_info.resize(80);
        EXPECT_STREQ(utils::trim(f_string).c_str(), utils::trim(cpp_info).c_str());
        std::free(f_string);
    }
};

TEST_F(LAMMPS_configuration, has_style)
{
    Info info(lmp);
    for (const auto &c : style_category) {
        std::vector<std::string> name = info.get_available_styles(c);
        for (const auto &s : name)
            EXPECT_EQ(f_lammps_has_style(c.c_str(), s.c_str()), info.has_style(c, s));
    }
    EXPECT_EQ(f_lammps_has_style("atom", "none"), 0);
};

TEST_F(LAMMPS_configuration, style_count)
{
    Info info(lmp);
    for (const auto &c : style_category)
        EXPECT_EQ(f_lammps_style_count(c.c_str()), info.get_available_styles(c.c_str()).size());
};

TEST_F(LAMMPS_configuration, style_name)
{
    char *buffer;
    Info info(lmp);
    for (const auto &c : style_category) {
        int nnames  = f_lammps_style_count(c.c_str());
        auto styles = info.get_available_styles(c);
        for (int i = 0; i < nnames; i++) {
            buffer = f_lammps_style_name(c.c_str(), i + 1);
            EXPECT_STREQ(buffer, styles[i].c_str());
            std::free(buffer);
        }
    }
};

TEST_F(LAMMPS_configuration, has_id)
{
    f_setup_has_id();
    EXPECT_EQ(f_lammps_has_id("compute", "com"), 0);
    EXPECT_EQ(f_lammps_has_id("compute", "COM"), 1);
    EXPECT_EQ(f_lammps_has_id("dump", "atom"), 0);
    EXPECT_EQ(f_lammps_has_id("dump", "1"), 1);
    EXPECT_EQ(f_lammps_has_id("fix", "nve"), 0);
    EXPECT_EQ(f_lammps_has_id("fix", "1"), 1);
    EXPECT_EQ(f_lammps_has_id("group", "one"), 1);
    EXPECT_EQ(f_lammps_has_id("group", "all"), 1);
    // Skip this one (we're testing thoroughly enough; molecules require files)
    // EXPECT_EQ(f_lammps_has_id("molecule")
    EXPECT_EQ(f_lammps_has_id("region", "simbox"), 1);
    EXPECT_EQ(f_lammps_has_id("region", "box"), 0);
    EXPECT_EQ(f_lammps_has_id("variable", "pi"), 1);
    EXPECT_EQ(f_lammps_has_id("variable", "PI"), 0);
};

TEST_F(LAMMPS_configuration, id_count)
{
    f_setup_has_id();
    // computes: thermo_temp, thermo_press, thermo_pe, COM
    EXPECT_EQ(f_lammps_id_count("compute"), 4);
    EXPECT_EQ(f_lammps_id_count("dump"), 1);  // only the one we defined
    EXPECT_EQ(f_lammps_id_count("fix"), 1);   // only the one we defined
    EXPECT_EQ(f_lammps_id_count("group"), 2); // "one" and "all"
    EXPECT_EQ(f_lammps_id_count("molecule"), 0);
    EXPECT_EQ(f_lammps_id_count("region"), 1);   // onle the one we created
    EXPECT_EQ(f_lammps_id_count("variable"), 3); // "zpos", "x", and "pi"
};

TEST_F(LAMMPS_configuration, id_name)
{
    f_setup_has_id();
    char *name;
    int nnames = f_lammps_id_count("compute");
    EXPECT_EQ(nnames, 4);

    name = f_lammps_id_name("compute", 1);
    EXPECT_STREQ(name, "thermo_temp");
    std::free(name);

    name = f_lammps_id_name("compute", 2);
    EXPECT_STREQ(name, "thermo_press");
    std::free(name);

    name = f_lammps_id_name("compute", 3);
    EXPECT_STREQ(name, "thermo_pe");
    std::free(name);

    name = f_lammps_id_name("compute", 4);
    EXPECT_STREQ(name, "COM");
    std::free(name);

    name = f_lammps_id_name("dump", 1);
    EXPECT_STREQ(name, "1");
    std::free(name);

    name = f_lammps_id_name("fix", 1);
    EXPECT_STREQ(name, "1");
    std::free(name);

    name = f_lammps_id_name("group", 1);
    EXPECT_STREQ(name, "all");
    std::free(name);

    name = f_lammps_id_name("group", 2);
    EXPECT_STREQ(name, "one");
    std::free(name);

    name = f_lammps_id_name("region", 1);
    EXPECT_STREQ(name, "simbox");
    std::free(name);

    name = f_lammps_id_name("variable", 1);
    EXPECT_STREQ(name, "zpos");
    std::free(name);

    name = f_lammps_id_name("variable", 2);
    EXPECT_STREQ(name, "x");
    std::free(name);

    name = f_lammps_id_name("variable", 3);
    EXPECT_STREQ(name, "pi");
    std::free(name);
};

TEST_F(LAMMPS_configuration, plugins)
{
#ifndef LMP_PLUGIN
    GTEST_SKIP();
#else
    int nplugins = f_lammps_plugin_count();
    for (int n = 0; n < nplugins; n++) {
        lammpsplugin_t *plugin = plugin_get_info(n);
        EXPECT_EQ(f_lammps_plugin_name(n + 1, plugin->style, plugin->name), 1);
    }
#endif
};
} // namespace LAMMPS_NS
