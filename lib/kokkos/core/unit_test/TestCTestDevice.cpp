#include <gtest/gtest.h>

namespace Kokkos {
namespace Impl {

int get_ctest_gpu(const char *local_rank_str);

}  // namespace Impl
}  // namespace Kokkos

#ifdef _WIN32
int setenv(const char *name, const char *value, int overwrite) {
  int errcode = 0;
  if (!overwrite) {
    size_t envsize = 0;
    errcode        = getenv_s(&envsize, NULL, 0, name);
    if (errcode || envsize) return errcode;
  }
  return _putenv_s(name, value);
}

int unsetenv(const char *name) { return _putenv_s(name, ""); }
#endif

// Needed because https://github.com/google/googletest/issues/952 has not been
// resolved
#define EXPECT_THROW_WITH_MESSAGE(stmt, etype, whatstring)      \
  EXPECT_THROW(                                                 \
      try { stmt; } catch (const etype &ex) {                   \
        EXPECT_EQ(std::string(ex.what()).find(whatstring), 0u); \
        throw;                                                  \
      },                                                        \
      etype)

class ctest_environment : public ::testing::Test {
 protected:
  void SetUp();
};

void ctest_environment::SetUp() {
  setenv("CTEST_KOKKOS_DEVICE_TYPE", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_COUNT", "10", 1);
  unsetenv("CTEST_RESOURCE_GROUP_0");
  setenv("CTEST_RESOURCE_GROUP_1", "threads", 1);
  setenv("CTEST_RESOURCE_GROUP_2", "threads,cores", 1);

  setenv("CTEST_RESOURCE_GROUP_3", "gpus", 1);
  unsetenv("CTEST_RESOURCE_GROUP_3_GPUS");

  setenv("CTEST_RESOURCE_GROUP_4", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_4_GPUS", "id:2", 1);

  setenv("CTEST_RESOURCE_GROUP_5", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_5_GPUS", "slots:1,id:2", 1);

  setenv("CTEST_RESOURCE_GROUP_6", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_6_GPUS", "id:2,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_7", "threads,gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_7_GPUS", "id:3,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_8", "gpus,threads", 1);
  setenv("CTEST_RESOURCE_GROUP_8_GPUS", "id:1,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_9", "cores,gpus,threads", 1);
  setenv("CTEST_RESOURCE_GROUP_9_GPUS", "id:4,slots:1", 1);
}

TEST_F(ctest_environment, no_device_type) {
  unsetenv("CTEST_KOKKOS_DEVICE_TYPE");
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("0"), 0);
}

TEST_F(ctest_environment, no_process_count) {
  unsetenv("CTEST_RESOURCE_GROUP_COUNT");
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("0"), 0);
}

TEST_F(ctest_environment, invalid_rank) {
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("10"), std::runtime_error,
      "Error: local rank 10 is outside the bounds of resource groups provided "
      "by CTest.");
}

TEST_F(ctest_environment, no_type_str) {
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("0"), std::runtime_error,
      "Error: CTEST_RESOURCE_GROUP_0 is not specified. Raised by "
      "Kokkos::Impl::get_ctest_gpu().");
}

TEST_F(ctest_environment, missing_type) {
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("1"), std::runtime_error,
      "Error: device type 'gpus' not included in CTEST_RESOURCE_GROUP_1. "
      "Raised by Kokkos::Impl::get_ctest_gpu().");
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("2"), std::runtime_error,
      "Error: device type 'gpus' not included in CTEST_RESOURCE_GROUP_2. "
      "Raised by Kokkos::Impl::get_ctest_gpu().");
}

TEST_F(ctest_environment, no_id_str) {
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("3"), std::runtime_error,
      "Error: CTEST_RESOURCE_GROUP_3_GPUS is not specified. Raised by "
      "Kokkos::Impl::get_ctest_gpu().");
}

TEST_F(ctest_environment, invalid_id_str) {
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("4"), std::runtime_error,
      "Error: invalid value of CTEST_RESOURCE_GROUP_4_GPUS: 'id:2'. Raised by "
      "Kokkos::Impl::get_ctest_gpu().");
  EXPECT_THROW_WITH_MESSAGE(
      Kokkos::Impl::get_ctest_gpu("5"), std::runtime_error,
      "Error: invalid value of CTEST_RESOURCE_GROUP_5_GPUS: 'slots:1,id:2'. "
      "Raised by Kokkos::Impl::get_ctest_gpu().");
}

TEST_F(ctest_environment, good) {
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("6"), 2);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("7"), 3);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("8"), 1);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu("9"), 4);
}
