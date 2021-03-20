
#include "gtest/gtest.h"
#include <sstream>
#include <string>

// include the implementation since ErrorStats is a standalone class
// this way we don't have to link against the style_tests and lammps libs
#include "error_stats.cpp"
#include "fmtlib_format.cpp"
#include "fmtlib_os.cpp"

TEST(ErrorStats, test)
{
    ErrorStats stats;

    stats.add(0.5);
    stats.add(0.1);
    stats.add(2.0);
    stats.add(0.3);
    stats.add(0);

    ASSERT_DOUBLE_EQ(stats.avg(), 0.58);
    ASSERT_DOUBLE_EQ(stats.dev(), 0.73047929470998685);
    ASSERT_DOUBLE_EQ(stats.max(), 2.0);
    ASSERT_EQ(stats.idx(), 3);

    std::stringstream out;
    out << stats;
    ASSERT_EQ(out.str(), "Average:  5.800e-01 StdDev:  7.305e-01 MaxErr:  2.000e+00 @ item: 3");

    stats.reset();
    ASSERT_EQ(stats.has_data(), false);
    ASSERT_DOUBLE_EQ(stats.avg(), 0.0);
    ASSERT_DOUBLE_EQ(stats.dev(), 0.0);
    ASSERT_DOUBLE_EQ(stats.max(), 0.0);
    ASSERT_EQ(stats.idx(), -1);

    stats.add(1.0);
    ASSERT_EQ(stats.has_data(), true);
    ASSERT_DOUBLE_EQ(stats.avg(), 1.0);
    ASSERT_DOUBLE_EQ(stats.dev(), 0.0);
    ASSERT_DOUBLE_EQ(stats.max(), 1.0);
    ASSERT_EQ(stats.idx(), 1);
}
