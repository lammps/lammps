// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_EXPERIMENTAL_ERROR_REPORTER_HPP
#define KOKKOS_TEST_EXPERIMENTAL_ERROR_REPORTER_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.error_reporter;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_ErrorReporter.hpp>
#endif

namespace Test {

// Just save the data in the report.  Informative text goes in the
// operator<<(..).
template <typename DataType1, typename DataType2, typename DataType3>
struct ThreeValReport {
  DataType1 m_data1;
  DataType2 m_data2;
  DataType3 m_data3;
};

template <typename DataType1, typename DataType2, typename DataType3>
std::ostream &operator<<(
    std::ostream &os,
    const ThreeValReport<DataType1, DataType2, DataType3> &val) {
  return os << "{" << val.m_data1 << " " << val.m_data2 << " " << val.m_data3
            << "}";
}

template <typename ReportType>
void checkReportersAndReportsAgree(const std::vector<int> &reporters,
                                   const std::vector<ReportType> &reports) {
  for (size_t i = 0; i < reports.size(); ++i) {
    EXPECT_EQ(1, reporters[i] % 2);
    EXPECT_EQ(reporters[i], reports[i].m_data1);
  }
}

template <typename DeviceType>
struct ErrorReporterDriverBase {
  using report_type = ThreeValReport<int, int, double>;
  using error_reporter_type =
      Kokkos::Experimental::ErrorReporter<report_type, DeviceType>;

  error_reporter_type m_errorReporter;

  ErrorReporterDriverBase(int reporter_capacity, int /*test_size*/)
      : m_errorReporter(reporter_capacity) {}

  KOKKOS_INLINE_FUNCTION bool error_condition(const int work_idx) const {
    return (work_idx % 2 != 0);
  }

  void check_expectations(int reporter_capacity, int test_size) {
    using namespace std;
    int num_reported = m_errorReporter.num_reports();
    int num_attempts = m_errorReporter.num_report_attempts();
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif
    EXPECT_EQ(num_reported, m_errorReporter.getNumReports());
    EXPECT_EQ(num_attempts, m_errorReporter.getNumReportAttempts());
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif
#endif

    int expected_num_reports = min(reporter_capacity, test_size / 2);
    EXPECT_EQ(expected_num_reports, num_reported);
    EXPECT_EQ(test_size / 2, num_attempts);

    bool expect_full   = (reporter_capacity <= (test_size / 2));
    bool reported_full = m_errorReporter.full();
    EXPECT_EQ(expect_full, reported_full);
  }
};

template <typename ErrorReporterDriverType>
void TestErrorReporter() {
  using tester_type = ErrorReporterDriverType;

  std::vector<int> reporters;
  std::vector<typename tester_type::report_type> reports;

  tester_type test1(100, 10);

  std::tie(reporters, reports) = test1.m_errorReporter.get_reports();
  checkReportersAndReportsAgree(reporters, reports);

  tester_type test2(10, 100);
  auto [reporters2, reports2] = test2.m_errorReporter.get_reports();
  checkReportersAndReportsAgree(reporters2, reports2);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
  KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif
  test2.m_errorReporter.getReports(reporters, reports);
  checkReportersAndReportsAgree(reporters, reports);

  typename Kokkos::View<int *,
                        typename ErrorReporterDriverType::execution_space>::
      host_mirror_type view_reporters;
  typename Kokkos::View<typename tester_type::report_type *,
                        typename ErrorReporterDriverType::execution_space>::
      host_mirror_type view_reports;
  test2.m_errorReporter.getReports(view_reporters, view_reports);

  int num_reports = view_reporters.extent(0);
  reporters.clear();
  reports.clear();
  reporters.reserve(num_reports);
  reports.reserve(num_reports);

  for (int i = 0; i < num_reports; ++i) {
    reporters.push_back(view_reporters(i));
    reports.push_back(view_reports(i));
  }
  checkReportersAndReportsAgree(reporters, reports);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
  KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif
#endif
}

template <typename DeviceType>
struct ErrorReporterDriver : public ErrorReporterDriverBase<DeviceType> {
  using driver_base = ErrorReporterDriverBase<DeviceType>;
  using execution_space =
      typename driver_base::error_reporter_type::execution_space;

  ErrorReporterDriver(int reporter_capacity, int test_size)
      : driver_base(reporter_capacity, test_size) {
    EXPECT_EQ(driver_base::m_errorReporter.capacity(), reporter_capacity);
    EXPECT_EQ(driver_base::m_errorReporter.num_reports(), 0);
    EXPECT_EQ(driver_base::m_errorReporter.num_report_attempts(), 0);

    execute(reporter_capacity, test_size);

    // Test that clear() and resize() work across memory spaces.
    if (reporter_capacity < test_size) {
      driver_base::m_errorReporter.clear();
      EXPECT_EQ(driver_base::m_errorReporter.capacity(), reporter_capacity);
      EXPECT_EQ(driver_base::m_errorReporter.num_reports(), 0);
      EXPECT_EQ(driver_base::m_errorReporter.num_report_attempts(), 0);

      driver_base::m_errorReporter.resize(test_size);
      EXPECT_EQ(driver_base::m_errorReporter.capacity(), test_size);
      EXPECT_EQ(driver_base::m_errorReporter.num_reports(), 0);
      EXPECT_EQ(driver_base::m_errorReporter.num_report_attempts(), 0);
      execute(test_size, test_size);
    }
  }

  void execute(int reporter_capacity, int test_size) {
    Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, test_size),
                         *this);
    Kokkos::fence();
    driver_base::check_expectations(reporter_capacity, test_size);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int work_idx) const {
    if (driver_base::error_condition(work_idx)) {
      double val = Kokkos::numbers::pi * static_cast<double>(work_idx);
      typename driver_base::report_type report = {work_idx, -2 * work_idx, val};
      driver_base::m_errorReporter.add_report(work_idx, report);
    }
  }
};

template <typename DeviceType>
struct ErrorReporterDriverUseLambda
    : public ErrorReporterDriverBase<DeviceType> {
  using driver_base = ErrorReporterDriverBase<DeviceType>;
  using execution_space =
      typename driver_base::error_reporter_type::execution_space;

  ErrorReporterDriverUseLambda(int reporter_capacity, int test_size)
      : driver_base(reporter_capacity, test_size) {
    execute(reporter_capacity, test_size);
  }

  void execute(int reporter_capacity, int test_size) {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space>(0, test_size),
        // NOLINTNEXTLINE(kokkos-implicit-this-capture)
        KOKKOS_CLASS_LAMBDA(const int work_idx) {
          if (driver_base::error_condition(work_idx)) {
            double val = Kokkos::numbers::pi * static_cast<double>(work_idx);
            typename driver_base::report_type report = {work_idx, -2 * work_idx,
                                                        val};
            driver_base::m_errorReporter.add_report(work_idx, report);
          }
        });
    Kokkos::fence();
    driver_base::check_expectations(reporter_capacity, test_size);
  }
};

#ifdef KOKKOS_ENABLE_OPENMP
struct ErrorReporterDriverNativeOpenMP
    : public ErrorReporterDriverBase<Kokkos::OpenMP> {
  using driver_base = ErrorReporterDriverBase<Kokkos::OpenMP>;
  using execution_space =
      typename driver_base::error_reporter_type::execution_space;

  ErrorReporterDriverNativeOpenMP(int reporter_capacity, int test_size)
      : driver_base(reporter_capacity, test_size) {
#pragma omp parallel for
    for (int work_idx = 0; work_idx < test_size; ++work_idx) {
      if (driver_base::error_condition(work_idx)) {
        double val = Kokkos::numbers::pi * static_cast<double>(work_idx);
        typename driver_base::report_type report = {work_idx, -2 * work_idx,
                                                    val};
        driver_base::m_errorReporter.add_report(work_idx, report);
      }
    };
    driver_base::check_expectations(reporter_capacity, test_size);
  }
};
#endif

// FIXME_MSVC MSVC just gets confused when using the base class in the
// KOKKOS_CLASS_LAMBDA
#ifndef KOKKOS_COMPILER_MSVC
TEST(TEST_CATEGORY, ErrorReporterViaLambda) {
  TestErrorReporter<ErrorReporterDriverUseLambda<TEST_EXECSPACE>>();
}
#endif

TEST(TEST_CATEGORY, ErrorReporter) {
  TestErrorReporter<ErrorReporterDriver<TEST_EXECSPACE>>();
}

TEST(TEST_CATEGORY, ErrorReporter_label_ctor) {
  Kokkos::Experimental::ErrorReporter<int, TEST_EXECSPACE> logger("Reporter",
                                                                  10);
}

void ErrorReporter_test_resize() {
  Kokkos::Experimental::ErrorReporter<int, TEST_EXECSPACE> logger("Reporter",
                                                                  10);

  // produce more errors when we can store
  Kokkos::parallel_for(
      "TestErrorReporter_resize", Kokkos::RangePolicy<TEST_EXECSPACE>(0, 20),
      KOKKOS_LAMBDA(int i) { logger.add_report(i, 0); });

  ASSERT_EQ(logger.num_reports(), 10);
  ASSERT_EQ(logger.num_report_attempts(), 20);

  logger.resize(15);
  ASSERT_EQ(logger.num_reports(), 10);
  ASSERT_EQ(logger.num_report_attempts(), 10);

  logger.resize(5);
  ASSERT_EQ(logger.num_reports(), 5);
  ASSERT_EQ(logger.num_report_attempts(), 10);
}

TEST(TEST_CATEGORY, ErrorReporter_resize) { ErrorReporter_test_resize(); }

}  // namespace Test
#endif  // #ifndef KOKKOS_TEST_ERROR_REPORTING_HPP
