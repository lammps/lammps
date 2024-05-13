//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_TEST_EXPERIMENTAL_ERROR_REPORTER_HPP
#define KOKKOS_TEST_EXPERIMENTAL_ERROR_REPORTER_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_ErrorReporter.hpp>

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
    int num_reported = m_errorReporter.getNumReports();
    int num_attempts = m_errorReporter.getNumReportAttempts();

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
  test1.m_errorReporter.getReports(reporters, reports);
  checkReportersAndReportsAgree(reporters, reports);

  tester_type test2(10, 100);
  test2.m_errorReporter.getReports(reporters, reports);
  checkReportersAndReportsAgree(reporters, reports);

  typename Kokkos::View<
      int *, typename ErrorReporterDriverType::execution_space>::HostMirror
      view_reporters;
  typename Kokkos::View<typename tester_type::report_type *,
                        typename ErrorReporterDriverType::execution_space>::
      HostMirror view_reports;
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
}

template <typename DeviceType>
struct ErrorReporterDriver : public ErrorReporterDriverBase<DeviceType> {
  using driver_base = ErrorReporterDriverBase<DeviceType>;
  using execution_space =
      typename driver_base::error_reporter_type::execution_space;

  ErrorReporterDriver(int reporter_capacity, int test_size)
      : driver_base(reporter_capacity, test_size) {
    execute(reporter_capacity, test_size);

    // Test that clear() and resize() work across memory spaces.
    if (reporter_capacity < test_size) {
      driver_base::m_errorReporter.clear();
      driver_base::m_errorReporter.resize(test_size);
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

#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)
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
#endif

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
#if !defined(KOKKOS_COMPILER_MSVC) && \
    (!defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA))
TEST(TEST_CATEGORY, ErrorReporterViaLambda) {
  TestErrorReporter<ErrorReporterDriverUseLambda<TEST_EXECSPACE>>();
}
#endif

TEST(TEST_CATEGORY, ErrorReporter) {
  TestErrorReporter<ErrorReporterDriver<TEST_EXECSPACE>>();
}

}  // namespace Test
#endif  // #ifndef KOKKOS_TEST_ERROR_REPORTING_HPP
