// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_EXPERIMENTAL_ERROR_REPORTER_HPP
#define KOKKOS_EXPERIMENTAL_ERROR_REPORTER_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#endif

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>

namespace Kokkos {
namespace Experimental {
template <typename ReportType,
          typename DeviceType = typename DefaultExecutionSpace::device_type>
class ErrorReporter {
 public:
  using report_type     = ReportType;
  using device_type     = DeviceType;
  using execution_space = typename device_type::execution_space;

  ErrorReporter(const std::string &label, int max_results)
      : m_numReportsAttempted(label + "::m_numReportsAttempted"),
        m_reports(label + "::m_reports", max_results),
        m_reporters(label + "::m_reporters", max_results) {
    clear();
  }

  ErrorReporter(int max_results)
      : ErrorReporter("ErrorReporter", max_results) {}

  int capacity() const { return m_reports.extent(0); }

  int num_reports() const {
    return std::clamp(num_report_attempts(), 0, capacity());
  }

  int num_report_attempts() const {
    int value;
    Kokkos::deep_copy(value, m_numReportsAttempted);
    return value;
  }

  auto get_reports() const {
    int num_reps = num_reports();
    std::vector<int> reporters_out(num_reps);
    std::vector<report_type> reports_out(num_reps);

    if (num_reps > 0) {
      Kokkos::View<int *, Kokkos::HostSpace> h_reporters(reporters_out.data(),
                                                         num_reps);
      Kokkos::View<report_type *, Kokkos::HostSpace> h_reports(
          reports_out.data(), num_reps);

      Kokkos::deep_copy(
          h_reporters, Kokkos::subview(m_reporters, Kokkos::pair{0, num_reps}));
      Kokkos::deep_copy(h_reports,
                        Kokkos::subview(m_reports, Kokkos::pair{0, num_reps}));
    }
    return std::pair{std::move(reporters_out), std::move(reports_out)};
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use capacity() instead")
  int getCapacity() const { return capacity(); }

  KOKKOS_DEPRECATED_WITH_COMMENT("Use num_reports() instead")
  int getNumReports() const { return num_reports(); }

  KOKKOS_DEPRECATED_WITH_COMMENT("Use num_report_attempts() instead")
  int getNumReportAttempts() const { return num_report_attempts(); }

  KOKKOS_DEPRECATED_WITH_COMMENT("Use get_reports() instead")
  void getReports(std::vector<int> &reporters_out,
                  std::vector<report_type> &reports_out);
  KOKKOS_DEPRECATED_WITH_COMMENT("Use get_reports() instead")
  void getReports(
      typename Kokkos::View<int *, typename DeviceType::execution_space>::
          host_mirror_type &reporters_out,
      typename Kokkos::View<
          report_type *, typename DeviceType::execution_space>::host_mirror_type
          &reports_out);
#endif

  bool full() const { return (num_report_attempts() >= capacity()); }

  void clear() const { Kokkos::deep_copy(m_numReportsAttempted, 0); }

  // This function keeps reports up to new_size alive
  // It may lose the information on attempted reports
  void resize(const size_t new_size) {
    // We have to reset the attempts so we don't accidently
    // report more stored reports than there actually are
    // after growing capacity.
    int num_reps = num_report_attempts();
    if (new_size > static_cast<size_t>(capacity()) && num_reps > capacity())
      Kokkos::deep_copy(m_numReportsAttempted, num_reports());

    Kokkos::resize(m_reports, new_size);
    Kokkos::resize(m_reporters, new_size);
  }

  KOKKOS_INLINE_FUNCTION
  bool add_report(int reporter_id, report_type report) const {
    int idx = Kokkos::atomic_fetch_inc(&m_numReportsAttempted());

    if (idx >= 0 && (idx < m_reports.extent_int(0))) {
      m_reporters(idx) = reporter_id;
      m_reports(idx)   = report;
      return true;
    } else {
      return false;
    }
  }

 private:
  Kokkos::View<int, device_type> m_numReportsAttempted;
  Kokkos::View<report_type *, device_type> m_reports;
  Kokkos::View<int *, device_type> m_reporters;
};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::getReports(
    std::vector<int> &reporters_out, std::vector<report_type> &reports_out) {
  reporters_out.clear();
  reports_out.clear();
  int num_reps = num_reports();

  if (num_reps > 0) {
    reporters_out.resize(num_reps);
    reports_out.resize(num_reps);

    Kokkos::View<int *, Kokkos::HostSpace> h_reporters(reporters_out.data(),
                                                       num_reps);
    Kokkos::View<report_type *, Kokkos::HostSpace> h_reports(reports_out.data(),
                                                             num_reps);

    Kokkos::deep_copy(h_reporters,
                      Kokkos::subview(m_reporters, Kokkos::pair{0, num_reps}));
    Kokkos::deep_copy(h_reports,
                      Kokkos::subview(m_reports, Kokkos::pair{0, num_reps}));
  }
}

template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::getReports(
    typename Kokkos::View<int *, typename DeviceType::execution_space>::
        host_mirror_type &reporters_out,
    typename Kokkos::View<report_type *, typename DeviceType::execution_space>::
        host_mirror_type &reports_out) {
  int num_reps  = num_reports();
  reporters_out = typename Kokkos::View<int *, DeviceType>::host_mirror_type(
      "ErrorReport::reporters_out", num_reps);
  reports_out =
      typename Kokkos::View<report_type *, DeviceType>::host_mirror_type(
          "ErrorReport::reports_out", num_reps);

  if (num_reps > 0) {
    Kokkos::deep_copy(reporters_out,
                      Kokkos::subview(m_reporters, Kokkos::pair{0, num_reps}));
    Kokkos::deep_copy(reports_out,
                      Kokkos::subview(m_reports, Kokkos::pair{0, num_reps}));
  }
}
#endif

}  // namespace Experimental
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#endif
#endif
