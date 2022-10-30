/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXPERIMENTAL_ERROR_REPORTER_HPP
#define KOKKOS_EXPERIMENTAL_ERROR_REPORTER_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#endif

#include <vector>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_DualView.hpp>

namespace Kokkos {
namespace Experimental {

template <typename ReportType, typename DeviceType>
class ErrorReporter {
 public:
  using report_type     = ReportType;
  using device_type     = DeviceType;
  using execution_space = typename device_type::execution_space;

  ErrorReporter(int max_results)
      : m_numReportsAttempted(""),
        m_reports("", max_results),
        m_reporters("", max_results) {
    clear();
  }

  int getCapacity() const { return m_reports.h_view.extent(0); }

  int getNumReports();

  int getNumReportAttempts();

  void getReports(std::vector<int> &reporters_out,
                  std::vector<report_type> &reports_out);
  void getReports(
      typename Kokkos::View<int *,
                            typename DeviceType::execution_space>::HostMirror
          &reporters_out,
      typename Kokkos::View<report_type *,
                            typename DeviceType::execution_space>::HostMirror
          &reports_out);

  void clear();

  void resize(const size_t new_size);

  bool full() { return (getNumReportAttempts() >= getCapacity()); }

  KOKKOS_INLINE_FUNCTION
  bool add_report(int reporter_id, report_type report) const {
    int idx = Kokkos::atomic_fetch_add(&m_numReportsAttempted(), 1);

    if (idx >= 0 && (idx < static_cast<int>(m_reports.d_view.extent(0)))) {
      m_reporters.d_view(idx) = reporter_id;
      m_reports.d_view(idx)   = report;
      return true;
    } else {
      return false;
    }
  }

 private:
  using reports_view_t     = Kokkos::View<report_type *, device_type>;
  using reports_dualview_t = Kokkos::DualView<report_type *, device_type>;

  using host_mirror_space = typename reports_dualview_t::host_mirror_space;
  Kokkos::View<int, device_type> m_numReportsAttempted;
  reports_dualview_t m_reports;
  Kokkos::DualView<int *, device_type> m_reporters;
};

template <typename ReportType, typename DeviceType>
inline int ErrorReporter<ReportType, DeviceType>::getNumReports() {
  int num_reports = 0;
  Kokkos::deep_copy(num_reports, m_numReportsAttempted);
  if (num_reports > static_cast<int>(m_reports.h_view.extent(0))) {
    num_reports = m_reports.h_view.extent(0);
  }
  return num_reports;
}

template <typename ReportType, typename DeviceType>
inline int ErrorReporter<ReportType, DeviceType>::getNumReportAttempts() {
  int num_reports = 0;
  Kokkos::deep_copy(num_reports, m_numReportsAttempted);
  return num_reports;
}

template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::getReports(
    std::vector<int> &reporters_out, std::vector<report_type> &reports_out) {
  int num_reports = getNumReports();
  reporters_out.clear();
  reporters_out.reserve(num_reports);
  reports_out.clear();
  reports_out.reserve(num_reports);

  if (num_reports > 0) {
    m_reports.template sync<host_mirror_space>();
    m_reporters.template sync<host_mirror_space>();

    for (int i = 0; i < num_reports; ++i) {
      reporters_out.push_back(m_reporters.h_view(i));
      reports_out.push_back(m_reports.h_view(i));
    }
  }
}

template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::getReports(
    typename Kokkos::View<
        int *, typename DeviceType::execution_space>::HostMirror &reporters_out,
    typename Kokkos::View<report_type *,
                          typename DeviceType::execution_space>::HostMirror
        &reports_out) {
  int num_reports = getNumReports();
  reporters_out   = typename Kokkos::View<int *, DeviceType>::HostMirror(
      "ErrorReport::reporters_out", num_reports);
  reports_out = typename Kokkos::View<report_type *, DeviceType>::HostMirror(
      "ErrorReport::reports_out", num_reports);

  if (num_reports > 0) {
    m_reports.template sync<host_mirror_space>();
    m_reporters.template sync<host_mirror_space>();

    for (int i = 0; i < num_reports; ++i) {
      reporters_out(i) = m_reporters.h_view(i);
      reports_out(i)   = m_reports.h_view(i);
    }
  }
}

template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::clear() {
  int num_reports = 0;
  Kokkos::deep_copy(m_numReportsAttempted, num_reports);
  m_reports.template modify<execution_space>();
  m_reporters.template modify<execution_space>();
}

template <typename ReportType, typename DeviceType>
void ErrorReporter<ReportType, DeviceType>::resize(const size_t new_size) {
  m_reports.resize(new_size);
  m_reporters.resize(new_size);
  typename DeviceType::execution_space().fence(
      "Kokkos::Experimental::ErrorReporter::resize: fence after resizing");
}

}  // namespace Experimental
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ERRORREPORTER
#endif
#endif
