// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_ErrorReporter.hpp>
#include <Kokkos_Random.hpp>

// Kokkos ErrorReporter can be used to record a certain
// number of errors up to a point for debugging purposes.
// The main benefit of ErrorReporter is that its thread safe
// and does not require storage that depends on the concurrency
// of the architecture you are running on.

// This little example assumes you want to sort particles
// based on their position into boxes, but it will report
// if any of the particles are outside of the boxes.
int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    Kokkos::View<double*> positions("Pos", 10000);
    Kokkos::View<int*> box_id("box_id");

    // Lets produce some random positions in the range of -5 to 105
    Kokkos::Random_XorShift64_Pool<> rand_pool(103201);
    Kokkos::fill_random(positions, rand_pool, -5., 105.);

    // Now create an error reporter that can store 10 reports
    // We will simply report the position, but it could be a user
    // defined type.
    Kokkos::Experimental::ErrorReporter<double> errors("MyErrors", 10);

    // Counting how many positions fall into the 0-50 and 50-100 range
    int num_lower_box = 0;
    int num_upper_box = 0;
    Kokkos::parallel_reduce(
        "ErrorReporter Example", positions.extent(0),
        KOKKOS_LAMBDA(int i, int& count_lower, int& count_upper) {
          double pos = positions(i);
          // Check for positions outside the range first
          if (pos < 0. || pos > 100.) {
            // add_report takes an id and a payload
            // Note that we don't have to check how many reports were already
            // submitted
            errors.add_report(i, pos);
          } else if (pos < 50.)
            count_lower++;
          else
            count_upper++;
        },
        num_lower_box, num_upper_box);

    // Lets report results
    printf(
        "Of %i particles %i fall into the lower box, and %i into the upper "
        "box\n",
        positions.extent_int(0), num_lower_box, num_upper_box);

    // Lets report errors
    printf(
        "There were %i particles outside of the valid domain (0 - 100). Here "
        "are the first %i:\n",
        errors.num_report_attempts(), errors.num_reports());

    // Using structured bindings to get the reporter ids and reports
    auto [reporter_ids, reports] = errors.get_reports();
    for (int e = 0; e < errors.num_reports(); e++)
      printf("%i %lf\n", reporter_ids[e], reports[e]);
  }
  Kokkos::finalize();
}
