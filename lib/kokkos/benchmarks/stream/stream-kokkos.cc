/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER
*/

#include "Kokkos_Core.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <sys/time.h>

#define STREAM_ARRAY_SIZE 100000000
#define STREAM_NTIMES     20

#define HLINE "-------------------------------------------------------------\n"

#if defined(KOKKOS_ENABLE_CUDA)
typedef Kokkos::View<double*, Kokkos::CudaSpace>::HostMirror StreamHostArray;
typedef Kokkos::View<double*, Kokkos::CudaSpace> StreamDeviceArray;
#else
typedef Kokkos::View<double*, Kokkos::HostSpace>::HostMirror StreamHostArray;
typedef Kokkos::View<double*, Kokkos::HostSpace> StreamDeviceArray;
#endif

typedef int StreamIndex;

double now() {
	struct timeval now;
	gettimeofday(&now, NULL);

	return (double) now.tv_sec + ((double) now.tv_usec * 1.0e-6);
}

void perform_copy(StreamDeviceArray& a, StreamDeviceArray& b, StreamDeviceArray& c) {

	Kokkos::parallel_for("copy", a.extent(0), KOKKOS_LAMBDA(const StreamIndex i) {
		c[i] = a[i];
	});

	Kokkos::fence();
}

void perform_scale(StreamDeviceArray& a, StreamDeviceArray& b, StreamDeviceArray& c,
       	const double scalar) {

	Kokkos::parallel_for("copy", a.extent(0), KOKKOS_LAMBDA(const StreamIndex i) {
		b[i] = scalar * c[i];
	});

	Kokkos::fence();
}

void perform_add(StreamDeviceArray& a, StreamDeviceArray& b, StreamDeviceArray& c) {
	Kokkos::parallel_for("add", a.extent(0), KOKKOS_LAMBDA(const StreamIndex i) {
                c[i] = a[i] + b[i];
        });

	Kokkos::fence();
}

void perform_triad(StreamDeviceArray& a, StreamDeviceArray& b, StreamDeviceArray& c,
	const double scalar) {

	Kokkos::parallel_for("triad", a.extent(0), KOKKOS_LAMBDA(const StreamIndex i) {
		a[i] = b[i] + scalar * c[i];
	});

	Kokkos::fence();
}

int perform_validation(StreamHostArray& a, StreamHostArray& b, StreamHostArray& c,
	const StreamIndex arraySize, const double scalar) {

	double ai = 1.0;
	double bi = 2.0;
	double ci = 0.0;

	for( StreamIndex i = 0; i < arraySize; ++i ) {
		ci = ai;
		bi = scalar * ci;
		ci = ai + bi;
		ai = bi + scalar * ci;
	};

	double aError = 0.0;
	double bError = 0.0;
	double cError = 0.0;

	for( StreamIndex i = 0; i < arraySize; ++i ) {
		aError = std::abs( a[i] - ai );
		bError = std::abs( b[i] - bi );
		cError = std::abs( c[i] - ci );
	}

	double aAvgError = aError / (double) arraySize;
	double bAvgError = bError / (double) arraySize;
	double cAvgError = cError / (double) arraySize;

	const double epsilon = 1.0e-13;
	int errorCount = 0;

	if( std::abs( aAvgError / ai ) > epsilon ) {
		fprintf(stderr, "Error: validation check on View a failed.\n");
		errorCount++;
	}

	if( std::abs( bAvgError / bi ) > epsilon ) {
		fprintf(stderr, "Error: validation check on View b failed.\n");
		errorCount++;
	}

	if( std::abs( cAvgError / ci ) > epsilon ) {
		fprintf(stderr, "Error: validation check on View c failed.\n");
		errorCount++;
	}

	if( errorCount == 0 ) {
		printf("All solutions checked and verified.\n");
	}

	return errorCount;
}

int run_benchmark() {

	printf("Reports fastest timing per kernel\n");
	printf("Creating Views...\n");

	printf("Memory Sizes:\n");
	printf("- Array Size:    %" PRIu64 "\n", static_cast<uint64_t>(STREAM_ARRAY_SIZE));
	printf("- Per Array:     %12.2f MB\n", 1.0e-6 * (double) STREAM_ARRAY_SIZE * (double) sizeof(double));
	printf("- Total:         %12.2f MB\n", 3.0e-6 * (double) STREAM_ARRAY_SIZE * (double) sizeof(double));

	printf("Benchmark kernels will be performed for %d iterations.\n", STREAM_NTIMES);

	printf(HLINE);

	StreamDeviceArray dev_a("a", STREAM_ARRAY_SIZE);
	StreamDeviceArray dev_b("b", STREAM_ARRAY_SIZE);
	StreamDeviceArray dev_c("c", STREAM_ARRAY_SIZE);

	StreamHostArray a = Kokkos::create_mirror_view(dev_a);
	StreamHostArray b = Kokkos::create_mirror_view(dev_b);
	StreamHostArray c = Kokkos::create_mirror_view(dev_c);

	const double scalar = 3.0;

	double copyTime  = std::numeric_limits<double>::max();
	double scaleTime = std::numeric_limits<double>::max();
	double addTime   = std::numeric_limits<double>::max();
	double triadTime = std::numeric_limits<double>::max();

	printf("Initializing Views...\n");

#if defined(KOKKOS_HAVE_OPENMP)
	Kokkos::parallel_for("init", Kokkos::RangePolicy<Kokkos::OpenMP>(0, STREAM_ARRAY_SIZE),
#else
	Kokkos::parallel_for("init", Kokkos::RangePolicy<Kokkos::Serial>(0, STREAM_ARRAY_SIZE),
#endif
		KOKKOS_LAMBDA(const int i) {

		a[i] = 1.0;
		b[i] = 2.0;
		c[i] = 0.0;
	});

	// Copy contents of a (from the host) to the dev_a (device)
	Kokkos::deep_copy(dev_a, a);
	Kokkos::deep_copy(dev_b, b);
	Kokkos::deep_copy(dev_c, c);

	double start;

	printf("Starting benchmarking...\n");

	for( StreamIndex k = 0; k < STREAM_NTIMES; ++k ) {
		start = now();
		perform_copy(dev_a, dev_b, dev_c);
		copyTime = std::min( copyTime, (now() - start) );

		start = now();
		perform_scale(dev_a, dev_b, dev_c, scalar);
		scaleTime = std::min( scaleTime, (now() - start) );

		start = now();
		perform_add(dev_a, dev_b, dev_c);
		addTime = std::min( addTime, (now() - start) );

		start = now();
		perform_triad(dev_a, dev_b, dev_c, scalar);
		triadTime = std::min( triadTime, (now() - start) );
	}

	Kokkos::deep_copy(a, dev_a);
	Kokkos::deep_copy(b, dev_b);
	Kokkos::deep_copy(c, dev_c);

	printf("Performing validation...\n");
	int rc = perform_validation(a, b, c, STREAM_ARRAY_SIZE, scalar);

	printf(HLINE);

	printf("Copy            %11.2f MB/s\n",
		( 1.0e-06 * 2.0 * (double) sizeof(double) * (double) STREAM_ARRAY_SIZE) / copyTime );
	printf("Scale           %11.2f MB/s\n",
		( 1.0e-06 * 2.0 * (double) sizeof(double) * (double) STREAM_ARRAY_SIZE) / scaleTime );
	printf("Add             %11.2f MB/s\n",
		( 1.0e-06 * 3.0 * (double) sizeof(double) * (double) STREAM_ARRAY_SIZE) / addTime );
	printf("Triad           %11.2f MB/s\n",
		( 1.0e-06 * 3.0 * (double) sizeof(double) * (double) STREAM_ARRAY_SIZE) / triadTime );

	printf(HLINE);

	return rc;
}

int main(int argc, char* argv[]) {

	printf(HLINE);
	printf("Kokkos STREAM Benchmark\n");
	printf(HLINE);

	Kokkos::initialize(argc, argv);
	const int rc = run_benchmark();
	Kokkos::finalize();

	return rc;
}
