#include <mpi.h>
#include <mutex>
#include <iostream>

#include "metatomic_timer.h"

using namespace LAMMPS_NS;

// lock to protect the other static from concurrent modification
static std::mutex METATOMIC_TIMER_MUTEX = {};
// depth of the timer, i.e. how many different timers are alive right now
static int64_t METATOMIC_TIMER_DEPTH = -1;
// strictly increasing timer counter, to know if a new one was created inside
// the scope of the current one.
static uint64_t METATOMIC_TIMER_COUNTER = 0;
// Is profiling enabled?
static bool METATOMIC_TIMER_ENABLED = false;


void MetatomicTimer::enable(bool toggle) {
    auto guard_ = std::lock_guard(METATOMIC_TIMER_MUTEX);

    METATOMIC_TIMER_ENABLED = toggle;
}


MetatomicTimer::MetatomicTimer(std::string name):
    enabled_(false),
    name_(std::move(name))
{
    auto guard_ = std::lock_guard(METATOMIC_TIMER_MUTEX);
    if (METATOMIC_TIMER_ENABLED) {
        METATOMIC_TIMER_DEPTH += 1;
        METATOMIC_TIMER_COUNTER += 1;

        this->enabled_ = true;
        this->starting_counter_ = METATOMIC_TIMER_COUNTER;
        this->start_ = std::chrono::high_resolution_clock::now();
        auto indent = std::string(METATOMIC_TIMER_DEPTH * 3, ' ');

        if (METATOMIC_TIMER_DEPTH == 0) {
            std::cerr << "\n";
        }

        std::cerr << "\n" << indent << this->name_ << " ...";
    }
}

MetatomicTimer::~MetatomicTimer() {
    auto guard_ = std::lock_guard(METATOMIC_TIMER_MUTEX);

    if (METATOMIC_TIMER_ENABLED && this->enabled_) {
        auto stop = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start_).count();

        if (METATOMIC_TIMER_COUNTER != starting_counter_) {
            auto indent = std::string(METATOMIC_TIMER_DEPTH * 3, ' ');
            std::cerr << "\n" << indent << this->name_;
        }

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::cerr << " took " << elapsed / 1e6 << "ms (rank " << rank << ")" << std::flush;
        METATOMIC_TIMER_DEPTH -= 1;
    }
}
