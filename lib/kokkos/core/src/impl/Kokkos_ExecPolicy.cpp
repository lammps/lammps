#include <Kokkos_Core.hpp>
namespace Kokkos {
namespace Impl {
    PerTeamValue::PerTeamValue(int arg):value(arg) {}

    PerThreadValue::PerThreadValue(int arg):value(arg) {}
}

Impl::PerTeamValue PerTeam(const int& arg)
{
  return Impl::PerTeamValue(arg);
}

Impl::PerThreadValue PerThread(const int& arg)
{
  return Impl::PerThreadValue(arg);
}

}
