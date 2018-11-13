#include<classes.hpp>

KOKKOS_FUNCTION
Foo::Foo() {
  val = 0;
} 

KOKKOS_FUNCTION
Foo_1::Foo_1() {
  val = 1;
}

KOKKOS_FUNCTION
int Foo_1::value() {
  return val;  
}

KOKKOS_FUNCTION
Foo_2::Foo_2() {
  val = 2;
}

KOKKOS_FUNCTION
int Foo_2::value() {
  return val;  
}
