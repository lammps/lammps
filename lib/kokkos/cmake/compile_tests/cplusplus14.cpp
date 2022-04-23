#include <type_traits>

int main() {
  // _t versions of type traits were added in C++14
  std::remove_cv_t<int> i = 0;

  return i;
}
