#include <string>

struct Up {};

int main(int argc, char* argv[]) {
  if (argc != 4 || std::string(argv[1]) != "one" ||
      std::string(argv[2]) != "2" || std::string(argv[3]) != "THREE") {
    throw Up{};
  }
  return 0;
}
