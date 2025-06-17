
#include "json.h"

#include <cerrno>
#include <cstdio>
#include <string>

using json = LAMMPS_NS::json;

int main(int argc, char **argv)
{
  if (argc < 3) {
    printf("Usage: %s <indent-width> <json-file-1> [<json-file-2> ...]\n", argv[0]);
    return 1;
  }

  int indent = std::stoi(argv[1]);

  // loop over files

  for (int i=2; i < argc; ++i) {
    std::string file = argv[i];
    std::string backup = file + ".bak";

    FILE *fp = fopen(file.c_str(), "r");
    if (!fp) {
      printf("Cannot open file %s for reading: %s\n", file.c_str(), strerror(errno));
      return 2;
    }

    auto jsondata = json::parse(fp, nullptr, true, true);
    fclose(fp);

    rename(file.c_str(), backup.c_str());

    fp = fopen(file.c_str(), "w");
    if (!fp) {
      printf("Cannot open file %s for writing: %s\n", file.c_str(), strerror(errno));
      return 3;
    }
    std::string data = jsondata.dump(indent);
    data += '\n';
    fputs(data.c_str(), fp);
    fclose(fp);
  }
  return 0;
}
