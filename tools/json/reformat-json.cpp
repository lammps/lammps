
#include "json.h"

#include <cerrno>
#include <cstdio>
#include <string>
#include <exception>

using json = LAMMPS_NS::json;

static void usage(const char *exe)
{
  printf("Usage: %s <indent-width> <json-file-1> [<json-file-2> ...]\n", exe);
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    usage(argv[0]);
    return 1;
  }

  int indent = 2;
  try {
    indent = std::stoi(argv[1]);
  } catch (std::exception &e) {
    printf("Error: %s: first argument not an integer\n", e.what());
    usage(argv[0]);
    return 1;
  }

  // loop over files

  for (int i=2; i < argc; ++i) {
    std::string file = argv[i];
    std::string backup = file + ".bak";

    FILE *fp = fopen(file.c_str(), "r");
    if (!fp) {
      printf("Cannot open file %s for reading: %s\n", file.c_str(), strerror(errno));
      usage(argv[0]);
      return 2;
    }

    try {
      auto jsondata = json::parse(fp, nullptr, true, true);
      fclose(fp);
      if (rename(file.c_str(), backup.c_str())) {
        printf("Cannot create backup for file %s: %s\n", file.c_str(), strerror(errno));
        return 3;
      }

      fp = fopen(file.c_str(), "w");
      if (!fp) {
        printf("Cannot open file %s for writing: %s\n", file.c_str(), strerror(errno));
        return 4;
      }
      std::string data = jsondata.dump(indent);
      data += '\n';
      fputs(data.c_str(), fp);
      fclose(fp);
    } catch (std::exception &e) {
      printf("%s: %s\nSkipping file...\n", argv[i], e.what());
      fclose(fp);
    }
  }
  return 0;
}
