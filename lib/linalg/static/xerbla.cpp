
#include <exception>
#include <string>

class LinalgException : public std::exception {
  std::string message;

 public:
  LinalgException() = delete;

  explicit LinalgException(const std::string &msg) { message = msg; }
  const char *what() const noexcept override { return message.c_str(); }
};

extern "C" {

#include "lmp_f2c.h"

integer xerbla_(const char *srname, integer *info)
{
  std::string mesg = " ** On entry to ";
  for (int i = 0; i < 1024; ++i) {
    if ((srname[i] == '\0') || (srname[i] == ' ')) break;
    mesg.push_back(srname[i]);
  }
  mesg += " parameter number " + std::to_string(*info) + " had an illegal value\n";
  throw LinalgException(mesg);
  return 0;
}
}
