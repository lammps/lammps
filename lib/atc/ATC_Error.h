// ATC_Error : a base class for atom-continuum errors

#ifndef ATC_ERROR
#define ATC_ERROR

#include <string>

// the following two convert __LINE__ to a string
#define STRING2(x) #x
#define STRING(x) STRING2(x)
// prints file and line number for error messages
#define ERROR(x) __FILE__":"STRING(__LINE__)" "x
//#define FILELINE __FILE__+to_string(__LINE__)
#define FILELINE __FILE__

#define ERROR_FOR_BACKTRACE
#define HACK(l,m)






namespace ATC {
  /**
   *  @class  ATC_Error
   *  @brief  Base class for throwing run-time errors with descriptions
   */

class ATC_Error {

 public:
  // constructor
  ATC_Error(std::string errorDescription)
  {
    errorDescription_ = "ERROR: " + errorDescription;
    ERROR_FOR_BACKTRACE
  };

  ATC_Error(std::string location, std::string errorDescription)
  {
    errorDescription_ = "ERROR: " + location + ": "+ errorDescription;
    ERROR_FOR_BACKTRACE
  };

  std::string error_description() {
    return errorDescription_;
  };

 private:
  // string describing the type of error
  std::string errorDescription_;
};

}
#endif
