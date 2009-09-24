// ATC_Error : a base class for atom-continuum errors

#ifndef ATC_ERROR
#define ATC_ERROR

#include <string>


namespace ATC {

class ATC_Error {

 public:

  // constructor
  ATC_Error(int procID, std::string errorDescription) {
    procID_ = procID;
    errorDescription_ = errorDescription;
  };

//ATC_Error(std::string errorDescription) {
//  procID_ = LammpsInterface::instance()->comm_rank(); 
//  errorDescription_ = errorDescription;
//};

  // data access
  virtual int get_id() {
    return procID_;
  };

  virtual std::string get_error_description() {
    return errorDescription_;
  };


 private:
  
  // id of the processor reporting the error
  int procID_;

  // string describing the type of error
  std::string errorDescription_;

};

}
#endif
