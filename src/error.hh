#ifndef __ROTD_ERROR__
#define __ROTD_ERROR__
#include <string>

using namespace std;

// Error classes
namespace Error {
  class Err {};

  // File handling errors
  class File_Err: public Err {};
  class Open_Err: public File_Err{}; // File opening errors
  class Form_Err: public File_Err{}; // Wrong format
  class EOF_Err:  public File_Err{}; // End-of-file encountered

  // Math errors
  class Math_Err: public Err {};
  
  // Indexing errors
  class Range_Err: public Err {};
  
  // Search errors
  class Find_Err: public Err {};
  
  // Initialization errors
  class Init_Err: public Err {};
  
  // Run errors
  class Run_Err: public Err {};

  // Floating-point exception error
  //
  class FPE {};
}
void error (const string&, int = 0);

#endif
