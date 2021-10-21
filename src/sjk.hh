#ifndef __SJK
#define __SJK

#include "tmatrix.hh"

#include <string>
#include <iostream>

namespace Sjk {
  using namespace std;
  // Error classes
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

  // initialization
  void init (const string&);
  bool isinit ();
  
  // potential
  int pot (const double*, double&);
}

#endif
