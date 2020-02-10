#ifndef LOG_HH
#define LOG_HH

#include <fstream>

namespace Log {
  enum {ERROR, WARN, NOTE, INFO, DEBUG};// logging levels

  extern std::ofstream out;
  extern std::ofstream debug;
  extern std::ofstream err; 
  extern std::ofstream global;

  bool more(int);
  int level();
  void set_level(int);
}
#endif
