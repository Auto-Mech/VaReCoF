#include "log.hh"

namespace Log {
  std::ofstream out;
  std::ofstream debug;
  std::ofstream err; 
  std::ofstream global;

  int _level = INFO;
  int level() { return _level; }
  void set_level(int l) { _level = l; }
  bool more(int l) { return _level >= l; }
}
