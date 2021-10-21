#ifndef __MULTIPOLE__
#define __MULTIPOLE__

#include <string>

namespace Multipole {
  using namespace std;
  void init (const string&);
  bool is_init ();
  const double* pos (int frag);
  double coupling ();
}

#endif
