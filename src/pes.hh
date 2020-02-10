#ifndef __PES
#define __PES

#include "tmatrix.hh"
typedef void (*pot_f) (int, Array<double>&);

namespace PES {// Potential energy surface
  extern pot_f pot;
  int size();
  void set_size (int);
}

#endif
