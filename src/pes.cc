#include "pes.hh"

namespace PES {// Potential Energy Surface namespace
  int _pot_size = 1;
  pot_f pot = 0;
}

int PES::size() { return _pot_size; }
void PES::set_size (int s) { _pot_size = s; }

