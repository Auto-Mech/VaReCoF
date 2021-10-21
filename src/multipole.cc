#include <fstream>
#include <iostream>
#include <cstdlib>

#include "multipole.hh"
#include "math.hh"

namespace Multipole {
  bool _is_init = false;
  double _coupling = 1.;    // interaction
  double _pos[2][3]; // direction of the dipole / quadrupole vector
}

bool Multipole::is_init() { return _is_init; }

double Multipole::coupling () { return _coupling; }

const double* Multipole::pos(int frag) 
{
  if(frag != 0 && frag != 1) {
    cout << "Multipole::pos: out of range\n";
    exit(1);
  }
  return _pos[frag]; 
}

void Multipole::init (const string& fname) 
{
  static const char funame[] = "Multipole::init: "; 


  string stemp;

  if(_is_init) {
    cout << funame << "double initialization\n";
    exit(1);
  }
  _is_init = true;

  for(int frag = 0; frag < 2; ++frag) {
    _pos[frag][0] = 1.;
    _pos[frag][1] = 0.;
    _pos[frag][2] = 0.;
  }

  std::ifstream from(fname.c_str());
  if(!from) {
    cout << funame << "cannot open " << fname << " file\n";
    exit(1);
  }

  while(from >> stemp) 
    if(stemp == "coupling") {
      from >> _coupling;
      if(!from) {
	cout << funame << "initialization file is corrupted\n";
	exit(1);
      }
    }
    else if(stemp == "position") {
      int frag;
      from >> frag;
      if(!from || frag != 0 && frag != 1) {
	cout << funame << "initialization file is corrupted\n";
	exit(1);
      }

      for(int i = 0; i < 3; ++i)
	from >> _pos[frag][i];
      if(!from) {
	cout << funame << "initialization file is corrupted\n";
	exit(1);
      }
      ::normalize(_pos[frag], 3);
    }
}
