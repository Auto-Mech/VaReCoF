#ifndef INTERFLUX_HH
#define INTERFLUX_HH

#include "tmatrix.hh"
#include "math.hh"

// flux interpolation
class InterFlux {
  Spline* _s;
  double _ezero;
  double _emin, _emax;
  double _nmin, _amin;
  double _nmax, _amax;

  bool _empty;

  // no copies
  InterFlux (const InterFlux&);
  InterFlux& operator= (const InterFlux&);

public:
  InterFlux (const Array<double>& ener, const Array<double>& flux, double vmin);

  ~InterFlux () { if(!_empty) delete _s; }
  double operator() (double) const; 
};

#endif
