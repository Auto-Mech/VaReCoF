#include "interflux.hh"
#include "units.hh"

#include <cmath>
#include <cstdlib>

InterFlux::InterFlux(const Array<double>& ener, const Array<double>& flux, double vmin)
{
  const char funame [] = "InterFlux::InterFlux:";

  if(ener.size() != flux.size()) {
    std::cout << funame << "flux array and energy grid dimensions mismatch\n";
    exit(1);
  }

  int i0;
  for(i0 = 0; i0 < ener.size(); ++i0)
    if(ener[i0] > vmin && flux[i0] > 0.)
      break;

  Array<double> x(ener.size() - i0);
  Array<double> y(ener.size() - i0);

  if(x.size() < 2) {
    _empty = true;
    return;
  }
  else
    _empty = false;

  if(!i0 || ener[i0 - 1] <= vmin)
    _ezero = vmin;
  else
    _ezero = ener[i0 - 1];

  _emin = ener[i0];
  _emax = ener[ener.size() - 1];

  for(int i = 0; i < x.size(); ++i) {
    x[i] = std::log(ener[i0 + i] - _ezero);
    if(flux[i0 + i] <= 0.) {
      std::cout << funame << "non-positive flux";
      exit(1);
    }
    y[i] = std::log(flux[i0 + i]);
  }

  _nmin = (y[1] - y[0]) / (x[1] - x[0]);
  _amin = std::exp(y[0] - x[0] * _nmin);

  _nmax = (y[x.size() - 1] - y[x.size() - 2]) / (x[x.size() - 1] - x[x.size() - 2]);
  _amax = std::exp(y[x.size() - 1] - x[x.size() - 1] * _nmax);

  _s = new Spline(x.data(), y.data(), x.size());
}

double InterFlux::operator() (double ener) const
{
  if(_empty)
    return 0.;

  if(ener <= _ezero)
    return 0.;

  if(ener <= _emin)
    return _amin * std::pow(ener - _ezero, _nmin);

  if(ener >= _emax)
    return _amax * std::pow(ener - _ezero, _nmax);

  return std::exp(_s->fit(std::log(ener - _ezero), 0));
}

