#ifndef ROTD_IO
#define ROTD_IO
#include <vector>

namespace IO {
  class Grid {
    std::vector<double> _data;
  public:
    int size () const { return _data.size(); }

    double& operator[] (int i) { return _data[i]; }
    double operator[] (int i) const { return _data[i]; }    
  };

  std::istream& operator>> (std::istream&, Grid&);
}
