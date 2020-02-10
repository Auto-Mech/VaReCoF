#ifndef RELAX_HH
#define RELAX_HH

namespace GeomRelax {

  //quantum chemistry method
  double (*qc_method) (const double* pos, double* force);
  void optimize(NormalModeCoor*);
}

#endif
