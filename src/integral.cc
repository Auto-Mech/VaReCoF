#include <cmath>
#include "math.hh"
#include "slatec.hh"
#include "error.hh"

namespace 
{
  const double eps = 1.0e-14;
  double tol = 1.e-4;

  double global_rot_en [3];
  int global_pow_num;

  double theta_integral (double a, double b, int n)
  {
    if (a <= 0.)
      return 0.;

    if (b <= eps * a)
      {
	double res = power(a, n/2);
	if (n%2)
	  res *= sqrt(a);
	return res;
      }

    if (n < -1)
      error("theta_integral: wrong n");

    if (n == 0)
      if (b <= a)
	return 1.0;
      else
	return sqrt(a/b);

    if (n == -1)
      {
	if (b < a)
	  return asin(sqrt(b/a)) / sqrt(b);
	else
	  return M_PI_2 / sqrt(b);
      }

    double res = a * static_cast<double>(n) * theta_integral(a, b, n-2);
    double ab = a - b;
    if (ab > 0.)
      {
	double dtemp = power(ab, n/2);
	if (n%2)
	  dtemp *= sqrt(ab);
	res += dtemp;
      }
    return res / static_cast<double>(n+1);
    
  }

  extern "C" double phi_integrand (const double& y)
  {

    double y_2 = y * y;
    double y_1 = 1.0 - y_2;

    double r[2];
    r[0] = global_rot_en[1] * y_1 + global_rot_en[2] * y_2;
    r[1] = global_rot_en[1] * y_2 + global_rot_en[2] * y_1;

    double res = 0.;
    for (int i = 0; i < 2; ++i)
      res += theta_integral(1.-r[i], global_rot_en[0]-r[i], global_pow_num);

    return res / sqrt(y_1);
  }
}

double mc_stat_weight (double kin_en, double ang_mom, 
		       const double* iner_mom, int dof_num)
{
  if (iner_mom[0] > iner_mom[1] || iner_mom[1] > iner_mom[2])
    error("mc_stat_weight: inertia moments are not monotonic");
  if (iner_mom[0] <= 0.)
    error("mc_stat_weight: inertia moments are not positive");

  if (kin_en <= 0.)
    return 0.;

  global_pow_num = dof_num - 4;

  for (int i = 0; i < 3; ++i)
    global_rot_en[i] = ang_mom * ang_mom / iner_mom[i] / 2. / kin_en;

  if (global_rot_en[2] >= 1. - eps)
    return 0.;

  double en_fac = power(kin_en, global_pow_num/2);
  if (global_pow_num%2)
    en_fac *= sqrt(kin_en);

  int err;
  double res;
  double dtemp;
  if ( global_rot_en[1] <= 1.)
    {
      dgaus8_(phi_integrand, 0., M_SQRT1_2, tol, res, err);
      return M_2_PI * en_fac * res;
    }


  double y_max;
  dtemp = (1. - global_rot_en[2]) / (global_rot_en[1] - global_rot_en[2]);
  if (dtemp <= 0.5)
    {
      y_max = sqrt(dtemp);
      dgaus8_(phi_integrand, 0., y_max, tol, res, err);
    }
  else
    {
      y_max = sqrt(1. - dtemp);      
      dgaus8_(phi_integrand, 0., y_max, tol, res, err);
      dgaus8_(phi_integrand, y_max, M_SQRT1_2, tol, dtemp, err);
      res += dtemp;
    }    
  return M_2_PI * en_fac * res;
}


