#include <cmath>
#include "math.hh"
#include "error.hh"
#include <fstream>
#include <vector>
#include "slatec.hh"
#include <iostream>

using namespace std;

namespace Limits {
  //
  // maximal value of exponent argument
  //
  double _exp_arg_max = 100.;
  
  double exp_arg_max () { return _exp_arg_max; }
}


/*************************************************************************
 *************************** Quaternions *********************************
 *************************************************************************/

void quat_product (const double* q1, const double* q2, double* q) {
  // q1*q2
  const double* v1 = q1 + 1;
  const double* v2 = q2 + 1;
  double*       v  = q  + 1;

  vector_product(v1, v2, v);
  for(int i = 0; i < 3; ++i)
    v[i] += q1[0] * v2[i] + q2[0] * v1[i];

  q[0] = q1[0] * q2[0] - vdot(v1, v2, 3);
}

bool quat2mat (const double* quat, TMatrix<double, 3>& mat)
{// quaternion-to-matrix transformation
  const double qmax = 10.0;
  const double qmin = 0.1;

  bool is_normalized = true; // is quaternion sufficiently normalized?

  register     double    q0, q1, q2, q3, qq0, qq1, qq2, qq3, qq;
  register     double    a01, a02, a03, a12, a13, a23;
  
  q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

  qq0 = q0 * q0;  qq1 = q1 * q1;  qq2 = q2 * q2;  qq3 = q3 * q3;
 
  qq = qq0 + qq1 + qq2 + qq3;

  if (qq < qmin || qq > qmax)
    {
      cout << "quat2mat: quaternion is not normalized; qq = " << qq << endl;
      is_normalized = false;
    }

  a01 = 2.0 * q0 * q1 / qq;  a02 = 2.0 * q0 * q2 / qq;  a03 = 2.0 * q0 * q3 / qq;
  a12 = 2.0 * q1 * q2 / qq;  a13 = 2.0 * q1 * q3 / qq;
  a23 = 2.0 * q2 * q3 / qq;

  mat (0, 1) = a12 + a03;
  mat (1, 0) = a12 - a03;
  mat (0, 2) = a13 - a02;
  mat (2, 0) = a13 + a02;
  mat (1, 2) = a23 + a01;
  mat (2, 1) = a23 - a01;
 
  mat (0, 0) = (qq0 + qq1 - qq2 - qq3)/qq;
  mat (1, 1) = (qq0 - qq1 + qq2 - qq3)/qq;
  mat (2, 2) = (qq0 - qq1 - qq2 + qq3)/qq;

  return is_normalized;
}

void euler2quat (const double* euler_ang, char conv, double* q)
{// transformation from euler angles in appropriate convention to quaternion

   const double& theta = euler_ang [0];
   double phi, psi;
   switch(conv)
     {
     case 'X':
       phi   = euler_ang [1];
       psi   = euler_ang [2];
       break;
     case 'Y':
       phi   = euler_ang [1] + M_PI_2;
       psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler2quat: wrong convention");
     }

   const double two = 2.0;

   const double sum = (phi+psi)/two, dif = (phi-psi)/two;

   const double cos_theta_2 = cos(theta/two), sin_theta_2 = sin(theta/two);

   q[0] = cos_theta_2 * cos(sum);
   q[1] = sin_theta_2 * cos(dif);
   q[2] = sin_theta_2 * sin(dif);
   q[3] = cos_theta_2 * sin(sum);
}

void euler2mat (const double* euler_ang, char conv, TMatrix<double,3>& mat)
{// rotational matrix in terms of euler angles in appropriate convention
 
   const double& theta = euler_ang [0];
   double phi, psi;
   switch(conv)
     {
     case 'X':
       phi   = euler_ang [1];
       psi   = euler_ang [2];
       break;
     case 'Y':
       phi   = euler_ang [1] + M_PI_2;
       psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler2mat: wrong convention");
     }

   const double cos_theta = cos(theta), cos_phi = cos(phi), cos_psi = cos(psi);
   const double sin_theta = sin(theta), sin_phi = sin(phi), sin_psi = sin(psi);

   mat(0,0) =  cos_phi * cos_psi - sin_phi * cos_theta * sin_psi;
   mat(0,1) =  sin_phi * cos_psi + cos_phi * cos_theta * sin_psi;
   mat(0,2) =  sin_theta * sin_psi;
   mat(1,0) = -cos_phi * sin_psi - sin_phi * cos_theta * cos_psi;
   mat(1,1) = -sin_phi * sin_psi + cos_phi * cos_theta * cos_psi;
   mat(1,2) =  sin_theta * cos_psi;
   mat(2,0) =  sin_phi * sin_theta;
   mat(2,1) = -cos_phi * sin_theta;
   mat(2,2) =  cos_theta;
}

void polar2cart (const double* polar_ang, double* n)
{// polar_angles-to-directional_cosines transformation

   const double& theta = polar_ang[0];
   const double& phi   = polar_ang[1];

   const double sin_theta = sin(theta);

   n[0] = sin_theta * cos(phi);
   n[1] = sin_theta * sin(phi);
   n[2] = cos(theta);
}


void euler2mf_av(const double* euler_ang, char conv, const double* euler_der, double* av)
{// transformation of euler angles derivatives to angular velocity in molecular frame

  const double& theta = euler_ang [0];
  double psi;
  switch(conv) {

  case 'X':

    psi   = euler_ang [2];
    break;
  case 'Y':

    psi   = euler_ang [2] - M_PI_2;
    break;
  default:
    error ("euler2mf_av: wrong convention");
  }

  const double& theta_d = euler_der [0];
  const double& phi_d   = euler_der [1];
  const double& psi_d   = euler_der [2];

  const double sin_theta = sin(theta), sin_psi = sin(psi), cos_psi = cos(psi);

  av[0] = phi_d * sin_theta * sin_psi + theta_d * cos_psi;
  av[1] = phi_d * sin_theta * cos_psi - theta_d * sin_psi;
  av[2] = phi_d * cos(theta) + psi_d;
}

void euler2lf_av(const double* euler_ang, char conv, const double* euler_der, double* av)
{// transformation of euler angles derivatives to angular velocity in laboratory frame

  const double& theta = euler_ang [0];
  double phi;
  switch(conv) {

  case 'X':
    phi   = euler_ang [1];

    break;
  case 'Y':
    phi   = euler_ang [1] + M_PI_2;

    break;
  default:
    error ("euler2lf_av: wrong convention");
  }

  const double& theta_d = euler_der [0];
  const double& phi_d   = euler_der [1];
  const double& psi_d   = euler_der [2];

  const double sin_theta = sin(theta), sin_phi = sin(phi), cos_phi = cos(phi);

  av[0] = theta_d * cos_phi + psi_d * sin_theta * sin_phi;
  av[1] = theta_d * sin_phi - psi_d * sin_theta * cos_phi;
  av[2] = phi_d + psi_d * cos(theta);
}

void polar2av(const double* polar_ang, const double* polar_der, double* av)
{// transformation from polar angles derivatives to angular velocity

   const double& theta = polar_ang [0];
   const double& phi   = polar_ang [1];

   const double& theta_d = polar_der [0];
   const double& phi_d   = polar_der [1];

   const double sin_phi = sin(phi), cos_phi = cos(phi), sin_theta = sin(theta);
   const double p = sin_theta * cos(theta);
   
   av[0] = - theta_d * sin_phi - phi_d * p * cos_phi;
   av[1] =   theta_d * cos_phi - phi_d * p * sin_phi;
   av[2] =                       phi_d * sin_theta * sin_theta;
}

void polar2lv(const double* polar_ang, const double* polar_der, double* lv)
{// transformation from polar angles derivatives to linear velocity on the unity sphere

   const double& theta = polar_ang [0];
   const double& phi   = polar_ang [1];

   const double& theta_d = polar_der [0];
   const double& phi_d   = polar_der [1];

   const double sin_phi = sin(phi), cos_phi = cos(phi), 
     sin_theta = sin(theta), cos_theta = cos(theta);
   
   lv[0] =   theta_d * cos_theta * cos_phi - phi_d * sin_theta * sin_phi;
   lv[1] =   theta_d * cos_theta * sin_phi + phi_d * sin_theta * cos_phi;
   lv[2] = - theta_d * sin_theta;
}

void euler_d2euler_m(const double* euler_ang, char conv, const double* iner_mom, 
		     const double* euler_d, double* euler_m)
{// transformation from euler angles derivatives to corresponding generalized momenta

  const double& theta = euler_ang [0];
  double psi;
  switch(conv) {

  case 'X':

    psi   = euler_ang [2];
    break;
  case 'Y':

    psi   = euler_ang [2] - M_PI_2;
    break;
  default:
    error ("euler_d2euler_m: wrong convention");
  }

  const double two = 2.0;

  const double cos_2psi = cos(two * psi), sin_2psi = sin(two * psi);
  const double cos_theta = cos(theta), sin_theta = sin(theta);
  const double im_sum = (iner_mom[0] + iner_mom[1]) / two;
  const double im_dif = (iner_mom[0] - iner_mom[1]) / two;

  TMatrix<double, 3> trans;

  trans(0,0) = im_sum + im_dif*cos_2psi;
  trans(1,1) = sin_theta * sin_theta * (im_sum - im_dif * cos_2psi) +  
    iner_mom[2] * cos_theta * cos_theta;
  trans(2,2) = iner_mom[2];

  double temp = im_dif * sin_theta * sin_2psi;
  trans(0,1) = temp;
  trans(1,0) = temp;

  temp = iner_mom[2] * cos_theta;
  trans(1,2) = temp;
  trans(2,1) = temp;

  trans(0,2) = 0.0;
  trans(2,0) = 0.0;

  matrix_vector_product (trans, euler_d, euler_m);
}

void euler_m2euler_d(const double* euler_ang, char conv, const double* iner_mom, 
		     const double* euler_m, double* euler_d)
{// transformation from euler angles momenta to derivatives

   const double& theta = euler_ang [0];
   double psi;
   switch(conv)   {

     case 'X':

        psi   = euler_ang [2];
       break;
     case 'Y':

        psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler_m2euler_d: wrong convention");
     }

   const double two = 2.0;

   const double cos_2psi = cos(two * psi), sin_2psi = sin(two * psi);
   const double cos_theta = cos(theta), sin_theta = sin(theta);
   const double im_sum = (iner_mom[0] + iner_mom[1]) / two;
   const double im_dif = (iner_mom[0] - iner_mom[1]) / two;

   TMatrix<double, 3> trans;

   trans(0,0) = im_sum + im_dif * cos_2psi;
   trans(1,1) = sin_theta * sin_theta * (im_sum - im_dif * cos_2psi) +  
      iner_mom[2] * cos_theta * cos_theta;
   trans(2,2) = iner_mom[2];

   double temp = im_dif * sin_theta * sin_2psi;
   trans(0,1) = temp;
   trans(1,0) = temp;

   temp = iner_mom[2] * cos_theta;
   trans(1,2) = temp;
   trans(2,1) = temp;

   trans(0,2) = 0.0;
   trans(2,0) = 0.0;

   trans.invert();

   matrix_vector_product (trans, euler_m, euler_d);
}

/*********************** general purpose functions ***********************/

double power (double x, int n)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;
  if (n > 0)
    return x * power(x, n-1);

  return power(x, n+1) / x;
}

double gamma_2 (int n) // gamma(n/2)
{
  static const double sqrt_pi = sqrt(M_PI);
  if (n == 1)
    return sqrt_pi;
  if (n == 2)
    return 1.0;
  if (n > 2)
    return double(n-2) / 2. * gamma_2(n-2);

  return 0.;
}

int find_min (const double* x, int dim)
{
  int min_ind = 0;
  double min_val = x[0];
  for (int i = 1; i < dim; ++i)
    if (x[i] < min_val)
      {
	min_ind = i;
	min_val = x[i];
      }
  return min_ind;
}

/*********************** Vector operations *******************************/

double normalize (double* vec, int dim)
{
  double norm = 0.0, temp;
  for (int i = 0; i < dim; ++i)
    {
      temp = vec[i];
      norm += temp * temp;
    }
  norm = sqrt (norm);

  if (norm == 0.)
    return norm;

  for (int i = 0; i < dim; ++i)
    {
      vec[i] /= norm;
    }

  return norm;
}

double orthogonalize (double* vec, const double* n, int dim)
{
  double norm = 0.0, proj = 0.0, norm1 = 0.0;
  for (int i = 0; i < dim; ++i) {
    norm += n[i] * n[i];
    norm1 += vec[i] * vec[i];
    proj += n[i] * vec[i];
  }

  if(norm == 0. || norm1 == 0.)
    return 0.;

  proj /= norm;
  for (int i = 0; i < dim; ++i) {
    vec[i] -= proj * n[i];
  }

  return proj * proj * norm / norm1;
}

double vdistance (const double* vec1, const double* vec2, int dim)
{
  double dist2 = 0.0, temp;
  for (int i = 0; i < dim; ++i)
    {
      temp = vec1[i] - vec2[i];
      dist2 += temp * temp;
    }
  return sqrt(dist2);
}

double vlength (const double* vec, int dim)
{
  double norm = 0.0, temp;
  for (int i = 0; i < dim; ++i)
    {
      temp = vec[i];
      norm += temp * temp;
    }
  return sqrt (norm);
}

double vdot (const double* v1, const double* v2, int n)
{
  double res = 0.;
  for (int i = 0; i < n; ++i)
    res += *v1++ * *v2++;
  return res;
}
/******************************************
 ************ 3-D Real Vector *************
 ******************************************/

D3::D3(const double* p1)
{
  double* p = _data;
  for(int i = 0; i < 3; ++i)
    *p++ = *p1++;
}
D3 D3::operator+ (const D3& a2) const
{
  D3 res;
  double* p = res._data;
  const double* p1 = _data;
  const double* p2 = a2._data;
  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ + *p2++;
  return res;
}
D3 D3::operator- (const D3& a2) const
{
  D3 res;
  double* p = res._data;
  const double* p1 = _data;
  const double* p2 = a2._data;
  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ - *p2++;
  return res;
}
D3 D3::operator* (double d) const
{
  D3 res;
  double* p = res._data;
  const double* p1 = _data;
  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ * d;
  return res;
}
D3 D3::operator/ (double d) const
{
  D3 res;
  double* p = res._data;
  const double* p1 = _data;
  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ / d;
  return res;
}
D3& D3::operator+= (const D3& a)
{
  double* p = _data;
  const double* p1 = a._data;
  for(int i = 0; i < 3; ++i)
    *p++ += *p1++;
  return *this;
}
D3& D3::operator-= (const D3& a)
{
  double* p = _data;
  const double* p1 = a._data;
  for(int i = 0; i < 3; ++i)
    *p++ -= *p1++;
  return *this;
}
D3& D3::operator*= (double d)
{
  double* p = _data;
  for(int i = 0; i < 3; ++i)
    *p++ *= d;
  return *this;
}
D3& D3::operator/= (double d)
{
  double* p = _data;
  for(int i = 0; i < 3; ++i)
    *p++ /= d;
  return *this;
}

D3 operator* (double d, const D3& a)
{
  D3 res;
  double* p = res._data;
  const double* p1 = a._data;
  for(int i = 0; i < 3; ++i)
    *p++ = d * *p1++;
  return res;
}
double vdot (const D3& a1, const D3& a2)
{
  double res = 0.;
  const double* p1 = a1._data;
  const double* p2 = a2._data;
  for(int i = 0; i < 3; ++i)
    res += *p1++ * *p2++;
  return res;
}
D3 vprod (const D3& a1, const D3& a2)
{
  D3 res;
  res[0] = a1[1] * a2[2] - a1[2] * a2[1];
  res[1] = a1[2] * a2[0] - a1[0] * a2[2];
  res[2] = a1[0] * a2[1] - a1[1] * a2[0];
  return res;
}
void vprod (const D3& a1, const D3& a2, D3& res)
{
  res[0] = a1[1] * a2[2] - a1[2] * a2[1];
  res[1] = a1[2] * a2[0] - a1[0] * a2[2];
  res[2] = a1[0] * a2[1] - a1[1] * a2[0];
}

/********************** Var_array *********************/

Var_array::Var_array (double s, double ds, double in, int n) 
  : _start(s), _step(ds), _incr(in), is_incr(true) , _size(n)
{
  if (_incr <= 1.00001) {
    _incr = 1.0;
    is_incr = false;
  }

  _thresh = int(-_start/_step);
  if(_thresh < 0) _thresh = 0;
}

/*
double Var_array::step (int i) const
{ 
  if (is_incr && i > _thresh)
    return _step * power(_incr, i-_thresh);
  else
    return _step; 
}
*/

double Var_array::operator[] (int i) const 
{ 
  if (is_incr && i > _thresh)
    return _start + _step * 
      (_thresh + (power(_incr, i - _thresh) - 1.0) / (_incr - 1.0));
  else
    return _start + _step * i;
}

/*****************************************
 ************ Spline fitting *************
 *****************************************/

Spline::Spline (const char* fname) : inbv(1)
{
  ifstream from (fname);
  if (!from)
    error("Spline::Spline: cannot open file");

  vector<double> vx, vy;
  double dtemp1, dtemp2;
  while (from >> dtemp1 >> dtemp2)
    {
      vx.push_back(dtemp1);
      vy.push_back(dtemp2);
    }

  Array<double> x(vx.size()), y(vy.size());
  for (int i = 0; i < vx.size(); ++i)
    {
      x[i] = vx[i];
      y[i] = vy[i];
    }

  dim = x.size() + 2;
  kn = new double [dim + 4];
  bc = new double [dim];
  Array<double> work(5*dim);

  int n, k;
  dbint4_ (x.data(), y.data(), x.size(), 2, 2, 0., 0., 1, 
	   kn, bc, n, k, work.data());

  if (n != dim || k != 4)
    error("Spline::Spline: dbint4 failed");

}

Spline::Spline (const double* x, const double* y, int data_dim) : inbv(1)
{

  dim = data_dim + 2;
  kn = new double [dim + 4];
  bc = new double [dim];
  Array<double> work(5*dim);

  int n, k;
  dbint4_ (x, y, data_dim, 2, 2, 0., 0., 1, 
	   kn, bc, n, k, work.data());

  if (n != dim || k != 4)
    error("Spline::Spline: dbint4 failed");

}

Spline::Spline (const Var_array& vx, const double* y) : inbv(1)
{
  Array<double> x(vx.size());
  for (int i = 0; i < vx.size(); ++i)
    x[i] = vx[i];
  
  dim = vx.size() + 2;
  kn = new double [dim + 4];
  bc = new double [dim];
  Array<double> work(5*dim);

  int n, k;
  dbint4_ (x.data(), y, vx.size(), 2, 2, 0., 0., 1, 
	   kn, bc, n, k, work.data());

  if (n != dim || k != 4)
    error("Spline::Spline: dbint4 failed");

}

Spline::~Spline()
{
  delete[] kn;
  delete[] bc;
}

double Spline::fit(double x, int deriv) const
{
  static double work[12];
  return dbvalu_ (kn, bc, dim, 4, deriv, x, inbv, work);
}
