#ifndef ROTD_MATH_HH
#define ROTD_MATH_HH

#include "tmatrix.hh"

namespace Limits {
  //
  // maximal value of exponent function argument
  //
  double exp_arg_max ();
}

/*************************************************************************
 *************************** Quaternions *********************************
 *************************************************************************/

void quat_product (const double*, const double*, double*);
bool quat2mat (const double*, TMatrix<double,3>&);
void euler2quat      (const double*, char,double*);
void euler2mat       (const double*, char, TMatrix<double,3>&);

void polar2cart (const double*, double*);

void polar2av   (const double*, const double*, double*);
void polar2lv   (const double*, const double*, double*);

void euler2mf_av     (const double*, char, const double*, double*);
void euler2lf_av     (const double*, char, const double*, double*);

void euler_d2euler_m (const double*, char, const double*, const double*, double*);
void euler_m2euler_d (const double*, char, const double*, const double*, double*);

/*********************** general purpose functions ***********************/

double power (double, int); // x**n
double gamma_2 (int); // gamma(n/2)
int find_min (const double*, int);

/********************** Vector Operations ******************************/

double normalize (double*, int);
double orthogonalize (double*, const double*, int);
double vdistance (const double*, const double*, int);
double vlength (const double*, int);
double vdot (const double*, const double*, int);

/******************************************
 ************ 3-D Real Vector *************
 ******************************************/

class D3
{
  double _data [3];

public:

  D3 () {}
  D3 (const double*);
  double& operator [] (int i) { return _data[i]; }
  double operator [] (int i) const { return _data[i]; }
  operator double* () { return _data; }
  operator const double* () const { return _data; }

  D3 operator+ (const D3&) const;
  D3 operator- (const D3&) const;
  D3 operator* (double) const;
  D3 operator/ (double) const;
  D3& operator+= (const D3&);
  D3& operator-= (const D3&);
  D3& operator*= (double);
  D3& operator/= (double);

  friend D3 operator* (double, const D3&);
  friend double vdot (const D3&, const D3&);
};

D3 operator* (double, const D3&);
double vdot (const D3&, const D3&);
D3 vprod (const D3&, const D3&);
void vprod (const D3&, const D3&, D3&);

/*************** Nonequidistant Set of  Values *************************/
 
class Var_array
{
  double _start;
  double _step;
  double _incr;
  bool is_incr;
  int _size;
  int _thresh;

  Var_array (const Var_array&);

public:

  Var_array () {}

  Var_array (double, double, double, int);
  int size () const { return _size; }
  //double start () const { return _start; }
  //double step (int) const;
  double operator[] (int i) const;
};

/*****************************************
 ************ Spline fitting *************
 *****************************************/

class Spline
{
  int dim;    // spline size
  double* kn; // spline knots
  double* bc; // spline coefficients
  mutable int inbv; // slatec internal parameter

  Spline (const Spline&);
  Spline& operator= (const Spline&);

public:

  Spline (const char*);
  Spline (const double*, const double*, int);
  Spline (const Var_array&, const double*);
  ~Spline ();

  double fit (double x, int i) const; // evaluate i-th derivative at x
};

class Grid {// arithmetical progression
  double _min;
  int    _num;
  double _step;
public:
  Grid() : _min(0.), _num(0), _step(0.) {}
  Grid(double m, int n, double s) : _min(m), _num(n), _step(s) {}
  void set(double m, int n, double s) { _min = m; _num = n; _step = s; }
  int size () const { return _num; }
  double step () const { return _step; }
  double start () const { return _min; }
  double operator[] (int i) const { return _min + (double)i * _step; }
};

#endif
