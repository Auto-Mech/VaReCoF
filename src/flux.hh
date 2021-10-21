#ifndef __FLUX_
#define __FLUX_

#include <iostream>
#include <vector>
#include <map>

#include <mpi.h>

#ifdef MPIPP_H
#include "comm.hh"
#endif


#include "rotd.hh"
#include "force.hh"
#include "error.hh"
#include "math.hh"
#include "divsur.hh"
#include "random.hh"
#include "comm.hh"

namespace Statistics {

  extern int smp_out_flag;

}

/****************************************************************
 **************************** Sampling **************************
 ****************************************************************/

class  Samp : public Dynvar
{
  double wfac;      // configurational weight
  Array<double> ener;
  double rval;      // random value between 0 and 1

  mutable bool need_update;
  mutable double _weight;

public:

  Samp () : wfac(0.0), ener(PES::size()), _weight(0.0), 
    rval(1.0), need_update(true) {}

  void set_val (double w, double e) 
  {wfac = w; ener[0] = e; rval = Random::flat(); need_update = true; }

  void set_val (double w, const Array<double>& e);

  double random_value () const { return rval; }
  double energy (int=0) const;
  double weight_factor () const { return wfac; }

  double weight (double, int=0) const;
  double weight () const;

  void send (int, int) const;
  void recv (int, int);

  friend ostream& operator<< (ostream&, const Samp&);
  friend istream& operator>> (istream&, Samp&);
};

ostream& operator<< (ostream&, const Samp&);
istream& operator>> (istream&, Samp&);

/*****************************************************************************
 ******************************* Flux Base  **********************************
 *****************************************************************************/

class FluxBase {

  static bool _is_stat_init; // are static variables initialized

  static int t_size, e_size, j_size; // dimensions of flux arrays
  static double *tmpr, *ener, *amom;   // grids

  int _acct_num;// successfull potential samplings
  int _fail_num;// failed potential samplings
  int _fake_num;// dummy potential sumplings
  int _face_num;// out-of-face samplings
  int _dist_num;// distance-between-atoms-close samplings

protected:
  // the number of samplings before it is decided 
  // that there are no potential samplings at all
  static const int vol_num_max = 1000; 
  // the number of failed potential samplings 
  // before it is decided that the node is dead
  static const int fail_num_max = 10;

public:
  static bool is_stat_init(); // are static variables initialized

  static void stat_init (const std::vector<double>&, 
			 const std::vector<double>&, 
			 const std::vector<double>&);

  static int    tm_size () { return t_size; }
  static int    en_size () { return e_size; }
  static int    am_size () { return j_size; }

  static double tm_grid (int i) { return tmpr[i]; }
  static double en_grid (int i) { return ener[i]; }
  static double am_grid (int i) { return amom[i]; }

  static int smp_out_flag; // raw sampling output flag

  int acct_smp () const { return _acct_num; }
  int fail_smp () const { return _fail_num; }
  int fake_smp () const { return _fake_num; }
  int face_smp () const { return _face_num; }
  int dist_smp () const { return _dist_num; }

  int pot_smp () const { return _acct_num + _fail_num + _fake_num; }
  int vol_smp () const { return _face_num + _dist_num; }
  int tot_smp () const { return vol_smp() + pot_smp(); }

  void add_acct_smp (int n) { _acct_num += n; }
  void add_fail_smp (int n) { _fail_num += n; }
  void add_fake_smp (int n) { _fake_num += n; }
  void add_face_smp (int n) { _face_num += n; }
  void add_dist_smp (int n) { _dist_num += n; }

  FluxBase () : _acct_num(0), _fail_num(0), _fake_num(0), 
		_face_num(0), _dist_num(0) {}
  FluxBase (const Div_surf&, int);
  FluxBase& operator+= (const FluxBase&);

  bool is_pot () const;

  // MPI communication
  void send (Comm::tag_t) const;
  void recv (int, Comm::tag_t);

  friend istream& operator>> (istream&, FluxBase&);
  friend ostream& operator<< (ostream&, const FluxBase&);
};

istream& operator>> (istream&, FluxBase&);
ostream& operator<< (ostream&, const FluxBase&);

/*****************************************************************
 *********************** Thermal Flux ****************************
 *****************************************************************/

class ThermalFlux : public FluxBase {

protected:
  Array<double> t_sum; // cannonical flux
  Array<double> t_var; // flux error

  double min_en;
  Dynvar min_dv;

public:
  ThermalFlux ();

  ThermalFlux& operator+= (const ThermalFlux&);
  ThermalFlux  operator+  (const ThermalFlux&) const;
  
  double  min_energy () const { return min_en; }
  Dynvar& min_dynvar () { return min_dv; }
  const Dynvar& min_dynvar () const { return min_dv; }

  double cn_flux (int i) const { return t_sum[i]; }
  double cn_dflx (int i) const { return t_var[i]; }

  double average (int) const;// average value of the cannonical flux
  double pot_var (int) const;// potential flux variation
  double vol_var (int) const;// volume flux variation

  void normalize ();// normalize fluxes and errors

  void slave_send  () const; // mpi send
  void master_recv (int);    // mpi receive

  friend istream& operator>> (istream&, ThermalFlux&);
  friend ostream& operator<< (ostream&, const ThermalFlux&);
};// Thermal Flux

ostream& operator<< (ostream&, const ThermalFlux&);// write the raw flux data
istream& operator>> (istream&, ThermalFlux&);// read the raw flux data

/*****************************************************
 ********************** Flux *************************
 *****************************************************/

class Flux : public ThermalFlux {
  Array<double> e_sum;  // energy resolved accumulated value
  Array<double> e_var;  // energy resolved variation

  Array_2<double> ej_sum; // ej_resolved accumulated value
  Array_2<double> ej_var; // ej_resolved variation

public:

  static double pr_exp;

  Flux();
  //real flux calculation HERE!
  Flux (const Div_surf&, int, std::vector<Samp>&);

  Flux& operator+= (const Flux&);
  Flux  operator+  (const Flux&) const;

  double mc_flux (int i) const { return e_sum[i]; }
  double mc_dflx (int i) const { return e_var[i]; }
  double ej_flux (int i) const { return ej_sum[i]; }
  double ej_dflx (int i) const { return ej_var[i]; }

  void normalize ();// normalize fluxes and errors

  void slave_send  () const; // mpi send
  void master_recv (int);    // mpi receive

  friend ostream& operator<< (ostream&, const Flux&);
  friend istream& operator>> (istream&, Flux&);

};

ostream& operator<< (ostream&, const Flux&);
istream& operator>> (istream&, Flux&);

/***************************************************************
 ************************ Flux Array ***************************
 ***************************************************************/

// flux throw all faces
class FluxArray {
  Array<Flux> _array;

public:
  explicit FluxArray(int n) : _array(n) {}
  int size() const { return _array.size(); }
  Flux& operator[] (int i) { return _array[i]; }
  const Flux& operator[] (int i) const { return _array[i]; }
  void normalize();
};

ostream& operator<< (ostream&, const FluxArray&);
istream& operator>> (istream&, FluxArray&);

/***********************************************************
 ***************** Distributed Calculation *****************
 ***********************************************************/

// Base for distributed flux calculation
class DistBase {
  std::set<int> _nodes;
public:
  void add_node(int);
  void remove_node(int);
  int node_size () const { return _nodes.size(); }
  int node (int i) const;
};

// Distributed Flux
class DistFlux : public Flux, public DistBase {
};

// Distributed Flux Array
class DistFluxArray {
  Array<DistFlux> _array;

public:
  explicit DistFluxArray (int n) : _array(n) {}

  DistFlux& operator[] (int i) { return _array[i]; }
  const DistFlux& operator[] (int i) const { return _array[i]; }

  int size() const { return _array.size(); }
  int node_size() const;
  void normalize();
};

ostream& operator<< (ostream&, const DistFluxArray&);
istream& operator>> (istream&, DistFluxArray&);

/***************************************************************
 ******************* Thermal Flux Array ************************
 ***************************************************************/

class ThermalFluxArray {
  Array<ThermalFlux> _array;

public:
  ThermalFluxArray(const FluxArray&);
  ThermalFluxArray(const DistFluxArray&);

  ThermalFluxArray& operator= (const FluxArray&);
  ThermalFluxArray& operator= (const DistFluxArray&);

  int size() const { return _array.size(); }

  ThermalFlux& operator[] (int i) { return _array[i]; }
  const ThermalFlux& operator[] (int i) const { return _array[i]; }

  void normalize();
};


/***************************************************************
 ******************* Thermal Flux Array Calculation ************
 ***************************************************************/

// Base for flux calculation
class Calc {
  bool _wait;

  static double _tol; // required accuracy
  static int _pmx;  // maximal number of potential samplings
  static int _pmn;  // minimal number of potential samplings
  static int _tmx;  // maximal number of total samplings
  static int _tmn;  // minimal number of total samplings
  static int _pln;  // number of potential samplings in one calculation

  static bool _isinit;  // are static vars initialized?
  
public:
  Calc () : _wait(false) {}
  void nowait () { _wait = false; }
  void wait () { _wait = true; }
  bool iswait () const { return _wait; }

  enum { 
    FLUX,    // flux sampling needed
    SURF,    // surface sampling needed
    WAIT,    // wait for sampling results
    STOP     // ready to go
  };

  static void stat_init(double, int, int, int, int, int);

  static double tol  () { return _tol; }
  static int pot_max () { return _pmx; }
  static int pot_min () { return _pmn; }
  static int tot_max () { return _tmx; }
  static int tot_min () { return _tmn; }
  static int pot_len () { return _pln; }
  static bool isinit () { return _isinit; }
};

// Main class for the calculation of flux through a single face
class DistFluxCalc : public DistFlux, public Calc {
public:
  int check_state(int&);
};

// Main class for the calculation of flux through multiple faces
class DistFluxArrayCalc : public DistFluxArray, public Calc {
  Array<int> _smp_num;
public:
  explicit DistFluxArrayCalc(int n): DistFluxArray(n), _smp_num(n) {}
  const int* surf_smp_num () const { return _smp_num.data(); }
  int check_state(int&);
};

/*****************************************************************
 ********************* Thermal MultiFlux *************************
 *****************************************************************/

class ThermalMultiFlux : public FluxBase {

protected:
  Array_2<double> t_sum; // cannonical flux
  Array_2<double> t_var; // flux error

  Array<double> min_en;
  Array<Dynvar> min_dv;

public:
  ThermalMultiFlux ();

  ThermalMultiFlux& operator+= (const ThermalMultiFlux&);
  ThermalMultiFlux  operator+  (const ThermalMultiFlux&) const;
  ThermalMultiFlux& operator*= (double);
  
  double  min_energy (int i) const { return min_en[i]; }
  Dynvar& min_dynvar (int i) { return min_dv[i]; }
  const Dynvar& min_dynvar (int i) const { return min_dv[i]; }

  double cn_flux (int i, int j) const { return t_sum(i, j); }
  double cn_dflx (int i, int j) const { return t_var(i, j); }

  double average (int, int) const;// average value of the cannonical flux
  double pot_var (int, int) const;// potential flux variation
  double vol_var (int, int) const;// volume flux variation

  void normalize ();// normalize fluxes and errors

  void slave_send  () const; // mpi send
  void master_recv (int);    // mpi receive

  friend istream& operator>> (istream&, ThermalMultiFlux&);
  friend ostream& operator<< (ostream&, const ThermalMultiFlux&);
};// Thermal MultiFlux

ostream& operator<< (ostream&, const ThermalMultiFlux&);
istream& operator>> (istream&, ThermalMultiFlux&);

/**********************************************************
 ********************** MultiFlux *************************
 **********************************************************/

class MultiFlux : public ThermalMultiFlux {
  Array_2<double> e_sum;  // energy resolved accumulated value
  Array_2<double> e_var;  // energy resolved variation

  Array_3<double> ej_sum; // ej_resolved accumulated value
  Array_3<double> ej_var; // ej_resolved variation

public:

  MultiFlux();
  //real flux calculation HERE!
  MultiFlux (const Div_surf&, int, std::vector<Samp>&);

  MultiFlux& operator+= (const MultiFlux&);
  MultiFlux  operator+  (const MultiFlux&) const;
  MultiFlux& operator*= (double);

  void normalize ();// normalize fluxes and errors

  void slave_send  () const; // mpi send
  void master_recv (int);    // mpi receive

  friend ostream& operator<< (ostream&, const MultiFlux&);
  friend istream& operator>> (istream&, MultiFlux&);

};

ostream& operator<< (ostream&, const MultiFlux&);
istream& operator>> (istream&, MultiFlux&);

/********************************************************************
 ************************ MultiFlux Array ***************************
 ********************************************************************/

// flux throw all faces
class MultiFluxArray {
  Array<MultiFlux> _array;

public:
  explicit MultiFluxArray(int n) : _array(n) {}
  int size() const { return _array.size(); }
  MultiFlux& operator[] (int i) { return _array[i]; }
  const MultiFlux& operator[] (int i) const { return _array[i]; }
  void normalize();
};

ostream& operator<< (ostream&, const MultiFluxArray&);
istream& operator>> (istream&, MultiFluxArray&);

/***********************************************************
 ***************** Distributed Calculation *****************
 ***********************************************************/

class DistMultiFlux : public MultiFlux, public DistBase {
};

class DistMultiFluxArray {
  Array<DistMultiFlux> _array;

public:
  explicit DistMultiFluxArray (int n) : _array(n) {}

  DistMultiFlux& operator[] (int i) { return _array[i]; }
  const DistMultiFlux& operator[] (int i) const { return _array[i]; }

  int size() const { return _array.size(); }
  int node_size() const;
  void normalize();
};

ostream& operator<< (ostream&, const DistMultiFluxArray&);
istream& operator>> (istream&, DistMultiFluxArray&);

/********************************************************************
 ******************* Thermal MultiFlux Array ************************
 ********************************************************************/

class ThermalMultiFluxArray {
  Array<ThermalMultiFlux> _array;

public:
  ThermalMultiFluxArray(const MultiFluxArray&);
  ThermalMultiFluxArray(const DistMultiFluxArray&);

  ThermalMultiFluxArray& operator= (const MultiFluxArray&);
  ThermalMultiFluxArray& operator= (const DistMultiFluxArray&);

  int size() const { return _array.size(); }

  ThermalMultiFlux& operator[] (int i) { return _array[i]; }
  const ThermalMultiFlux& operator[] (int i) const { return _array[i]; }

  void normalize();
};


/******************************************************
 ******************* MultiFlux Calculation ************
 ******************************************************/

// Main class for the calculation of flux through a single face
class DistMultiFluxCalc : public DistMultiFlux, public Calc {
public:
  int check_state(int&);
};

// Main class for the calculation of flux through multiple faces
class DistMultiFluxArrayCalc : public DistMultiFluxArray, public Calc {
  Array<int> _smp_num;
public:
  explicit DistMultiFluxArrayCalc(int n): DistMultiFluxArray(n), _smp_num(n) {}
  const int* surf_smp_num () const { return _smp_num.data(); }
  int check_state(int&);
};

/*****************************************************************************
 *********************** Dividing Surface Identifier *************************
 *****************************************************************************/

class Sid
{
  int sid_[2]; // dividing surface number and face
  //void check () const;

public:

  int num() const { return sid_[0]; }
  int face() const { return sid_[1]; }
  void set_num(int n) { sid_[0] = n; }
  void set_face(int f) { sid_[1] = f; }

  Sid (int num_ =0, int face_ =0);
  Sid (const Sid& s) 
    { sid_[0] = s.num(); sid_[1] = s.face();}
  Sid& operator= (const Sid& s) 
    { sid_[0] = s.num(); sid_[1] = s.face(); return *this; }

  // comparisons
  bool operator< (const Sid& s) const
  { return num() < s.num() || num() == s.num() && face() < s.face(); } 
  bool operator> (const Sid& s) const
  { return num() > s.num() || num() == s.num() && face() > s.face(); } 
  bool operator== (const Sid& s) const
  { return num() == s.num() && face() == s.face(); }
  bool operator!= (const Sid& s) const
  { return num() != s.num() || face() != s.face(); }
 
  // mpi send and receive
  void send (int dest, int tag) const;
  int slave_recv ();
  MPI::Status master_recv ();

  friend ostream& operator<< (ostream&, const Sid&);
  friend istream& operator>> (istream&, Sid&);

};

inline Sid::Sid (int n, int f) { sid_[0] = n; sid_[1] = f; }// check(); }

ostream& operator<< (ostream&, const Sid&);
istream& operator>> (istream&, Sid&);

#endif
