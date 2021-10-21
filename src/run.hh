#ifndef __RUN__
#define __RUN__

#include "rotd.hh"
#include "tmatrix.hh"
#include "force.hh"

namespace Traj {

  extern double rel_tol;       // integration relative tolerance
  extern double abs_tol;       // integration absolute tolerance
  extern double time_step;     // integration time step

  extern int dist_out_flag;    // distance output
  extern int aux_out_flag;     // auxiliary output
  extern int anim_out_flag;    // animated trajectory output
  extern int dvd_out_flag;     // internal dvd output
  extern int ang_out_flag;     // angular momentum info output
  extern int dyn_out_flag;     // dynamic variables output
  extern int ener_out_flag;    // energy output

  extern Div_surf reac_surf;   // reactive surface
  extern Div_surf tran_surf;   // transition state surface
  extern Div_surf diss_surf;   // dissociative surface

  extern int reac_ener_flag;   // using reactive energy condition
  extern vector<double> reac_ener;     // array of reactive energies
  extern int reac_surf_flag;   // using reactive surface mechanism

  //extern pot_f pot;            // potential function
  extern string base_name;     // trajectory base name

  extern int save_time_step; // time interval between savings traj_info
  extern int debug;             // debugging level
  extern string debug_mesg;
  enum { INFO = 1, DEBUG };     // debug levels

  class Error {};
  class Range_Err : public Error {};
  class Pot_Err   : public Error {};
  class Form_Err  : public Error {};

  class Res {// result of trajectory propagation

    int _lf [2];

  public:

    Res () { _lf[0] = 0; _lf[1] = -1; } // level = 0, face = -1

    void set_ener_level (int l) { _lf[0]=l; }
    int ener_level () const { return _lf[0]; }
    int add_level () { return ++_lf[0]; } 
    int reac_face  () const { return _lf[1]; }
    void set_face (int f) { _lf[1] = f; }
    bool is_reac () { return ener_level() > 0 || reac_face() >= 0; }

    void send (int, int) const; // should be sent and received in one chank
    void recv (int*, int*);
    void recv (int, int);

    friend istream& operator>> (istream&, Res&);
    friend ostream& operator<< (ostream&, const Res&);

  };

  istream& operator>> (istream&, Res&);
  ostream& operator<< (ostream&, const Res&);

  class State : public Dynvar, public Res {

    double time;

  public:

    State () : time(0.) {}
    State (const Dynvar& dynvar) : Dynvar(dynvar), time(0.) {}

    void send (int, int) const;
    void recv (int, int);

    void run (int);

    friend istream& operator>> (istream&, State&);
    friend ostream& operator<< (ostream&, const State&);
  };

  istream& operator>> (istream&, State&);
  ostream& operator<< (ostream&, const State&);
  
}

#endif

