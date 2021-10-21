#ifndef ROTD_MOLPRO_HH
#define ROTD_MOLPRO_HH

#include <string>
#include <iostream>
#include <vector>

#include "tmatrix.hh"
#include "math.hh"

namespace Molpro {
  using namespace std;

  extern int debug;

  // Error classes
  class Err {};

  // File handling errors
  class File_Err: public Err {};
  class Open_Err: public File_Err{}; // File opening errors
  class Form_Err: public File_Err{}; // Wrong format
  class EOF_Err:  public File_Err{}; // End-of-file encountered

  // Math errors
  class Math_Err: public Err {};
  // Indexing errors
  class Range_Err: public Err {};
  // Search errors
  class Find_Err: public Err {};
  // Initialization errors
  class Init_Err: public Err {};
  // Run errors
  class Run_Err: public Err {}; 

  struct Atom {
    string name;                 // atomic number
    double pos [3];          // atom's position
    Atom () { pos[0] = pos[1] = pos[2] = 0.;}
  };
  extern vector<Atom> atoms;
  extern vector<D3> forces;
  extern vector<double> energies;

  // different flags
  enum {WF = 1, FORCE = 2};

  // wave function file name
  extern std::string wf_name;

  // initialization
  void init (const string&);
  bool isinit ();

  // molpro potential energy
  void pot (int method = 0, int flags = 0);

  void print_dist ();
  ostream& print_geom (ostream&);

}

#endif
