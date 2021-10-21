#ifndef __ROT_DYN__
#define __ROT_DYN__

#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <iostream>

#include "tmatrix.hh"
#include "surface.hh"
#include "math.hh"

using namespace std;

struct Atom;
class Molecule;

extern vector<Molecule*> mol_array;
typedef vector<Molecule*>::iterator Moter;
typedef vector<Molecule*>::const_iterator Const_moter;

void update_atoms();
void init_dv();
void mol_init(const string&); // read molecular structure from the file

bool are_atoms_close(double);
void print_dist(ostream&);
void print_dist(ostream&, const std::vector<D3>&);
void print_geom(ostream&);
void print_geom(ostream&, const std::vector<D3>&);
void print_ang_mom_info(ostream&);
void print_cm_info(ostream&);

int  get_tot_iner_mom (double*);
void get_cm_pos (double*);
void get_cm_vel (double*);
void get_orb_ang_mom (double*);
double get_ref_dist (const Div_surf&, int* =0);

/*************************************************************
 ********************** Dynamical Variables ******************
 *************************************************************/

class Dynvar :  public Array<double> {
public:
  Dynvar();
  void attach();// attach to the molecular array
};

/************************* Atom *************************************/

extern vector<Atom> atoms; // array of all atoms of the system
typedef vector<Atom>::iterator Ater;
typedef vector<Atom>::const_iterator Const_ater;

enum {
  Dummy        = 0  ,
  Hydrogen     = 1  ,  	
  Helium       = 2  ,  	
  Lithium      = 3  ,  	
  Beryllium    = 4  ,	
  Boron        = 5  ,	
  Carbon       = 6  ,	
  Nitrogen     = 7  ,	
  Oxygen       = 8  ,	
  Fluorine     = 9  ,	
  Neon         = 10 ,	
  Sodium       = 11 ,	
  Magnesium    = 12 ,	
  Aluminum     = 13 ,	
  Silicon      = 14 ,	
  Phosphorus   = 15 ,	
  Sulfur       = 16 ,	
  Chlorine     = 17 ,	
  Argon        = 18 ,	
  Potassium    = 19 ,	
  Calcium      = 20 ,	
  Scandium     = 21 ,	
  Titanium     = 22 ,	
  Vanadium     = 23 ,	
  Chromium     = 24 ,	
  Manganese    = 25 ,	
  Iron         = 26 ,	
  Cobalt       = 27 ,	
  Nickel       = 28 ,	
  Copper       = 29 ,	
  Zinc         = 30 ,	
  Gallium      = 31 ,	
  Germanium    = 32 ,	
  Arsenic      = 33 ,	
  Selenium     = 34 ,	
  Bromine      = 35 ,	
  Krypton      = 36 ,	
  Rubidium     = 37 ,	
  Strontium    = 38 ,	
  Yttrium      = 39 ,	
  Zirconium    = 40 ,	
  Niobium      = 41 ,	
  Molybdenum   = 42 ,	
  Technetium   = 43 ,	
  Ruthenium    = 44 ,	
  Rhodium      = 45 ,	
  Palladium    = 46 ,	
  Silver       = 47 ,	
  Cadmium      = 48 ,	
  Indium       = 49 ,	
  Tin          = 50 ,	
  Antimony     = 51 ,	
  Tellurium    = 52 ,	
  Iodine       = 53 ,	
  Xenon        = 54 ,	
  Cesium       = 55 ,	
  Barium       = 56 ,	
  Lanthanum    = 57 ,	
  Cerium       = 58 ,	
  Praseodymium = 59 ,	
  Neodymium    = 60 ,	
  Promethium   = 61 ,	
  Samarium     = 62 ,	
  Europium     = 63 ,	
  Gadolinium   = 64 ,	
  Terbium      = 65 ,	
  Dysprosium   = 66 ,	
  Holmium      = 67 ,	
  Erbium       = 68 ,	
  Thulium      = 69 ,	
  Ytterbium    = 70 ,	
  Lutetium     = 71 ,	
  Hafnium      = 72 ,	
  Tantalum     = 73 ,	
  Tungsten     = 74 ,	
  Rhenium      = 75 ,	
  Osmium       = 76 ,	
  Iridium      = 77 ,	
  Platinum     = 78 ,	
  Gold         = 79 ,	
  Mercury      = 80 ,	
  Thallium     = 81 ,	
  Lead         = 82 ,	
  Bismuth      = 83 ,	
  Polonium     = 84 ,	
  Astatine     = 85 ,	
  Radon        = 86 ,	
  Francium     = 87 ,	
  Radium       = 88 ,	
  Actinium     = 89 ,	
  Thorium      = 90 ,	
  Protactinium = 91 ,	
  Uranium      = 92 ,	
  Neptunium    = 93 ,	
  Plutonium    = 94 ,	
  Americium    = 95 ,	
  Curium       = 96 ,	
  Berkelium    = 97 ,	
  Californium  = 98 ,	
  Einsteinium  = 99 ,	
  Fermium      = 100,	
  Mendelevium  = 101,	
  Nobelium     = 102,	
  Lawrencium   = 103,	
  Rutherfordium= 104,	
  Dubnium      = 105,	
  Seaborgium   = 106,	
  Bohrium      = 107,	
  Hassium      = 108,	
  Meitnerium   = 109,	
  Darmstadtium = 110,	
  Roentgenium  = 111,	
  Copernicium  = 112,	
  Nihonium     = 113,	
  Flerovium    = 114,	
  Moscovium    = 115,	
  Livermorium  = 116,	
  Tennessine   = 117,	
  Oganesson    = 118	
};		     	  
		     	  
enum Mol_type {ATOM, LINEAR, NONLINEAR};

class Atom {
  static const char* _name (int);

public:
  Atom (const string&, int = 0, const double* = 0);

  int type;                    // type of the atom
  int isotope;                 // isotope
  double lf_pos  [3];          // position of the atom in the laboratory frame
  double mf_pos  [3];          // position of the atom in the molecular frame
  double rel_pos [3];          // position of the atom in LF relative to CM
  double force   [3];          // force in the laboratory frame

  static bool is_initialized;  // is the trajectory starting point chosen?

  double      mass () const;
  const char* name () const { return _name(type); }
};

ostream& operator<< (ostream&, const Atom&);

/********************** Molecule **************************/

class Molecule    // Monoatomic molecule
{
  int start;         // index of the first atom of the atoms array

protected:

  bool _ismulti, _isdipole, _isquad, _ispolar, _isip, _ischarge;
  double _charge; // charge
  double _polar;// polarizability
  double _ip; // ionization potential

  // dynamical variables
  double*  dyn_var;       // pointer to the  dynamical variables
  double* dyn_var_der;    // pointer to the dynamical variables derivatives

  Molecule () : start(atoms.size()), _ismulti(false), _isip(false),
		_isdipole(false), _isquad(false), _ispolar(false), 
		_ischarge(false), dyn_var(0), dyn_var_der(0) {}

public:
  bool is_charge () const { return _ischarge; }
  bool is_dipole () const { return _isdipole; }
  bool is_quad   () const { return _isquad; }
  bool is_polar ()  const { return _ispolar; }
  bool is_multi ()  const { return _ismulti; }
  bool is_ip () const { return _isip; }

  double charge  () const { return _charge; }
  virtual double dipole () const = 0;
  virtual double quad   () const = 0;
  double polar  () const { return _polar; }
  double ip () const { return _ip; }

  Ater begin () { return atoms.begin() + start; }
  Const_ater begin () const { return atoms.begin() + start; }
  Atom& atom (int i) { return *(atoms.begin() + start + i); }
  const Atom& atom (int i) const { return *(atoms.begin() + start + i); }
  virtual Ater end () = 0;
  virtual Const_ater end () const = 0;

  virtual int size    () const = 0;       // number of atoms
  virtual Mol_type type () const = 0;
    
  virtual int dv_size () const = 0;
  virtual void set_dv  (double*&) = 0;
  virtual void set_dvd (double*&) = 0;
  virtual void init_dv ();

  const double* read_cm_pos () const { return dyn_var; }      // CM position
  const double* read_cm_vel () const { return dyn_var + 3; }  // CM velocity
  virtual const double* read_ang_pos () const = 0;            // orientation
  virtual const double* read_ang_vel () const  = 0;        // angular velocity

  double* write_cm_pos () { return dyn_var; }
  double* write_cm_vel () { return dyn_var + 3; }
  virtual double* write_ang_pos () = 0;
  virtual double* write_ang_vel ()  = 0;

  virtual int ang_pos_size () const = 0;
  virtual int ang_vel_size () const = 0;

  virtual void update_dvd   () = 0;  // derivatives of dynamic variables
  virtual void update_atoms () = 0;  // updates the coordinates of atoms

  virtual bool is_ang_normalized () const = 0;
  virtual void normalize () = 0; /* normalize angular coordinates represen-
				    tation (for polyatomic molecules) */

  virtual bool is_ang_vel_orthogonal () const { return true; }
  virtual void orthogonalize () {} /* make angular velocity orthogonal to
				      the linear  molecule axis */

  virtual void lf2mf (const double*, double*) const = 0; //LF to MF
  virtual void mf2lf (const double*, double*) const = 0; // MF to LF

  virtual int nm_size () const =0;
  virtual double nm_force (int mode) const =0;
  virtual void nm2lf(const double*, std::vector<D3>&) const =0;  // convert normal 

  // mode coordinates to laboratory frame coordinates
  virtual double mass () const = 0;       // mass of the molecule
  virtual double iner_mom (int) const = 0; // inertia moment
  virtual double imm (int, int) const = 0; // inertia moment matrix in LF
  virtual double kin_energy () const = 0;
  virtual void get_ang_mom (double*) const = 0;
  virtual double stat_sum () const =0;//partition function without temp.factor

  // rotate molecule
  // virtual void rotate (double angle, int axis) = 0;
};

class Monoatomic : public Molecule
{
  void read_multipole_moments (istream&);

public:


  double dipole () const {
    std::cout << "Monoatomic::dipole: should not be here\n";
    return 0.;
  }
  double quad   () const {
    std::cout << "Monoatomic::quad: should not be here\n";
    return 0.;
  }

  Monoatomic (istream&);

  Ater end () { return begin() + 1; }
  Const_ater end () const { return begin() + 1; }

  int size ()    const { return 1; } // number of atoms
  Mol_type type () const { return ATOM; }

  int dv_size () const { return 6; }
  void set_dv  (double*& dv)  { dyn_var = dv; dv += dv_size(); }
  void set_dvd (double*& dvd) { dyn_var_der = dvd; dvd += dv_size(); }

  const double* read_ang_pos () const { return 0; }
  const double* read_ang_vel () const { return 0; }

  double* write_ang_pos () { return 0; }
  double* write_ang_vel () { return 0; }

  int ang_pos_size () const { return 0; }
  int ang_vel_size () const { return 0; }

  void update_atoms ();
  void update_dvd   ();

  bool is_ang_normalized () const { return true; }
  void normalize () {}

  void lf2mf (const double* lf_vec, double* mf_vec) const
  {
    for(int i = 0; i < 3; ++i)
      mf_vec[i] = lf_vec[i];
  }
  void mf2lf (const double* mf_vec, double* lf_vec) const
  {
    for(int i = 0; i < 3; ++i)
      lf_vec[i] = mf_vec[i];
  }

  int nm_size () const { return 0; }
  double nm_force (int mode) const { return 0.; }
  void nm2lf (const double* q, std::vector<D3>& r) const {
    for(int i = 0; i < 3; ++i)
      r[0][i] = begin()->lf_pos[i];
  }

  double mass () const { return begin()->mass(); } // mass of the molecule
  double iner_mom (int i) const { return 0.0; }
  double imm (int, int) const { return 0.0; }
  double kin_energy () const;
  void get_ang_mom (double*) const;
  double stat_sum () const { return 1.; }
};

class Normal_mode {
  double _freq;
  std::vector<D3> _pos;
public:
  Normal_mode (int s) : _pos(s) {}
  double frequency () const { return _freq; }
  void set_frequency (double f) { _freq = f; }
  D3& operator[] (int i) { return _pos[i]; }
  const D3& operator[] (int i) const { return _pos[i]; }
  int size () const { return _pos.size(); }
};

class NormalModeCoor {// normal mode coordinates
  Array<double> q0;
  Array<double> q1;
public:
  NormalModeCoor () : q0(mol_array[0]->nm_size()), q1(mol_array[1]->nm_size()) 
  {q0.init(); q1.init(); }
  int size(int frag) const 
  { if(frag) return q1.size(); else return q0.size(); }
  double& operator() (int frag, int i) 
  { if(frag) return q1[i]; else return q0[i]; }
  double operator() (int frag, int i) const 
  { if(frag) return q1[i]; else return q0[i]; }
  double* data(int frag) 
  {if(frag) return q1.data(); else return q0.data(); }
  const double* data(int frag) const 
  {if(frag) return q1.data(); else return q0.data(); }
  NormalModeCoor& operator+= (const NormalModeCoor& n) { q0 += n.q0; q1 += n.q1; return *this;}
  NormalModeCoor& operator= (const NormalModeCoor& n);
  void init () {q0.init(); q1.init();}
};

void nm2pos(const NormalModeCoor&, std::vector<D3>&);

istream& operator>> (istream&, Normal_mode&);

class Polyatomic : public Molecule
{
  double mm;               // mass of the molecule
  int nn;                  // the number of atoms
  void set_iner_mom ();
  void set_mass ();
  double _torque [3];
  double _force [3];

  double _dipole; //  z-value
  double _quad;   // zz-value

protected:

  bool _isnm;
  void set_normal_modes (istream&);

  TMatrix<double, 3> orig_mfo;  // initial molecular frame orientation matrix
  TMatrix<double, 3> mfo;  // molecular frame orientation matrix

  double ang_mass [3];     // inertia moment
  bool ang_normalized;     // check if angular variables need normalization

  Polyatomic (istream&);

  double force (int) const;
  double torque (int) const;

  vector<Normal_mode> nm_array;
  
  void read_multipole_moments (istream&);

public:

  double dipole () const { return _dipole; }
  double quad   () const { return _quad; }

  void set_torque (double* t) { 
    for(int i = 0; i < 3; ++i)
      _torque[i] = t[i];
  }
  void set_forces (double* f) { 
    for(int i = 0; i < 3; ++i)
      _force[i] = f[i];
  }

  Ater end () { return begin() + nn; }
  Const_ater end () const { return begin() + nn; }

  int size ()    const { return nn; } // number of atoms

  const double* read_ang_pos () const { return dyn_var + 6; }
  double* write_ang_pos () { return dyn_var + 6; }
  int ang_vel_size () const { return 3; }

  bool is_ang_normalized () const { return ang_normalized; }

  void lf2mf (const double* lf_vec, double* mf_vec) const
  { matrix_vector_product (mfo, lf_vec, mf_vec); }
  void mf2lf (const double* mf_vec, double* lf_vec) const
  { vector_matrix_product (mf_vec, mfo, lf_vec); }

  int nm_size () const { return nm_array.size(); }
  double nm_force (int) const;
  void nm2lf (const double*, std::vector<D3>&) const;

  double mass () const { return mm; } // mass of the molecule
};

class Linear : public Polyatomic
{
  mutable bool ang_vel_orthogonal;  // check if the ang. vel. is orthogonal
  void set_rot_mat ();

public:

  Linear (istream& from);

  Mol_type type () const { return LINEAR; }

  int dv_size () const { return 12; }
  void set_dv  (double*& dv)  { dyn_var = dv; dv += dv_size(); }
  void set_dvd (double*& dvd) { dyn_var_der = dvd; dvd += dv_size(); }
  void init_dv ();

  const double* read_ang_vel () const { return dyn_var + 9; }
  double* write_ang_vel () { return dyn_var + 9; }
  int ang_pos_size () const { return 3; }

  void update_atoms ();
  void update_dvd   ();

  void normalize ();
  bool is_ang_vel_orthogonal () const { return ang_vel_orthogonal; }
  void orthogonalize ();

  double iner_mom (int i) const { return ang_mass [2]; }
  double imm (int, int) const;
  double kin_energy () const;
  void get_ang_mom (double*) const;
  double stat_sum () const { return 2. * ang_mass[2]; }
};

class Nonlinear : public Polyatomic
{
  double euler_coef [3];

public:

  Nonlinear (istream& from);

  Mol_type type () const { return NONLINEAR; }

  int dv_size () const { return 13; }
  void set_dv  (double*& dv)  { dyn_var = dv; dv += dv_size(); }
  void set_dvd (double*& dvd) { dyn_var_der = dvd; dvd += dv_size(); }
  void init_dv ();

  const double* read_ang_vel () const { return dyn_var + 10; }
  double* write_ang_vel () { return dyn_var + 10; }
  int ang_pos_size () const { return 4; }

  void update_atoms ();
  void update_dvd   ();

  void normalize ();

  double iner_mom (int i) const { return ang_mass [i]; }
  double imm (int, int) const;
  double kin_energy () const;
  void get_ang_mom (double*) const;
  double stat_sum () const { return 2. * sqrt(2.* M_PI * ang_mass[0]
					      * ang_mass[1] * ang_mass[2]); }
};

#endif
