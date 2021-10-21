#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>

#include "rotd.hh"
#include "error.hh"
#include "units.hh"
#include "gauss.hh"
#include "molpro.hh"

using namespace std;

/*************************************************************
 ********************** Dynamical Variables ******************
 *************************************************************/

Dynvar::Dynvar () 
  : Array<double>(mol_array[0]->dv_size() + mol_array[1]->dv_size()) 
{init();}

void Dynvar::attach ()
{
  double* dp = data();
  for (int i = 0; i < 2; ++i)
    mol_array[i]->set_dv(dp);
}

/************************* Global Definitions ************************/

vector<Molecule*> mol_array;

bool are_atoms_close(double d)
{
  for (Ater at0 = mol_array[0]->begin(); at0 != mol_array[0]->end(); ++at0)
    for (Ater at1 = mol_array[1]->begin(); at1 != mol_array[1]->end(); ++at1)
      if(at0->type != Dummy && at1->type != Dummy && 
	 vdistance(at0->lf_pos, at1->lf_pos, 3) < d)
	return true;
  return false;
}

void update_atoms()
{
  for(int i = 0; i < 2; ++i)
    mol_array[i]->update_atoms();
}

void init_dv ()
{
  for(int i = 0; i < 2; ++i)
    mol_array[i]->init_dv();
}

void get_cm_pos (double* cm_pos)
{
  double mass = 0.;
  for(int frag = 0; frag < 2; ++frag)
    mass += mol_array[frag]->mass();
 
  for(int i = 0; i < 3; ++i) {
    cm_pos[i] = 0.;
    for(int frag = 0; frag < 2; ++frag)
      cm_pos[i] += mol_array[frag]->mass() 
	* mol_array[frag]->read_cm_pos()[i];
    cm_pos[i] /= mass;
  }
}

void get_cm_vel (double* cm_vel)
{
  double mass = 0.;
  for(int frag = 0; frag < 2; ++frag)
    mass += mol_array[frag]->mass();
 
  for(int i = 0; i < 3; ++i) {
    cm_vel[i] = 0.;
    for(int frag = 0; frag < 2; ++frag)
      cm_vel[i] += mol_array[frag]->mass() 
	* mol_array[frag]->read_cm_vel()[i];
    cm_vel[i] /= mass;
  }
}

void get_orb_ang_mom (double* orb_am)
{
  double cm_pos [3];
  get_cm_pos(cm_pos);

  double pos [3];
  double vprod [3];

  for(int i = 0; i < 3; ++i)
      orb_am[i] = 0.;

  for(int frag = 0; frag < 2; ++frag) {
    for(int i = 0; i < 3; ++i)
      pos[i] = mol_array[frag]->read_cm_pos()[i] - cm_pos[i];

    vector_product(pos,  mol_array[frag]->read_cm_vel(), vprod);

    for(int i = 0; i < 3; ++i)
      orb_am[i] += mol_array[frag]->mass() * vprod[i];
  }
}

int get_tot_iner_mom (double* tim)
{
  const char funame [] = "get_tot_iner_mom: ";

  if(!tim) {
    cout << funame <<  "null pointer\n";
    return -1;
  }

  const double emass = mol_array[0]->mass() * mol_array[1]->mass() /
    (mol_array[0]->mass() + mol_array[1]->mass());

  double pos_0 [3];
  for(int i = 0; i < 3; ++i)
    pos_0[i] = mol_array[1]->read_cm_pos()[i] 
      - mol_array[0]->read_cm_pos()[i];
 
  TMatrix<double, 3> tim_mat;
  double r0_2 = 0.;
  for (int i = 0; i < 3; ++i)
    r0_2 += pos_0[i] * pos_0[i];
    
  for (int i = 0; i < 3; ++i) {
    tim_mat(i, i) = mol_array[0]->imm(i, i) + mol_array[1]->imm(i, i) +
      emass * (r0_2 - pos_0[i] * pos_0[i]);
      
    for (int i1 = i+1; i1 < 3; ++i1)
      tim_mat(i, i1) = 
	mol_array[0]->imm(i, i1) + mol_array[1]->imm(i, i1) -
	emass * pos_0[i] * pos_0[i1];
    }

  tim_mat.eigenv(tim);
  if (tim[0] <= 0.) {
    cout << funame << "inertia moment = " 
	 << tim[0] << " is negative\n";
    return -1;
  }
  
  return 0;
}

void mol_init(const string& spec_file)
{
  string stemp;
  string comment;

  ifstream from(spec_file.c_str());
  if (!from)
    error("mol_init: specification file not found");
  
  for(int frag = 0; frag < 2; ++frag) {

    while (from.peek() == '#')
      getline(from, comment);

    from >> stemp; // molecule type

    if(stemp == "Monoatomic")
      mol_array.push_back(new Monoatomic(from));
    else if(stemp == "Linear")
      mol_array.push_back(new Linear(from)); 
    else if(stemp == "Nonlinear")
      mol_array.push_back(new Nonlinear(from));
    else
      error ("mol_init: wrong molecular type");
  }

  // g98 description initialization
  Gauss::cluster.resize(atoms.size());
  for (int i = 0; i < atoms.size(); ++i)
    Gauss::cluster[i].num = atoms[i].type;

  // molpro initialization
  Molpro::atoms.resize(atoms.size());
  for (int i = 0; i < atoms.size(); ++i)
    Molpro::atoms[i].name = atoms[i].name();

  // dynamical variables initialization
  static Dynvar dynvar;
}

// minimal  distance for a dividing surface
double get_ref_dist (const Div_surf& ds, int* facep)
{
  const char funame [] = "get_ref_dist: ";

  double dtemp;

  double min_dist = 0.;
  int min_face = -1;
  double lf_ref_pos [2] [3];
  for(int face = 0; face < ds.end(); ++face) {
    if(ds.dist(face) <= 0.)
      continue;
    // pivot point positions in laboratory frame
    for (int frag = 0; frag < 2; ++frag) {
      switch (mol_array[frag]->type()) {
      case ATOM:
	for (int i = 0; i < 3; ++i)
	  lf_ref_pos[frag][i] = 0.;
	break;

      case LINEAR:
	dtemp = vlength(mol_array[frag]->read_ang_pos(),3);
	for (int i = 0; i < 3; ++i)
	  lf_ref_pos[frag][i] = mol_array[frag]->read_ang_pos()[i]
	    * ds.ref_pos(frag, ds.ref_ind(frag, face))[0] / dtemp;
	break;

      case NONLINEAR:
	mol_array[frag]-> mf2lf(ds.ref_pos(frag, ds.ref_ind(frag, face)), 
				lf_ref_pos[frag]);
	break;

      default:
	cout << funame << "wrong type\n";
	exit(1);
      }
      for (int i = 0; i < 3; ++i)
	lf_ref_pos[frag][i] += mol_array[frag]->read_cm_pos()[i];
    }

    dtemp = vdistance(lf_ref_pos[0], lf_ref_pos[1], 3) / ds.dist(face);
    if(min_face < 0 || dtemp < min_dist) {
      min_face = face;
      min_dist = dtemp;
    }
  }
  if(facep)
    *facep = min_face;
  return min_dist;
}

// print interatomic distances
void print_dist (ostream& out)
{
  int old_precision = out.precision(4);
  out << "interatomic distance matrix (Angstrom):\n\n";
  out << " F1\\F2";
  for (int i1 = 0; i1 < mol_array[1]->size(); ++i1)
    if(mol_array[1]->begin()[i1].type != Dummy)
      out << setw(6) << mol_array[1]->begin()[i1].name() << i1+1;
  out << "\n";
  
  for (int i0 = 0; i0 < mol_array[0]->size(); ++i0)
    if(mol_array[0]->begin()[i0].type != Dummy){
      out << setw(2) << mol_array[0]->begin()[i0].name() << i0+1 << "   ";
      for (int i1 = 0; i1 < mol_array[1]->size(); ++i1)
	if(mol_array[1]->begin()[i1].type != Dummy)
	  out << setw(7) << Phys_const::bohr *
	    vdistance(mol_array[0]->begin()[i0].lf_pos, 
		      mol_array[1]->begin()[i1].lf_pos, 3);
      out << "\n";
    }
  out << "\n";
  out.precision(old_precision);
}

void print_dist(ostream& out, const std::vector<D3>& pos) {

  int old_precision = out.precision(4);
  out << "interatomic distance matrix (Angstrom):\n\n";
  out << " \\  ";
  for(int i1 = 0; i1 < atoms.size()-1; ++i1)
    if(atoms[i1].type != Dummy)
      out << setw(6) << atoms[i1].name() << i1+1;
  out << "\n";
  
  for(int i0 = 1; i0 < atoms.size(); ++i0)
    if(atoms[i0].type != Dummy){
      out << setw(2) << atoms[i0].name() << i0+1 << " ";
      for(int i1 = 0; i1 < i0; ++i1)
	if(atoms[i1].type != Dummy)
	  out << setw(7) << Phys_const::bohr *
	    vdistance((const double*)pos[i0], (const double*)pos[i1], 3);
      out << "\n";
    }
  out << "\n";
  out.precision(old_precision);
}
	
void print_geom (ostream& out, const std::vector<D3>& pos) {
  int old_precision = out.precision(4);
  out << "geometry (angstrom):\n\n";
  for(int at = 0; at < atoms.size(); ++at) {
    out << setw(2) <<  atoms[at].name();
    for(int i = 0; i < 3; ++i)
      out << " " << setw(9) << Phys_const::bohr * pos[at][i];
    out << "\n";
  }
  out << "\n";
  out.precision(old_precision);
}

void print_geom (ostream& out) {
  int old_precision = out.precision(6);
  ios_base::fmtflags old_flags = out.flags();
  out.setf(ios_base::fixed, ios_base::floatfield);

  for(int at = 0; at < atoms.size(); ++at) {
    out << setw(2) <<  atoms[at].name();
    for(int i = 0; i < 3; ++i)
      out << " " << setw(11) << Phys_const::bohr * atoms[at].lf_pos[i];
    out << "\n";
  }

  out.precision(old_precision);
  out.flags(old_flags);
}

void print_ang_mom_info (ostream& aux_out)
{
  const int tab = 13;
  double tot_ang_mom [3];
  get_orb_ang_mom(tot_ang_mom);
  aux_out << "orbital angular momentum:\n";
  aux_out << setw(tab) << "X, au"
	  << setw(tab) << "Y, au"
	  << setw(tab) << "Z, au"
	  << setw(tab+1) << "Length\n";
  for(int i = 0; i < 3; ++i)
    aux_out << setw(tab) << tot_ang_mom[i];
  aux_out << setw(tab) << vlength(tot_ang_mom, 3);
  aux_out << "\n\n";

  double ang_mom [3];
  for(int frag = 0; frag < 2; ++frag) {
    mol_array[frag]->get_ang_mom(ang_mom);
    aux_out << frag << "-th fragment angular momentum:\n";
    aux_out << setw(tab) << "X, au"
	    << setw(tab) << "Y, au"
	    << setw(tab) << "Z, au"
	    << setw(tab+1) << "Length\n";
    for(int i = 0; i < 3; ++i)
      aux_out << setw(tab) << ang_mom[i];
    aux_out << setw(tab) << vlength(ang_mom, 3);
    aux_out << "\n\n";
    for(int i = 0; i < 3; ++i)
      tot_ang_mom[i] += ang_mom[i]; 
  }
  aux_out << "total angular momentum:\n";
  aux_out << setw(tab) << "X, au"
	  << setw(tab) << "Y, au"
	  << setw(tab) << "Z, au"
	  << setw(tab+1) << "Length\n";
  for(int i = 0; i < 3; ++i)
    aux_out << setw(tab) << tot_ang_mom[i];
  aux_out << setw(tab) << vlength(tot_ang_mom, 3);
  aux_out << "\n\n";
  
  double pos[3];
  for(int i = 0; i < 3; ++i)
    pos[i] = mol_array[1]->read_cm_pos()[i] - mol_array[0]->read_cm_pos()[i];
  ::normalize(pos, 3);
  aux_out << "angular momentum projection onto interfragment axis = "
	  << ::vdot(tot_ang_mom, pos, 3) << " au\nn";
  
}

void print_cm_info (ostream& aux_out)
{
  const int tab = 13;

  double cm_vel [3];
  get_cm_vel(cm_vel);

  aux_out << "center-of-mass velocity:\n";
  aux_out << setw(tab) << "X, au"
	  << setw(tab) << "Y, au"
	  << setw(tab) << "Z, au\n";
  for(int i = 0; i < 3; ++i)
    aux_out << setw(tab) << cm_vel[i];
  aux_out << "\n\n";
}

/************************* Atom *************************************/

vector<Atom> atoms;
bool Atom::is_initialized = false;

ostream& operator<< (ostream& to, const Atom& at)
{
  int old_prec = to.precision(14);
  to << at.name();
  for (int i = 0; i < 3; ++i)
    to << " " << setw(20) << at.lf_pos[i];
  to << "\n";
  to.precision(old_prec);
  return to;
}

const char* Atom::_name (int t)
{
  switch(t) {
  case Dummy:           return "X" ;
  case Hydrogen:        return "H" ; 
  case Helium:          return "He";
  case Lithium:         return "Li";
  case Beryllium:       return "Be";
  case Boron:           return "B" ; 
  case Carbon:          return "C" ;
  case Nitrogen:        return "N" ; 
  case Oxygen:          return "O" ; 
  case Fluorine:        return "F" ; 
  case Neon:            return "Ne";
  case Sodium:          return "Na";
  case Magnesium:       return "Mg";
  case Aluminum:        return "Al";
  case Silicon:         return "Si";
  case Phosphorus:      return "P" ;
  case Sulfur:          return "S" ;
  case Chlorine:        return "Cl";
  case Argon:           return "Ar";
  case Potassium:       return "K" ;
  case Calcium:         return "Ca";
  case Scandium:        return "Sc";
  case Titanium:        return "Ti";
  case Vanadium:        return "V" ;
  case Chromium:        return "Cr";
  case Manganese:       return "Mn";
  case Iron:            return "Fe";
  case Cobalt:          return "Co";
  case Nickel:          return "Ni";
  case Copper:          return "Cu";
  case Zinc:            return "Zn";
  case Gallium:         return "Ga";
  case Germanium:       return "Ge";
  case Arsenic:         return "As";
  case Selenium:        return "Se";
  case Bromine:         return "Br";
  case Krypton:         return "Kr";
  case Rubidium:        return "Rb";
  case Strontium:       return "Sr";
  case Yttrium:         return "Y" ;
  case Zirconium:       return "Zr";
  case Niobium:         return "Nb";
  case Molybdenum:      return "Mo";
  case Technetium:      return "Tc";
  case Ruthenium:       return "Ru";
  case Rhodium:         return "Rh";
  case Palladium:       return "Pd";
  case Silver:          return "Ag";
  case Cadmium:         return "Cd";
  case Indium:          return "In";
  case Tin:             return "Sn";
  case Antimony:        return "Sb";
  case Tellurium:       return "Te";
  case Iodine:          return "I" ;
  case Xenon:           return "Xe";
  case Cesium:          return "Cs";
  case Barium:          return "Ba";
  case Lanthanum:       return "La";
  case Cerium:          return "Ce";
  case Praseodymium:    return "Pr";
  case Neodymium:       return "Nd";
  case Promethium:      return "Pm";
  case Samarium:        return "Sm";
  case Europium:        return "Eu";
  case Gadolinium:      return "Gd";
  case Terbium:         return "Tb";
  case Dysprosium:      return "Dy";
  case Holmium:         return "Ho";
  case Erbium:          return "Er";
  case Thulium:         return "Tm";
  case Ytterbium:       return "Yb";
  case Lutetium:        return "Lu";
  case Hafnium:         return "Hf";
  case Tantalum:        return "Ta";
  case Tungsten:        return "W" ;
  case Rhenium:         return "Re";
  case Osmium:          return "Os";
  case Iridium:         return "Ir";
  case Platinum:        return "Pt";
  case Gold:            return "Au";
  case Mercury:         return "Hg";
  case Thallium:        return "Tl";
  case Lead:            return "Pb";
  case Bismuth:         return "Bi";
  case Polonium:        return "Po";
  case Astatine:        return "At";
  case Radon:           return "Rn";
  case Francium:        return "Fr";
  case Radium:          return "Ra";
  case Actinium:        return "Ac";
  case Thorium:         return "Th";
  case Protactinium:    return "Pa";
  case Uranium:         return "U" ;
  case Neptunium:       return "Np";
  case Plutonium:       return "Pu";
  case Americium:       return "Am";
  case Curium:          return "Cm";
  case Berkelium:       return "Bk";
  case Californium:     return "Cf";
  case Einsteinium:     return "Es";
  case Fermium:         return "Fm";
  case Mendelevium:     return "Md";
  case Nobelium:        return "No";
  case Lawrencium:      return "Lr";
  case Rutherfordium:   return "Rf";
  case Dubnium:         return "Db";
  case Seaborgium:      return "Sg";
  case Bohrium:         return "Bh";
  case Hassium:         return "Hs";
  case Meitnerium:      return "Mt";
  case Darmstadtium:    return "Ds";
  case Roentgenium:     return "Rg";
  case Copernicium:     return "Cn";
  case Nihonium:        return "Nh";
  case Flerovium:       return "Fl";
  case Moscovium:       return "Mc";
  case Livermorium:     return "Lv";
  case Tennessine:      return "Ts";
  case Oganesson:       return "Og";    
  default:
    std::cerr << "Atom::_name: unknown atomic number: " << t << "\n";
    throw Error::Range_Err();
  }
}

Atom::Atom (const string& nm, int isot, const double* xyz)
{
  static map<string, int> _n2n;
  
  if(!_n2n.size())
    for(int n = 0; n < 119; ++n)
      _n2n[_name(n)] = n;
  
  if(_n2n.find(nm) == _n2n.end()) {
    std::cerr << "Atom::Atom: unknown atom name: " << nm << "\n";
    throw Error::Range_Err();
  }

  isotope = isot;
  type = _n2n[nm];
  
  // angstroms to Bohr conversion
  if(xyz)
    for (int i = 0; i < 3; ++i)
      lf_pos [i] = Phys_const::ang2bohr(xyz[i]);
}

double Atom::mass () const
{
  switch (type) {
  case Dummy:
    return Phys_const::amu * isotope;
    
  case Hydrogen:
    switch (isotope) {
    case 0:   return Phys_const::amu * 1.00784;
    case 1:   return Phys_const::amu * 1.007825;
    case 2:   return Phys_const::amu * 2.014;
    case 3:   return Phys_const::amu * 3.01605;
    }
    break;
    
  case Helium:
    switch (isotope) {
    case 0:   return   Phys_const::amu * 4.002602;
    case 3:   return   Phys_const::amu * 3.0160293201;
    case 4:   return   Phys_const::amu * 4.00260325413;
    }
    break;
    
  case Lithium:
    switch (isotope) {
    case 0:   return   Phys_const::amu * 6.997;
    case 6:   return   Phys_const::amu * 6.0151228874;
    case 7:   return   Phys_const::amu * 7.0160034366;
    }
    break;
    
  case Beryllium:
    switch (isotope) {
    case 0:   return   Phys_const::amu * 9.0121831;
    case 9:   return   Phys_const::amu * 9.012183065;
    }
    break;
    
  case Boron:
    switch (isotope) {
    case 0:    return   Phys_const::amu * 10.806;
    case 10:   return   Phys_const::amu * 10.01293695;
    case 11:   return   Phys_const::amu * 11.00930536;
    }
    break;
    
  case Carbon:
    switch (isotope) {
    case 0: return Phys_const::amu * 12.0096;
    case 12: return Phys_const::amu * 12.0;
    case 13: return Phys_const::amu * 13.00335;
    }
    break;
    
  case Nitrogen:
    switch (isotope) {
    case 0:  return Phys_const::amu * 14.00643;
    case 14: return Phys_const::amu * 14.00307;
    case 15: return Phys_const::amu * 15.00011;
    }
    break;
    
  case Oxygen:
    switch (isotope) {
    case 0:  return Phys_const::amu * 15.99903;
    case 16: return Phys_const::amu * 15.99491;
    case 17: return Phys_const::amu * 16.99913175650;
    case 18: return Phys_const::amu * 17.99915961286;
    }
    break;
    
  case Fluorine:
    switch(isotope) {
    case 0:  return Phys_const::amu * 18.9984;
    case 19: return Phys_const::amu * 18.9984;
    }
    break;
    
  case Neon:
    switch(isotope) {
    case 0:    return Phys_const::amu * 20.1797;
    case 20:   return Phys_const::amu * 19.9924401762;
    case 21:   return Phys_const::amu * 20.993846685;
    case 22:   return Phys_const::amu * 21.991385114;
    }
    break;
    
  case Sodium:
    switch(isotope) {
    case 0:  return Phys_const::amu * 22.9898;
    case 23: return Phys_const::amu * 22.9898;
    }
    break;
    
  case Magnesium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 24.304;
    case 24:   return Phys_const::amu * 23.985041697;
    case 25:   return Phys_const::amu * 24.985836976;
    case 26:   return Phys_const::amu * 25.982592968;
    }
    break;

  case Aluminum:
    switch(isotope) {
    case 0:    return Phys_const::amu * 26.98153853;
    case 27:   return Phys_const::amu * 26.98153853;
    }
    break;
    
  case Silicon:
    switch(isotope) {
    case 0:  return Phys_const::amu * 28.084;
    case 28: return Phys_const::amu * 27.97693;
    case 29: return Phys_const::amu * 28.97649;
    case 30: return Phys_const::amu * 29.97376;
    }
    break;
    
  case Phosphorus:
    switch(isotope) {
    case 0:  return Phys_const::amu * 30.97376;
    case 31: return Phys_const::amu * 30.97376;
    }
    break;
    
  case Sulfur:
    switch(isotope) {
    case 0:  return Phys_const::amu * 32.059;
    case 32: return Phys_const::amu * 31.97207;
    case 33: return Phys_const::amu * 32.97146;
    case 34: return Phys_const::amu * 33.96786;
    case 36: return Phys_const::amu * 35.96709;
    }
    break;
    
  case Chlorine:
    switch(isotope) {
    case 0:  return Phys_const::amu * 35.446;
    case 35: return Phys_const::amu * 34.96885;
    case 37: return Phys_const::amu * 36.965902602;
    }
    break;
    
  case Argon:
    switch(isotope) {
    case 0:    return Phys_const::amu * 39.948;
    case 36:   return Phys_const::amu * 35.967545105;
    case 38:   return Phys_const::amu * 37.96273211;
    case 40:   return Phys_const::amu * 39.9623831237;
    }
    break;
    
  case Potassium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 39.0983;
    case 39:   return Phys_const::amu * 38.9637064864;
    case 40:   return Phys_const::amu * 39.963998166;
    case 41:   return Phys_const::amu * 40.9618252579;
    }
    break;
    
  case Calcium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 40.078;
    case 40:   return Phys_const::amu * 39.962590863;
    case 42:   return Phys_const::amu * 41.95861783;
    case 43:   return Phys_const::amu * 42.95876644;
    case 44:   return Phys_const::amu * 43.95548156;
    case 46:   return Phys_const::amu * 45.9536890;
    case 48:   return Phys_const::amu * 47.95252276;
    }
    break;
    
  case Scandium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 44.95590828;
    case 45:   return Phys_const::amu * 44.95590828;
    }
    break;
    
  case Titanium:
    switch(isotope) {
    case 0:  return Phys_const::amu * 47.867;
    case 39: return Phys_const::amu * 39.0013;
    case 40: return Phys_const::amu * 39.9905;
    case 41: return Phys_const::amu * 40.98313;
    case 42: return Phys_const::amu * 41.97303;
    case 43: return Phys_const::amu * 42.96852;
    case 44: return Phys_const::amu * 43.95969;
    case 45: return Phys_const::amu * 44.958124;
    case 46: return Phys_const::amu * 45.95263;
    case 47: return Phys_const::amu * 46.951764;
    case 48: return Phys_const::amu * 47.947947;
    case 49: return Phys_const::amu * 48.947871;
    case 50: return Phys_const::amu * 49.944792;
    case 51: return Phys_const::amu * 50.946616;
    case 52: return Phys_const::amu * 51.9469;
    case 53: return Phys_const::amu * 52.9497;
    case 54: return Phys_const::amu * 53.951;
    case 55: return Phys_const::amu * 54.9552;
    case 56: return Phys_const::amu * 55.958;
    case 57: return Phys_const::amu * 56.964;
    }
    break;
      
  case Vanadium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 50.9415;
    case 50:   return Phys_const::amu * 49.94715601;
    case 51:   return Phys_const::amu * 50.94395704;
    }	       
    break;

  case Chromium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 51.9961;
    case 50:   return Phys_const::amu * 49.94604183;
    case 52:   return Phys_const::amu * 51.94050623;
    case 53:   return Phys_const::amu * 52.94064815;
    case 54:   return Phys_const::amu * 53.93887916;
    }
    break;

  case Manganese:
    switch(isotope) {
    case 0:    return Phys_const::amu * 54.938044;
    case 55:   return Phys_const::amu * 54.93804391;
    }
    break;

  case Iron:
    switch(isotope) {
    case 0:    return Phys_const::amu * 55.845;
    case 54:   return Phys_const::amu * 53.93960899;
    case 56:   return Phys_const::amu * 55.93493633;
    case 57:   return Phys_const::amu * 56.93539284;
    case 58:   return Phys_const::amu * 57.93327443;
    }
    break;

  case Cobalt:
    switch(isotope) {
    case 0:    return Phys_const::amu * 58.933194;
    case 59:   return Phys_const::amu * 58.93319429;
    }
    break;

  case Nickel:
    switch(isotope) {
    case 0:    return Phys_const::amu * 58.6934;
    case 58:   return Phys_const::amu * 57.93534241;
    case 60:   return Phys_const::amu * 59.93078588;
    case 61:   return Phys_const::amu * 60.93105557;
    case 62:   return Phys_const::amu * 61.92834537;
    case 64:   return Phys_const::amu * 63.92796682;
    }
    break;

  case Copper:
    switch(isotope) {
    case 0:    return Phys_const::amu * 63.546;
    case 63:   return Phys_const::amu * 62.92959772;
    case 65:   return Phys_const::amu * 64.92778970;
    }
    break;

  case Zinc:
    switch(isotope) {
    case 0:    return Phys_const::amu * 65.38;
    case 64:   return Phys_const::amu * 63.92914201;
    case 66:   return Phys_const::amu * 65.92603381;
    case 67:   return Phys_const::amu * 66.92712775;
    case 68:   return Phys_const::amu * 67.92484455;
    case 70:   return Phys_const::amu * 69.9253192;
    }
    break;

  case Gallium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 69.723;
    case 69:   return Phys_const::amu * 68.9255735;
    case 71:   return Phys_const::amu * 70.92470258;
    }
    break;

  case Germanium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 72.630;
    case 70:   return Phys_const::amu * 69.92424875;
    case 72:   return Phys_const::amu * 71.922075826;
    case 73:   return Phys_const::amu * 72.923458956;
    case 74:   return Phys_const::amu * 73.921177761;
    case 76:   return Phys_const::amu * 75.921402726;
    }
    break;

  case Arsenic:
    switch(isotope) {
    case 0:    return Phys_const::amu * 74.921595;
    case 75:   return Phys_const::amu * 74.92159457;
    }
    break;

  case Selenium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 78.971;
    case 74:   return Phys_const::amu * 73.922475934;
    case 76:   return Phys_const::amu * 75.919213704;
    case 77:   return Phys_const::amu * 76.919914154;
    case 78:   return Phys_const::amu * 77.91730928;
    case 80:   return Phys_const::amu * 79.9165218;
    case 82:   return Phys_const::amu * 81.9166995;
    }
    break;

  case Bromine:
    switch(isotope) {
    case 0:    return Phys_const::amu * 79.901;
    case 79:   return Phys_const::amu * 78.9183;
    case 81:   return Phys_const::amu * 80.9163;
    }
    break;
    
  case Krypton:
    switch(isotope) {
    case 0:    return Phys_const::amu * 83.798;
    case 78:   return Phys_const::amu * 77.92036494;
    case 80:   return Phys_const::amu * 79.91637808;
    case 82:   return Phys_const::amu * 81.91348273;
    case 83:   return Phys_const::amu * 82.91412716;
    case 84:   return Phys_const::amu * 83.9114977282;
    case 86:   return Phys_const::amu * 85.9106106269;
    }
    break;

  case Rubidium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 85.4678;
    case 85:   return Phys_const::amu * 84.9117897379;
    case 87:   return Phys_const::amu * 86.9091805310;
	}
    break;

  case Strontium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 87.62;
    case 84:   return Phys_const::amu * 83.9134191;
    case 86:   return Phys_const::amu * 85.9092606;
    case 87:   return Phys_const::amu * 86.9088775;
    case 88:   return Phys_const::amu * 87.9056125;
    }
    break;

  case Yttrium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 88.90584;
    case 89:   return Phys_const::amu * 88.9058403;
    }
    break;

  case Zirconium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 91.224;
    case 90:   return Phys_const::amu * 89.9047026;
    case 91:   return Phys_const::amu * 90.9056439;
    case 92:   return Phys_const::amu * 91.9050386;
    case 94:   return Phys_const::amu * 93.9063148;
    case 96:   return Phys_const::amu * 95.908275;
    }
    break;

  case Niobium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 92.90637;
    case 93:   return Phys_const::amu * 92.9063730;
    }
    break;

  case Molybdenum:
    switch(isotope) {
    case 0:    return Phys_const::amu * 95.95;
    case 92:   return Phys_const::amu * 91.90680796;
    case 94:   return Phys_const::amu * 93.90508490;
    case 95:   return Phys_const::amu * 94.90583877;
    case 96:   return Phys_const::amu * 95.90467612;
    case 97:   return Phys_const::amu * 96.90601812;
    case 98:   return Phys_const::amu * 97.90540482;
    case 100:  return Phys_const::amu * 99.9074718;
    }
    break;

  case Technetium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 98.;
    case 97:   return Phys_const::amu * 96.9063667;
    case 98:   return Phys_const::amu * 97.9072124;
    case 99:   return Phys_const::amu * 98.9062508;
    }
    break;

  case Ruthenium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 101.07;
    case 96:   return Phys_const::amu * 95.90759025;
    case 98:   return Phys_const::amu * 97.9052868;
    case 99:   return Phys_const::amu * 98.9059341;
    case 100:  return Phys_const::amu * 99.9042143;
    case 101:  return Phys_const::amu * 100.9055769;
    case 102:  return Phys_const::amu * 101.9043441;
    case 104:  return Phys_const::amu * 103.9054275;
    }
    break;

  case Rhodium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 102.90550;
    case 103:  return Phys_const::amu * 102.9054980;
    }
    break;

  case Palladium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 106.42;
    case 102:  return Phys_const::amu * 101.9056022;
    case 104:  return Phys_const::amu * 103.9040305;
    case 105:  return Phys_const::amu * 104.9050796;
    case 106:  return Phys_const::amu * 105.9034804;
    case 108:  return Phys_const::amu * 107.9038916;
    case 110:  return Phys_const::amu * 109.90517220;
    }
    break;

  case Silver:
    switch(isotope) {
    case 0:    return Phys_const::amu * 107.8682;
    case 107:  return Phys_const::amu * 106.9050916;
    case 109:  return Phys_const::amu * 108.9047553;
    }	       
    break;

  case Cadmium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 112.414;
    case 106:  return Phys_const::amu * 105.9064599;
    case 108:  return Phys_const::amu * 107.9041834;
    case 110:  return Phys_const::amu * 109.90300661;
    case 111:  return Phys_const::amu * 110.90418287;
    case 112:  return Phys_const::amu * 111.90276287;
    case 113:  return Phys_const::amu * 112.90440813;
    case 114:  return Phys_const::amu * 113.90336509;
    case 116:  return Phys_const::amu * 115.90476315;
    }	       
    break;

  case Indium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 114.818;
    case 113:  return Phys_const::amu * 112.90406184;
    case 115:  return Phys_const::amu * 114.903878776;
    }	       
    break;

  case Tin:
    switch(isotope) {
    case 0:    return Phys_const::amu * 118.710;
    case 112:  return Phys_const::amu * 111.90482387;
    case 114:  return Phys_const::amu * 113.9027827;
    case 115:  return Phys_const::amu * 114.903344699;
    case 116:  return Phys_const::amu * 115.90174280;
    case 117:  return Phys_const::amu * 116.90295398;
    case 118:  return Phys_const::amu * 117.90160657;
    case 119:  return Phys_const::amu * 118.90331117;
    case 120:  return Phys_const::amu * 119.90220163;
    case 122:  return Phys_const::amu * 121.9034438;
    case 124:  return Phys_const::amu * 123.9052766;
    }
    break;

  case Antimony:
    switch(isotope) {
    case 0:    return Phys_const::amu * 121.760;
    case 121:  return Phys_const::amu * 120.9038120;
    case 123:  return Phys_const::amu * 122.9042132;
    }	       
    break;

  case Tellurium:
    switch(isotope) {
    case 0:    return Phys_const::amu * 127.60;
    case 120:  return Phys_const::amu * 119.9040593;
    case 122:  return Phys_const::amu * 121.9030435;
    case 123:  return Phys_const::amu * 122.9042698;
    case 124:  return Phys_const::amu * 123.9028171;
    case 125:  return Phys_const::amu * 124.9044299;
    case 126:  return Phys_const::amu * 125.9033109;
    case 128:  return Phys_const::amu * 127.90446128;
    case 130:  return Phys_const::amu * 129.906222748;
    }
    break;

  case Iodine:
    switch(isotope) {
    case 0:    return Phys_const::amu * 126.90447;
    case 127:  return Phys_const::amu * 126.904468;
    case 129:  return Phys_const::amu * 128.904988;
    }
    break;
    
  case Xenon: //54  Xe
    switch(isotope) {
    case 0:    return Phys_const::amu * 131.293;
    case 124:  return Phys_const::amu * 123.9058920;
    case 126:  return Phys_const::amu * 125.9042983;
    case 128:  return Phys_const::amu * 127.9035310;
    case 129:  return Phys_const::amu * 128.9047808611;
    case 130:  return Phys_const::amu * 129.903509349;
    case 131:  return Phys_const::amu * 130.90508406;
    case 132:  return Phys_const::amu * 131.9041550856;
    case 134:  return Phys_const::amu * 133.90539466;
    case 136:  return Phys_const::amu * 135.907214484;
    }
    break;

  case Cesium: //55  Cs
    switch(isotope) {
    case 0:    return Phys_const::amu * 132.90545196;
    case 133:  return Phys_const::amu * 132.9054519610;
    }
    break;

  case Barium: //56  Ba
    switch(isotope) {
    case 0:    return Phys_const::amu * 137.327;
    case 130:  return Phys_const::amu * 129.9063207;
    case 132:  return Phys_const::amu * 131.9050611;
    case 134:  return Phys_const::amu * 133.90450818;
    case 135:  return Phys_const::amu * 134.90568838;
    case 136:  return Phys_const::amu * 135.90457573;
    case 137:  return Phys_const::amu * 136.90582714;
    case 138:  return Phys_const::amu * 137.90524700;
    }
    break;

  case Lanthanum: //57  La
    switch(isotope) {
    case 0:    return Phys_const::amu * 138.90547;
    case 138:  return Phys_const::amu * 137.9071149;
    case 139:  return Phys_const::amu * 138.9063563;
    }	       
    break;

  case Cerium: //58  Ce
    switch(isotope) {
    case 0:    return Phys_const::amu * 140.116;
    case 136:  return Phys_const::amu * 135.90712921;
    case 138:  return Phys_const::amu * 137.905991;
    case 140:  return Phys_const::amu * 139.9054431;
    case 142:  return Phys_const::amu * 141.9092504;
    }
    break;

  case Praseodymium: //59  Pr
    switch(isotope) {
    case 0:    return Phys_const::amu * 140.90766;
    case 141:  return Phys_const::amu * 140.9076576;
    }
    break;

  case Hafnium: //72  Hf
    switch(isotope) {
    case 0:    return Phys_const::amu * 178.49;
    case 174:  return Phys_const::amu * 173.9400461;
    case 176:  return Phys_const::amu * 175.9414076;
    case 177:  return Phys_const::amu * 176.9432277;
    case 178:  return Phys_const::amu * 177.9437058;
    case 179:  return Phys_const::amu * 178.9458232;
    case 180:  return Phys_const::amu * 179.9465570;
    }
    break;

  case Tantalum: //73  Ta
    switch(isotope) {
    case 0:    return Phys_const::amu * 180.94788;
    case 180:  return Phys_const::amu * 179.9474648;
    case 181:  return Phys_const::amu * 180.9479958;
    }	       
    break;

  case Tungsten: //74 W
    switch(isotope) {
    case 0:    return Phys_const::amu * 183.84;
    case 180:  return Phys_const::amu * 179.9467108;
    case 182:  return Phys_const::amu * 181.94820394;
    case 183:  return Phys_const::amu * 182.95022275;
    case 184:  return Phys_const::amu * 183.95093092;
    case 186:  return Phys_const::amu * 185.9543628;
    }	       
    break;

  case Rhenium: //75  Re
    switch(isotope) {
    case 0:    return Phys_const::amu * 186.207;
    case 185:  return Phys_const::amu * 184.9529545;
    case 187:  return Phys_const::amu * 186.9557501;
    }	       
    break;

  case Osmium: //76  Os
    switch(isotope) {
    case 0:    return Phys_const::amu * 190.23;
    case 184:  return Phys_const::amu * 183.9524885;
    case 186:  return Phys_const::amu * 185.9538350;
    case 187:  return Phys_const::amu * 186.9557474;
    case 188:  return Phys_const::amu * 187.9558352;
    case 189:  return Phys_const::amu * 188.9581442;
    case 190:  return Phys_const::amu * 189.9584437;
    case 192:  return Phys_const::amu * 191.9614770;
    }
    break;

  case Iridium: //77  Ir
    switch(isotope) {
    case 0:    return Phys_const::amu * 192.217;
    case 191:  return Phys_const::amu * 190.9605893;
    case 193:  return Phys_const::amu * 192.9629216;
    }
    break;

  case Platinum: //78  Pt
    switch(isotope) {
    case 0:    return Phys_const::amu * 195.084;
    case 190:  return Phys_const::amu * 189.9599297;
    case 192:  return Phys_const::amu * 191.9610387;
    case 194:  return Phys_const::amu * 193.9626809;
    case 195:  return Phys_const::amu * 194.9647917;
    case 196:  return Phys_const::amu * 195.96495209;
    case 198:  return Phys_const::amu * 197.9678949;
    }
    break;

  case Gold: //79  Au
    switch(isotope) {
    case 0:    return Phys_const::amu * 196.966569;
    case 197:  return Phys_const::amu * 196.96656879;
    }
    break;

  case Mercury: //80  Hg
    switch(isotope) {
    case 0:    return Phys_const::amu * 200.592;
    case 196:  return Phys_const::amu * 195.9658326;
    case 198:  return Phys_const::amu * 197.96676860;
    case 199:  return Phys_const::amu * 198.96828064;
    case 200:  return Phys_const::amu * 199.96832659;
    case 201:  return Phys_const::amu * 200.97030284;
    case 202:  return Phys_const::amu * 201.97064340;
    case 204:  return Phys_const::amu * 203.97349398;
    }
    break;

  case Thallium: //81  Tl
    switch(isotope) {
    case 0:    return Phys_const::amu * 204.382;
    case 203:  return Phys_const::amu * 202.9723446;
    case 205:  return Phys_const::amu * 204.9744278;
    }	       
    break;

  case Lead: //82  Pb
    switch(isotope) {
    case 0:    return Phys_const::amu * 207.2;
    case 204:  return Phys_const::amu * 203.9730440;
    case 206:  return Phys_const::amu * 205.9744657;
    case 207:  return Phys_const::amu * 206.9758973;
    case 208:  return Phys_const::amu * 207.9766525;
    }
    break;

  case Bismuth: //83  Bi
    switch(isotope) {
    case 0:    return Phys_const::amu * 208.98040;
    case 209:  return Phys_const::amu * 208.9803991;
    }
    break;

  case Polonium: //84  Po
    switch(isotope) {
    case 0:    return Phys_const::amu * 209.;
    case 209:  return Phys_const::amu * 208.9824308;
    case 210:  return Phys_const::amu * 209.9828741;
    }	       
    break;

  case Astatine: //85  At
    switch(isotope) {
    case 0:    return Phys_const::amu * 210.;
    case 210:  return Phys_const::amu * 209.9871479;
    case 211:  return Phys_const::amu * 210.9874966;
    }
    break;

  case Radon: //86  Rn
    switch(isotope) {
    case 0:    return Phys_const::amu * 222.;
    case 211:  return Phys_const::amu * 210.9906011;
    case 220:  return Phys_const::amu * 220.0113941;
    case 222:  return Phys_const::amu * 222.0175782;
    }
    break;

  default:
    std::cerr << "Atom::mass: unknown atomic number: " << type << "\n";
    throw Error::Range_Err();
  }
  
  std::cerr << "Atom::mass: " << isotope << " isotope is not available for " << name() << "\n";
  throw Error::Range_Err();
}

/**************************** MOLECULE ******************************/

void Molecule::init_dv () {
  for (int i = 0; i < 3; ++i)
    write_cm_pos()[i] = write_cm_vel()[i] = 0.0;
}

/************************ Monoatomic ********************************/

void Monoatomic::get_ang_mom (double* ang_mom) const
{
  for(int i = 0; i < 3; ++i)
    ang_mom[i] = 0.;
}

Monoatomic::Monoatomic (istream& from) {

  const char funame [] = "Monoatomic::Monoatomic: ";
  string comment;
  getline(from, comment);

  int isot;
  string at_type;
  from >> at_type >> isot;
  if(!from) {
    cout << funame <<  "format error\n";
    exit(1);
  }
  atoms.push_back(Atom(at_type, isot));
  getline(from, comment);
	  
  for(int i = 0; i < 3; ++i)
    atoms[atoms.size()-1].mf_pos[i] = 0.;
  
  string stemp;
  while(from >> stemp)
    if(stemp == "End") 
      return;
    else if(stemp == "MultipoleMoments")
      read_multipole_moments(from);
    else {
      cout << funame << "unknown keyword " << stemp << endl;
      exit(1);
    }
  cout << funame << "format error\n";
  exit(1);
}

void Monoatomic::read_multipole_moments (istream& from) 
{
  const char funame [] = "Monoatomic::read_multipole_moments: ";
  string stemp;

  _ismulti = true;
  while(from >> stemp)
    if(stemp == "End")
      return;
    else if(stemp == "Charge") {
      _ischarge = true;
      from >> _charge;
    }
    else if(stemp == "Polarizability") {
      _ispolar = true;
      from >> _polar;
    }
    else if(stemp == "IP") {
      _isip = true;
      from >> _ip;
    }
    else {
      cout << funame << "unknown key word " << stemp << endl;
      exit(1);
    }

  cout << funame << "format error\n";
  exit(1);
}

void Monoatomic::update_atoms ()
{
  for (int i = 0; i < 3; i++)
  {
    begin()->lf_pos[i] = read_cm_pos()[i];
  }
}

void Monoatomic::update_dvd ()
{ 
  // set dynamical variables derivatives pointers
  double* cm_pos_d = dyn_var_der;        // time derivative of position
  double* cm_vel_d = dyn_var_der + 3;    // time derivative of momentum

  for (int i = 0; i < 3; i++)
  {
    // CM coordinates derivatives
    cm_pos_d [i] = read_cm_vel() [i];
    // CM momenta derivatives
    cm_vel_d [i] = begin()->force[i]/mass();
  }
}

double Monoatomic::kin_energy () const
{
  return mass() * ::vdot(read_cm_vel(), read_cm_vel(), 3) / 2.;
}


/************************* Normal mode *************************/

NormalModeCoor& NormalModeCoor::operator= (const NormalModeCoor& n) {
  for(int mode = 0; mode < q0.size(); ++mode)
    q0[mode] = n.q0[mode];
  for(int mode = 0; mode < q1.size(); ++mode)
    q1[mode] = n.q1[mode];
  return *this;
}

void nm2pos(const NormalModeCoor& q, std::vector<D3>& pos)
{
  int all_at = 0;
  for(int frag = 0; frag < 2; ++frag) {
    std::vector<D3> r(mol_array[frag]->size());
    mol_array[frag]->nm2lf(q.data(frag), r);
    for (int at = 0; at < mol_array[frag]->size(); ++at) {
      for (int i = 0; i < 3; ++i)
	pos[all_at][i] = r[at][i];
      ++all_at;
    }
  }
}

/************************* Polyatomic *************************/

Polyatomic::Polyatomic (istream& from) : _isnm(false)
{
  string comment;

  from >> nn;
  if (nn < 2)
    error ("Polyatomic::Polyatomic: wrong number of atoms");
  getline(from, comment);

  int isot;
  string at_type;
  double xyz[3];

  for(int i = 0; i < nn; ++i) {
    from >> at_type >> isot >> xyz [0] >> xyz [1] >> xyz [2];
    getline(from, comment);
    if(!from)
      error("Polyatomic::Polyatomic: input error");
    atoms.push_back(Atom(at_type, isot, xyz));
  }  

  ang_normalized = true;
  set_mass(); 
  set_iner_mom(); 

}

void Polyatomic::set_normal_modes (istream& from)
{
  const char funame [] = "Polyatomic::set_normal_modes: ";

  string comment;
  double dtemp;
  int itemp;

  _isnm = true;

  from >> itemp;
  const int tot_mod_num = itemp;
  
  getline(from, comment);
  nm_array.resize(tot_mod_num, Normal_mode(nn));
  
  int read_mod_num = 0;
  while(read_mod_num < tot_mod_num) {
    std::getline(from, comment);
    if(!from) {
      std::cout << funame << "cannot read frequencies, mode number = " 
		<< read_mod_num << "\n";
      exit(1);
    }
    std::istringstream iss(comment);
    int new_mod_num = 0;
    while(iss >> dtemp) {
      itemp = read_mod_num + new_mod_num;
      if(itemp >= tot_mod_num) {
	std::cout << funame 
		  << "current number of the normal mode frequencies, " 
		  << itemp 
		  << ", is larger than declared normal modes number\n";
	exit(1);
      }
      nm_array[itemp].set_frequency(dtemp * Phys_const::incm);
      ++new_mod_num;
    }
    for(int at = 0; at < nn; ++at)
      for(int i = 0; i < 3; ++i)
	for(int m = 0; m < new_mod_num; ++m) {
	  if(!(from >> dtemp)) {
	    std::cout << funame << "cannot read normal mode coordinate:\n"
		      << ", atom number = " << at << ", axis = " << i 
		      << ", read modes number = " << read_mod_num
		      << ", relative mode number = " << m << "\n";
	    exit(1);
	  }
	  nm_array[read_mod_num+m][at][i] = dtemp;
      }
    read_mod_num += new_mod_num;
  }

  double temp_pos [3];
  //orient normal modes;
  for(int mode = 0; mode < nm_array.size(); ++mode)
    for(int at = 0; at < nn; ++at) {
      vector_matrix_product((const double*)nm_array[mode][at], orig_mfo, temp_pos);
      for(int i = 0; i < 3; ++i)
	nm_array[mode][at][i] = temp_pos[i];
    }

  // normalize normal modes, so M=omega^-2
  double norm;
  for(int mode = 0; mode < tot_mod_num; ++mode) {
    norm = 0.;
    for(int at = 0; at < nn; ++at) {
      dtemp = vdot(nm_array[mode][at], nm_array[mode][at], 3);
      norm += (begin() + at)->mass() * dtemp;
    }   
    norm = sqrt(norm);
    for(int at = 0; at < nn; ++at)
      for(int i = 0; i < 3; ++i)
	nm_array[mode][at][i] /= norm * nm_array[mode].frequency();
  }

  //check on normal modes
  // normal mode cross products
  double max_val = 0.;
  for(int mode = 0; mode < nm_array.size(); ++mode)
    for(int mode1 = mode + 1; mode1 < nm_array.size(); ++mode1) {
      double prod = 0.;
      for(int at = 0; at < nn; ++at)
	prod += (begin()+at)->mass() * 
	  vdot(nm_array[mode][at], nm_array[mode1][at], 3);

      prod *= nm_array[mode].frequency() * nm_array[mode1].frequency();
      if(!mode && mode1 == mode + 1 || max_val < prod)
	max_val = prod;
      if(abs(prod) > .001) {
	cout << funame << "WARNING!!! cross product of the normal modes "
	     << mode << " and " << mode1 << " is " << prod << "\n";
	//exit(1);
      }
    }
  //cout << funame << "maximum cross product of the normal modes is "
  //     << max_val << "\n";

  // normal modes and rotations products
  norm = 0.;
  for(Ater at = begin(); at != end(); ++at)
    norm += vdot(at->mf_pos, at->mf_pos, 3) * at->mass();
  norm = sqrt(norm);
  for(int mode = 0; mode < nm_array.size(); ++mode) {
    Array<double> prod(3);
    prod.init();
    Array<double> pos(3);
    for(int at = 0; at < nn; ++at) {
      vector_product((begin()+at)->mf_pos,(const double*)nm_array[mode][at], pos.data());
      pos *= (begin() + at)->mass();
      prod += pos;
    }
    prod *= nm_array[mode].frequency() / norm;
    for(int i = 0; i < 3; ++i) {
      if(!mode && !i || max_val < abs(prod[i]))
	max_val = abs(prod[i]);
      if(abs(prod[i]) > 0.001) {
	cout << funame << "WARNING!!! " <<  mode <<  "-th mode and " << i
	     << "-th rotation cross product is "
	     << prod[i] << "\n";
	//exit(1);
      }
    }
  }
  //cout << funame << "maximum normal mode and rotation cross product is "
  //     << max_val << "\n";

  // translations
  for(int mode = 0; mode < nm_array.size(); ++mode) {
    norm = 0.;
    for(int at = 0; at < nn; ++at)
      norm += vdot(nm_array[mode][at], nm_array[mode][at], 3) 
	* (begin() + at)->mass() * (begin() + at)->mass();
    norm = sqrt(norm);
    Array<double> pos(3);
    pos.init();
    for(int at = 0; at < nn; ++at)
      for(int i = 0; i < 3; ++i)
	pos[i] += (begin() + at)->mass() * nm_array[mode][at][i];
    pos /= norm;
    for(int i = 0; i < 3; ++i) {
      if(!mode && !i || max_val < abs(pos[i]))
	max_val = abs(pos[i]);
      if(abs(pos[i]) > 0.001) {
	cout << funame << "WARNING!!! " << mode << "-th mode " 
	     << i <<"-th axis displacement " << i << " is "
	     << pos[i] << "\n";
	//exit(1);
      }
    }
  }
  //cout << funame << "maximum normal mode displacement is "
  //     << max_val << "\n";
}

void Polyatomic::set_mass ()
{
  mm = 0.0;
  for (Ater at = begin(); at != end(); ++at)
    mm += at->mass();
}

void Polyatomic::set_iner_mom ()
{
  double dtemp;
  double temp_pos [3];
  // center of mass position
  for (int i = 0; i < 3; ++i)
    temp_pos [i] = 0.0;

  for (Ater at = begin(); at != end(); ++at)
    for (int i = 0; i < 3; ++i)
      temp_pos [i] += at->mass() * at->lf_pos [i];
 
  for (int i = 0; i < 3; ++i)
    temp_pos [i] /= mass ();

  // atom positions relative to CM
  for (Ater at = begin(); at != end(); ++at)
    for (int i = 0; i < 3; ++i)
      at->rel_pos [i] = at->lf_pos [i] - temp_pos [i];
   
  // inertia moment matrix
  orig_mfo.base().init();
  for (Ater at = begin(); at != end(); ++at)
      for (int i1 = 0; i1 < 3; ++i1)
	for (int i2 = i1; i2 < 3; ++i2)
	  orig_mfo (i1, i2) += - at->mass() * at->rel_pos [i1] * at->rel_pos [i2];
  orig_mfo.diagonal(0) -= orig_mfo.trace();

  // principal inertia moments and rotational matrix
  orig_mfo.eigenv(ang_mass);

  // check the chirality
  dtemp = 0.;
  for(int i0 = 0; i0 < 3; ++i0) {
    int i1 = (i0 + 1) % 3;
    int i2 = (i0 + 2) % 3;
    dtemp += orig_mfo(0, i0) * (orig_mfo(1, i1) * orig_mfo(2, i2) - orig_mfo(1, i2) * orig_mfo(2, i1));
  }

  if(dtemp < 0.)
    for(int i = 0; i < 3; ++i)
      orig_mfo(i, 0) = -orig_mfo(i, 0);

  // molecular frame coordinates
  for (Ater at = begin(); at != end(); ++at)
    vector_matrix_product(at->rel_pos, orig_mfo, at->mf_pos);
}

double Polyatomic::force (int i) const
{
    double res = 0.0;
    for (Const_ater at = begin(); at != end(); ++at)
      res  += at->force[i];
    return res;
}

double Polyatomic::torque (int i) const
{
    double res = 0.0;
    for (Const_ater at = begin(); at != end(); ++at)
	res += vector_product (i, at->rel_pos, at->force);
    return res;
}

double Polyatomic::nm_force (int mode) const
{
  double res = 0.;
  double mf_force [3];
  for(int at = 0; at < size(); ++at) {
    lf2mf((begin()+at)->force, mf_force);
    res += vdot(nm_array[mode][at], mf_force, 3);
  }
  return res;
}

void Polyatomic::nm2lf (const double* q, std::vector<D3>& r) const
{
  for(int at = 0; at < size(); ++at)
    for(int i = 0; i < 3; ++i) {
      r[at][i] = 0; 
      for(int mode = 0; mode < nm_array.size(); ++mode)
	r[at][i] += q[mode] * nm_array[mode][at][i];
    }

  double pos [3];
  for(int at = 0; at < size(); ++at) {
    mf2lf(r[at], pos);
    for(int i = 0; i < 3; ++i)
      r[at][i] = (begin() + at)->lf_pos[i] + pos[i];
  }
}
void Polyatomic::read_multipole_moments (istream& from) 
{
  const char funame [] = "Polyatomic::read_multipole_moments: ";
  string stemp;

  _ismulti = true;
  while(from >> stemp)
    if(stemp == "End")
      return;
    else if(stemp == "Dipole") {
      _isdipole = true;
      from >> _dipole;
      getline(from, stemp);
    }
    else if(stemp == "Quadrupole") {
      _isquad = true;
      from >> _quad;
      getline(from, stemp);
    }
    else if(stemp == "Polarizability") {
      _ispolar = true;
      from >> _polar;
    }
    else if(stemp == "Charge") {
      _ischarge = true;
      from >> _charge;
    }
    else if(stemp == "IP") {
      _isip = true;
      from >> _ip;
    }
    else {
      cout << funame << "unknown key word " << stemp << endl;
      exit(1);
    }

  cout << funame << "format error\n";
  exit(1);
}

/************************** Linear *********************************/

Linear::Linear (istream& from) : Polyatomic (from), ang_vel_orthogonal(true) 
{
  const char funame [] = "Linear::Linear: ";

  string stemp;
  while(from >> stemp)
    if(stemp == "End") 
      return;
    else if(stemp == "NormalModes")
      set_normal_modes(from);
    else if(stemp == "MultipoleMoments")
      read_multipole_moments(from);
    else {
      cout << funame << "unknown keyword " << stemp << endl;
      exit(1);
    }
  cout << funame << "format error\n";
  exit(1);
}

void Linear::get_ang_mom (double* ang_mom) const
{
  for(int i = 0; i < 3; ++i)
    ang_mom[i] = ang_mass[2] * read_ang_vel()[i];
  ::orthogonalize(ang_mom, read_ang_pos(), 3);
}

void Linear::init_dv () 
{
  Molecule::init_dv();

  write_ang_pos() [0] = 1.;
  for (int i = 1; i < ang_pos_size(); ++i)
    write_ang_pos()[i] = 0.;
  for (int i = 0; i < ang_vel_size(); ++i)
    write_ang_vel()[i] = 0.;
}

void Linear::update_atoms ()
{
  // normalization for the axis director
  const double norm_min = 0.1;
  const double norm_max = 10.0;

  double temp;
  double norm = 0.0;
  for (int i = 0; i < 3; i++)
  {
    temp = read_ang_pos() [i];
    norm += temp*temp;
  }
  if (norm < norm_min || norm > norm_max)
    {
      cout << "Linear::update_atoms: unity vector length = " << norm << endl;
      ang_normalized = false;
    }

  norm = sqrt(norm);
  // update atoms coordinates
  double factor;
  for (Ater at = begin (); at != end (); ++at)
    {
      factor = *(at->mf_pos)/norm;
      for (int i = 0; i < 3; i++)
	{
	  temp = factor*read_ang_pos()[i];
	  at->rel_pos [i] = temp;
	  at->lf_pos  [i]  = read_cm_pos() [i] + temp; 
	}
    }

  set_rot_mat();
}

void Linear::set_rot_mat()
{
  double rot_mat[3][3];
  for(int i = 0; i < 3; ++i)
    rot_mat[0][i] = read_ang_pos()[i];
  ::normalize(rot_mat[0], 3);

  rot_mat[1][0] = 0.;
  rot_mat[1][1] =-rot_mat[0][2];
  rot_mat[1][2] = rot_mat[0][1];
  if(::normalize(rot_mat[1], 3) < 1.e-14) {
    rot_mat[1][1] = 1.;
    rot_mat[1][2] = 0.;
  }
  vector_product(rot_mat[0], rot_mat[1], rot_mat[2]);

  for(int i = 0; i < 3; ++i)
    mfo.row(i) = rot_mat[i];
}

void Linear::update_dvd ()
{ 

  double*  cm_pos_d  = dyn_var_der;
  double*  cm_vel_d  = dyn_var_der + 3;
  double*  ang_pos_d = dyn_var_der + 6;   // derivatives of orientation coordinates
  double*  ang_vel_d = dyn_var_der + 9;   // derivatives of rotational frequencies

  for (int i = 0; i < 3; i++)
  {
    // Center of mass coordinates derivatives
    cm_pos_d [i] = read_cm_vel() [i];
    // Molecular axis director derivative
    ang_pos_d [i] = vector_product (i, read_ang_vel(), read_ang_pos());
    // Center of mass velocity derivative
    cm_vel_d  [i] = force(i)/mass();
    // Angular velocity derivative
    ang_vel_d [i] = torque(i)/iner_mom(0);
  } 
}

void Linear::normalize ()
{
  double temp, norm;
  for (int i = 0; i < ang_pos_size(); ++i)
    {
      temp = read_ang_pos() [i];
      norm += temp*temp;
    }
  norm = sqrt(norm);
  for (int i = 0; i < ang_pos_size(); ++i)
    write_ang_pos() [i] /= norm;
}

double Linear::kin_energy () const
{
  const char funame [] = "Linear::kin_energy: ";
  const double max_ratio = 0.01;

  // rotational part
  double av [3];
  for(int i = 0; i < 3; ++i)
    av[i] = read_ang_vel()[i];
  if(::orthogonalize(av, read_ang_pos(), 3) > max_ratio) {
    cout << funame << "angular velocity needs orthogonalization\n";
    ang_vel_orthogonal = false;
  }

  return (iner_mom(2) * vdot(av, av, 3) 
	 // translational part
	 + mass() * vdot(read_cm_vel(),read_cm_vel(), 3)) / 2.;
}

void Linear::orthogonalize ()
{
   double temp;

   //orthogonal part of angular velocity
   double norm = 0.0, proj = 0.0;
   for (int i = 0; i < 3; ++i)
   {
      temp = read_ang_pos()[i];
      norm += temp*temp;
      proj += temp*read_ang_vel()[i];
   }
   proj /= norm;
   
   for (int i = 0; i < 3; ++i)
      write_ang_vel()[i] -= read_ang_pos()[i]*proj;
}

double Linear::imm (int i, int j) const
{
  double val = 0.0;
  for (int k = 0; k < 3; ++k)
    val += read_ang_pos()[k] * read_ang_pos()[k];

  if (i == j)
    return iner_mom(2) * (1.0 - read_ang_pos()[i] * read_ang_pos()[j] / val);
  else
    return -iner_mom(2) * read_ang_pos()[i] * read_ang_pos()[j] / val;
}

/************************** Nonlinear *********************************/

Nonlinear::Nonlinear (istream& from) : Polyatomic (from)
{ 
  const char funame [] = "Nonlinear::Nonlinear: ";

  euler_coef [0] = (iner_mom(1) - iner_mom(2)) / iner_mom(0);
  euler_coef [1] = (iner_mom(2) - iner_mom(0)) / iner_mom(1);
  euler_coef [2] = (iner_mom(0) - iner_mom(1)) / iner_mom(2);

  string stemp;
  while(from >> stemp)
    if(stemp == "End") 
      return;
    else if(stemp == "NormalModes")
      set_normal_modes(from);
    else if(stemp == "MultipoleMoments")
      read_multipole_moments(from);
    else {
      cout << funame << "unknown keyword " << stemp << endl;
      exit(1);
    }
  cout << funame << "format error\n";
  exit(1);
}

void Nonlinear::get_ang_mom (double* ang_mom) const
{
  double mam [3];
  for(int i = 0; i < 3; ++i)
    mam[i] = ang_mass[i] * read_ang_vel()[i];
  mf2lf(mam, ang_mom);
}

void Nonlinear::init_dv () 
{
  Molecule::init_dv();

  write_ang_pos() [0] = 1.;
  for (int i = 1; i < ang_pos_size(); ++i)
    write_ang_pos()[i] = 0.;
  for (int i = 0; i < ang_vel_size(); ++i)
    write_ang_vel()[i] = 0.;
}

void Nonlinear::update_atoms ()
{
  // update molecular frame orientation matrix
  if (ang_normalized != quat2mat(read_ang_pos(), mfo))
    {
      ang_normalized = !ang_normalized;
      cout << "Nonlinear::update_atoms: ang_normalized has changed to "
	   << ang_normalized << endl;
    }

  // update atoms coordinates
  for (Ater at = begin (); at != end (); ++at)
    {
      // relative coordinates
      vector_matrix_product(at->mf_pos, mfo, at->rel_pos);
      // laboratory frame coordinates
      for (int i = 0; i < 3; ++i)
	at->lf_pos [i] = at->rel_pos [i] + read_cm_pos() [i];
    }
}

void Nonlinear::update_dvd ()
{ 
  double* cm_pos_d  = dyn_var_der;
  double* cm_vel_d  = dyn_var_der + 3;
  double* ang_pos_d = dyn_var_der + 6;
  double* ang_vel_d = dyn_var_der + 10;

  double lf_torque [3];
  double mf_torque [3];

  for (int i = 0; i < 3; i++)
  {
    cm_pos_d  [i] = read_cm_vel() [i];          // CM coordinates derivatives
    cm_vel_d  [i] = force(i)/mass();     // CM velocity derivatives
    lf_torque [i] = torque(i);           // LF torque
  }

  matrix_vector_product (mfo, lf_torque, mf_torque); // transform to MF

  register double ox = read_ang_vel() [0], 
                  oy = read_ang_vel() [1], 
                  oz = read_ang_vel() [2];

  register double q0 = read_ang_pos() [0], 
                  q1 = read_ang_pos() [1], 
                  q2 = read_ang_pos() [2], 
                  q3 = read_ang_pos() [3];

  // quaternion derivatives

  ang_pos_d [0] = (- q1 * ox - q2 * oy - q3 * oz) / 2.0;
  ang_pos_d [1] = (  q0 * ox - q3 * oy + q2 * oz) / 2.0;
  ang_pos_d [2] = (  q3 * ox + q0 * oy - q1 * oz) / 2.0;
  ang_pos_d [3] = (- q2 * ox + q1 * oy + q0 * oz) / 2.0;


  // rotational frequencies derivatives

  ang_vel_d [0] = mf_torque [0] / iner_mom (0) + euler_coef [0] * oy * oz;
  ang_vel_d [1] = mf_torque [1] / iner_mom (1) + euler_coef [1] * oz * ox;
  ang_vel_d [2] = mf_torque [2] / iner_mom (2) + euler_coef [2] * ox * oy;
}

void Nonlinear::normalize ()
{
  double temp, norm;
  for (int i = 0; i < ang_pos_size(); ++i)
    {
      temp = read_ang_pos() [i];
      norm += temp*temp;
    }
  norm = sqrt(norm);
  for (int i = 0; i < ang_pos_size(); ++i)
    write_ang_pos() [i] /= norm;

  cout << "Nonlinear::normalize: quaternion is normalized" << endl;

}

double Nonlinear::kin_energy () const
{
  // rotational part
   double res = 0.;
   for (int i = 0; i < 3; ++i)
     res += iner_mom(i) * read_ang_vel()[i] * read_ang_vel()[i];

   // translational part
   res += mass() * ::vdot(read_cm_vel(), read_cm_vel(), 3);
   res /= 2.0;

   return res;
}

double Nonlinear::imm (int i, int j) const
{
  double res = 0;
  for (int k = 0; k < 3; ++k)
    res += iner_mom(k) * mfo(k,i) * mfo(k,j);
  return res;
}

