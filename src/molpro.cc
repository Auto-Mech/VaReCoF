#include "molpro.hh"
#include "system.hh"
#include "log.hh"
#include "units.hh"
#include "math.hh"

#include <cctype>
#include <unistd.h>
#include <cstdio>
#include <iomanip>
#include <dlfcn.h>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <cstring>

// available potential corrections
namespace Molpro {

  extern "C" {
    typedef double  (*corr_t) (const int&, const double [], 
			       const double [], const int []);
  }
  struct Corr {
    corr_t sym;// potential correction symbol pointer
    std::vector<int>    mpe_index;// molpro energy indices
    std::vector<double> rpar_data;
    std::vector<int>    ipar_data;
    double* rpar;// array of real parameters
    int* ipar;// array of integer parameters

    Corr(corr_t =0);
    ~Corr();
    void set_par ();
  };
  
  vector<Corr> corr_table;

  bool is_init = false;
  std::vector<Atom> atoms;
  std::vector<D3> forces;
  std::vector<double> energies;

  int debug = Log::INFO;          // debug output

  string wf_name; // wave function file name
  string exe_name = "/usr/local/bin/molpro";

  string geom_patt = "geometry";// geometry pattern in the template file
  string ener_patt = "molpro_energy";// energy pattern in the output
  string forc_patt = "molpro_forces";
  string fail_patt = "molpro_failure";
  string symmetry = "nosym";

  vector<string> base_name;// input-output file names
  vector<string> prolog;
  vector<string> epilog;
  vector<int> fail_counter;

  const string tml_sfx = ".tml"; // template file name extension
  const string inp_sfx = ".inp"; // input file name extension
  const string bac_sfx = ".bac"; // backup file name extension

  const char symm_key [] = "Symmetry";
  const char exec_key [] = "MolproExecFile";
  const char geom_key [] = "GeometryPattern";
  const char ener_key [] = "EnergyPattern";
  const char fail_key [] = "FailurePattern";
  const char base_key [] = "BaseFileName";
  const char debug_key[] = "DebugLevel";
  const char csym_key [] = "PotentialCorrectionTable";
  const char rpar_key [] = "ParameterReal";
  const char ipar_key [] = "ParameterInteger";
  const char  end_key [] = "End";
  const char forc_key [] = "ForcePattern";

  const char* all_keys [] = {
    symm_key, exec_key, geom_key, ener_key, fail_key, base_key, debug_key, 
    csym_key, rpar_key, ipar_key, end_key,  forc_key, 0
  };

  // backup molpro output file
  //
  void backup (int = 0);
  ostream& operator<< (ostream&, const Atom&);
}

Molpro::Corr::Corr(corr_t pc) : sym(pc), rpar(0), ipar(0) {}
Molpro::Corr::~Corr() { 
  if(rpar)
    delete[] rpar; 
  if(ipar)
    delete[] ipar;
}
void Molpro::Corr::set_par () {
  if(rpar_data.size()) {
    if(rpar)
      delete[] rpar;
    rpar = new double[rpar_data.size()];
    for(int i = 0; i < rpar_data.size(); ++i)
      rpar[i] = rpar_data[i];
  }
  if(ipar_data.size()) {
    if(ipar)
      delete[] ipar;
    ipar = new int[ipar_data.size()];
    for(int i = 0; i < ipar_data.size(); ++i)
      ipar[i] = ipar_data[i];
  }
}


void Molpro::backup (int method)
{
  const char funame [] = "Molpro::backup: ";

  std::ostringstream back_name;
  back_name << base_name[method] << bac_sfx << "_" << fail_counter[method];
  ++fail_counter[method];
  const string out_name  = base_name[method] + ".out";
  if(rename(out_name.c_str(), back_name.str().c_str())) {
    Log::out << funame << "failed" << std::endl;
    throw Init_Err();
  }
  Log::out << funame << "check " << back_name.str() << " file" << std::endl;
}

ostream& Molpro::print_geom(ostream& to) {
  for(int i = 0; i < atoms.size(); ++i)
    to << atoms[i];
  return to;
}

void Molpro::print_dist ()
{
  int old = Log::out.precision(3);

  Log::out << "interatomic distances matrix (angstrom):\n\n";
  Log::out << "\\ ";
  for(int i = 0; i < atoms.size()-1; ++i)
    Log::out << setw(5) << atoms[i].name.c_str() << i+1;
  Log::out << "\n";

  for(int j = 1; j < atoms.size(); ++j) {
    Log::out << atoms[j].name << j+1;
    for(int i = 0; i < j; ++i)
      Log::out << setw(6) << Phys_const::bohr * 
	vdistance(atoms[i].pos, atoms[j].pos, 3);
    Log::out << "\n";
  }
  Log::out << "\n";
  Log::out.precision(old);
}

ostream& Molpro::operator<< (ostream& to, const Atom& a)
{
  to << a.name;
  for (int i = 0; i < 3; ++i)
    to << ", " << Phys_const::bohr * a.pos[i];
  to << "\n";
  return to;
}

bool Molpro::isinit () { return is_init; }

void Molpro::init(const string& fname)
{
  const char funame [] = "Molpro::init: ";

  if(is_init) {
    std::cout << funame << "molpro environment is allready initialized\n";
    throw Init_Err();
  }
  is_init = true;

  string token, stemp;
  int itemp;
  double dtemp;

  ifstream from(fname.c_str());
  if(!from) {
    std::cout << funame << "cannot open input file " << fname << "\n";
    throw Open_Err();
  }

  bool is_geom = false;
  bool is_ener = false;
  bool is_fail = false;
  bool is_exec = false;
  bool is_debug = false;
  bool is_symm = false;
  while(from >> token)
    if(token == geom_key) {
      if(is_geom) {
	std::cout << funame << "another " << token << " declaration\n";
	throw Form_Err();
      }
      is_geom = true;

      from >> geom_patt;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
    }
    else if(token == symm_key) {
      is_symm = true;
      from >> symmetry;
    }
    else if(token == forc_key) {
      from >> forc_patt;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
      for(int i = 0; i < forc_patt.size(); ++i)
	forc_patt[i] = toupper(forc_patt[i]); // convert to upper case
    }
    else if(token == ener_key) {
      if(is_ener) {
	std::cout << funame << "another " << token << " declaration\n";
	throw Form_Err();
      }
      is_ener = true;

      from >> ener_patt;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
      for(int i = 0; i < ener_patt.size(); ++i)
	ener_patt[i] = toupper(ener_patt[i]); // convert to upper case
    }
    else if(token == fail_key) {
      if(is_fail) {
	std::cout << funame << "another " << token << " declaration\n";
	throw Form_Err();
      }
      is_fail = true;

      from >> fail_patt;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
      for(int i = 0; i < fail_patt.size(); ++i)
	fail_patt[i] = toupper(fail_patt[i]); // convert to upper case
    }
    else if(token == base_key) {
      getline(from, token);
      istringstream iss(token);

      while(iss >> stemp)
	base_name.push_back(stemp);
    }
    else if(token == exec_key) {
      if(is_exec) {
	std::cout << funame << "another " << token << " declaration\n";
	throw Form_Err();
      }
      is_exec = true;

      from >> exe_name;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
      if(access(exe_name.c_str(), X_OK)) {
	std::cout << funame << "no access to " << exe_name << "\n";
	throw Init_Err();
      }
    }
    else if(token == debug_key) {
      if(is_debug) {
	std::cout << funame << "another " << token << " declaration\n";
	throw Form_Err();
      }
      is_debug = true;

      from >> debug;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }
    }
    else if(token == csym_key) {
      from >> token; // read library name
      void* clib = dlopen(token.c_str(), RTLD_NOW);
      const char* messg = dlerror();
      if(messg) {
	cout << funame << messg << "\n";
	throw Run_Err();
      }
      std::string symbol;
      while(from >> token) {// read library symbols
	if(token == end_key)
	  break;
	else if(token == rpar_key) {//real params for the previous symbol
	  if(!corr_table.size()) {
	    std::cout << funame 
		      << "parameters definition should follow the symbol\n";
	    throw Form_Err();
	  }
	  getline(from, token);
	  std::istringstream iss(token);
	  while(iss >> dtemp)
	    corr_table.rbegin()->rpar_data.push_back(dtemp);
	  if(!iss.eof() && iss.fail()) {
	    std::cout << funame << "the list of real parameters for " 
		      << symbol 
		      << " correction symbol is corrupted\n";
	    throw Form_Err();
	  }
	}
	else if(token == ipar_key) {//integer params for the previous symbol
	  if(!corr_table.size()) {
	    std::cout << funame 
		      << "parameters definition should follow the symbol\n";
	    throw Form_Err();
	  }
	  getline(from, token);
	  std::istringstream iss(token);
	  while(iss >> itemp)
	    corr_table.rbegin()->ipar_data.push_back(itemp);
	  if(!iss.eof() && iss.fail()) {
	    std::cout << funame << "the list of integer parameters for " 
		      << symbol 
		      << " correction symbol is corrupted\n";
	    throw Form_Err();
	  }
	}
	else {//assume that it is a library symbol
	  symbol = token;
	  corr_table.push_back(Corr((corr_t)dlsym(clib, symbol.c_str())));
	  messg = dlerror();
	  if(messg) {
	    cout << funame << messg << "\n";
	    throw Run_Err();
	  }
	  getline(from, token);
	  std::istringstream iss(token);
	  while(iss >> itemp) {// read molpro energies indices 
	    if(itemp < 1) {
	      std::cout << funame << "the index value " << itemp 
			<< " for the " << symbol 
			<< " symbol in the correction table is out of range\n";
	      throw Form_Err();
	    }
	    corr_table.rbegin()->mpe_index.push_back(itemp-1);
	  }
	  if(!iss.eof() && iss.fail()) {
	    std::cout << funame << "the list of indices for the "
		      << symbol 
		      << " symbol in the correction table is corrupted\n";
	    throw Form_Err();
	  }
	  if(!corr_table.rbegin()->mpe_index.size())
	    corr_table.rbegin()->mpe_index.push_back(0);
	}
      }
      if(!from) {
	std::cout << funame << "cannot read " << csym_key << "\n";
	throw Form_Err();
      }
    }
  // unknown key
    else {
      std::cout << funame << "unknown key: " << token
		<< "\n\navailable keys:\n";
      for(int i = 0; all_keys[i] != 0; ++i)
	std::cout << all_keys[i] << "\n";
      throw Form_Err();
    }
  from.close();

  for(int corr = 0; corr < corr_table.size(); ++corr)
    corr_table[corr].set_par();

  if(debug >= Log::WARN) {
    if(!corr_table.size())
      std::cout << funame << csym_key << " is not found in "
		<< fname << ". Potential correction will not be used\n";
    if(!is_geom) 
      std::cout << funame << geom_key << " is not found in "
		<< fname << "; using the default: " 
		<< geom_patt << "\n";

    if(!is_ener) 
      std::cout << funame << ener_key << " is not found in "
		<< fname << "; using the default: " 
		<< ener_patt << "\n";

    if(!is_fail) 
      std::cout << funame << fail_key << " is not found in "
		<< fname << "; using the default: " 
		<< fail_patt << "\n";

    if(!is_debug) 
      std::cout << funame << debug_key << " is not found in "
		<< fname << "; using the default: " 
		<< debug << "\n";
  
    if(!is_exec) 
      std::cout << funame << exec_key << " is not found in "
		<< fname << "; using the default: " 
		<< exe_name << "\n";
    
  }

  if(!base_name.size()) {
    std::cout << funame << "no templates are defined\n";
    throw Form_Err();
  }

  fail_counter.resize(base_name.size());
  for(int i = 0; i < fail_counter.size(); ++i)
    fail_counter[i] = 0;

  wf_name = base_name[0] + ".wf";

  // read template files
  prolog.resize(base_name.size());
  epilog.resize(base_name.size());
  for(int method = 0; method < base_name.size(); ++method) {
    string tname = base_name[method] + tml_sfx;
    from.clear();
    from.open(tname.c_str());
    if(!from) {
      std::cout << funame << "cannot open input file " << tname << "\n";
      throw Open_Err();
    }

    bool is_prolog = true;
    while(getline(from,token)) {
      if(is_prolog) {
	istringstream iss(token);
	iss >> stemp;
	if(stemp == geom_patt) {
	  is_prolog = false;
	  continue;
	}
	prolog[method] += token + "\n";
      }
      else
	epilog[method] += token + "\n";
    }
    if(is_prolog) {
      std::cout << funame << "geometry pattern " << geom_patt 
		<< " is not found in " << tname << "\n";
      throw Open_Err();
    }
    from.close();
  }
}

//void Molpro::pot (Array<double>& resultp, int method, int flags) throw(Err)
void Molpro::pot (int method, int flags)
{
  const char funame [] = "Molpro::pot: ";

  if(method < 0 || method >= base_name.size()) {
    std::cout << funame  << "template index (" << method <<") out of range\n";
    throw Range_Err();
  }
  energies.clear();

  double dtemp;
  int itemp;
  string stemp;

  // create an input file for molpro
  const string inp_name = base_name[method] + inp_sfx;
  ofstream to(inp_name.c_str());
  if(!to) {
    Log::out << funame << "cannot open " << inp_name << "\n";
    throw Open_Err();
  }
  to << prolog[method];

  if(flags & WF) {
    to << "file, 2, " << wf_name << "\n";
  }
  if(symmetry == "nosym")
    to << "symmetry, nosym; orient, noorient\n";
  else if(symmetry != "default")
    to << "symmetry, " << symmetry << "; orient, mass\n";
  else
    to << "orient, mass\n";
  to << "angstrom\ngeomtyp=xyz\ngeometry={\n";
  to << atoms.size() << "\nrotational dynamics\n";
  print_geom(to) << "}\n" << epilog[method];
  to.close();

  // run molpro executable
  remove((base_name[method] + ".out_1").c_str());
  if(call_exe(exe_name.c_str(), inp_name.c_str(), (char*) 0)) {
    Log::out << funame << "molpro failed\n";
    backup(method);
    throw Run_Err();
  }

  // read the molpro results
  const string out_name  = base_name[method] + ".out";
  ifstream from(out_name.c_str());
  if(!from) {
    Log::out << funame << "cannot open " << out_name << "\n";
    throw Open_Err();
  }

  string token, line, comment;
  vector<double> ener_array;

  bool is_force = false;

  const char* nptr;
  char* endptr;

  while(from >> token) {// scan molpro out
    if(token == ener_patt) {
      std::getline(from, line);
      std::istringstream iss(line);

      // read "=" sign
      iss >> stemp;

      // read energy
      if(!(iss >> stemp)) {
	Log::out << funame << "cannot read the " << ener_array.size() 
		 << "-th energy from " << out_name << "\n";
	backup(method);
	throw Form_Err();
      }

      
      // remove optional "D" sign
      for(int s = 0; s < stemp.size(); ++s)
	if(stemp[s] == 'D' || stemp[s] == 'd')
	  stemp[s] = 'e';

      nptr = stemp.c_str();
      dtemp = std::strtod(nptr, &endptr);
      if(std::strlen(endptr)) {
	Log::out << funame << "energy: cannot convert properly string to double: the residue is " << endptr << "\n";
	backup(method);
	throw Form_Err();
      }

      ener_array.push_back(dtemp);

      if(debug>=Log::INFO)
	Log::out << funame << ener_array.size() << "-th molpro energy = "
		 << dtemp / Phys_const::kcal << " kcal/mole\n";

    }
    if(token == fail_patt) {
      Log::out << funame << "internal molpro failure" << "\n";
      backup(method);
      throw Run_Err();
    }
    // energy gradient for mcscf and rhf using keyword forces
    else if(!is_force && flags & FORCE && token == forc_patt) {//read gradients
      is_force = true;

      for(int i = 0; i < 4; ++i)
	getline(from, token);
      
      forces.resize(atoms.size());
      for(int at = 0; at < atoms.size(); ++at)
	for(int i = 0; i < 3; ++i) {

	  from >> stemp;
	  for(int s = 0; s < stemp.size(); ++s)
	    if(stemp[s] == 'D' || stemp[s] == 'd')
	      stemp[s] = 'e';

	  nptr = stemp.c_str();
	  dtemp = std::strtod(nptr, &endptr);
	  if(std::strlen(endptr)) {
	    Log::out << funame << "forces: cannot convert properly string to double: the residue is " << endptr << "\n";
	    backup(method);
	    throw Form_Err();
	  }

	  forces[at][i] = -dtemp;
	}
      
      if(!from) {
	Log::out << funame 
		   << "cannot read gradients from the molpro output file "
		   << out_name << "\n";
	  backup(method);
	  throw Form_Err();
      }
    }// read gradients
    
  }// scan molpro out
  if(debug>=Log::INFO)
    Log::out << "\n";

  if(flags & FORCE && !is_force) {
    Log::out << funame << "energy gradients were not found in " << out_name
	     << "\n";
    backup(method);
    throw Form_Err();
  }

  // energy correction
  Array_2<double> coord(3, atoms.size());
  if(corr_table.size()) {
    for(int at = 0; at < atoms.size(); ++at)
      for(int i = 0; i < 3; ++i)
	coord(i, at) = atoms[at].pos[i];
  }
  
  int ener_count = 0;
  if(corr_table.size())
    for(int ci = 0; ci < corr_table.size(); ++ci) {
      double encorr = corr_table[ci].sym(atoms.size(), coord.data(), 
					 corr_table[ci].rpar, 
					 corr_table[ci].ipar);
      if(debug>=Log::INFO)
	Log::out << funame << ci << "-th energy correction = " 
		 << encorr / Phys_const::kcal << " kcal/mole\n";
      for(int ei = 0; ei < corr_table[ci].mpe_index.size(); ++ei) {
	/*if(ener_count == resultp.size()) {
	  Log::out << funame 
	  << "the number of generated energies exceed "
	  "the number of allocated energies\n";
	  backup(method);
	  throw Form_Err();
	  }
	*/
	if(corr_table[ci].mpe_index[ei] >= ener_array.size()) {
	  Log::out << funame << ei << "-th index of the " << ci 
		   << "-th correction " << corr_table[ci].mpe_index[ei] 
		   << " exceeds the generated number of energies "
		   << ener_array.size() << "\n";
	  throw Form_Err();
	}
	//resultp[ener_count++] = ener_array[corr_table[ci].mpe_index[ei]] 
	//  + encorr;
	energies.push_back(ener_array[corr_table[ci].mpe_index[ei]] 
			   + encorr);
      }
    }
  else
    for(;ener_count < ener_array.size(); ++ener_count) {
      /*if(ener_count == resultp.size()) {
	Log::out << funame 
	<< "the number of generated energies exceed "
	"the number of allocated energies\n";
	backup(method);
	throw Form_Err();
	}
	resultp[ener_count] = ener_array[ener_count];
      */
      energies.push_back(ener_array[ener_count]);
    }

  /*
    if(ener_count != resultp.size()) {
    Log::out << funame << "the number of generated  energies, " 
    << ener_count  << ", is not consistent with the number of"
    " allocated energies, " << resultp.size() << "\n";
    backup(method);
    throw Form_Err();
    }
  */
}
