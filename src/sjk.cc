#include "sjk.hh"
#include "tmatrix.hh"
#include "pes.hh"
#include "rotd.hh"

#include <fstream>
#include <dlfcn.h>
#include <cstdlib>
#include <vector>
#include <sstream>

namespace Sjk {

  enum {NEW_FORM, SJK_FORM};

  extern "C" {
    typedef double (*ener_t) (const double* r, const double* rpar, const int* ipar, int& status);
    typedef double  (*sjk_t) (const double* r, const double* rpar, const int* ipar);
    typedef void   (*init_t) (const char*);
  }

  int _pot_format = NEW_FORM;
  int _cor_format = NEW_FORM;

  void* _pot_ener = 0;
  void* _cor_ener = 0;

  double* _pot_rpar = 0;
  int*    _pot_ipar = 0;
  double* _cor_rpar = 0;
  int*    _cor_ipar = 0;

  bool       _init = false;

  const char pot_form_key [] = "Format";
  const char pot_libr_key [] = "Library";
  const char pot_ener_key [] = "EnergyMethod";
  const char pot_init_key [] = "InitMethod";
  const char pot_data_key [] = "InitData";
  const char pot_rpar_key [] = "ParameterReal";
  const char pot_ipar_key [] = "ParameterInteger";

  const char cor_form_key [] = "CorrectionFormat";
  const char cor_libr_key [] = "CorrectionLibrary";
  const char cor_ener_key [] = "CorrectionEnergyMethod";
  const char cor_init_key [] = "CorrectionInitMethod";
  const char cor_data_key [] = "CorrectionInitData";
  const char cor_rpar_key [] = "CorrectionParameterReal";
  const char cor_ipar_key [] = "CorrectionParameterInteger";

  const char* all_keys [] = {
    pot_form_key,
    pot_libr_key,
    pot_ener_key,
    pot_init_key,
    pot_data_key,
    pot_rpar_key,
    pot_ipar_key,
    cor_form_key,
    cor_libr_key,
    cor_ener_key,
    cor_init_key,
    cor_data_key,
    cor_rpar_key,
    cor_ipar_key,
    0
  };
}

int Sjk::pot (const double* r, double& ener) 
{ 
  const char funame [] = "Sjk::pot: ";

  int itemp;

  int status = 0;

  Array_3<double> sjk_r(2, 100, 3);
  if(_pot_format == SJK_FORM || _cor_format == SJK_FORM) {
    for(int f = 0, aa = 0; f < 2; ++f)
      for(int a = 0; a < mol_array[f]->size(); ++a, ++aa) {
	for(int i = 0; i < 3; ++i)
	  sjk_r(f, a, i) = r[aa * 3 + i];
	}
  }

  switch(_pot_format) {
  case NEW_FORM:
    ener = ((ener_t)_pot_ener)(r, _pot_rpar, _pot_ipar, status);
    if(status)
      return status;
    break;

  case SJK_FORM:
    ener = ((sjk_t)_pot_ener)(sjk_r.data(), _pot_rpar, _pot_ipar);
    break;
  }
    
  // correction energy
  if(_cor_ener)
    switch(_cor_format) {
    case NEW_FORM:
      ener += ((ener_t)_cor_ener)(r, _pot_rpar, _pot_ipar, status);
      if(status)
	return status;
      break;

    case SJK_FORM:
      ener += ((sjk_t)_cor_ener)(sjk_r.data(), _pot_rpar, _pot_ipar);
      break;
    }

  return 0;
}

bool Sjk::isinit () { return _init; }

void Sjk::init(const string& fname)
{
  const char funame [] = "Sjk::init: ";

  if(_init) {
    cerr << funame << "sjk environment is allready initialized\n";
    throw Init_Err();
  }
  _init = true;

  int    itemp;
  double dtemp;
  string stemp;

  string token, line;

  ifstream from(fname.c_str());
  if(!from) {
    cerr << funame << "cannot open input file " << fname << "\n";
    throw Open_Err();
  }

  const char* messg;

  void*  pot_lib = 0;
  void*  cor_lib = 0;
  init_t pot_init = 0, cor_init = 0;

  std::string pot_data, cor_data;

  while(from >> token) {
    // potential format
    if(token == pot_form_key) {
      if(!(from >> stemp)) {
	cerr << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }

      if(stemp == "sjk")
	_pot_format = SJK_FORM;
      else if(stemp == "new")
	_pot_format = NEW_FORM;
      else {
	cerr << funame << "unknown format : " << stemp << ": available formats: sjk, new\n";
	throw Form_Err();
      }
    }
    // potential dynamical library
    else if(token == pot_libr_key) {
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      pot_lib = dlopen(stemp.c_str(), RTLD_NOW);
      if(messg = dlerror()) {
	cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential energy method symbol
    else if(token == pot_ener_key) {
      if(!pot_lib) {
	std::cerr << funame << token << ": library should be initialized first\n";
	throw Run_Err();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      _pot_ener = dlsym(pot_lib, stemp.c_str());
      if(messg = dlerror()) {
	std::cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential initialization method symbol
    else if(token == pot_init_key) {
      if(!pot_lib) {
	std::cerr << funame << token << ": library should be initialized first\n";
	throw Run_Err();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      pot_init = (init_t)dlsym(pot_lib, stemp.c_str());
      if(messg = dlerror()) {
	std::cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential initialization data
    else if(token == pot_data_key) {
      if(!(from >> pot_data)) {
	std:: cerr << funame << token << " : corrupted\n";
	throw Form_Err();
      }
    }
    // potential real parameters
    else if(token == pot_rpar_key) {

      if(_pot_rpar){
	cerr << funame << token << ": already defined\n";
	throw Form_Err();
      }

      std::vector<double> vtemp;

      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> dtemp)
	vtemp.push_back(dtemp);

      if(!vtemp.size()) {
	std::cerr << funame << token << ": no data found\n";
	throw Form_Err();
      }

      _pot_rpar = new double[vtemp.size()];
      for(int i = 0; i < vtemp.size(); ++i)
	_pot_rpar[i] = vtemp[i];
    }
    // potential integer parameters
    else if(token == pot_ipar_key) {

      if(_pot_ipar){
	cerr << funame << token << ": already defined\n";
	throw Form_Err();
      }

      std::vector<int> vtemp;

      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> itemp)
	vtemp.push_back(itemp);

      if(!vtemp.size()) {
	std::cerr << funame << token << ": no data found\n";
	throw Form_Err();
      }

      _pot_ipar = new int[vtemp.size()];
      for(int i = 0; i < vtemp.size(); ++i)
	_pot_ipar[i] = vtemp[i];
    }
    // potential format
    else if(token == cor_form_key) {
      if(!(from >> stemp)) {
	cerr << funame << "cannot read " << token << "\n";
	throw Form_Err();
      }

      if(stemp == "sjk")
	_cor_format = SJK_FORM;
      else if(stemp == "new")
	_cor_format = NEW_FORM;
      else {
	cerr << funame << "unknown format : " << stemp << ": available formats: sjk, new\n";
	throw Form_Err();
      }
    }
    // potential dynamical library
    else if(token == cor_libr_key) {
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      cor_lib = dlopen(stemp.c_str(), RTLD_NOW);
      if(messg = dlerror()) {
	cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential energy method symbol
    else if(token == cor_ener_key) {
      if(!cor_lib) {
	std::cerr << funame << token << ": library should be initialized first\n";
	throw Run_Err();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      _cor_ener = dlsym(cor_lib, stemp.c_str());
      if(messg = dlerror()) {
	std::cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential initialization method symbol
    else if(token == cor_init_key) {
      if(!cor_lib) {
	std::cerr << funame << token << ": library should be initialized first\n";
	throw Run_Err();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Form_Err();
      }

      cor_init = (init_t)dlsym(cor_lib, stemp.c_str());
      if(messg = dlerror()) {
	std::cerr << funame << messg << "\n";
	throw Run_Err();
      }
    }
    // potential initialization data
    else if(token == cor_data_key) {
      if(!(from >> cor_data)) {
	std:: cerr << funame << token << " : corrupted\n";
	throw Form_Err();
      }
    }
    // potential real parameters
    else if(token == cor_rpar_key) {

      if(_cor_rpar){
	cerr << funame << token << ": already defined\n";
	throw Form_Err();
      }

      std::vector<double> vtemp;

      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> dtemp)
	vtemp.push_back(dtemp);

      if(!vtemp.size()) {
	std::cerr << funame << token << ": no data found\n";
	throw Form_Err();
      }

      _cor_rpar = new double[vtemp.size()];
      for(int i = 0; i < vtemp.size(); ++i)
	_cor_rpar[i] = vtemp[i];
    }
    // potential integer parameters
    else if(token == cor_ipar_key) {

      if(_cor_ipar){
	cerr << funame << token << ": already defined\n";
	throw Form_Err();
      }

      std::vector<int> vtemp;

      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> itemp)
	vtemp.push_back(itemp);

      if(!vtemp.size()) {
	std::cerr << funame << token << ": no data found\n";
	throw Form_Err();
      }

      _cor_ipar = new int[vtemp.size()];
      for(int i = 0; i < vtemp.size(); ++i)
	_cor_ipar[i] = vtemp[i];
    }
    else {
      cerr << funame << "unknown key: " << token
	   << "\n\navailable keys:\n";
      for(int i = 0; all_keys[i] != 0; ++i)
	cerr << all_keys[i] << "\n";
      throw Form_Err();
    }
  }

  if(pot_init)
    pot_init(pot_data.c_str());

  if(cor_init)
    cor_init(cor_data.c_str());

  if(!_pot_ener) {
    std::cerr << funame << "potential energy method not initialized\n";
    throw Init_Err();
  }
}
