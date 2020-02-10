#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <dlfcn.h>

#include "gauss.hh"
#include "system.hh"
#include "log.hh"
#include "units.hh"

namespace Gauss {

  // debug level
  //
  int debug = Log::DEBUG;

  // gausssian root directory
  //
  std::string exe_dir; 

  // gaussian executable name
  //
  std::string exe_name;

  // failure counter
  //
  int fail_count = 0;

  // geometry pattern
  //
  std::string geom_str = "GEOMETRY_HERE";

  // scratch directory
  //
  std::string scratch_dir;
  
  // potential correction
  //
  extern "C" {
    typedef double  (*pcorr_t) (const int&, const double [], 
				const double [], const int [], double []);
  }
  pcorr_t _pcorr = 0;

  bool  _ispcorr = false;

  // array of real parameters for potential correction
  //
  Array<double> _rpar;

  // array of integer parameters for potential correction
  //
  Array<int>    _ipar;

  // atom specification
  //
  std::vector<Atom> cluster;

  // all methods
  //
  Method::_mpool_t Method::_method_pool;

  // interpolation specification
  //
  Interpolation interpol;

  // minimization specification
  //
  mid_t min_pool;

  // optimization method
  //
  MetId opt_method;

  double low_ener_limit = -1.; // minimal possible relative energy
  
  // test stream
  //
  void test_stream (const std::istream&, const std::string&);

  // convert integer to string
  //
  std::string num2str (int);

  // scan file for a key
  //
  void scan_file (std::ifstream&, const std::string&);
}

double Gauss::pcorr(const int& n, const double r [], double f [])
{
  if(_ispcorr)
    return _pcorr(n, r, _rpar.data(), _ipar.data(), f);
  cout << "Gauss::pcorr: potential is not defined\n";
  return 0.;    
}

bool Gauss::ispcorr () { return _ispcorr; }

void Gauss::test_stream (const std::istream& from, const std::string& funame) 
{
  if(!from) {
    //
    std::cerr << funame << "is corrupted" << std::endl;

    throw Init_Err();
  }
}

std::string Gauss::num2str (int num)
{
  std::ostringstream to;

  to << num;

  return to.str();
}

void Gauss::scan_file (std::ifstream& from, const std::string& key)
{
  std::string token;

  while(from >> token)
    //
    if(key == token)
      //
      return;
}

/*****************************************************
 ******************* Atom class **********************
 *****************************************************/

ostream& Gauss::operator<< (ostream& to, const Atom& a)
{
  to << a.num;
  for (int i = 0; i < 3; ++i)
    to << "   " << a.pos[i] * Phys_const::bohr;
  return to;
}

/*********************************************************************************************
 **************************************** METHOD ID ******************************************
 *********************************************************************************************/

Gauss::Method* Gauss::MetId::operator-> () const 
{
  const char funame [] = "Gauss::MetId::operator->: ";
  
  Method* res = Method::find_method(*this); 

  if(res)
    //
    return res;

  std::cerr << funame << *this << " method is not defined. Available methods:\n";

  Method::print_available_methods(std::cerr);

  throw Init_Err();
}

/*****************************************************************************************
 ************************************** PATTERNS *****************************************
 *****************************************************************************************/

// search pattern base
//
Gauss::Pattern::Pattern (std::istream& from)
{
  const char* funame = "Gauss::Pattern::Pattern: ";

  std::string stemp, separator;

  from >> separator;

  while(from >> stemp) {
    //
    if(stemp == separator)
      //
      break;

    push_back(stemp);
  }

  test_stream(from, funame);
  
  if (!size()) {
    //
    std::cerr << funame << "there should be at least one search key in the pattern" << std::endl;

    throw Init_Err();
  }
}

bool Gauss::Pattern::scan (std::istream& from) const
{
  const char* funame = "Gauss::Pattern::scan: ";

  test_stream(from, funame);
  
  std::string stemp;

  std::list<std::string> test;
  
  while(from >> stemp) {
    //
    if(test.size() < size()) {
      //
      test.push_back(stemp);

      continue;
    }

    test.push_back(stemp);

    test.pop_front();

    if(test == *this)
      //
      return true;
  }

  from.seekg(0);

  from.clear();

  return false;
}

void Gauss::Pattern::print (std::ostream& to) const
{
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    if(it != begin())
      //
      to << " ";
	    
    to << *it;
  }
}

// energy search pattern
//
Gauss::EnergyPattern::EnergyPattern (std::istream& from) : Pattern(from)
{
  const char* funame = "Gauss::EnergyPattern::EnergyPattern: ";
  
  from >> _weight;
  
  test_stream(from, funame);
}

Gauss::MessagePattern::MessagePattern (std::istream& from) : Pattern(from)
{
  const char* funame = "Gauss::MessagePattern::MessagePattern: ";

  std::string stemp, separator;
  
  // read message
  //
  from >> separator;
  
  while(from >> stemp) {
    //
    if (stemp == separator)
      //
      break;
    
    if(!_message.size()) {
      //
      _message = stemp;
    }
    else {
      //
      _message += " " + stemp;
    }//
    //
  }//
  
  test_stream(from, funame);

  if(!_message.size()) {
    //
    std::cerr << funame << "there is no message" << std::endl;

    throw Init_Err();
  }
}

Gauss::FailurePattern::FailurePattern (std::istream& from) : MessagePattern(from)
{
  const char* funame = "Gauss::FailurePattern::FailurePattern: ";

  std::string stemp, separator;
  
  // read failsafe methods
  //
  from >> separator;
  
  while(from >> stemp) {
    //
    if (stemp == separator)
      //
      break;
    
    mid_t::push_back(stemp);
  }

  test_stream(from, funame);

  if(!mid_t::size()) {
    //
    std::cerr << funame << "there should be at least one failsafe method provided" << std::endl;

    throw Init_Err();
  }
}

/*************************************************************************
 ******************************** METHOD *********************************
 *************************************************************************/

void Gauss::Method::init_all ()
{
  for(_mpool_t::const_iterator pit = _method_pool.begin(); pit != _method_pool.end(); ++pit)
    //
    pit->second->init();
}

Gauss::Method* Gauss::Method::find_method (const std::string& mid)
{
  _mpool_t::const_iterator pit = _method_pool.find(mid);

  if(pit != _method_pool.end())
    //
    return pit->second;

  return 0;
}

void Gauss::Method::print_available_methods (std::ostream& to)
{
  for(_mpool_t::const_iterator pit = _method_pool.begin(); pit != _method_pool.end(); ++pit)
    //
    to << pit->first << "\n";
}

void Gauss::Method::_check_for_loops (std::set<std::string>& mstack) const
{
  const char* funame = "Gauss::Method::_check_for_loops: ";

  if (_checked)
    //
    return;  

  std::string stemp;
  int         itemp;

  // check if the current method falls into an infinite loop
  //
  if(!mstack.insert(id()).second) {
    //
    std::cerr << funame << id() << "method is in the infinite loop" << std::endl;

    throw Init_Err();
  }
	
  // check initial guess methods
  //
  for (mid_t::const_iterator mit = guess_methods.begin(); mit != guess_methods.end(); ++mit)
    //
    (*mit)->_check_for_loops(mstack);    

  // check default failsafe methods
  //
  for (mid_t::const_iterator mit = default_failsafe_methods.begin(); mit != default_failsafe_methods.end(); ++mit)
    //
    (*mit)->_check_for_loops(mstack);
   

  // check failure-specific failsafe methods
  //
  for (fpat_t::const_iterator fit = fail_patt.begin(); fit != fail_patt.end(); ++fit)
    //
    for (mid_t::const_iterator mit = fit->mid_t::begin(); mit != fit->mid_t::end(); ++mit)
      //
      (*mit)->_check_for_loops(mstack);

  _checked = true;
  
  mstack.erase(id());
}

// check all methods for infinite loops
//
void Gauss::Method::check_all ()
{
  std::set<std::string> mstack;
  
  for (_mpool_t::const_iterator mit = _method_pool.begin(); mit != _method_pool.end(); ++mit)
    //
    mit->second->_check_for_loops(mstack);
}

Gauss::Method::Method (istream& from)
  : guess_required(false), guess_min(false), 
    _ref_ener(0.), _abs_ener(0.), _checked(false), _status(INIT), _failsafe(0) 
{
  std::string funame = "Gauss::Method::Method: ";

  static const std::string init_guess_key = "InitialGuessMethods";
  static const std::string  guess_req_key = "InitialGuessRequired";
  static const std::string  guess_min_key = "InitialGuessMinimum";
  static const std::string  ener_patt_key = "EnergySearchPattern";
  static const std::string  fail_patt_key = "FailSearchPattern";
  static const std::string   run_patt_key = "RerunPattern";
  static const std::string    default_key = "DefaultFailSafeMethods";
  static const std::string      refen_key = "ReferenceEnergy";
  static const std::string        end_key = "EndMethod";              

  static const char* all_keys [] = {
    init_guess_key.c_str(),
    guess_req_key.c_str(),  
    guess_min_key.c_str(),  
    ener_patt_key.c_str(), 
    fail_patt_key.c_str(), 
    run_patt_key.c_str(),  
    default_key.c_str(),   
    refen_key.c_str(),     
    end_key.c_str(),    
    0
  };

  std::string stemp;
  
  int         itemp;

  std::string separator, comment, token, line, prefix;

  // method name
  //
  from >> _id;

  std::getline(from, comment);

  // check for duplicated methods
  //
  if(find_method(id())) {
    //
    std::cerr << funame << "duplicated method id: " << _id << std::endl;

    throw Init_Err();
  }

  // redefine function name
  //
  funame += id() + " method: ";
  
  // input cycle
  //
  while (from >> token) {
    //
    prefix = funame + token + ": ";
    
    // end key
    //
    if (token == end_key) {
      //
      break;
    }
    // initial guess methods
    //
    else if (token == init_guess_key) {
      //
      if (guess_methods.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;
	
	throw Init_Err();
      }
      
      from >> separator;
      
      while (from >> stemp) {
	//
	if (stemp == separator)
	  //
	  break;
	
	guess_methods.push_back(stemp);
      }

      test_stream(from, prefix);
    }
    // is initial guess required
    //
    else if (token == guess_req_key) {
      //
      from >> guess_required;

      test_stream(from, prefix);
    }
    // use minimal initial guess
    //
    else if (token == guess_min_key) {
      //
      from >> guess_min;
      
      test_stream(from, prefix);
    }
    // energy search pattern
    //
    else if (token == ener_patt_key) {
      //
      EnergyPattern ptemp(from);
      
      ener_patt.push_back(ptemp);
    }
    // method failure pattern and associated failsafe methods
    //
    else if (token == fail_patt_key) {
      //
      FailurePattern ptemp(from);
      
      fail_patt.push_back(ptemp);
    }
    // rerun required pattern
    //
    else if(token == run_patt_key) {
      
      MessagePattern ptemp(from);
      
      run_patt.push_back(ptemp);
    }
    // default failsafe method
    //
    else if(token == default_key) {
      //
      if (default_failsafe_methods.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;
	
	throw Init_Err();
      }
      
      from >> separator;
      
      while(from >> stemp) {
	//
	if (stemp == separator)
	  //
	  break;
	
	default_failsafe_methods.push_back(stemp);
      }

      test_stream(from, prefix);
    }
    // reference energy
    //
    else if(token == refen_key) {
      //
      from >> _ref_ener;

      test_stream(from, prefix);
    }
    // unknown keyword
    //
    else {
      //
      std::cerr << funame << "unknown keyword: " << token << ": available keywords:\n";
      
      for(int i = 0; all_keys[i] != 0; ++i)
	//
	std::cerr << all_keys[i] << "\n";

      throw Init_Err();
    }
  }// input cycle

  test_stream(from, funame);

  if(!guess_methods.size() && guess_required) {
    //
    std::cerr << funame << "there should be at least one initial guess method" << std::endl;

    throw Init_Err();
  }

  if(!ener_patt.size()) {
    //
    std::cerr << funame << "there should be at least one energy pattern" << std::endl;

    throw Init_Err();
  }

  // set template string
  //
  std::string tml_file = id() +".tml";
  
  std::ifstream tin(tml_file.c_str());
  
  if(!tin) {
    std::cerr << funame << "cannot open template file " << tml_file << std::endl;
    
    throw Init_Err();
  }

  while(getline(tin, stemp))
    //
    _template += stemp + "\n";

  _geom_pos = _template.find(geom_str);
  
  if(_geom_pos < 0 || _geom_pos >= _template.size()) {
    //
    std::cerr << funame << "cannot find geometry pattern " << geom_str
	      << " in template file " << tml_file << std::endl;
    
    throw Init_Err();
  }

  _template.erase(_geom_pos, geom_str.size());
  
  //add the current method to the methods pool
  //
  _method_pool[id()] = this;
}

void Gauss::Method::backup()
{
  static const char funame[] = "Gauss::Method::backup: ";

  _check_status();
  
  file_copy((id() + ".chk").c_str(), (id() + ".chk.back").c_str());

  file_copy((id() + ".log").c_str(), (id() + ".log.back").c_str());
  
  file_copy((id() + ".fchk").c_str(),(id() + ".fchk.back").c_str());
  
  _save_ener = _abs_ener;
}

void Gauss::Method::restore()
{
  static const char funame[] = "Gauss::Method::restore: ";

  if(rename((id() + ".chk.back").c_str(), (id() + ".chk").c_str()) ||
     rename((id() + ".log.back").c_str(), (id() + ".log").c_str()) ||
     rename((id() + ".fchk.back").c_str(),(id() + ".fchk").c_str())) {
    //
    std::cerr << funame << "failed" << std::endl;
				       
    throw Init_Err();
  }
  
  _abs_ener = _save_ener;
  
  _status = SUCCESS;
}

void Gauss::Method::save_log (int count) const
{
  // gaussian file extensions
  //
  const char* fext [3] = {".com", ".log", ".chk"};

  const std::string suffix = "." + num2str(count);

  std::string from, to;
  
  for(int i = 0; i < 3; ++i) {
    //
    from = id() + fext[i];

    to = from + suffix;

    rename(from.c_str(), to.c_str());
  }

  // rename standard output and standard error files
  //
  const char* extra [2] = {"std.out", "std.err"};

  for(int i = 0; i < 2; ++i) {
    //
    to = extra[i] + suffix;

    rename(extra[i], to.c_str());
  }
}

void Gauss::Method::_check_status () const
{
  const char* funame = "Gauss::Method::_check_status: ";
  
  if (status() != SUCCESS) {
    //
    std::cerr << funame << id() << " method: wrong status: " << status() << std::endl;
											 
    throw Logic_Err();
  }
}

void Gauss::Method::read_forces (std::vector<D3>& force_array) const
{
  static const char* funame = "Gauss::Method::read_forces: ";

  _check_status();

  double      dtemp;
  std::string stemp;

  std::string comment, token, line;

  if (force_array.size() != cluster.size()) {
    //
    std::cerr << funame << "numbers of forces, " << force_array.size()
	      << ", and atoms, " << cluster.size() << " mismatch" << std::endl;

    throw Logic_Err();
  }

  // formatted checkpoint file
  //
  const std::string file = id() + ".fchk";

  std::ifstream from(file.c_str());

  if(!from) {
    //
    std::cerr << funame <<  "cannot open " << file << std::endl;

    throw Open_Err();
  }

  scan_file(from, "Gradient");

  if(!from) {
    //
    std::cerr << funame << file << ": energy gradient is not found: check files *." 
	      << fail_count << std::endl;

    Log::out << funame << file << ": energy gradient is not found: check files *." 
	      << fail_count << std::endl;

    save_log(fail_count++);

    throw Find_Err();
  }

  std::getline(from, line);

  for (int at = 0; at < cluster.size(); ++at) {
    //
    for (int i = 0; i < 3; ++i) {
      //
      if(!(from >> dtemp)) {
	//
	Log::out << funame << file <<  ": reading energy gradient failed: check files *."
		 << fail_count << std::endl;

	std::cerr << funame << file <<  ": reading energy gradient failed: check files *."
		  << fail_count << std::endl;

	save_log(fail_count++);

	throw Form_Err();
      }
      
      force_array[at][i] = -dtemp;
    }//
    //
  }//
  //
}// read forces

Gauss::Method* Gauss::Method::apply (int flags)
{

  static const char* funame = "Gauss::Method::apply: ";

  std::string stemp;
  double      dtemp;

  // check if the method was allready applied
  //
  if(status() == SUCCESS)
    //
    return this;

  if(status() == FAIL) 
    //
    if(flags & DIRECT) {
      //
      return 0;
    }
    else
      //
      return _failsafe;

  _status = FAIL;
  
  _failsafe = 0;

  // save checkpoint file
  //
  if(flags & READ) {
    //
    file_copy((id() + ".chk").c_str(), (id() + ".chk.read").c_str());
  }
  // initial guess calculation
  //
  else {
    //
    Method* guess = 0;
    
    for (mid_t::const_iterator mit = guess_methods.begin(); mit != guess_methods.end(); ++mit) {
      //
      Method* mp = (*mit)->apply();
	
      if(guess_min) {
	//
	if(mp && (!guess || mp->abs_energy() < guess->abs_energy()))
	  //
	  guess = mp;
      }
      else if(mp) {
	//
	guess = mp;

	break;
      }
    }

    // no initial guess
    //
    if(!guess && guess_required) {
      //
      std::cerr << funame << "initial guess calculation for "
		<< id() << " method failed" << std::endl;

      if(debug >= Log::WARN)
	//
	Log::out << funame << "Oops!!! initial guess calculation for "
		 << id() << " method failed" << std::endl;

      // no failsafe methods are allowed
      //
      if(flags & DIRECT)
	//
	return 0;
	
      // default failsafe methods
      //
      for(mid_t::const_iterator mit = default_failsafe_methods.begin(); mit != default_failsafe_methods.end(); ++mit) {
	//
	if (debug >= Log::WARN)
	  //
	  Log::out << funame <<"trying default failsafe method " << (*mit)->id() << std::endl;
	  
	if(_failsafe = (*mit)->apply(flags))
	  //
	  return _failsafe;
      }
	
      return 0;
    }

    if(guess && file_copy((guess->id()+".chk").c_str(),(id()+".chk").c_str())) {
      //
      std::cerr << funame << "copy of initial guess " << guess->id()
		<< " method checkpoint file failed" << std::endl;

      throw File_Err();
    }//
    //
  }//initial guess

  // molecular geometry
  //
  std::ostringstream geom;

  for (int i = 0; i < cluster.size(); ++i) {
    //
    if(i != 0)
      //
      geom << "\n";
    
    geom << cluster[i];
  }
    
  std::string input = _template;

  input.insert(_geom_pos, geom.str());

  stemp = id() + ".com";
  
  std::ofstream gout(stemp.c_str());

  if (!gout) {
    //
    std::cerr << funame << "cannot open  " << stemp<< " file for writing" << std::endl;
    
    throw Open_Err();
  }

  gout << input;
  
  gout.close();

  // gaussian call
  //
  int exit_code;
 
  // gaussian call failure
  //
  if(exit_code = call_exe(exe_name.c_str(), id().c_str(), (char*) 0)) {
    //
    if(debug >= Log::ERROR)
      //
      Log::out << funame << "Oops!!! method " << id() << " failed with exit code " << exit_code << "\n"
	       << "Check files *." << fail_count << "\n\n";
    
    std::cerr << funame << "Oops!!! method " << id() << " failed with exit code " << exit_code << "\n"
	      << "Check files *." << fail_count << "\n\n";

    save_log(fail_count++);

    std::ifstream from((id() + ".log").c_str());
    
    // rerun pattern search
    //
    for (rpat_t::const_iterator rit = run_patt.begin(); rit != run_patt.end(); ++rit) {
      //
      if(rit->scan(from)) {
	//
	from.close();

	std::cerr << funame << rit->message() << std::endl;

	if (debug >= Log::ERROR)
	  //
	  Log::out << funame << rit->message() << std::endl;
	    
	_status = INIT;
	    
	if(flags & READ)
	  //
	  rename((id() + ".chk.read").c_str(), (id() + ".chk").c_str());
	    
	return apply(flags);
      }
    }

    if(flags & DIRECT)
      //
      return 0;
    
    // specific failsafe methods
    //
    for(fpat_t::const_iterator fit = fail_patt.begin(); fit != fail_patt.end(); ++fit) {
      //
      if(fit->scan(from)) {
	//
	from.close();

	std::cerr << funame << fit->message() << std::endl;
	    
	if(debug >= Log::ERROR)
	  //
	  Log::out << funame << fit->message() << std::endl;
	    
	for(mid_t::const_iterator mit = fit->mid_t::begin(); mit != fit->mid_t::end(); ++mit) {
	  //
	  if(debug >= Log::INFO)
	    //
	    Log::out << funame <<"trying failsafe method " << (*mit)->id() << std::endl;
	      
	  if(_failsafe = (*mit)->apply(flags))
	    //
	    return _failsafe;
	}

	return 0;
      }
    }

    // default failsafe methods 
    //
    for(mid_t::const_iterator mit = default_failsafe_methods.begin(); mit != default_failsafe_methods.end(); ++mit) {

      if (debug >= Log::INFO)
	//
	Log::out << funame <<"trying default failsafe method " << (*mit)->id() << std::endl;

      if(_failsafe = (*mit)->apply(flags))
	//
	return _failsafe;
    }

    return 0;
    //
  }// gaussian failed
  //
  // gaussian succeeded
  //
  else {
    //
    // convert checkpoint file
    //
    if (call_exe("formchk", (id() + ".chk").c_str(), (id() + ".fchk").c_str(), (char*) 0)) {
      //
      std::cerr << funame << "checkpoint file conversion failure" << std::endl;

      return 0;
    }

    // find energies in checkpoint file
    //
    std::ifstream from((id() + ".fchk").c_str());

    _abs_ener = 0.;

    for (epat_t::const_iterator pit = ener_patt.begin(); pit != ener_patt.end(); ++pit) {
      //
      if(pit->scan(from)) {
	//
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << "cannot read energy value for <<" << *pit << ">> energy search pattern\n"
		    << funame << "Check files *." << fail_count << std::endl;

	  save_log(fail_count++);

	  throw Form_Err();
	}

	_abs_ener += dtemp * pit->weight();

	// move to the beginning of the file
	//
	from.seekg(0);
      }
      // pattern not found
      //
      else {
	//
	std::cerr << funame << "<<" << *pit << ">> energy search pattern not found\n"
		  << funame << "Check files *." << fail_count << std::endl;

	save_log(fail_count++);

	throw Form_Err();
      }
    }
 
    _status = SUCCESS;
      
    return this;
    //
  }// g98 succeeded
  //
}// apply()

void Gauss::Method::print_all_info (std::ostream& to)
{
  int old_precision = to.precision(9);
  
  for (_mpool_t::const_iterator mit = _method_pool.begin(); mit != _method_pool.end(); ++mit)
    //
    mit->second->print_info(to);
  
  to.precision(old_precision);
  
  to << std::endl;
}

void Gauss::Method::print_info (std::ostream& to) const
{
  const char funame [] = "Gauss::Method::print_info: ";
  
  to << std::setw(10) << id() << " method: ";
  
  switch(status()) {
  case INIT:
    //
    to << "not initialized yet";
    
    break;
  case SUCCESS:
    //
    to << "Energy[a.u] = " << abs_energy();

    break;
  case FAIL:
    //
    to << "failed";

    break;
  default:
    //
    std::cerr << funame << "wrong case: " << status() << std::endl;

    throw Logic_Err();
  }

  to << std::endl;
}

/************************************************************
 ***************** Interpolation class **********************
 ************************************************************/

void Gauss::Interpolation::init (ifstream& from)
{
  static const char* funame = "Gauss::Interpolation::init: ";

  std::string stemp;
  double      dtemp;

  std::string separator;

  from >> separator;
  
  while (from >> stemp) {
    //
    if (stemp == separator)
      //
      break;
    
    if(!(from >> dtemp)) {
      //
      std::cerr << funame << "there should be a weight factor after each method name" << std::endl;

      throw Init_Err();
    }

    if(dtemp <= 0.)
      //
      std::cerr << funame << "WARNING: " << stemp << " method: negative weight: " << dtemp << std::endl;

    push_back(std::make_pair(MetId(stemp), dtemp));
  }

  test_stream(from, funame);

  if (!size()) {
    //
    std::cerr << funame << "there should be at least one method to apply" << std::endl;

    throw Init_Err();
  }

  // normalize pool weights
  //
  dtemp = 0.;

  for(const_iterator pit = begin(); pit != end(); ++pit)
    //
    dtemp += pit->second;

  if(dtemp <= 0.) {
    //
    std::cerr << funame << "cumulative weight should be positive: " << std::endl;

    throw Init_Err();
  }

  for(iterator pit = begin(); pit != end(); ++pit)
    //
    pit->second /= dtemp;
}	

void Gauss::Interpolation::check () const
{
  static const char* funame = "Gauss::Interpolation::check: ";

  for (const_iterator pit = begin(); pit != end(); ++pit)
    //
    if(!Method::find_method(pit->first)) {
      //
      std::cerr << funame << pit->first << " interpolation method is not defined" << std::endl;

      throw Init_Err();
    }
}

double Gauss::Interpolation::apply(int flags) const
{
  static const char* funame = "Gauss::Interpolation::apply: ";

  if (flags & FORCE)
    //
    for (int at = 0; at < cluster.size(); ++at)
      //
      for (int j = 0; j < 3; ++j)
	//
	cluster[at].force[j] = 0.;

  double res = 0.;
  //
  for(const_iterator pit = begin(); pit != end(); ++pit) {
    //
    Method* mp = pit->first->apply(flags);

    if(!mp) {
      //
      std::cerr << funame << pit->first << " method failed" << std::endl;

      Log::out << funame << pit->first << " method failed" << std::endl;
      
      throw Run_Err();
    }
    
    res += mp->rel_energy() * pit->second;
    
    if (flags & FORCE) {
      //
      std::vector<D3> force_array(cluster.size());
  
      mp->read_forces(force_array);
      
      for(int at = 0; at < cluster.size(); ++at)
	//
	for (int j = 0; j < 3; ++j)
	  //
	  cluster[at].force[j] += force_array[at][j] * pit->second;
    }
  }
  
  return res;
}

// main initialization
//
void Gauss::init (const std::string& fname)
{
  static const char* funame = "Gauss::init: ";

  static const std::string    end_key = "EndInit";
  static const std::string  inter_key = "Interpolation";
  static const std::string method_key = "Method";
  static const std::string    min_key = "MinMethods";
  static const std::string  debug_key = "DebugLevel";
  static const std::string   exec_key = "ExecDir";
  static const std::string    g98_key = "ExecName";
  static const std::string   csym_key = "PotentialCorrectionSymbol";
  static const std::string   clib_key = "PotentialCorrectionLibrary";
  static const std::string   rpar_key = "PotentialParameterReal";
  static const std::string   ipar_key = "PotentialParameterInteger";
  static const std::string    opt_key = "OptimizationMethod";
  static const std::string   read_key = "OptimizationReadMethod";
  static const std::string   geom_key = "GeometryPattern";
  static const std::string    scr_key = "ScratchDirectory";

  static const char* all_keys [] = {
    scr_key.c_str(),
    geom_key.c_str(),
    end_key.c_str(),
    inter_key.c_str(),
    method_key.c_str(),
    min_key.c_str(),
    debug_key.c_str(),
    exec_key.c_str(),
    g98_key.c_str(),
    csym_key.c_str(),
    clib_key.c_str(),
    rpar_key.c_str(),
    ipar_key.c_str(),
    opt_key.c_str(),
    read_key.c_str(),
    0
  };

  double      dtemp;
  int         itemp;
  std::string stemp;
  
  std::string token, separator, comment, line, prefix;

  std::ifstream from(fname.c_str());
  
  if(!from) {
    //
    std::cerr << funame << "cannot open " << fname << " input file" << std::endl;
    
    throw Open_Err();
  }

  // potential correction library and symbol
  //
  std::string clib, csym;

  // input cycle
  //
  while (from >> token) {
    //
    prefix = funame + token + ": ";
    
    // end key
    //
    if (token == end_key) {
      //
      break;
    }
    // template geometry keyword
    //
    else if (token == geom_key) {
      //
      from >> geom_str;

      getline(from, comment);

      test_stream(from, prefix);
    }
    // scratch directory
    //
    else if (token == scr_key) {
      //
      from >> scratch_dir;
      
      getline(from, comment);

      test_stream(from, prefix);
    }
    // debug level
    //
    else if (token == debug_key) {
      //
      from >> debug;
      
      test_stream(from, prefix);
    }
    // gaussian directory
    //
    else if(token == exec_key) {
      //
      if(exe_dir.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      if(!(from >> exe_dir)) {
	//
	std::cerr << prefix << "corrupted" << std::endl;

	throw Init_Err();
      }
    }
    // gaussian executable name
    //
    else if(token == g98_key) {
      //
      if(exe_name.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      if(!(from >> exe_name)) {
	//
	std::cerr << prefix << "corrupted" << std::endl;

	throw Init_Err();
      }
    }
    // interpolation initialization
    //
    else if (token == inter_key) {
      //
      if(interpol.isinit()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      interpol.init(from);
    }
    // new method initialization
    //
    else if (token == method_key) {
      //
      new Method(from);
    }
    // minimization initialization
    //
    else if (token == min_key) {
      //
      from >> separator;

      while (from >> token) {
	//
	if (token == separator)
	  //
	  break;
	
	min_pool.push_back(token);
      }

      test_stream(from, prefix);
    }
    // potential correction library
    //
    else if(token == clib_key) {
      //
      if(clib.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      if(!(from >> clib)) {
	//
	std::cerr << prefix << "corrupted" << std::endl;
	
	throw Init_Err();
      }

      std::getline(from, comment);
    }
    // potential correction method
    //
    else if(token == csym_key) {
      //
      if(csym.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      if(!(from >> csym)) {
	//
	std::cerr << prefix << "corrupted" << std::endl;
	
	throw Init_Err();
      }

      std::getline(from, comment);
    }
    // real parameters
    //
    else if(token == rpar_key) {
      //
      if(_rpar.size()){
	//
	std::cerr << prefix << "already defined" << std::endl;
	
	throw Init_Err();
      }

      std::getline(from, line);

      std::istringstream iss(line);
      
      std::list<double> data;

      while (iss >> dtemp) {
	//
	data.push_back(dtemp);
      }

      if(!data.size()) {
	//
	std::cerr << prefix << "no readable data" << std::endl;

	throw Init_Err();
      }
      
      _rpar.resize(data.size());

      itemp = 0;
      
      for(std::list<double>::const_iterator dit = data.begin(); dit != data.end(); ++dit, ++itemp)
	//
	_rpar[itemp] = *dit;
    }
    // integer parameters
    //
    else if(token == ipar_key) {
      //
      if(_ipar.size()){
	//
	std::cerr << prefix << "already defined" << std::endl;

	throw Init_Err();
      }

      std::getline(from, line);

      std::istringstream iss(line);
      
      std::list<int> data;

      while (iss >> itemp) {
	//
	data.push_back(itemp);
      }

      if(!data.size()) {
	//
	std::cerr << prefix << "no readable data" << std::endl;

	throw Init_Err();
      }
      
      _ipar.resize(data.size());
      
      itemp = 0;
      
      for(std::list<int>::const_iterator dit = data.begin(); dit != data.end(); ++dit, ++itemp)
	//
	_ipar[itemp] = *dit;
    }
    // optimization method
    //
    else if(token == opt_key) {
      //
      if(opt_method.size()) {
	//
	std::cerr << prefix << "already initialized" << std::endl;

	throw Init_Err();
      }
      
      if(!(from >> opt_method)) {
	//
	std::cerr << prefix << "corrupted" << std::endl;

	throw Init_Err();
      }
    }
    // unknown keyword
    //
    else {
      std::cerr << funame << "unknown key: " << token
		<< "\navailable keys:\n";
      
      for(int i = 0; all_keys[i] != 0; ++i)
	//
	std::cerr << all_keys[i] << "\n";

      throw Init_Err();
    }
  }

  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted" << std::endl;

    throw Init_Err();
  }

  if(!exe_dir.size()) {
    //
    std::cerr << funame << "gaussian exec directory is not defined" << std::endl;

    throw Init_Err();
  }
  
  if(!exe_name.size()) {
    //
    std::cerr << funame << "gaussian executable name is not defined" << std::endl;

    throw Init_Err();
  }

  if (!Method::pool_size()) {
    //
    std::cerr << funame << "there should be at least one method" << std::endl;
    
    throw Init_Err();
  }

  if (!min_pool.size() && !interpol.isinit()) {
    //
    std::cerr << funame << "there should be either minimization or "
      "interpolation specifications available" << std::endl;

    throw Init_Err();
  }

  // check the methods tree for infinite loops
  //
  Method::check_all();

  // check interpolation specification
  //
  if(interpol.isinit())
    //
    interpol.check();

  // check minimization methods
  //
  for (mid_t::const_iterator mit = min_pool.begin(); mit != min_pool.end(); ++mit)
    //
    if(!Method::find_method(*mit)) {
      //
      std::cerr << funame << "minimization pool method " << *mit << " is not defined" << std::endl;

      throw Init_Err();
    }

  // check optimization method
  //
  if(opt_method.size() && !Method::find_method(opt_method)) {
    //
    std::cerr << funame << opt_method << " optimization method is not defined" << std::endl;

    throw Init_Err();
  }

  // gaussian directory check
  //
  if (exe_dir[0] != '/') {
    //
    std::cerr << funame << "the gaussian directory " << exe_dir
	      << " should be an absolute path" << std::endl;

    throw Init_Err();
  }

  // scratch directory check
  //
  if(!scratch_dir.size()) {
    //
    std::cerr << funame << "scratch directory is not defined" << std::endl;

    throw Init_Err();
  }
  
  // set environment
  //
  setenv("GAUSS_SCRDIR", scratch_dir.c_str(), 1);
  
  setenv("GAUSS_EXEDIR", exe_dir.c_str(), 1);
  
  setenv("PATH", exe_dir.c_str(), 1);

  // set potential correction
  //
  if(clib.size() && csym.size()) {
    //
    _ispcorr = true;
    
    void* handle = dlopen(clib.c_str(), RTLD_NOW);
    
    const char* messg = dlerror();
    //
    if(messg) {
      //
      std::cerr << funame << messg << std::endl;
      
      throw Init_Err();
    }
    
    _pcorr = (pcorr_t)dlsym(handle, csym.c_str());
    
    messg = dlerror();
    //
    if(messg) {
      //
      std::cerr << funame << messg << std::endl;
      
      throw Init_Err();
    }
  }

  if(_ipar.size()) {
    //
    std::ostringstream to;
    
    to << funame << "correction potential integer parameters:";
    //
    for(int i = 0; i < _ipar.size(); ++i)
      //
      to << "   " << _ipar[i];

    to << "\n";

    std::cout << to.str();
  }
    
  if(_rpar.size()) {
    //
    std::ostringstream to;
    
    to << funame << "correction potential real parameters:";
    //
    for(int i = 0; i < _rpar.size(); ++i)
      //
      to << "   " << _rpar[i];

    to << "\n";

    std::cout << to.str();
  }
}

void Gauss::print_geom ()
{
  Log::out << "Geometry(Angstrom):\n\n";

  for (int i = 0; i < cluster.size(); ++i)
    //
    Log::out << cluster[i];

  Log::out << "\n";
}

Gauss::Method* Gauss::min_pot (int flags)
{
  static const char* funame = "Gauss::min_pot: ";
  
  Method::init_all();

  if(debug >= Log::INFO)
    //
    print_geom();

  Method* res = 0;
  
  for (mid_t::const_iterator mit = min_pool.begin(); mit != min_pool.end(); ++mit) {
    //
    Method* mp = (*mit)->apply(flags);
    
    if(mp && (!res || mp->abs_energy() < res->abs_energy()))
      //
      res = mp;
  }

  if (debug >= Log::INFO)
    //
    Method::print_all_info(Log::out);

  if (!res) {
    //
    std::cerr << funame << "energy calculation failed" << std::endl;
   
    if(debug >= Log::ERROR)
      //
      Log::out << funame << "energy calculation failed" << std::endl;

    throw Run_Err();
  }
  
  if(res->rel_energy() < low_ener_limit) {
    //
    std::cerr << funame << "energy[kcal/mol], " << res->rel_energy() / Phys_const::kcal
	      << ", is too negative, check *." << fail_count << " files" << std::endl;
    
    if(debug >= Log::ERROR) {
      //
      Log::out << funame << "energy[kcal/mol], " << res->rel_energy() / Phys_const::kcal
	       << ", is too negative, check *." << fail_count << " files" << std::endl;
    }

      
    res->save_log(fail_count++);
    
    throw Range_Err();
  }

  if (flags & FORCE) {
    //
    std::vector<D3> force_array(cluster.size());
  
    res->read_forces(force_array);
      
    for(int at = 0; at < cluster.size(); ++at)
      //
      for (int j = 0; j < 3; ++j)
	//
	cluster[at].force[j] = force_array[at][j];
  }

  return res;
}

double Gauss::inter_pot (int flags)
{
  Method::init_all();

  return interpol.apply(flags);
}
