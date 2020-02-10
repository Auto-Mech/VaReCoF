#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <iomanip>
#include <sstream>

#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <cerrno>
#include <sys/stat.h>
#include <cstdio>
#include <cstring>
#include <csignal>
#include <csetjmp>
#include <cfenv>

#include <mpi.h>

#include "rotd.hh"    
#include "tmatrix.hh"
#include "raninit.hh"
#include "divsur.hh"
#include "integral.hh"
#include "error.hh"
#include "math.hh"
#include "units.hh"
#include "input.hh"

#include "comm.hh"
#include "random.hh"
#include "force.hh"
#include "ref.hh"
#include "flux.hh"
#include "gauss.hh"
#include "multipole.hh"
#include "molpro.hh"
#include "sjk.hh"
#include "log.hh"

void test_except (const std::string& prefix)
{
  const char funame [] = "floating-point exception: ";
    
  int res = fetestexcept(FE_ALL_EXCEPT);

  int stat = 1;

  std::ostringstream to;
  
  if(res & FE_INEXACT) {
    //
    to << prefix << funame << "inexact\n";

    stat = 0;
  }

  if(res & FE_INVALID) {
    //
    to << prefix << funame << "invalid\n";

    stat = 0;
  }

  if(res & FE_DIVBYZERO) {
    //
    to << prefix << funame << "divide by zero\n";

    stat = 0;
  }
  
  if(res & FE_OVERFLOW) {
    //
    to << prefix << funame << "overflow\n";

    stat = 0;
  }

  if(res & FE_UNDERFLOW) {
    //
    to << prefix << funame << "underflow\n";

    stat = 0;
  }

  if(stat) {
    //
    to << prefix << funame << "OK\n";
  }

  std::cerr << to.str();
  
  // clear all exceptions
  //
  feclearexcept(FE_ALL_EXCEPT);
}
 
// thermal flux array;
std::map<int, Ref<ThermalMultiFluxArray> > thermal_flux;

void print_thermal_flux (const string& fname)
{
  const double conv_fac = 612.6;
  const int out_size = 6;

  static int tm_size = 0;
  static bool is_init = false;
  if(!is_init) {
    is_init = true;
    for(; tm_size < FluxBase::tm_size(); ++tm_size)
      if(FluxBase::tm_grid(tm_size) / Phys_const::kelv > 5000.)
	break;
  }
  
  double dtemp;
  int itemp;
  int count;

  std::ofstream to(fname.c_str());
  std::map<int, Ref<ThermalMultiFluxArray> >::iterator it, begin;

  int max_face = 0;
  for(it = thermal_flux.begin(); it != thermal_flux.end(); ++it)
    if(max_face < it->second.value().size())
      max_face = it->second.value().size();
	
  for(int face = 0; face < max_face; ++face) {// face cycle
    to << "FACE = " << face << "\n\n";

    for(begin = thermal_flux.begin(); begin != thermal_flux.end(); 
	begin = it) {// flux output cycle
      
      to << std::setw(5) << "Surf:";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->first;
	}
      to << "\n\n";

      to << "Minimal energy (kcal/mol):\n";
      for(int pes = 0; pes < PES::size(); ++pes) {
	to << std::setw(5) << " ";
	for(it = begin, count = 0; 
	    it != thermal_flux.end() && count != out_size; ++it)
	  if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	    ++count;
	    to << std::setw(14) 
	       << it->second.value()[face].min_energy(pes)/Phys_const::kcal;
	  }
	to << "\n";
      }
      to << "\nTemperature (K) - Flux(10^11 cm^3/sec)\n";
      for(int tm = 0; tm < tm_size; ++tm) {// temperature cycle
	for(int pes = 0; pes < PES::size(); ++pes) {
	  if(!pes) 
	    to << std::setprecision(4) << std::setw(5) 
	       << FluxBase::tm_grid(tm) / Phys_const::kelv;
	  else
	    to << std::setw(5) << " ";
	  for(it = begin, count = 0; 
	      it != thermal_flux.end() && count != out_size; ++it)
	    if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	      ++count;
	      // flux
	      if(it->second.value()[face].cn_flux(tm, pes) > 1.e99)
		to << " " << std::setw(8) << "****";
	      else
		to << " " << std::setprecision(3) << std::setw(8) 
		   << it->second.value()[face].cn_flux(tm, pes) * conv_fac;
	      // error
	      to << "(";
	      if(it->second.value()[face].cn_dflx(tm, pes) < 1.)
		to << "< 1";
	      else if(it->second.value()[face].cn_dflx(tm, pes) > 99.)
		to << "***";
	      else
		to << std::setprecision(2) << std::setw(3) 
		   << it->second.value()[face].cn_dflx(tm, pes);
	      to << ")";
	    }
	  to << "\n";
	}
      }// temperature cycle
      to << "\n";

      to << std::setw(5) << "Surf:";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->first;
	}
      to << "\n\n";

      to << "Successful samplings:\n";
      to << std::setw(5) << " ";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->second.value()[face].acct_smp();
	}
      to << "\n";

      to << "Failed samplings:\n";
      to << std::setw(5) << " ";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->second.value()[face].fail_smp();
	}
      to << "\n";

      to << "Dummy samplings:\n";
      to << std::setw(5) << " ";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && 
	   it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->second.value()[face].fake_smp();
	}
      to << "\n";

      to << "Out-of-face samplings:\n";
      to << std::setw(5) << " ";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->second.value()[face].face_smp();
	}
      to << "\n";

      to << "Close-atoms samplings:\n";
      to << std::setw(5) << " ";
      for(it = begin, count = 0; 
	  it != thermal_flux.end() && count != out_size; ++it)
	if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	  ++count;
	  to << std::setw(14) << it->second.value()[face].dist_smp();
	}
      to << "\n\n";

      // minimal thermal flux
      to << "Minimal thermal flux (10^11 cm^3/sec):\n";
      to.precision(4);
      for(int pes = 0; pes < PES::size(); ++pes) {// PES cycle
	vector<int> min_tm;
	for(it = begin, count = 0; 
	    it != thermal_flux.end() && count != out_size; ++it)
	  if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	    ++count;
	    for(int tm = 0; tm < FluxBase::tm_size(); ++tm)
	      if(!tm || it->second.value()[face].cn_flux(tm, pes) < dtemp) {
		dtemp = it->second.value()[face].cn_flux(tm, pes);
		itemp = tm;
	      }
	    min_tm.push_back(itemp);
	  }

	to << std::setw(5) << "T, K:";
	for(it = begin, count = 0; 
	    it != thermal_flux.end() && count != out_size; ++it)
	  if(face < it->second.value().size() && it->second.value()[face].acct_smp())
	    to << std::setw(14) 
	       << FluxBase::tm_grid(min_tm[count++])/ Phys_const::kelv;
	to << "\n";

	to << std::setw(5) << "Flux:";
	for(it = begin, count = 0; 
	    it != thermal_flux.end() && count != out_size; ++it)
	  if(face < it->second.value().size() && it->second.value()[face].acct_smp())
	    to << std::setw(14) 
	       << it->second.value()[face].cn_flux(min_tm[count++], pes) * conv_fac;
	to << "\n";
      }// PES cycle
      to << "\n\n*******************************************************\n\n";
    }// flux output cycle

    to.precision(5);
    for(it = thermal_flux.begin(); it != thermal_flux.end(); ++it)
      if(face < it->second.value().size() && it->second.value()[face].acct_smp()) {
	to << "Surface = " << it->first << "\n\n";
	for(int pes = 0; pes < PES::size(); ++pes) {
	  it->second.value()[face].min_dynvar(pes).attach();
	  update_atoms();
	  to.setf(ios::fixed, ios::floatfield);
	  to << "Minimal energy configuration (Angstrom):\n\n";
	  for(int at = 0; at < atoms.size(); ++at)
	    if(atoms[at].type != Dummy) { 
	      to << std::setw(2) << atoms[at].name();
	      for(int i = 0; i < 3; ++i)
		to << std::setw(10) << atoms[at].lf_pos[i] * Phys_const::bohr;
	      to << "\n";
	    }
	  to << "\n";
	  to.setf(ios::fmtflags(0), ios::floatfield);
	  print_dist(to);
	}
      }
    to << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
  }// face cycle
}

std::map<int, Ref<DistMultiFluxArrayCalc> > state_map;
typedef map<int, Ref<DistMultiFluxArrayCalc> >::iterator state_map_it;

std::map<Sid, std::vector<Samp> > smp_map;
typedef std::map<Sid, std::vector<Samp> >::iterator smp_map_it;

string save_file, smp_save_file;

extern "C" void signal_handler (int sig)
{
  static const char funame [] = "signal_handler: ";

  switch(sig) {
    //
  case SIGUSR1:
    //
    return;

  case SIGFPE:
    //
    std::cerr << Comm::mesg() << funame << "Oops!!! Floating-point exception" << std::endl; 

    throw Error::FPE();
    
  default:
      //
      std::cerr << Comm::mesg() << funame << "Oops!!! signal " 
		<< sig << " has been received" << std::endl;
						  
    // save flux
    if(!Comm::rank()) {
      ofstream to(save_file.c_str());
      if(to)
	for (state_map_it it = state_map.begin(); it != state_map.end(); ++it)
	  to << it->first << " " << it->second.value().size() << "\n" << it->second;
      to.close();
      to.clear();
      // save samplings
      if(FluxBase::smp_out_flag) {
	to.open(smp_save_file.c_str());
	for(smp_map_it it=smp_map.begin(); it != smp_map.end(); ++it) {
	  to << it->first << it->second.size() << "\n";
	  for(int i = 0; i < it->second.size(); ++i)
	    to << it->second[i];
	}
      }
      // report
      cout << Comm::mesg() << funame 
	   << "intermediate flux data saved in " 
	   << save_file;
      if(FluxBase::smp_out_flag)
	cout << " file, intermediate samplings saved in "
	     << smp_save_file << " file \n";
      else
	cout << " file\n";
    }

    exit(1);
  }
}

int main (int argc, char* argv[])
{// main

  const char funame [] = "main: ";
  const char debug_mesg [] = "DEBUG: ";

  // temporary variables
  int itemp;
  Sid sid;

  /*************************************************************
   *********************  MPI STARTS HERE!  ********************
   *************************************************************/
  // test for floating-point exception
  //
  Comm::init(argc, argv); // working directory changed here!

  test_except(Comm::mesg());

  //
  // floating point exception trap
  //
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    
  // Random numbers initialization
  //
  Random::init();

  // Parameters input
  //
  map <string, Read> input_data;
  typedef map<string, Read>::iterator input_data_it;
 
  // grids
  std::vector<double> tmpr_grid, ener_grid, amom_grid;
  input_data ["tmpr_grid"] = Read(tmpr_grid);
  input_data ["ener_grid"] = Read(ener_grid);
  input_data ["amom_grid"] = Read(amom_grid);

  int pot_smp_min, pot_smp_max, tot_smp_min, tot_smp_max, pot_smp_len;
  input_data ["pot_smp_min"] = Read(pot_smp_min);
  input_data ["pot_smp_max"] = Read(pot_smp_max);
  input_data ["tot_smp_min"] = Read(tot_smp_min);
  input_data ["tot_smp_max"] = Read(tot_smp_max);
  input_data ["pot_smp_len"] = Read(pot_smp_len);
  
  double flux_rel_err; // maximal flux relative error
  input_data ["flux_rel_err"]   = Read(flux_rel_err);
  input_data ["min_atom_dist"]  = Read(MIN_ATOM_DIST, 1.5);

  // potential specification
  string pot_name, pot_type,  opt_method, opt_type;
  int pes_size;
  input_data ["pot_name"]  = Read(pot_name);
  input_data ["pot_type"]  = Read(pot_type);
  input_data ["pes_size"]  = Read(pes_size, 1);
  input_data ["opt_method"]  = Read(opt_method);
  input_data ["opt_type"]  = Read(opt_type, "none");

  // specific methods
  std::string gauss_file, molpro_file, multipole_file, sjk_file;
  input_data["gauss_inp_file"] = Read(gauss_file, "gauss.inp");
  input_data["molpro_inp_file"] = Read(molpro_file, "molpro.inp");
  input_data["sjk_inp_file"] = Read(sjk_file, "sjk.inp");
  input_data["multipole_inp_file"] = Read(multipole_file, "multipole.inp");

  // sampling specification
  string smp_type, ds_inp_file;
  vector<int> ds_face, ds_symm;
  input_data ["sampling"]  = Read(smp_type, "multifacet");
  input_data ["ds_inp_file"]  = Read(ds_inp_file, "divsur.inp");
  input_data ["face"]  = Read(ds_face, vector<int>());
  input_data ["face_symm"]  = Read(ds_symm, vector<int>());

  // Input/Output
  string mol_spec_file; // molecular specifications file name
  string flux_file; // flux output file
  input_data ["mol_spec_file"]  = Read(mol_spec_file, "structure.inp");
  input_data ["save_file"]  = Read(save_file, "flux.save");      
  input_data ["flux_file"]  = Read(flux_file, "flux.dat");      

  int therm_flux_out;
  string therm_flux_name;
  input_data ["therm_flux_out"]  = Read(therm_flux_out, 1);      
  input_data ["therm_flux_name"]  = Read(therm_flux_name, "flux.out");      

  // raw sampling output
  string raw_smp_file;
  input_data ["raw_smp_file"] = Read(raw_smp_file, "smp.raw");
  input_data ["smp_out_flag"] = Read(FluxBase::smp_out_flag, 0);
  input_data ["smp_save_file"] = Read(smp_save_file, "smp.save");

  // environment
  int run_time, time_step, flux_debug;
  string work_dir; // slave working directory
  string log_file; // local log file
  string log_level;
  input_data ["log_file"]   = Read(log_file);
  input_data ["log_level"]  = Read(log_level, "info");
  input_data ["work_dir"]   = Read(work_dir);  
  input_data ["run_time"]   = Read(run_time, 3000000);      
  input_data ["save_time"]  = Read(time_step, 300);      
  input_data ["flux_debug"]  = Read(flux_debug, 0);      
  input_data ["comm_debug"] = Read(Comm::debug, 0);      
  int nice_val;     // nice increment
  input_data ["nice"]        = Read(nice_val, 0);      
  
  if (argc < 2) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame
	   << "usage: multi input_file\n";
    exit(1);
  }

  ifstream from(argv[1]);
  if(!from) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame 
	   << "input file " << argv[1] << " is not found\n";
    exit(1);
  }

  string key, comment;
  while(from >> key) {
    input_data_it it = input_data.find(key);
    if (it == input_data.end())
      getline(from, comment);
    else
      from >> it->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  bool is_init = true;
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(!it->second.is_init()) {
      if(!Comm::rank())
	cout << Comm::mesg() << funame 
	     << it->first << " is not initialized\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  if(!Comm::rank()) {
    std::cout << "Default parameters:\n";
    for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
      if(it->second.is_default())
	cout << it->first << "  " << it->second << "\n";
    std::cout << "\n";
  }

  //setting logging level
  if(log_level == "error")
    Log::set_level(Log::ERROR);
  else if(log_level == "warning")
    Log::set_level(Log::WARN);
  else if(log_level == "note")
    Log::set_level(Log::NOTE);
  else if(log_level == "info")
    Log::set_level(Log::INFO);
  else if(log_level == "debug")
    Log::set_level(Log::DEBUG);
  else {
    if(!Comm::rank())
      cout << Comm::mesg() << funame << "unknown logging level: " << log_level
	   << ". Available levels: error, warning, note, info, debug\n";
    exit(1);
  }

  // check faces
  if(ds_face.size() == 1 && !ds_symm.size())
    ds_symm.resize(1, 1);
  else if(ds_face.size() != ds_symm.size()) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame 
	   << "number of symmetry factors should be the same as the number of faces\n"; 
    exit(1);
  }

  // initialize molecular structure
  mol_init(mol_spec_file);

  // dividing surfaces specification
  surf_init(ds_inp_file.c_str());

  // sampling function specification
  name2ranp(smp_type);

  // flux static variables initialization
  FluxBase::stat_init(tmpr_grid, ener_grid, amom_grid);
  
  // potential
  PES::pot = name2pot(pot_name);
  PES::set_size(pes_size);

  // Initializing flux calculation static variables
  Calc::stat_init(flux_rel_err, pot_smp_max, pot_smp_min, 
		  tot_smp_max, tot_smp_min, pot_smp_len);

  // set the process' file creation mask to 022
  umask(S_IWGRP | S_IWOTH);

  // setting nice
  if(!Comm::rank() && nice_val < 0)
    cout << Comm::mesg() << funame 
	 << "WARNING!!! negative nice increment, ignoring" 
	 << endl;
  if(Comm::rank() && nice_val > 0)
    nice(nice_val);

  // list nodes
  char hostname [100];
  gethostname(hostname, 99);
  cout << Comm::mesg() << funame 
       << "host = " << hostname 
       << ", process =  " << getpid() << endl;

  MPI::COMM_WORLD.Barrier();

  // need to block signals while writing
  sigset_t block_sig;
  sigset_t old_sig;
  sigfillset(&block_sig); // block all signals

  struct sigaction sigact;
  sigact.sa_handler = &signal_handler;
  sigact.sa_mask = block_sig;
  sigact.sa_flags = 0;

  // signals to save the state map and continue
  sigaction(SIGUSR1, &sigact, 0);
  sigaction(SIGUSR2, &sigact, 0);
  
  // signals to save the step map and exit
  sigaction(SIGINT,  &sigact, 0);
  sigaction(SIGTERM, &sigact, 0);
  sigaction(SIGHUP,  &sigact, 0);
  sigaction(SIGABRT, &sigact, 0);
  sigaction(SIGQUIT, &sigact, 0);
  sigaction(SIGTRAP, &sigact, 0);
  sigaction(SIGALRM, &sigact, 0);
  
  sigaction(SIGFPE,  &sigact, 0);
  sigaction(SIGILL,  &sigact, 0);
  sigaction(SIGBUS,  &sigact, 0);
  sigaction(SIGSEGV, &sigact, 0);
  sigaction(SIGPIPE, &sigact, 0);

  /***********************************************************************
   ****************************** MASTER *********************************
   ***********************************************************************/
  //
  if (Comm::rank() == 0) {
    //
    // set run-time
    //
    alarm(run_time);

    // set save time
    //
    time_t save_time = time(0) + (time_t)time_step;

    cout << "\n"<< Comm::mesg() << funame << "Starting calculation ...\n";
    
    // initialize available nodes
    //
    set<int> avail_nodes;
    for (int i = 1; i < Comm::size(); ++i)
      avail_nodes.insert(i);

    cout << endl;

    // initialize Sid-to-State map
    int curr_surf_ind = 0; // first unsampled surface
    MultiFlux temp_fl;

    int ds, ds_size;
    from.open(save_file.c_str());
    while (from >> ds >> ds_size) {// read the dividing surface index & the number of faces
      if(state_map.find(ds) != state_map.end()) {
	std::cout << debug_mesg << Comm::mesg() << funame 
		  << "duplicate surface index " << ds
		  << " in the temporary flux data file " 
		  << save_file << std::endl;
	exit(1);
      }

      state_map[ds] = Ref<DistMultiFluxArrayCalc>(DistMultiFluxArrayCalc(ds_size));
      if(!(from >> state_map[ds])) {
	std::cout << debug_mesg << Comm::mesg() << funame 
		  << "temporary flux file " << save_file 
		  << " is corrupted\n";
	exit(1);
      }

      if(therm_flux_out) {
	thermal_flux[ds] = state_map[ds].value();
	thermal_flux[ds].value().normalize();
      }

      if(ds >= curr_surf_ind)
	curr_surf_ind = ds + 1;
    }
    from.close();
    from.clear();

    if(flux_debug)
      cout << debug_mesg << Comm::mesg() << funame 
	   << "save file scanned, current surface index = "
	   << curr_surf_ind << endl;

    // search for the last processed surface index in the flux file
    std::set<int> finished_set;
    from.open(flux_file.c_str());
    while (from >> ds >> ds_size) {
      if(!finished_set.insert(ds).second) {
	std::cout << debug_mesg << Comm::mesg() << funame 
		  << "duplicate surface index " << ds
		  << " in the flux data file " << flux_file 
		  << std::endl;
	exit(1);
      }
      if(ds >= curr_surf_ind)
	curr_surf_ind = ds + 1;

      // read flux data
      MultiFluxArray srfl(ds_size);
      from >> srfl;

      if(therm_flux_out)
	thermal_flux[ds] = srfl;

      if(!from) {
	std::cout << debug_mesg << Comm::mesg() << funame 
		  << "flux file " << flux_file << " is corrupted\n";
	exit(1);
      }
    }
    from.close();
    from.clear();

    if(flux_debug)
      cout << debug_mesg << Comm::mesg() << funame 
	   << "flux file scanned, current surface index = "
	   << curr_surf_ind << endl;

    // check if the flux.dat and flux.save are consistent
    for(int i = 0; i < curr_surf_ind; ++i) {
      bool is_sav = (state_map.find(i) != state_map.end());
      bool is_fin = (finished_set.find(i) != finished_set.end());
      if(is_sav && is_fin) {
	std::cout << "master: duplicate surface index " << i
		  << " in the  files " << flux_file 
		  << " and " << save_file << std::endl;
	exit(1);
      }
      if(!is_sav && !is_fin) {
	std::cout << "master: surface index " << i
		  << " is not found in " << flux_file 
		  << " and " << save_file << std::endl;
	exit(1);
      }
    }

    // output
    ofstream flux_out(flux_file.c_str(), ios::app);

    Samp samp;

    // raw sampling output
    ofstream smp_out;
    if(FluxBase::smp_out_flag && curr_surf_ind) {
      from.open(smp_save_file.c_str());
      while(from >> sid) {
	from >> itemp;
	smp_map[sid].resize(itemp);
	for(int i = 0; i < smp_map[sid].size(); ++i)
	  from >> smp_map[sid][i];
      }
      from.close();
      from.clear();

      if(curr_surf_ind)
	smp_out.open(raw_smp_file.c_str(), ios::app);
      else
	smp_out.open(raw_smp_file.c_str());
    }
	
    while (1) {// dividing surface cycle
      //
      while (avail_nodes.size()) {//send
	//
	/**********************************************************
	 ************    SEND REQUEST TO A SLAVE    ***************
	 **********************************************************/
	int node = *avail_nodes.begin();
	bool is_node_used = false;
	state_map_it it = state_map.begin();
	int sm_ind = 0;
	int work_face;
	while (it != state_map.end()) {// state
	  switch (it->second.value().check_state(work_face)) {// sampling	    
	  case Calc::WAIT:// skip sampling
	    it->second.value().wait();
	    ++it;
	    ++sm_ind;
	    continue;
	  case Calc::FLUX:// sent for flux sampling
	    Sid(it->first, work_face).send(node, Comm::FLUX_TAG);

	    if(flux_debug)
	      cout << debug_mesg << Comm::mesg() << it->first 
		   << "-th surface, " << work_face 
		   << "-the face: flux sampling requested" << endl;

	    avail_nodes.erase(node);
	    it->second.value()[work_face].add_node(node);
	    is_node_used = true;
	    break;
	  case Calc::SURF:// sent for surface sampling
	    it->second.value().wait();
	    Sid(it->first, 0).send(node, Comm::SURF_TAG);

	    if(flux_debug)
	      cout << debug_mesg << Comm::mesg() << it->first 
		   << "-th surface: volume sampling requested" << endl;
	    try {
	      MPI::COMM_WORLD.Send(it->second.value().surf_smp_num(), 
				   it->second.value().size(), 
				   MPI::INT, node, Comm::SURF_TAG);
	      if (Comm::debug)
		cout << Comm::mesg() << funame 
		     << Comm::SURF_TAG << "-th tag send to "
		     << node << "-th node" << endl;
	    }
	    catch (MPI::Exception) {
	      cout << Comm::mesg() << funame 
		   << "error in sending smp_num to " 
		   << node << "-th node" << endl;
	      throw;
	    }
	    avail_nodes.erase(node);
	    it->second.value()[0].add_node(node);
	    is_node_used = true;
	    break;
	  case Calc::STOP:// finish sampling
	    // message
	    cout << Comm::mesg() << funame << "flux through the "
		 << it->first << "-th surface is done" << endl;

	    sigprocmask(SIG_SETMASK, &block_sig, &old_sig);// block signals
	    // save flux
	    flux_out << it->first << " " << it->second.value().size() << "\n";
	    it->second.value().normalize();
	    flux_out << it->second;
	    flux_out.flush();
	    // save sampling data
	    for(int face=0; face < it->second.value().size(); ++face) {
	      smp_out << "Surface = " << it->first << " Facet = " << face << "\n";
	      for(int i = 0; i< smp_map[Sid(it->first, face)].size(); ++i) {
		// print energies
		smp_out << "Energies(kcal/mol): ";
		for(int en = 0; en < PES::size(); ++en)
		  smp_out << smp_map[Sid(it->first, face)][i].energy(en)/Phys_const::kcal
			  << " ";
		smp_out << "\n";
		// print geometry
		smp_map[Sid(it->first, face)][i].attach();
		update_atoms();
		smp_out << "Geometry(Angstrom):\n";
		print_geom(smp_out);
		smp_out << "\n";
	      }
	      // remove from samplings map
	      smp_map.erase(Sid(it->first, face));
	    }
	    smp_out.flush();

	    // remove the state from the state map
	    itemp = distance(state_map.begin(), it);
	    state_map.erase(it);
	    it = state_map.begin();
	    advance(it, itemp);

	    // update save file
	    {
	      ofstream save_out(save_file.c_str());
	      if(save_out)
		for (state_map_it i = state_map.begin(); i != state_map.end(); ++i)
		  save_out << i->first << " " << i->second.value().size() << "\n" << i->second;
	    }

	    sigprocmask(SIG_SETMASK, &old_sig, 0);// unblock signals
	    continue;
	  default:
	    error("master: wrong case, exitting");
	  }// sampling
	  break;
	}// state

	if (is_node_used)
	  continue;

	// ... did not find in the state map => then start new surface
	if (curr_surf_ind < ds_array.size()) {
	  sigprocmask(SIG_SETMASK, &block_sig, &old_sig);// block signals

	  if(ds_face.size())
	    state_map[curr_surf_ind] = 
	      Ref<DistMultiFluxArrayCalc>(DistMultiFluxArrayCalc(ds_face.size()));
	  else
	    state_map[curr_surf_ind] = 
	      Ref<DistMultiFluxArrayCalc>(DistMultiFluxArrayCalc(ds_array[curr_surf_ind]->size()));

	  sigprocmask(SIG_SETMASK, &old_sig, 0);// unblock signals

	  Sid(curr_surf_ind, 0).send(node, Comm::FLUX_TAG); // send sid
	  if(flux_debug)
	    cout << debug_mesg << Comm::mesg() << curr_surf_ind 
		 << "-th surface, 0-th face: "
	      "first flux sampling requested" << endl;

	  avail_nodes.erase(node); // update available nodes database
	  state_map[curr_surf_ind].value()[0].add_node(node);
	  is_node_used = true;
	  ++curr_surf_ind;
	}
	else
	  break;

      }// send 

      /*************************************************************
       ********************   STOP CONDITIONS   ********************
       *************************************************************/
      //
      // no work to be done
      //
      if(curr_surf_ind == ds_array.size() && state_map.size() == 0) {
	//
	cout << Comm::mesg() << funame << "no work to be done => exitting ..." << endl;
	
	break;
      }

      //no nodes available
      //
      int run_nodes = 0;
      
      for (state_map_it it = state_map.begin(); it != state_map.end(); ++it)
	//
	run_nodes += it->second.value().node_size();
      
      if (!avail_nodes.size() && !run_nodes)  {
	//
	cout << Comm::mesg() << funame << "no nodes are available => exitting..." << endl;
	
	break;
      }

      /*************************************************************
       **********     RECEIVE FLUX FROM A SLAVE    *****************
       *************************************************************/
      //
      MPI::Status stat = sid.master_recv(); // <- wait here!

      switch(stat.Get_tag()) {// receive

      case Comm::RAND_TAG:// Random Number Generator Seed
	Random::send_seed(stat.Get_source());

	break;
      case Comm::FAIL_TAG:// Node Failure
	cout << Comm::mesg() << funame << "node " << stat.Get_source()
	     << " failed; EXCLUDED!!!" << endl;
	//Comm::set_state(stat.Get_source(), Comm::FAIL_STATE);

	// update states map
	state_map[sid.num()].value()[sid.face()]
	  .remove_node(stat.Get_source());
	state_map[sid.num()].value().nowait();
	break;
	

      case Comm::SAMP_TAG: // sampling
	MPI::COMM_WORLD.Recv(&itemp, 1, MPI::INT, stat.Get_source(), 
			     Comm::SAMP_TAG);
	
	for(int i = 0; i < itemp; ++i) {
	  samp.recv(stat.Get_source(), Comm::SAMP_TAG);
	  sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
	  smp_map[sid].push_back(samp);
	  sigprocmask(SIG_SETMASK, &old_sig, 0);
	}
	break;
	

      case Comm::FLUX_TAG:// Flux Sampling
	if (state_map.find(sid.num()) == state_map.end())
	  error("master: the sid is not the state map, exitting");

	// receive the  flux data
	temp_fl.master_recv(stat.Get_source());
	if(flux_debug)
	  cout << debug_mesg << Comm::mesg() << sid.num() 
	       << "-th surface, "  << sid.face()
	       << "-th face: flux sampling acknowledged" << endl;

	// update the flux data
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
	state_map[sid.num()].value()[sid.face()] += temp_fl;
	sigprocmask(SIG_SETMASK, &old_sig, 0);

	if(therm_flux_out) {
	  thermal_flux[sid.num()] = state_map[sid.num()].value();
	  thermal_flux[sid.num()].value().normalize();
	  print_thermal_flux(therm_flux_name);
	}

	if (avail_nodes.find(stat.Get_source()) 
	    != avail_nodes.end())
	  error("master: the set of avail. nodes is corrupted, exitting");

	// update available nodes
	avail_nodes.insert(stat.Get_source());
	state_map[sid.num()].value()[sid.face()].remove_node(stat.Get_source());
	state_map[sid.num()].value().nowait();
	break;

      case Comm::SURF_TAG:// Surface Sampling
	if (state_map.find(sid.num()) == state_map.end())
	  error("master: the sid is not the state map, exitting");

	// receive the surface sampling results
	Array<FluxBase> smp_res(state_map[sid.num()].value().size());
	try {
	  for(int i = 0; i < smp_res.size(); ++i)
	    smp_res[i].recv(stat.Get_source() , Comm::SURF_TAG);

	  if (Comm::debug)
	    cout << Comm::mesg() << funame 
		 << Comm::SURF_TAG << "-th tag received from "
		 << stat.Get_source() << "-th node" << endl;
	}
	catch (MPI::Exception) {
	  cout << Comm::mesg() << funame 
	       << "error in receiving surface sampling results" 
	       << endl;
	  throw;
	}
	if(flux_debug)
	  cout << debug_mesg << Comm::mesg() << sid.num() 
	       << "-th surface: volume sampling acknowledged" << endl;

	// update the flux data
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	for (int i = 0; i < state_map[sid.num()].value().size(); ++i)
	  static_cast<FluxBase&>(state_map[sid.num()].value()[i]) += 
	    smp_res[i];

	if(therm_flux_out) {
	  thermal_flux[sid.num()] = state_map[sid.num()].value();
	  thermal_flux[sid.num()].value().normalize();
	  print_thermal_flux(therm_flux_name);
	}

	//  write total flux to a file
	flux_out << sid.num() << " " << state_map[sid.num()].value().size() << "\n";
	state_map[sid.num()].value().normalize();
	flux_out << state_map[sid.num()];
	flux_out.flush();
	// save sampling data
	for(int face=0; face < state_map[sid.num()].value().size(); ++face) {
	  smp_out << "Surface = " << sid.num() << " Facet = " << face << "\n";
	  for(int i = 0; i< smp_map[Sid(sid.num(), face)].size(); ++i) {
	    // print energies
	    smp_out << "Energies(kcal/mol): ";
	    for(int en = 0; en < PES::size(); ++en)
	      smp_out << smp_map[Sid(sid.num(), face)][i].energy(en)/Phys_const::kcal
		      << " ";
	    smp_out << "\n";
	    // print geometry
	    smp_map[Sid(sid.num(), face)][i].attach();
	    update_atoms();
	    smp_out << "Geometry(Angstrom):\n";
	    print_geom(smp_out);
	    smp_out << "\n";
	  }
	  // remove from samplings map
	  smp_map.erase(Sid(sid.num(), face));
	}
	smp_out.flush();

	//  remove the state from the state map
	state_map.erase(sid.num());

	// update save file
	{
	  ofstream save_out(save_file.c_str());
	  if(save_out)
	    for (state_map_it i = state_map.begin(); i != state_map.end(); ++i)
	      save_out << i->first << " " << i->second.value().size() << "\n" << i->second;
	}

	// unblock the signals
	sigprocmask(SIG_SETMASK, &old_sig, 0);

	// update set of available nodes
	if (avail_nodes.find(stat.Get_source()) != avail_nodes.end())
	  error("master: the set of avail. nodes is corrupted, exitting");
	avail_nodes.insert(stat.Get_source());

	// message
	cout << Comm::mesg() << funame 
	     << "flux through the " << sid.num() 
	     << "-th surface is done" << endl;
	break;

      }// receive

      // save flux data at time_step intervals
      //
      time_t new_time = time(0);
      if(new_time > save_time) {// save
	save_time = new_time + (time_t)time_step;
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	ofstream to(save_file.c_str());
	if(to)
	  for (state_map_it it = state_map.begin(); it != state_map.end(); ++it)
	    to << it->first << " " << it->second.value().size() << "\n" << it->second;
	to.close();
	to.clear();
	if(FluxBase::smp_out_flag) {
	  to.open(smp_save_file.c_str());
	  for(smp_map_it it=smp_map.begin(); it != smp_map.end(); ++it) {
	    to << it->first  << it->second.size() << "\n";
	    for(int i = 0; i < it->second.size(); ++i)
	      to << it->second[i];
	  }
	}

	sigprocmask(SIG_SETMASK, &old_sig, 0);
      
	cout << Comm::mesg() << funame << "intermediate flux data saved in " << save_file << " file\n";
	
	if(FluxBase::smp_out_flag)
	  //
	  cout << Comm::mesg() << funame << "intermediate samplings saved in " << smp_save_file << " file\n";
	//
      }// save
      //
    }// dividing surface cycle
   
    // send signal for slaves to finish
    //
    for (int i = 1; i < Comm::size(); ++i)
      //
      Sid().send(i, Comm::STOP_TAG);

    // remove save file
    //
    remove(save_file.c_str());
    
    remove(smp_save_file.c_str());

  }// master
  //
  /********************************************************************
   **************************    SLAVE    *****************************
   ********************************************************************/    
  //
  else {
    //
    alarm(run_time + 5);

    // initialize potential environment
    //
    // gaussian
    //
    if (pot_type == "g98" || opt_method == "gauss")
      //
      Gauss::init(gauss_file);

    // molpro
    //
    if(pot_type == "molpro" || opt_method == "molpro") {
      //
      Molpro::init(molpro_file);
      
      std::ostringstream oss;
      
      oss << Molpro::wf_name << Comm::rank();
      
      Molpro::wf_name = oss.str();
    }
    
    // multipole
    //
    if(pot_type == "multipole")
      //
      Multipole::init(multipole_file);

    // user supplied analytical potential
    //
    if(pot_type == "sjk")
      //
      Sjk::init(sjk_file);

    // optimaztion method
    //
    set_opt(opt_type, opt_method);

    // all common files should be open before changing directory!!!
    //
    Comm::change_work_dir(work_dir);

    // set local logging
    //
    Log::out.open(log_file.c_str(), ios::app);

    bool is_run = true;
    
    // loop
    //
    while(is_run) {
      //
      Div_surf ds;
      
      MultiFlux fl;

      // sampling
      //
      switch(sid.slave_recv()) {
	//
	// stop sampling
	//
      case Comm::STOP_TAG:
	//
	cout << Comm::mesg() << funame << ": exitting..." << endl;
	
	is_run = false;
	
	break;

	// sample surface
	//
      case Comm::SURF_TAG: {
	//
	ds = *ds_array[sid.num()];
	  
	Array<int> smp_num(ds_face.size() ? ds_face.size(): ds.size());

	try {
	  //
	  MPI::COMM_WORLD.Recv(smp_num.data(), smp_num.size(),  
			       MPI::INT, 0, Comm::SURF_TAG);
	  if (Comm::debug)
	    //
	    cout << Comm::mesg() << funame << Comm::SURF_TAG << "-th tag received from master" << endl;
	}
	catch (MPI::Exception) {
	  //
	  cout << Comm::mesg() << funame << "error in receiving sampling numbers" << endl;
	  
	  throw;
	}

	// sampling surface
	//
	std::vector<FluxBase> smp_res;
	
	for (int face = 0; face < smp_num.size(); ++face) {
	  //
	  if(ds_face.size()) {
	    //
	    ds.set_face(ds_face[face]);
	  }
	  else
	    //
	    ds.set_face(face);
	  
	  smp_res.push_back(FluxBase(ds, smp_num[face]));
	}

	// send results
	//
	sid.send(0, Comm::SURF_TAG);
	
	try {
	  //
	  for(int i = 0; i < smp_res.size(); ++i)
	    //
	    smp_res[i].send(Comm::SURF_TAG);
	  
	  if (Comm::debug)
	    //
	    cout << Comm::mesg() << funame << Comm::SURF_TAG << "-th tag sent to master" << endl;
	}
	catch (MPI::Exception) {
	  //
	  cout << Comm::mesg() << funame << "error in sending surface sampling results to master" << endl;
	  
	  throw;
	}
      }
	
	break;

	// sample flux
	//
      case Comm::FLUX_TAG:
	//
	try {
	  //
	  Div_surf ds = *ds_array[sid.num()];
	  
	  if(ds_face.size()) {
	    //
	    ds.set_face(ds_face[sid.face()]);
	  }
	  else
	    //
	    ds.set_face(sid.face());
	  
	  std::vector<Samp> smp_array(pot_smp_len);

	  Log::out << "\n("<< sid.num() << "," << sid.face() << ")-th SURFACE:\n\n";

	  MultiFlux fl(ds, pot_smp_len, smp_array);

	  // symmetry factors
	  //
	  if(ds_face.size())
	    //
	    fl *= (double)ds_symm[sid.face()];

	  // send flux
	  //
	  sid.send(0, Comm::FLUX_TAG);
	  
	  fl.slave_send();

	  // send successful samplings
	  //
	  if(FluxBase::smp_out_flag) {
	    //
	    sid.send(0, Comm::SAMP_TAG);
	    
	    MPI::COMM_WORLD.Send(&pot_smp_len, 1, MPI::INT, 0, Comm::SAMP_TAG);
	    
	    for(int i = 0; i < pot_smp_len; ++i)
	      //
	      smp_array[i].send(0, Comm::SAMP_TAG);
	  }
	}
	catch(Error::Err) {
	  //
	  cout << Comm::mesg() << funame << "flux sampling failed\n";
	  
	  sid.send(0, Comm::FAIL_TAG);
	}

	break;
	
      default:
	//
	cout << Comm::mesg() << funame << "slave: wrong tag, exitting\n";
	
	exit(1);
	//
      }// sampling
      //
    }// loop
    //
  }// slave

  MPI::COMM_WORLD.Barrier();
  
  MPI::Finalize();
  
  sleep(1);

  return 0;
}
