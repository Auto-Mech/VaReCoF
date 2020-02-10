#include<vector>
#include<string>
#include<fstream>
#include<map>
#include<iomanip>
#include<set>
#include <exception>

#include<cstdio>
#include<cmath>
#include<unistd.h>
#include<sys/types.h>
#include <sys/stat.h>
#include <csignal>
#include <ctime>
#include <sstream>
#include <cstdlib>

#include <mpi.h>

#include "input.hh"
#include "units.hh"
#include "rotd.hh"    
#include "force.hh"
#include "gauss.hh"
#include "run.hh"
#include "error.hh"
#include "comm.hh"
#include "ref.hh"
#include "flux.hh"
#include "multipole.hh"
#include "sjk.hh"
#include "molpro.hh"
#include "log.hh"

map<int, Ref<Traj::State> > start_set; // trajectories started
map<int, Ref<Traj::State> > run_set;   // running trajectories
set<int> fail_set;                     // failed trajectories
set<int> init_set;                     // trajectories not yet started
map<int, Traj::Res> res_set;           // trajectories finished

vector<Samp> smp_array;                // initial configurations

typedef map<int, Ref<Traj::State> >::const_iterator Run_iter;
typedef set<int>::const_iterator Fail_iter;
typedef map<int, Traj::Res>::const_iterator Res_iter;

sigset_t block_sig;
sigset_t old_sig;

int start_new_trajectory (int node) {

  static const char funame [] = "start_new_trajectory: ";

  int traj = -1;
  if(start_set.size()) {// start set
    // trajectory to run
    traj = start_set.begin()->first;
    if(Traj::debug >= Traj::DEBUG &&
       (run_set.find(traj)  != run_set.end() ||
	fail_set.find(traj) != fail_set.end() ||
	res_set.find(traj)  != res_set.end() ||
	init_set.find(traj) != init_set.end())) {
      cout << Traj::debug_mesg << Comm::mesg() << funame << traj 
	   << "-th trajectory  in a wrong set\n";
      exit(1);
    }

    // update control sets
    sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
    run_set[traj] = start_set.begin()->second;
    start_set.erase(start_set.begin());
    sigprocmask(SIG_SETMASK, &old_sig, 0);

    // send the data to the working node
    MPI::COMM_WORLD.Send(&traj, 1, MPI::INT, node, Comm::RUN_TAG);
    if(Comm::debug)
      cout << Comm::mesg() << funame
	   << Comm::RUN_TAG << "-th tag sent to " 
	   << node << "-th node" << endl;
    run_set[traj].value().send(node, Comm::RUN_TAG);

    // debug output
    if(Traj::debug >= Traj::DEBUG)
      cout << Traj::debug_mesg << Comm::mesg() << funame 
	   << traj << "-th trajectory restarted on "
	   << node <<"-th node" << endl;

    return traj;
  }// start set

  if(init_set.size()) {// init set
    // trajectory to run
    traj = *init_set.begin();
    if(Traj::debug >= Traj::DEBUG &&
       (run_set.find(traj)   != run_set.end() ||
	fail_set.find(traj)  != fail_set.end() ||
	res_set.find(traj)   != res_set.end() ||
	start_set.find(traj) != start_set.end())) {
      cout << Traj::debug_mesg << Comm::mesg() << funame << traj 
	   << "-th trajectory  in a wrong set\n";
      exit(1);
    }

    // update control sets
    sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
    run_set[traj] = smp_array[traj / 2];
    sigprocmask(SIG_SETMASK, &old_sig, 0);
    init_set.erase(init_set.begin());

    // send data to the working node
    MPI::COMM_WORLD.Send(&traj, 1, MPI::INT, node, Comm::RUN_TAG);
    if(Comm::debug)
      cout << Comm::mesg() << funame
	   << Comm::RUN_TAG << "-th tag sent to " << node 
	   << "-th node" << endl;
    run_set[traj].value().send(node, Comm::RUN_TAG);

    // debug output
    if(Traj::debug >= Traj::DEBUG)
      cout << Traj::debug_mesg << Comm::mesg() << funame
	   << "new " << traj << "-th trajectory started on "
	   << node <<"-th node" << endl;

    return traj;
  }// init set

  // debug output
  if(Traj::debug >= Traj::DEBUG)
    cout << Traj::debug_mesg << Comm::mesg() << funame
	 << "no trajectory started, " << node 
	 << "-th node is idle" << endl;

  return traj;
}

void save_traj_data (const string& sf)
{
  static const char funame [] = "save_traj_data: ";

  ofstream to(sf.c_str());
  // failed trajectories
  to << fail_set.size() << "\n";
  for(Fail_iter it = fail_set.begin(); it != fail_set.end(); ++it)
    to << *it << "  ";
  to << "\n\n";

  // finished trajectories
  to << res_set.size() << "\n";
  int count = 0;
  for(Res_iter it = res_set.begin(); it != res_set.end(); ++it) {
    to << setw(4) << it->first << " " 
       << it->second << " ";
    if(++count >= 7) {
      count = 0;
      to << "\n";
    }
  }
  if(count)
    to << "\n";
  to << "\n";
  
  // running trajectories
  to << start_set.size() + run_set.size() << "\n";
  for(Run_iter it = run_set.begin(); it != run_set.end(); ++it)
    to << it->first << "\t" << it->second << "\n";
  for(Run_iter it = start_set.begin(); it != start_set.end(); ++it)
    to << it->first << "\t" << it->second << "\n";

  if(to) {
  if(Traj::debug >= Traj::INFO)
    cout << Comm::mesg() << funame
	 << "trajectory data have been saved to "
	 << sf << " file" << endl;
  }
  else
    cout << Comm::mesg() << funame << "failed\n";
}

string traj_save_file; // backup file name

extern "C" void signal_handler (int sig)
{
  static const char funame [] = "signal_handler: ";

  cout << Comm::mesg() << funame << "Oops!!! signal " 
       << sig << " has been received\n";

  if(!Comm::rank())
    save_traj_data(traj_save_file);
  
  if(sig == SIGUSR1)
    return;

  exit(1);
}

int main (int argc, char* argv[])
{// main

  const char funame [] = "main: ";

  /*************************************************************
   *********************  MPI STARTS HERE!  ********************
   *************************************************************/
  Comm::init(argc, argv); // working directory changed here!

  /*************** Input of parameters from rotd.inp ********************/

  map<string, Read> input_data;
  typedef map<string, Read>::iterator Inter;

  // The potential 
  string pot_name, pot_type, pot_file;
  input_data ["pot_name"]  = Read(pot_name);
  input_data ["pot_type"]  = Read(pot_type);
  input_data ["pot_file"]  = Read(pot_file);

  // molecular specification
  string mol_spec_file;
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");

  // Input/Output
  string wght_smp_file;  // random samplings file
  input_data ["wght_smp_file"] = Read(wght_smp_file);
  input_data ["traj_save_file"] = Read(traj_save_file, "rfactor.save");

  // dividing surfaces definitions
  // transition state surface
  string ds_inp_file;
  int ds_ind;
  input_data ["ds_inp_file"]  = Read(ds_inp_file, "divsur.inp");
  input_data ["ds_ind"]       = Read(ds_ind);
  // reactive surface
  string reac_surf_file;
  double reac_ener_max, reac_ener_step;
  int reac_ener_num;
  input_data ["is_reac_ener"]   = Read(Traj::reac_ener_flag);
  input_data ["reac_ener_max"]  = Read(reac_ener_max);
  input_data ["reac_ener_step"] = Read(reac_ener_step);
  input_data ["reac_ener_num"]   = Read(reac_ener_num);
  input_data ["is_reac_surf"]   = Read(Traj::reac_surf_flag);
  input_data ["reac_surf_file"]  = Read(reac_surf_file);
  // dissociative surface
  double traj_diss_dist;
  input_data ["traj_diss_dist"]   = Read(traj_diss_dist);

  // Trajectory propagation parameters
  input_data ["traj_time_step"] = Read(Traj::time_step, 10.);
  input_data ["traj_rel_tol"]   = Read(Traj::rel_tol, 1.e-5);
  input_data ["traj_abs_tol"]   = Read(Traj::abs_tol,  1.e-5);

  input_data ["traj_aux_out"]   = Read(Traj::aux_out_flag);
  input_data ["traj_ener_out"]  = Read(Traj::ener_out_flag);
  input_data ["traj_dist_out"]  = Read(Traj::dist_out_flag);
  input_data ["traj_ang_out"]   = Read(Traj::ang_out_flag);
  input_data ["traj_dyn_out"]   = Read(Traj::dyn_out_flag);
  input_data ["traj_anim_out"]  = Read(Traj::anim_out_flag);

  input_data ["traj_base_name"] = Read(Traj::base_name);
  input_data ["traj_debug"]     = Read(Traj::debug);      

  // environment
  int run_time;
  string work_dir;  // slave working directory
  int nice_val;     // nice increment
  string log_file;  // local log file
  input_data ["log_file"]   = Read(log_file);      
  input_data ["work_dir"]      = Read(work_dir);  
  input_data ["run_time"]    = Read(run_time, 3000000);      
  input_data ["save_time"]   = Read(Traj::save_time_step, 300);      
  input_data ["comm_debug"]  = Read(Comm::debug, 0);      
  input_data ["nice"]        = Read(nice_val, 0);      
 
  if (argc < 2) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame
	   << "usage: rfactor input_file\n";
    exit(1);
  }

  string key, comment;

  ifstream from(argv[1]);
  if(!from) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame 
	   << "input file " << argv[1] << " is not found\n";
    exit(1);
  }
  while(from >> key) {
    Inter i = input_data.find(key);
    if (i == input_data.end())
      getline(from, comment);
    else
      from >> i->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  bool is_init = true;
  for(Inter it = input_data.begin(); it != input_data.end(); ++it)
    if(!it->second.is_init()) {
      if(!Comm::rank())
	cout << Comm::mesg() << funame 
	     << it->first << " is not initialized, exitting\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  if(!Comm::rank()) {
    std::cout << "Default parameters:\n";
    for(Inter it = input_data.begin(); it != input_data.end(); ++it)
      if(!it->second.is_default())
	cout << it->first << "  " << it->second << "\n";
    std::cout << "\n";
  }

  // read molecular structure from the file  
  mol_init(mol_spec_file);

  // dividing surfaces
  surf_init(ds_inp_file.c_str());  // dividing surfaces specification
  Traj::tran_surf = *ds_array[ds_ind];  // transition state surface

  if(Traj::reac_surf_flag) {// reactive surface
    from.open(reac_surf_file.c_str());
    if(from)
      from >> Traj::reac_surf;
    else {
      cout << Comm::mesg() << funame 
	   << "WARNING!!! cannot read reactive surface specification,"
	" assuming none\n";
      Traj::reac_surf_flag = 0;
    }
    from.close();
    from.clear();
  }

  if(Traj::reac_ener_flag) {// reactive energies
    if(reac_ener_step > 0.) {
      if(!Comm::rank())
	cout << Comm::mesg() << funame
	     << "reactive energy step should be negative, exitting\n";
      exit(1);
    }
    if(reac_ener_num < 1) {
      if(!Comm::rank())
	cout << Comm::mesg() << funame
	     << "number of reactive energies should be positive, exitting\n";
      exit(1);
    }

    Traj::reac_ener.resize(reac_ener_num);
    for(int i = 0; i < Traj::reac_ener.size(); ++i)
      Traj::reac_ener[i] = reac_ener_max + i * reac_ener_step;
  }

  if(!Traj::reac_surf_flag && !Traj::reac_ener_flag) {
    cout << Comm::mesg() << funame 
	 << "at least one reaction condition mechanism " 
      "should be provided, exitting\n";
    exit(1);
  }

  Traj::diss_surf.write_dist(0) = traj_diss_dist;// dissociative surface
  for(int frag = 0; frag < 2; ++frag)
    for(int i = 0; i < 3; ++i)
      Traj::diss_surf.write_ref_pos(frag, 0)[i] = 0.;

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
  cout << Comm::mesg() << "host = " << hostname 
       << ", process =  " << getpid() << endl;

  MPI::COMM_WORLD.Barrier();

  int traj;
  Traj::Res traj_res;
  if (Comm::rank() == 0) {// master
    /***********************************************************************
     **************************    MASTER    *******************************
     ***********************************************************************/

    // need to block signals while writing
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

    // set step time
    alarm(run_time);
    time_t save_time = std::time(0) + (time_t)Traj::save_time_step;

    cout << "\nStarting calculation ..." << endl;
    
    sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
    from.open(traj_save_file.c_str());
    if(from) {
      int set_size;
      // failed trajectories
      from >> set_size;
      for(int i = 0; i < set_size; ++i) {
	from >> traj;
	fail_set.insert(traj);
      }
      // finished trajectories
      from >> set_size;
      for(int i = 0; i < set_size; ++i) {
	from >> traj >> traj_res;
	res_set[traj] = traj_res;
      }
      // running trajectories
      Traj::State traj_state;
      from >> set_size;
      for(int i = 0; i < set_size; ++i) {
	from >> traj >> traj_state;
	start_set[traj] = traj_state;
      }
      if(!from) {
	cout << Comm::mesg() << funame << "save file is corrupted, exitting\n";
	exit(1);
      }
    }
    from.close();
    from.clear();
    sigprocmask(SIG_SETMASK, &old_sig, 0);

    from.open(wght_smp_file.c_str());
    if (!from) {
      cout << Comm::mesg() << funame << "sampling file " << wght_smp_file 
	   << " is not found, exitting\n";
      exit(1);
    }
    traj = 0;
    Samp samp;
    while(from >> samp) {
      smp_array.push_back(samp);
      // backward trajectories
      if(fail_set.find(traj)  == fail_set.end()  &&
	 start_set.find(traj) == start_set.end() &&
	 res_set.find(traj)   == res_set.end())
	init_set.insert(traj);
      ++traj;
      // forward trajectories
      if(fail_set.find(traj)  == fail_set.end()  &&
	 start_set.find(traj) == start_set.end() &&
	 res_set.find(traj)   == res_set.end()   &&
	 res_set.find(traj - 1) != res_set.end() &&
	 res_set[traj - 1].is_reac())
	init_set.insert(traj);
      ++traj;
    }
    from.close();
    from.clear();

    if(Traj::debug >= Traj::INFO)
      cout << Comm::mesg() << funame 
	   << traj << " samplings were processed"
	   << "\n\t>>> finished trajectories:  " << res_set.size()
	   << "\n\t>>> failed trajectories:  "   << fail_set.size()
	   << "\n\t>>> started trajectories:  "  << start_set.size()
	   << "\n\t>>> new trajectories:  "      << init_set.size()
	   << endl;

    // initial load
    vector<int> node_map(Comm::size(), -1);
    for(int node = 1; node < Comm::size(); ++node) {
      traj = start_new_trajectory(node);
      if(traj < 0)
	break;
      node_map[node] = traj;
    }

    while(run_set.size()) {// receive-send loop

      // save flux data at time_step intervals
      time_t curr_time = std::time(0);
      if(curr_time > save_time) {// save
	save_time = curr_time + (time_t)Traj::save_time_step;
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	save_traj_data(traj_save_file);
	sigprocmask(SIG_SETMASK, &old_sig, 0);
      }// save

      // receive result of trajectory propagation
      int tag, node;
      traj_res.recv(&node, &tag);
      traj = node_map[node];
      if(Traj::debug >= Traj::DEBUG &&
	 (run_set.find(traj)   == run_set.end()   ||
	  start_set.find(traj) != start_set.end() ||
	  res_set.find(traj)   != res_set.end()   ||
	  fail_set.find(traj)  != fail_set.end()  ||
	  init_set.find(traj)  != init_set.end())) {
	cout << Traj::debug_mesg << Comm::mesg() << funame << traj 
	     << "-th trajectory in a wrong set, node = " << node 
	     << ", exitting\n";
	exit(1);
      }

      switch(tag) {// receive
      case Comm::FAIL_TAG:// trajectory calculation failure
	// update control sets
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	fail_set.insert(traj);
	run_set.erase(traj);
	sigprocmask(SIG_SETMASK, &old_sig, 0);

	node_map[node] = -1;
	if(Traj::debug >= Traj::INFO)
	  cout << Comm::mesg() << funame 
	       << traj << "-th trajectory failed on "
	       << node << "-th node" << endl;
	break;
	  
      case Comm::STOP_TAG:// trajectory calculation successfully finished
	// update control sets
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	res_set[traj] = traj_res;
	run_set.erase(traj);
	sigprocmask(SIG_SETMASK, &old_sig, 0);

	if(traj % 2 == 0 && traj_res.is_reac())
	  init_set.insert(traj + 1);

	node_map[node] = -1;
	if(Traj::debug >= Traj::INFO)
	  cout << Comm::mesg() << funame 
	       << traj << "-th trajectory finished on "
	       << node << "-th node" << endl;
	break;
	
      case Comm::RUN_TAG:// current position
	// receive intermediate trajectory info
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
	run_set[traj].value().recv(node, Comm::RUN_TAG);
	sigprocmask(SIG_SETMASK, &old_sig, 0);

	if(Traj::debug >= Traj::DEBUG)
	  cout << Traj::debug_mesg << Comm::mesg() << funame << traj 
	       << "-th trajectory info on "
	       << node << "-th node has been updated" << endl;
	continue;

      default:
	cout << Comm::mesg() << funame << "wrong tag received, exitting\n";
	exit(1);
      }// receive

      // start new trajectory
      if(tag != Comm::FAIL_TAG) {
	traj = start_new_trajectory(node);
	if(traj >= 0)
	  node_map[node] = traj;
      }
    }// receive-send loop

    // save the results
    sigprocmask(SIG_SETMASK, &block_sig, &old_sig); 
    save_traj_data(traj_save_file);
    sigprocmask(SIG_SETMASK, &old_sig, 0);

    // summary results
    vector<int> traj_num_reac(Traj::reac_ener.size(), 0);
    std::vector<std::vector<int> > react_table(Traj::reac_surf.size() + 1, 
					       std::vector<int>(Traj::reac_surf.size()));
    int traj_num_tot = 0;
    int recross = 0;
    for(int smp = 0; smp < smp_array.size(); ++smp) {// sampling cycle
      int back = 2 * smp;
      int forw = back + 1;
      if(fail_set.find(back) != fail_set.end() ||
	 fail_set.find(forw) != fail_set.end())
	continue;
      if(res_set.find(back) == res_set.end()) {
	cout << Comm::mesg() << funame
	     << back << "-th trajectory should be in res_set, exitting\n";
	exit(1);
      }
      ++traj_num_tot;

      if(res_set[back].reac_face() >= 0 && res_set.find(forw) == res_set.end()) {
	cout << Comm::mesg() << funame
	     << forw << "-th trajectory should be in res_set, exitting\n";
	exit(1);
      }

      if(Traj::reac_ener_flag) {// energy level 
	for(int level = 0; level < Traj::reac_ener.size(); ++level)
	  if((res_set[back].ener_level() > level ||  
	      res_set[back].reac_face() >= 0) &&
	     (res_set[forw].ener_level() <= level && 
	      res_set[forw].reac_face() < 0))
	    ++traj_num_reac[level];
      }
      else if(res_set[back].reac_face() < 0) // trajectory recrossed the dividing surface
	++recross;
      else if(res_set[forw].reac_face() < 0) // reactants dissociated
	++react_table[0][res_set[back].reac_face()];
      else // reactants rearranged
	++react_table[res_set[forw].reac_face() + 1][res_set[back].reac_face()];
    }// sampling cycle
    
    cout.precision(4);
    std::cout.setf(ios_base::left, ios_base::adjustfield);
    cout << "\nTotal number of trajectories:  " << traj_num_tot << "\n\n";
    cout << "Trajectories recrossed TS surface: " << recross << "\n\n";
    if(Traj::reac_ener_flag) {// energy levels
      cout << "Reactive Energy (kcal/mol): "; 
      for(int i = 0; i < Traj::reac_ener.size(); ++i)
	cout << setw(7) << Traj::reac_ener[i] / Phys_const::kcal << " ";
      cout   << "\n\nRecrossing Factor:          ";
      for (int i = 0; i < traj_num_reac.size(); ++i)
	cout << setw(7) << static_cast<double>(traj_num_reac[i]) 
	  / static_cast<double>(traj_num_tot) << " ";
      cout << "\n\n";
    }
    else {// reactive surface
      std::cout << "Reaction table:\n\n";
      std::cout << std::setw(4) << "P\\R";
      for(int r = 0; r < Traj::reac_surf.size(); ++r)
	std::cout << std::setw(6) << r + 1;
      std::cout << "\n";
      for(int p = 0; p < react_table.size(); ++p) {// products cycle
	std::cout << std::setw(4);
	if(!p)
	  std::cout << "BI";
	else
	  std::cout << p;
	for(int r = 0; r < Traj::reac_surf.size(); ++r) // reactants cycle
	  std::cout << std::setw(6) << react_table[p][r];
	std::cout << "\n";
      } //products cycle
      std::cout << "\n";
    }// reactive surface
   
    // stop slaves
    for(int node = 1; node < Comm::size(); ++node) {
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, Comm::STOP_TAG);
      if(Comm::debug)
	cout << Comm::mesg() << funame
	     << Comm::STOP_TAG << "-th tag sent to " 
	     << node << "-th node" << endl;
    }

  }// master
  else {// slave
    /***********************************************************************
     **************************    SLAVE   *********************************
     ***********************************************************************/

    // initialize potential environment
    if (pot_type == "g98")
      Gauss::init(pot_file);
    else if(pot_type == "molpro")
      Molpro::init(pot_file);
    else if(pot_type == "multipole")
      Multipole::init(pot_file);
    else if(pot_type == "sjk")
      Sjk::init(pot_file);
      
    // all common files should be open before changing directory!!!
    Comm::change_work_dir(work_dir);

    // set unique wf name
    std::ostringstream oss;
    oss << "n" << Comm::rank() << ".wf";
    Molpro::wf_name = oss.str();

    // set local logging
    Log::out.open(log_file.c_str(), ios::app);

    // the potential function to use
    PES::pot = name2pot(pot_name);

    bool is_run = true;
    Traj::State traj_state;
    while(is_run) {// run loop
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(&traj, 1, MPI::INT, 0, MPI::ANY_TAG, stat);
      if(Comm::debug)
	cout << Comm::mesg() << funame
	     << stat.Get_tag() << "-th tag received from master" << endl;

      switch(stat.Get_tag()) {// receive
      case Comm::STOP_TAG:
	is_run = false;
	break;
	
      case Comm::RUN_TAG:
	traj_state.recv(0, Comm::RUN_TAG);
	// run trajectory
	try {
	  traj_state.run(traj);
	}
	catch(Traj::Error) {
	  static_cast<Traj::Res&>(traj_state).send(0, Comm::FAIL_TAG);
	  break;
	}
	// send trajectory propagation result
	static_cast<Traj::Res&>(traj_state).send(0, Comm::STOP_TAG);
	break;

      default:
	cout <<  Comm::mesg() << funame << "Oops!!! wrong tag, exitting\n";
	exit(1);
      }// receive
    }// run loop
  }// slave

  MPI::COMM_WORLD.Barrier();
  sleep(1);
  return 0;
}// main

