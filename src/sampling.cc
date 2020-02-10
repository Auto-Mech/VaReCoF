#include<vector>
#include<string>
#include<fstream>
#include<map>
#include<iomanip>

#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <cerrno>
#include <sys/stat.h>
#include <cstdio>
#include <cstring>
#include <csignal>

#include <mpi.h>

#include "input.hh"
#include "rotd.hh"    
#include "force.hh"
#include "gauss.hh"
#include "raninit.hh"
#include "random.hh"
#include "error.hh"
#include "comm.hh"
#include "flux.hh"
#include "multipole.hh"
#include "sjk.hh"
#include "molpro.hh"
#include "log.hh"

vector<Samp> smp_array;
double temperature;

void update_smp_array (const Samp& new_samp)
{
  static int smp_num_accept = 0;
  static double max_weight = 0.;       // maximum weight

  vector<Samp>::iterator samp_iter; 
  if(new_samp.weight(temperature) > max_weight) {
    max_weight = new_samp.weight();
    samp_iter = smp_array.begin();
    while(samp_iter != smp_array.end())
      if(samp_iter->random_value() * max_weight > samp_iter->weight())
	samp_iter = smp_array.erase(samp_iter);
      else
	++samp_iter;

    smp_array.push_back(new_samp);
  }
  else if(new_samp.random_value() * max_weight < new_samp.weight())
    smp_array.push_back(new_samp);

  cout << "sampling = "     << ++smp_num_accept 
       << "  traj. num. = " << smp_array.size()
       << "  max. weight = " << max_weight << "\n";
}
     
int main (int argc, char* argv[])
{// main

  const char funame [] = "main: ";

  // temporary variables
  string stemp;

  /*************************************************************
   *********************  MPI STARTS HERE!  ********************
   *************************************************************/
  Comm::init(argc, argv); // working directory changed here!

  // Random numbers initialization
  Random::init();

  /*************** Input of parameters ********************/

  map <string, Read> input_data;
  typedef map<string, Read>::iterator Inter;

  // The potential 
  string pot_name, pot_type, pot_file;
  input_data ["pot_name"]  = Read(pot_name);
  input_data ["pot_type"]  = Read(pot_type);
  input_data ["pot_file"]  = Read(pot_file);

  string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");

  // sampling specification
  int face, ds_ind;
  string ds_inp_file, smp_type; // dividing surface array initializer
  input_data ["sampling"]  = Read(smp_type, "multifacet");
  input_data ["ds_inp_file"]  = Read(ds_inp_file, "divsur.inp");
  input_data ["face"] = Read(face);
  input_data ["ds_ind"] = Read(ds_ind);

  string raw_smp_file, wght_smp_file;
  input_data ["is_smp_out"]  = Read(Flux::smp_out_flag);
  input_data ["raw_smp_file"]  = Read(raw_smp_file, "sampling.raw");
  input_data ["wght_smp_file"]  = Read(wght_smp_file, "sampling.wght");

  int smp_num_max; // target weighted samplings number
  input_data ["smp_num_max"] = Read(smp_num_max);
  input_data ["temperature"] = Read(temperature);

  // minimal distance between atoms
  input_data ["min_atom_dist"]  = Read(MIN_ATOM_DIST, 1.5);

  string work_dir; // slave working directory
  int nice_val; // nice increment
  string log_file; // local log file
  input_data ["log_file"]    = Read(log_file);      
  input_data ["work_dir"]    = Read(work_dir);  
  input_data ["comm_debug"]  = Read(Comm::debug, 0);      
  input_data ["nice"]        = Read(nice_val, 0);

  if (argc < 2) {
    if(!Comm::rank())
      cout << Comm::mesg() << funame
	   << "usage: sampling input_file\n";
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

  // dividing surfaces specification
  surf_init(ds_inp_file.c_str());

  // sampling function specification
  name2ranp(smp_type);

  // set the process' file creation mask to 022
  umask(S_IWGRP | S_IWOTH);

  // setting nice
  if(!Comm::rank() && nice_val < 0)
    cout << Comm::mesg() << funame 
	 << "WARNING!!! negative nice increment, ignoring" 
	 << endl;
  if(Comm::rank() && nice_val > 0)
    nice(nice_val);

  MPI::COMM_WORLD.Barrier();

  // list nodes
  char hostname [100];
  gethostname(hostname, 99);
  cout << Comm::mesg() << funame
       << "host = " << hostname 
       << ", process =  " << getpid() << endl;

  MPI::COMM_WORLD.Barrier();

  Samp samp;
  if (Comm::rank() == 0) {// master
    /***********************************************************************
     **************************    MASTER    *******************************
     ***********************************************************************/

    sigset_t block_sig;
    sigset_t old_sig;
    sigfillset(&block_sig); // block all signals

    // Surface id
    const Sid tran_sid(ds_ind, face);

    if(Flux::smp_out_flag) {
      from.open(raw_smp_file.c_str());
      if (from) {
	cout << "\nReading raw samplings from " 
	     << raw_smp_file << " file:\n\n";
	Sid curr_sid;
	while (from >> curr_sid >> samp) {
	  if(tran_sid == curr_sid)
	    update_smp_array(samp);
	  if(smp_array.size() >= smp_num_max)
	    break;
	}
      }
      from.close();
      from.clear();
    }

    ofstream smp_out;
    if(Flux::smp_out_flag)
      smp_out.open(raw_smp_file.c_str(), ios::app);

    int smp_num_fail = 0;
    int run_node_num = 0;
    if(smp_array.size() < smp_num_max) {
      cout << "\nStarting new samplings:\n\n";
      for(int node = 1; node < Comm::size(); ++node) {
	MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, Comm::SAMP_TAG);
	if(Comm::debug)
	  cout << Comm::mesg() << funame 
	       << Comm::SAMP_TAG
	       << "-th tag sent to "
	       << node << "-th node" << endl;
	++run_node_num;
      }
    }
    while(run_node_num) {// sampling cycle

      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, 
			   MPI::ANY_TAG, stat);
      int node = stat.Get_source();
      if(Comm::debug)
	cout << Comm::mesg() << funame 
	     << stat.Get_tag()
	     << "-th tag received from "
	     << node << "-th node" << endl;

      switch(stat.Get_tag()) {// tag menu

      case Comm::RAND_TAG:// Random Number Generator Seed
	Random::send_seed(stat.Get_source());
	continue;

      case Comm::FAIL_TAG:// Node Failure
	cout << "master: node " << stat.Get_source()
	     << " failed\n" << endl;

	++smp_num_fail;
	--run_node_num;
	break;

      case Comm::SAMP_TAG:// Flux Sampling

	samp.recv(node, Comm::SAMP_TAG);
	if(Flux::smp_out_flag) {// save sampling
	  sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
	  smp_out << tran_sid << samp;
	  sigprocmask(SIG_SETMASK, &old_sig, 0);
	}
	update_smp_array(samp);
	--run_node_num;
	break;
      
      default:
	cout << Comm::mesg() << funame << "wrong tag, exitting\n";
	exit(1);

      }// tag menu

      if(smp_array.size() < smp_num_max) {// request new sampling
	MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, Comm::SAMP_TAG);
	if(Comm::debug)
	  cout << Comm::mesg() << funame 
	       << Comm::SAMP_TAG
	       << "-th tag sent to "
	       << node << "-th node" << endl;
	++run_node_num;
      }

    }// sampling cycle

    if(Flux::smp_out_flag) {
      smp_out.close();
      smp_out.clear();
    }

    // sampling output
    const double temper2 = sqrt(temperature);
    smp_out.open(wght_smp_file.c_str());
    for (int smp = 0; smp < smp_array.size(); ++smp) {
      smp_array[smp].attach();
      for(int frag = 0; frag < 2; ++frag) {
	for(int i = 0; i < 3; ++i)
	  mol_array[frag]->write_cm_vel()[i] *= temper2;
	for(int i = 0; i < mol_array[frag]->ang_vel_size(); ++i)
	  mol_array[frag]->write_ang_vel()[i] *= temper2;
      }
      sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
      smp_out << smp_array[smp];
      sigprocmask(SIG_SETMASK, &old_sig, 0);
    }
    smp_out.close();
    smp_out.clear();

    cout << "\nnecessary number of trajectories to run has been reached\n"
	 << "\nNumber of failed samplings = " 
	 << smp_num_fail << endl;

    for(int node = 1; node < Comm::size(); ++node) {
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, Comm::STOP_TAG);
      if(Comm::debug)
	cout << Comm::mesg() << funame 
	     << Comm::STOP_TAG
	     << "-th tag sent to "
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

    // set local logging
    Log::out.open(log_file.c_str(), ios::app);

    // the potential function to use
    PES::pot = name2pot(pot_name);

    // transition state surface
    Div_surf tran_surf = *ds_array[ds_ind];
    tran_surf.set_face(face);

    bool is_run = true;
    while(is_run) {// run

      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, 0, MPI::ANY_TAG, stat);
      if(Comm::debug)
	cout << Comm::mesg() << funame 
	     << stat.Get_tag()
	     << "-th tag received from master" << endl;

      double energy, weight;
      switch(stat.Get_tag()) {// tag menu

      case Comm::STOP_TAG:// stop calculation
	cout << Comm::mesg() << funame 
	     << "finishing calculation...\n";
	is_run = false;
	break;

      case Comm::SAMP_TAG:// sample SurFace
	while(rand_pos(tran_surf, samp, &weight)) {}
	try {
	  Array<double> ener_arr(1);
	  PES::pot(0, ener_arr);
	  energy = ener_arr[0];
	}
	catch (Pot_error perr) {
	  MPI::COMM_WORLD.Send(0, 0, MPI::INT, 0, Comm::FAIL_TAG);
	  if(Comm::debug)
	    cout << Comm::mesg() << funame 
		 << Comm::FAIL_TAG
		 << "-th tag sent to master" << endl;
	  break;
	}
	samp.set_val(weight, energy);

	MPI::COMM_WORLD.Send(0, 0, MPI::INT, 0, Comm::SAMP_TAG);
	if(Comm::debug)
	  cout << Comm::mesg() << funame 
	       << Comm::SAMP_TAG
	       << "-th tag sent to master" << endl;
	samp.send(0, Comm::SAMP_TAG);
	break;

      default:
	cout << Comm::mesg() << funame << "wrong tag, exitting\n";
	exit(1);
      }// tag menu
    }// run
  }// slave
  cout.flush();
  MPI::COMM_WORLD.Barrier();
  sleep(1);
  return 0;
}
