#include <cmath>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>

#include "random.hh"
#include "error.hh"
#include "comm.hh"

#include <mpi.h>

#include "random_common.cc"

namespace Random {

  const int count_max = 100000;
  bool is_init = false;
  unsigned short server_seed [3];
  unsigned short seed [3];
  MPI::Request rqst;

  void generate (); // generate random seed (master)
  void request ();  // request random seed (slave)
  void update ();   // update random seed
  void cleanup ();  // finish


}

void Random::init () 
{
  const char funame [] = "Random::init: ";

  if(is_init) {
    cout << Comm::mesg() << funame 
	 << "has been allready initialized" << endl;
    return;
  }
  is_init = true;

  if(!Comm::is_initialized()) {
    cout << Comm::mesg() << funame
	 << "communication environment has not been initialized, exitting\n";
    exit(1);
  }

  // cleanup on exit
  if(atexit(&cleanup))
    cout << Comm::mesg() << funame
	 << "failed to register cleaning procedure" 
	 << endl;

  if(Comm::rank() == 0) { // master
    // initialize
    server_seed[0] = getpid();
    server_seed[1] = time(0);
    server_seed[2] = time(0) << 16;
    
    for (int i = 0; i < 3; ++i)
      seed[i] = server_seed[i];
    generate();
    
    for (int rank = 1; rank < Comm::size(); ++rank)
      send_seed(rank);
  }
  else { // slave
    MPI::COMM_WORLD.Recv(server_seed, 3, MPI::UNSIGNED_SHORT, 
			 0, Comm::RAND_TAG);
      if(Comm::debug)
	cout << Comm::mesg() << funame
	     << Comm::RAND_TAG << "-th tag received from master" 
	     << endl;
    
    for (int i = 0; i < 3; ++i)
      seed[i] = server_seed[i];
    
    request();	
  }
}

void Random::cleanup ()
{
  //const char funame [] = "Random::cleanup: ";

  //  if (Comm::rank() == 0)
  // cout << Comm::mesg() << funame 
  // << "deactivating random generator ..." 
  // << endl;

  //  if (Comm::rank() != 0)
  // rqst.Wait();
}

void Random::generate ()
{
  const char funame [] = "Random::generate: ";

  if (Comm::rank() == 0)
    for (int i = 0; i < count_max; ++i)
      nrand48(server_seed);
  else {
    cout << Comm::mesg() << funame 
	 << "wrong host, exitting\n";
    exit(1);
  }
}

void Random::request ()
{
  const char funame [] = "Random::request: ";

  if (Comm::rank() != 0) {
    MPI::COMM_WORLD.Send(0, 0, MPI::CHAR, 0, Comm::RAND_TAG);
    if(Comm::debug)
      cout << Comm::mesg() << funame
	   << Comm::RAND_TAG << "-th tag sent to master" 
	   << endl;
    
    rqst = MPI::COMM_WORLD.Irecv(server_seed, 3, MPI::UNSIGNED_SHORT, 
				 0, Comm::RAND_TAG);
  }
  else {
    cout << Comm::mesg() << funame 
	 << "wrong host, exitting\n";
    exit(1);
  }
}

void Random::update ()
{
  const char funame [] = "Random::update: ";

  if (Comm::rank() == 0) {
    for (int i = 0; i < 3; ++i)
      seed[i] = server_seed[i];
    generate();
  }
  else {
    rqst.Wait();
    if(Comm::debug)
      cout << Comm::mesg() << funame
	   << Comm::RAND_TAG << "-th tag received from master" 
	   << endl;
    for (int i = 0; i < 3; ++i)
      seed[i] = server_seed[i];
    request();
  }
}
  
void Random::send_seed (int rank)
{
  const char funame [] = "Random::send_seed: ";

  if (!is_init) {
    cout << Comm::mesg() << funame
	 << "not initialized, exitting\n";
    exit(1);
  }

  if (rank == 0) {
      cout << Comm::mesg() << funame
	   << "trying to send a seed to the server, exitting\n";
      exit(1);
  }

  if (Comm::rank() == 0) {
    MPI::COMM_WORLD.Send(server_seed, 3, MPI::UNSIGNED_SHORT, 
			 rank, Comm::RAND_TAG);
    if(Comm::debug)
      cout << Comm::mesg() << funame
	   << Comm::RAND_TAG << "-th tag sent to "
	   << rank << "-th node" << endl;
    
    generate();
  }
  else {
    cout << Comm::mesg() << funame 
	 << "wrong host, exitting\n";
    exit(1);
  }
}

double Random::flat ()
{
  const char funame [] = "Random::flat: ";

  static int count = 0;

  if (!is_init) {
    cout << Comm::mesg() << funame
	 << "not initialized, exitting\n";
    exit(1);
  }

  if (++count > count_max) {
    update();
    count = 0;
  }

  return erand48(seed);
}

