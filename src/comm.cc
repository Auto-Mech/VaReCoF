#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
#include <cerrno>
#include <unistd.h>
#include <cstdio>
#include <cstring>

#include "comm.hh"
#include "error.hh"

namespace Comm {

  int debug = 0;
  bool is_init = false;
  int _rank;
  int _size;
  int* nstat; // states of the nodes
  
  void cleanup ();
}

bool Comm::is_initialized ()
{
  return is_init;
}

void Comm::init (int& argc, char**& argv)
{
  const char funame [] = "Comm::init: ";

  if (is_init) {
    cout << funame << "has been allready initialized" << endl;
    return;
  }
  is_init = true;

  // mpi initialization
  //
  MPI::Init(argc, argv);
  
  _rank = MPI::COMM_WORLD.Get_rank();
  
  _size = MPI::COMM_WORLD.Get_size();

  // cleanup on exit
  //
  if (atexit(&cleanup))
    //
    std::cerr << funame << "failed to register at-exit clean-up procedure" << endl;
 
  // node states array
  //
  if (!rank())
    //
    nstat = new int[size()];
}

void Comm::change_work_dir(string work_dir)
{// set slave working directory
  const char funame [] = "Comm::change_work_dir: ";

  if (!rank())
    return;

  char my_rank [7];
  sprintf(my_rank, "%i", rank());
  work_dir += my_rank;
  struct stat wstat;
  if(mkdir(work_dir.c_str(), S_IRWXU | S_IRGRP | 
	   S_IXGRP | S_IROTH | S_IXOTH )) {
    if(errno != EEXIST)
      cout << funame << strerror(errno) << endl;
    switch(errno) {
    case EEXIST:
      if(stat(work_dir.c_str(), &wstat)) {
	cout << funame << strerror(errno) << endl;
	exit(1);
      }
      if(!S_ISDIR(wstat.st_mode)) {
	cout << funame << "file " << work_dir <<" is not a directory\n";
	exit(1);
      }
      if(!(wstat.st_mode & S_IRUSR) || !(wstat.st_mode & S_IWUSR) ||
	 !(wstat.st_mode & S_IXUSR)) {
	cout << funame << "directory " << work_dir 
	     << " has wrong permissions\n";
	exit(1);
      }
      break;
    default:
      exit(1);
    }   
  }
  if(chdir(work_dir.c_str())) {
    cout << funame << strerror(errno) << endl;
    exit(1);
  }
}

void Comm::gather (node_t ns)
{
  int itemp = ns;
  MPI::COMM_WORLD.Gather(&itemp, 1, MPI::INT, nstat, 1, MPI::INT, 0);
}

int Comm::rank () 
{ 
  return _rank;
}

const string& Comm::mesg ()
{
  static string _mesg;

  static bool is_first = true;
  if(is_init && is_first) {
    is_first = false;
    if(rank()) {
    char c_str [33];
    sprintf(c_str, "node %i: ", rank());
    _mesg = c_str;
    }
    else
      _mesg = "master: ";
  }

  return _mesg;
}

int Comm::size ()
{
  return _size;
}

Comm::node_t Comm::state (int i)
{
  const char func [] = "state";
 
  if (rank())
    throw Host_err(func);

  if (i <= 0 || i >= size())
    throw Range_err(func);

  return (node_t) nstat[i];
}

void Comm::set_state (int i, node_t ns)
{
  const char func [] = "set_state";
 
  if (rank())
    throw Host_err(func);

  if (i <= 0 || i >= size())
    throw Range_err(func);

  nstat[i] = ns;
}

int Comm::find ()
{
  const char func [] = "find";
 
  if (rank())
    throw Host_err(func);
  for (int i = 1; i < size(); ++i)
    if(nstat[i] == IDLE_STATE)
      return i;
  return 0;
}
  
void Comm::cleanup ()
{
  if (Comm::rank() == 0)
    cout << "Comm::cleanup: MPI round up ..." << endl;
  MPI::Finalize();
}
