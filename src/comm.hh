#ifndef __COMM__
#define __COMM__

#include <string>
#include <iostream>

#include <mpi.h>

namespace Comm
{
  using namespace std;

  enum tag_t  { RAND_TAG, FLUX_TAG, SURF_TAG, SAMP_TAG, TRAJ_TAG,
		RUN_TAG, STOP_TAG, FAIL_TAG};
  enum node_t { BUSY_STATE, IDLE_STATE, FAIL_STATE };

  // error classes
  class Error {
  public:
    Error (const char* func) { cout << "Comm::" << func; }
  }; 
  class Host_err : public Error {
  public:
    Host_err (const char* func) : Error(func) 
    { cout << ": wrong host" << endl; }
  };
  class Range_err : public Error {
  public:
    Range_err (const char* func) : Error(func) 
    { cout << ": range error" << endl; }
  };

  extern int debug;
  void init (int&, char**&);      // initialize communication environment
  void change_work_dir(string);
  bool is_initialized ();         // is Comm initialized ?
  int rank ();                    // rank of the node
  const string& mesg ();          // node message
  int size ();                    // number of the nodes
  void gather (node_t);           // gather node states
  node_t state (int i);           // state of the i-th node (master)
  void set_state (int i, node_t); // set state of the i-th node (master)
  int find();                     // find available node
}

#endif
