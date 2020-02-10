#ifndef ROTD_SYSTEM_HH
#define ROTD_SYSTEM_HH

#include <sys/types.h>
#include <fstream>
#include <iostream>

using namespace std;

// delete all files in the directory
void clean_dir (const char*);

// file copy
int file_copy(const char*, const char*);

/**********************************************************************
 * call to external executable (substitute for system()); 
 * list of arguments in call_exe() should be terminated by 0 pointer
 **********************************************************************/

int call_exe (const char* ...);

/*************************************************************************
 *                      Pipes
 *************************************************************************/

class Pipe_base
{
  int fd [2];

  Pipe_base (const Pipe_base&);
  Pipe_base& operator= (const Pipe_base&);

  Pipe_base ();

  friend class Pipe;
};

class Pipe : public Pipe_base
{
  Pipe (const Pipe&);
  Pipe& operator= (const Pipe&);

public:

  ifstream pin;
  ofstream pout;
  Pipe();
};

/*************************** Semaphores **********************************/

class Semaphore
{
  key_t key;   // semaphore key
  int   id;    // semaphore id
  int   num;   // number of semaphores in the set
  pid_t creator; // creator of semaphore

  Semaphore(const Semaphore&);
  Semaphore& operator= (const Semaphore&);

public:

  explicit Semaphore (int); // creates set of n semaphores & initilizes them 
  Semaphore (key_t, int); // initializes existing semaphore set
  ~Semaphore();
  
  key_t get_key () const {return key;}
  void busy (int) const; // raise n-th semaphore (P, wait)
  void free (int) const; // free n-th semaphore  (V, signal)
};

#endif

