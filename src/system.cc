#include <iostream>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/types.h>
#include <wait.h>
#include <cstdarg>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdlib>

#include <string>
#include <vector>

#include "system.hh"
#include "error.hh"
#include "tmatrix.hh"

// remove files from directory

#include <ftw.h>

int rm_file (const char* fname, const struct stat* fstat, int flag)
{
  remove(fname);
  return 0;
}

void clean_dir (const char* dname)
{
  ftw(dname, rm_file, 100); // 100 descriptors is allowed
}

// copy files
int file_copy(const char* old_fname, const char* new_fname)
{
  ifstream from(old_fname);
  if(!from)
    return 1;
  ofstream   to(new_fname);
  if (!to)
    return 2;

  char ch;
  while(from.get(ch))
    if (!to.put(ch))
      return -1;
  return 0;
}

/************************************************************************
 **************         External call      ******************************
 ************************************************************************/

int call_exe (const char* exename ...) 
// list of arguments should be terminated by 0 pointer
{
  char* argv [7]; //maximum number of parameters = 5
  Array<char> xn(strlen(exename) + 1);

  va_list ap;
  va_start(ap, exename);
  char** argp = argv;
  *argp++ = basename(strcpy(xn.data(), exename));
  while (*argp++ = va_arg(ap, char*)) {}
  va_end(ap);

  cout.flush();
  int fd, exe_stat;
  switch(fork()) {      
  case -1: //error

    cout << "call_exe: fork failed\n";
    return 1000;

  case 0: //child

    // redirect standard output
    fd = open("std.out", O_WRONLY | O_CREAT | O_TRUNC, 
	      S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
    dup2(fd, 1);
    close(fd);
    // redirect standard error
    fd = open("std.err", O_WRONLY | O_CREAT | O_TRUNC,
	      S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
    dup2(fd, 2);
    close(fd);
    // substitute another image
    exe_stat = execvp(exename, argv);
    std::cout << "running <" << exename << "> with arguments ";
    argp = argv;
    while(*argp)
      std::cout << "<" << *argp++ << ">";
    std::cout << "failed with the exit status " << exe_stat << std::endl;
    exit(1);

  default: //parent

    int child_stat;
    wait(&child_stat);
    return WEXITSTATUS(child_stat);

  }
}

/*************************************************************************
 *                      Pipes
 *************************************************************************/

Pipe_base::Pipe_base()
{
  if (pipe(fd)) //create a pipe
    switch (errno)
      {
      case EMFILE:
	error("pipe: too many file descriptors are in use");
      case ENFILE:
	error("pipe: system file table is full");
      case EFAULT:
	error("pipe: file descriptor is not valid");
      default:
	error("pipe: something wrong");
      }
}
// it seems that a new gcc standard does not support attaching streams to file descriptors
//Pipe::Pipe () : pin(fd[0]), pout(fd[1])
Pipe::Pipe ()
{
  pin.precision (14);
  pout.precision(14);
}

/*********************************************************
 *                   Semaphores 
 *********************************************************/
union semun
{
  int val;                          // value for SETVAL
  struct semid_ds *buf;             // buffer for IPC_STAT & IPC_SET
  unsigned short int *array;        // array for GETALL & SETALL
  struct seminfo *__buf;            // buffer for IPC_INFO
};

Semaphore::Semaphore(int n)
   : key(IPC_PRIVATE), num(n), creator(getpid())
{
   // create
   while((id = semget(++key, num, 0666 | IPC_CREAT | IPC_EXCL)) == -1)
   {} 

   // initialize
   semun su;
   Array<unsigned short> init_val (n);
   for (int i = 0; i < n; ++i)
      init_val[i] = 1;
   su.array = init_val.data();

   if (semctl(id, 0, SETALL, su) == -1)
      {
	 semctl(id, 0, IPC_RMID);
	 error("Semaphore::Semaphore: couldn't initialize");
      }
}

Semaphore::Semaphore(key_t k, int n)
   : key(k), num(n), creator(0)
{
   if ((id = semget(k, n, 0666)) == -1)
      error("semaphore::semaphore: couldn't open existing semaphore set");
}

Semaphore::~Semaphore()
{
   if (getpid() == creator)
      semctl(id, 0, IPC_RMID);
} 

void Semaphore::busy (int n) const
{
   if (n >= num)
      error("Semaphore::busy: wrong semaphore number");

   struct sembuf sb;
   sb.sem_num = n;
   sb.sem_op = -1;
   sb.sem_flg = SEM_UNDO;

   if (semop(id, &sb, 1) == -1)
      error("Semaphore::busy: failed");
}

void Semaphore::free (int n) const
{
   if (n >= num)
      error("Semaphore::busy: wrong semaphore number");

   struct sembuf sb;
   sb.sem_num = n;
   sb.sem_op = 1;
   sb.sem_flg = SEM_UNDO;

   if (semop(id, &sb, 1) == -1)
      error("Semaphore::free: failed");
}
