#include <string>
#include <iostream>
#include <cstdlib>
#include "error.hh"

/***************** Error handling ***********************/
void error (const string& message, int exit_code)
{
  cout << message << endl;
  exit (exit_code);
}

