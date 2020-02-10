#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <iomanip>
#include <cstdlib>

#include "rotd.hh"    
#include "divsur.hh"
#include "gauss.hh"
#include "error.hh"
#include "multipole.hh"
#include "molpro.hh"
#include "sjk.hh"
#include "input.hh"

int main (int argc, char* argv[])
{// main

  // Parameters input
  map <string, Read> input_data;
  typedef map<string, Read>::iterator Inter;
 
  string ds_inp, ds_out; // dividing surface array initializer
  input_data ["ds_inp_file"]  = Read(ds_inp, "divsur.inp");
  input_data ["ds_out_file"]  = Read(ds_out, "divsur.out");

  string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");
  
  ifstream from(argv[1]);
  if (!from)
    error("main: input file is not found");

  string key, comment;
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
      std::cout << it->first << " is not initialized, exitting\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // read molecular structure from the file  
  // dividing surfaces specification
  mol_init(mol_spec_file);
  surf_init(ds_inp.c_str(), ds_out.c_str());

  // molecular frame atomic coordinates
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

  double rm = 0.;
  for(int frag = 0; frag < 2; ++frag) {
    cout << "Fragment " << frag << ": \n";
    for(int at = 0; at < mol_array[frag]->size(); ++at) {
      cout << setw(2) << mol_array[frag]->begin()[at].name();
      for(int i = 0; i < 3; ++i)
	cout << setw(6) <<  mol_array[frag]->begin()[at].mf_pos[i];
      cout << "\n";
    }
  }
  return 0;
}
