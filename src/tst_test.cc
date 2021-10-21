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
 
  string pot_type, opt_method;
  input_data ["pot_type"]  = Read(pot_type);
  input_data ["opt_method"]  = Read(opt_method, "none");

  string ds_inp, ds_out; // dividing surface array initializer
  input_data ["ds_inp_file"]  = Read(ds_inp, "divsur.inp");
  input_data ["ds_out_file"]  = Read(ds_out, "divsur.out");

  string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");
 
  // specific methods
  std::string gauss_file, molpro_file, multipole_file, sjk_file;
  input_data["gauss_inp_file"] = Read(gauss_file, "gauss.inp");
  input_data["molpro_inp_file"] = Read(molpro_file, "molpro.inp");
  input_data["sjk_inp_file"] = Read(sjk_file, "sjk.inp");
  input_data["multipole_inp_file"] = Read(multipole_file, "multipole.inp");


  if (argc != 2)
    error("usage: tst_test input_file");

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

  // checking for default values
  std::cout << "Default parameters:\n";
  for(Inter it = input_data.begin(); it != input_data.end(); ++it)
    if(it->second.is_default())
      cout << it->first << "  " << it->second << "\n";
  std::cout << "\n";

  std::cout << "Initializing potential ..." << std::endl;
  if (pot_type == "g98" || opt_method == "gauss")
    Gauss::init(gauss_file);
  if(pot_type == "molpro" || opt_method == "molpro")
    Molpro::init(molpro_file);
  if(pot_type == "multipole")
    Multipole::init(multipole_file);
  if(pot_type == "sjk")
    Sjk::init(sjk_file);
  std::cout << "Initializing potential done" << std::endl;


  // read molecular structure from the file  
  std::cout << "Reading molecular structure ... ";
  std::cout.flush();
  mol_init(mol_spec_file);
  std::cout << "done" << std::endl;

  // dividing surfaces specification
  std::cout << "Reading dividing surfaces ... ";
  std::cout.flush();
  surf_init(ds_inp.c_str(), ds_out.c_str());
  std::cout << "done" << std::endl;

  // molecular frame atomic coordinates
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

  double rm = 0.;
  for(int frag = 0; frag < 2; ++frag) {
    cout << "Fragment " << frag << ": \n";
    for(int at = 0; at < mol_array[frag]->size(); ++at) {
      cout << setw(2) << mol_array[frag]->begin()[at].name() << at+1;
      for(int i = 0; i < 3; ++i)
	cout << setw(6) <<  mol_array[frag]->begin()[at].mf_pos[i];
      cout << "\n";
    }
    cout << "Fragment mass: " << mol_array[frag]->mass() << endl;
    cout << "Inertia moments: "; 
    for(int i = 0; i < 3; ++i)
      cout << "  " << mol_array[frag]->iner_mom(i);
    cout << endl;
    rm += 1. / mol_array[frag]->mass();
  }
  cout << "Reduced mass = " << 1./rm << endl;
  return 0;
}
