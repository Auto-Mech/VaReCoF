#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "input.hh"
#include "units.hh"
#include "rotd.hh"
#include "math.hh"
#include "tmatrix.hh"
#include "slatec.hh"
#include "error.hh"
#include "divsur.hh"
#include "pes.hh"

int main (int argc, char* argv[])
{
  const char funame [] = "cut_multi: ";

  if (argc != 2)
    error("usage: cut_multi input_file");

  // temporary variables
  int itemp;
  double dtemp, dtemp1;
  std::string key, comment, token, line;

  // input
  const InputKey multi_key = "MultiInputFile";
  std::string multi_input_file;

  const InputKey ener_key  = "Energy[kcal/mol]";
  double ener;

  const InputKey amom_key  = "AngularMomentum[au]";
  double amom;

  const InputKey temp_key  = "Temperature[K]";
  double tempr;

  const InputKey  out_key  = "OutputFile";
  std::string out_file;

  // reading from input file
  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
    exit(1);
  }
  while(from >> token) {
    if(token == multi_key) {
      from >> multi_input_file;
      std::getline(from, comment);
      if(!from) {
	std::cout << funame << token << " reading failed\n";
	exit(1);
      }
    }
    else if(token == ener_key) {
      from >> ener;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
      ener *=Phys_const::kcal;
    }
    else if(token == amom_key) {
      from >> amom;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
    }
    else if(token == temp_key) {
      from >> tempr;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
      tempr *=Phys_const::kelv;
    }
    else if(token == out_key) {// output file name
      from >> out_file;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
    }
    else {
      std::cout << funame << "unknown keyword " << token << ". Known keywords are:\n";
      InputKey::show_all();
      exit(1);
    }
  }
  from.close();
  from.clear();

  // reading multi input file
  std::map <string, Read> input_data;
  typedef map<string, Read>::iterator input_data_it;
 
  std::vector<double> temperature;
  input_data ["tmpr_grid"] = Read(temperature);
  std::vector<double> ener_grid;
  input_data ["ener_grid"] = Read(ener_grid);
  std::vector<double> amom_grid;
  input_data ["amom_grid"] = Read(amom_grid);

  string flux_file; // flux output file
  input_data ["flux_file"]  = Read(flux_file, "flux.dat");      
  string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");

  int pes_size;
  input_data ["pes_size"] = Read(pes_size, 1);

  from.open(multi_input_file.c_str());
  if(!from) {
    std::cout << funame << "cannot open " << multi_input_file << " file\n";
    exit(1);
  }
  while(from >> key) {
    input_data_it it = input_data.find(key);
    if (it == input_data.end())
      getline(from, comment);
    else {
      std::cout << key << "\n";
      from >> it->second;
    }
  }
  from.close();
  from.clear();
  
  // check if all parameters were initialized
  bool is_init = true;
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(!it->second.is_init()) {
      std::cout << funame << it->first << " is not initialized\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  std::cout << "Default parameters:\n";
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(it->second.is_default())
      std::cout << it->first << "  " << it->second << "\n";
  std::cout << "\n";

  // read molecular structure from the file  
  mol_init(mol_spec_file);


  // number of degrees of freedom with cm fixed
  int dof_num = 3;
  for (int i = 0; i < mol_array.size(); ++i)
    switch(mol_array[i]->type())
      {
      case ATOM:
	break;
      case LINEAR:
	dof_num += 2;
	break;
      case NONLINEAR:
	dof_num += 3;
      }

  // indices
  int ener_index = -1;
  for(int i = 0; i < ener_grid.size(); ++i)
    if(ener < ener_grid[i]) {
      ener_index = i;
      break;
    }
  if(ener_index < 0)
    ener_index = ener_grid.size() - 1;

  int temp_index = -1;
  for(int i = 0; i < temperature.size(); ++i)
    if(tempr < temperature[i]) {
      temp_index = i;
      break;
    }
  if(temp_index < 0)
    temp_index = temperature.size() - 1;

  int amom_index = -1;
  for(int i = 0; i < amom_grid.size(); ++i)
    if(amom < amom_grid[i]) {
      amom_index = i;
      break;
    }
  if(amom_index < 0)
    amom_index = amom_grid.size() - 1;

  std::map<int, std::vector<double> > flux, e_flux, ej_flux;

  int ds, ds_size;
  Dynvar temp_dv;

  ifstream flux_in(flux_file.c_str());

  std::cout << "input file read\n";
  while(flux_in >> ds >> ds_size) {//ds cycle
    flux[ds].resize(pes_size, 0.);
    e_flux[ds].resize(pes_size, 0.);
    ej_flux[ds].resize(pes_size, 0.);
    /*
    for(int pes = 0; pes < pes_size; ++pes) {
      flux   [ds][pes] = 0.;
      e_flux [ds][pes] = 0.;
      ej_flux[ds][pes] = 0.;

    }
    */
    for(int face = 0; face < ds_size; ++face) {//face cycle
      for(int i = 0; i < 4; ++i)// sampling numbers
	flux_in >> itemp;
      getline(flux_in, comment);

      for(int pes = 0; pes < pes_size; ++pes)// minimal energies
	flux_in >> dtemp;

      for(int pes = 0; pes < pes_size; ++pes)// minimal configurations
	flux_in >> temp_dv;

      for(int pes = 0; pes < pes_size; ++pes)
	for(int i = 0; i < temperature.size(); ++i) {// thermal flux
	  flux_in >> dtemp >> dtemp1;
	  if(i == temp_index)
	    flux[ds][pes] += dtemp;
	}

      for(int pes = 0; pes < pes_size; ++pes)
	for(int en = 0; en < ener_grid.size(); ++en) {// microcannonical flux
	  flux_in >> dtemp >> dtemp1;
	  if(en == ener_index)
	    e_flux[ds][pes] += dtemp;
	}

      for(int pes = 0; pes < pes_size; ++pes)
	for(int am = 0; am < amom_grid.size(); ++am)// ej-resolved flux
	  for(int en = 0; en < ener_grid.size(); ++en) {// ej-resolved flux
	    flux_in >> dtemp >> dtemp1;
	    if(en == ener_index && am == amom_index)
	      ej_flux[ds][pes] += dtemp;
	  }
    }// face cycle
    if(!flux_in)
      error("main: flux file is corrupted");
  }//ds cycle

  std::cout << "flux read\n";

  // normalization factors
  double mc_fac = 2. * sqrt(2.*M_PI);
  dtemp =sqrt(mol_array[0]->mass() * mol_array[1]->mass()
	      /(mol_array[0]->mass() + mol_array[1]->mass()));
  mc_fac *= dtemp*dtemp*dtemp;
  mc_fac *= mol_array[0]->stat_sum() * mol_array[1]->stat_sum();
  if(ener_index >= 0) {
    /*
    for(ds = 0; ds < flux.size(); ++ds)
      for(int pes = 0; pes < pes_size; ++pes)
	e_flux[ds][pes] *= mc_fac;
    */
  }

  double ej_fac = mc_fac / 2.;
  if(amom_index >= 0) {
    ej_fac *= amom_grid[amom_index];

    /*
    for(ds = 0; ds < flux.size(); ++ds)
      for(int pes = 0; pes < pes_size; ++pes)
	ej_flux[ds][pes] *= ej_fac;
    */
  }

  // output
  ofstream main_out (out_file.c_str());
  if(!main_out) {
    std::cout << funame << "cannot open " << out_file << "\n";
    exit(1);
  }
  main_out.precision(4);

  if(temp_index >= 0)
    main_out << "grid-point temperature = " 
	     << temperature[temp_index] / Phys_const::kelv 
	     << " K\n";
  if(ener_index >= 0)
    main_out << "grid-point energy = " 
	     << ener_grid[ener_index] / Phys_const::kcal 
	     << " kcal/mol\n";
  if(amom_index >= 0)
    main_out << "grid-point angular momentum = " 
	     << amom_grid[amom_index]
	     << " a.u.\n\n";

  main_out << "DS#       PES =  1";
  for(int pes = 1; pes < pes_size; ++pes)
    main_out << std::setw(15) << pes + 1;
  main_out << "\n\n";
  if(temp_index >= 0) {
    main_out << "canonical flux:\n";
    for(std::map<int, std::vector<double> >::iterator it = flux.begin(); it != flux.end(); ++it) {
      main_out << std::setw(3) << it->first;
      for(int pes = 0; pes < pes_size; ++pes)
	main_out << std::setw(15) << it->second[pes];
      main_out << "\n";
    }
    main_out << "\n";
  }
  if(ener_index >= 0) {
    main_out << "microcanonical flux:\n";
    for(std::map<int, std::vector<double> >::iterator it = e_flux.begin(); it != e_flux.end(); ++it) {
      main_out << std::setw(3) << it->first;
      for(int pes = 0; pes < pes_size; ++pes)
	main_out << std::setw(15) << it->second[pes] * mc_fac;
      main_out << "\n";
    }
    main_out << "\n";
  }

  if(ener_index >= 0 && amom_index >= 0) {
    main_out << "E,J-resolved flux:\n";
    for(std::map<int, std::vector<double> >::iterator it = ej_flux.begin(); it != ej_flux.end(); ++it) {
      main_out << std::setw(3) << it->first;
      for(int pes = 0; pes < pes_size; ++pes)
	main_out << std::setw(15) << it->second[pes] * ej_fac;
      main_out << "\n";
    }
    main_out << "\n";
  }  
  return 0;
}
