#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
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
#include "interflux.hh"

int main (int argc, char* argv[])
{
  const char funame [] = "convert_multi: ";

  // usage
  if(argc != 2)
    error("usage: convert_multi input_file");

  // temporary variables
  int itemp, size;
  double dtemp, dtemp1;
  std::string key, comment, token, line;

  // input
  std::string multi_input_file;
  const InputKey  multi_key = "MultiInputFile"; // Multi input file

  int pes_val = 0;
  const InputKey pes_key = "ElectronicSurface"; // electronic state weights

  std::map<std::pair<int, int>, double> ener_shift;
  const InputKey corr_key = "EnergyCorrection";

  // reading input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << "cannot open input file is not found";
  }

  while(from >> key)
    if(key == multi_key) {
      from >> multi_input_file;
      std::getline(from, comment);
      if(!from) {
	std::cout << funame << key << " reading failed\n";
	exit(1);
      }
    }
    else if(key == pes_key) {
      from >> pes_val;
      if(!from) {
	std::cout << funame << key << " reading failed\n";
	exit(1);
      }
    }
    else if(key == corr_key) {
      from >> size;
      if(!from) {
	std::cout << funame << key << "failed to read number of entries\n";
	exit(1);
      }
      std::getline(from, comment);
      for(int i = 0; i < size; ++i) {
	std::getline(from, line);
	std::istringstream ist(line);
	int sur, face;
	ist >> sur >> face >> dtemp;
	if(!ist) {
	  std::cout << funame << key 
		    << ": entry format: sur# face# shift\n";
	  exit(1);
	}
	ener_shift[std::make_pair(sur, face)] = dtemp * Phys_const::kcal;
      }
      if(!from) {
	std::cout << funame << key << " wrong number of entries\n";
	exit(1);
      }
    }
    else {
      std::cout << funame << "unknown keyword " << key << ". Known keywords are:\n";
      InputKey::show_all();
      exit(1);
    }
  from.close();
  from.clear();

  // reading multi input file
  std::map<string, Read> input_data;
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
  if (!from) {
    std::cout << "cannot open " << multi_input_file << " file\n";
    exit(1);
  }
  
  while(from >> key) {
    input_data_it i = input_data.find(key);
    if (i == input_data.end())
      getline(from, comment);
    else
      from >> i->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  bool is_init = true;
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(!it->second.is_init()) {
      std::cout << it->first << " is not initialized, exitting\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  std::cout << "Default parameters:\n";
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(it->second.is_default())
      cout << it->first << "  " << it->second << "\n";
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

  // fluxes
  Array<std::vector<double> > flux(temperature.size());
  Array<std::vector<double> > fluc(temperature.size());
  Array<std::vector<double> > e_flux(ener_grid.size());
  Array_2<std::vector<double> > ej_flux(ener_grid.size(), amom_grid.size());

  Array<std::vector<double> > mlux(temperature.size());
  Array<std::vector<double> > mluc(temperature.size());
  Array<std::vector<double> > e_mlux(ener_grid.size());
  Array_2<std::vector<double> > ej_mlux(ener_grid.size(), amom_grid.size());

  // minimum flux dividing surface
  Array<int> min_ds(temperature.size());

  // shifted grid and flux
  Array<double> shift_grid(ener_grid.size());
  Array<double> shift_flux(ener_grid.size());

  // minimum energy
  double vmin;
  
  int ds, ds_size;
  Dynvar temp_dv;

  ifstream flux_in(flux_file.c_str());
  bool is_first_ds = true;
  int count = 0;
  int ds_size_max = 0;
  while(flux_in >> ds >> ds_size) {//ds cycle
    ds_size_max = ds_size > ds_size_max ? ds_size : ds_size_max;

    ++count;

    for(int i = 0; i < flux.size(); ++i) {
      flux[i].resize(ds_size);
      fluc[i].resize(ds_size);
    }
    for(int i = 0; i < e_flux.size(); ++i)
      e_flux[i].resize(ds_size);
    for(int i = 0; i < ej_flux.size(); ++i)
      ej_flux[i].resize(ds_size);

    for(int face = 0; face < ds_size; ++face) {//face cycle
      std::map<std::pair<int, int>, double>::const_iterator esp=
	ener_shift.find(std::make_pair(ds, face));

      flux_in >> itemp;// sampling numbers
      std::getline(flux_in, token);

      for(int pes = 0; pes < pes_size; ++pes)// minimal energies
	if(pes != pes_val)
	  flux_in >> dtemp;
	else
	  flux_in >> vmin;
      
      
      for(int pes = 0; pes < pes_size; ++pes)// minimal configurations
	flux_in >> temp_dv;

      if(esp != ener_shift.end()) {// energy shift

	vmin += esp->second;
	for(int en = 0; en < ener_grid.size(); ++en)
	  shift_grid[en] = ener_grid[en] + esp->second;

	for(int pes = 0; pes < pes_size; ++pes) // thermal flux
	  if(pes != pes_val)
	    for(int t = 0; t < temperature.size(); ++t)
	      flux_in >> dtemp >> dtemp;
	  else {
	    for(int t = 0; t < temperature.size(); ++t) {
	      flux_in >> flux[t][face] >> dtemp;
	      flux[t][face] /= exp(esp->second / temperature[t]);
	      dtemp *= flux[t][face] / 100.;
	      fluc[t][face] = dtemp * dtemp;
	    }
	  }

	for(int pes = 0; pes < pes_size; ++pes)// microcannonical flux
	  if(pes != pes_val)
	    for(int en = 0; en < ener_grid.size(); ++en)
	      flux_in >> dtemp >> dtemp;
	  else {
	    for(int en = 0; en < ener_grid.size(); ++en)
	      flux_in >> shift_flux[en] >> dtemp;
	    InterFlux inter_flux(shift_grid, shift_flux, vmin);
	    for(int en = 0; en < ener_grid.size(); ++en)
	      e_flux[en][face] = inter_flux(ener_grid[en]);
	  }

	for(int pes = 0; pes < pes_size; ++pes)// ej_resolved flux
	  for(int am = 0; am < amom_grid.size(); ++am)
	    if(pes != pes_val)
	      for(int en = 0; en < ener_grid.size(); ++en)
		flux_in >> dtemp >> dtemp;
	    else {
	      for(int en = 0; en < ener_grid.size(); ++en)
		flux_in >> shift_flux[en] >> dtemp;
	      InterFlux inter_flux(shift_grid, shift_flux, vmin);
	      for(int en = 0; en < ener_grid.size(); ++en)
		ej_flux(en, am)[face] = inter_flux(ener_grid[en]);
	    }
      }// energy shift
      else {// no energy shift
	for(int pes = 0; pes < pes_size; ++pes) // thermal flux
	  if(pes != pes_val)
	    for(int t = 0; t < temperature.size(); ++t)
	      flux_in >> dtemp >> dtemp;
	  else {
	    for(int t = 0; t < temperature.size(); ++t) {
	      flux_in >> flux[t][face] >> dtemp;
	      dtemp *= flux[t][face] / 100.;
	      fluc[t][face] = dtemp * dtemp;
	    }
	  }

	for(int pes = 0; pes < pes_size; ++pes)// microcannonical flux
	  if(pes != pes_val)
	    for(int en = 0; en < ener_grid.size(); ++en)
	      flux_in >> dtemp >> dtemp;
	  else
	    for(int en = 0; en < ener_grid.size(); ++en)
	      flux_in >> e_flux[en][face] >> dtemp;

	for(int pes = 0; pes < pes_size; ++pes)// ej_resolved flux
	  for(int am = 0; am < amom_grid.size(); ++am)
	    for(int en = 0; en < ener_grid.size(); ++en)
	      if(pes != pes_val)
		flux_in >> dtemp >> dtemp;
	      else
		flux_in >> ej_flux(en, am)[face] >> dtemp;

      }// no energy shift
    }// face cycle

    if(!flux_in)
      error("main: flux file is corrupted");

    // find minima
    for(int i = 0; i < flux.size(); ++i)
      if(is_first_ds) {
	mlux[i] = flux[i];
	mluc[i] = fluc[i];
        min_ds[i] = ds;
      }
      else {
	dtemp = 0.;
	for(int face = 0; face < flux[i].size(); ++face)
	  dtemp += flux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < mlux[i].size(); ++face)
	  dtemp1 += mlux[i][face];
	
	if(dtemp < dtemp1) {
	  mlux[i] = flux[i];
	  mluc[i] = fluc[i];
	  min_ds[i] = ds;
	}
      }
    
    for (int i = 0; i < e_flux.size(); ++i)
      if(is_first_ds) 
	e_mlux[i] = e_flux[i];
      else {
	dtemp = 0.;
	for(int face = 0; face < e_flux[i].size(); ++face)
	  dtemp += e_flux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < e_mlux[i].size(); ++face)
	  dtemp1 += e_mlux[i][face];
	if(dtemp < dtemp1)
	  e_mlux[i] = e_flux[i];
      }

    for (int i = 0; i < ej_flux.size(); ++i)
      if(is_first_ds) 
	ej_mlux[i] = ej_flux[i];
      else {
	dtemp = 0.;
	for(int face = 0; face < ej_flux[i].size(); ++face)
	  dtemp += ej_flux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < ej_mlux[i].size(); ++face)
	  dtemp1 += ej_mlux[i][face];
	if(dtemp < dtemp1)
	  ej_mlux[i] = ej_flux[i];
      }

    is_first_ds = false;
  }//ds cycle

  std::cout << "flux through " << count << " surfaces was read\n";

  int en_start = 0;
  for (; en_start < ener_grid.size(); ++en_start)
    if(Phys_const::hart2kelv(ener_grid[en_start]) > 0.1)
      break;

  Array<double> ener_arr(ener_grid.size() - en_start + 1);
  ener_arr[0] = 0.;
  for (int en_ind = 1; en_ind < ener_arr.size(); ++en_ind)
    ener_arr[en_ind] = ener_grid[en_ind + en_start - 1];

  //Array_3<double>   cn_rate(temperature.size(), wsize);
  Array_2<double>   mc_rate(temperature.size(), ds_size_max);
  Array_2<double>   ej_rate(temperature.size(), ds_size_max);


  Array<double> enint(ener_arr.size());
  enint[0] = 0.;

  // microcanonical integration
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {// t_cycle
    double t_val = temperature[t_ind];
    double temper_fac = power(t_val, dof_num/2); 
    if (dof_num%2)
      temper_fac *= sqrt(t_val); 
    for(int face = 0; face < ds_size_max; ++face) {//facet cycle
      for(int en_ind = 1; en_ind < ener_arr.size(); ++en_ind) {
	itemp = en_ind + en_start - 1;
	if(face < e_mlux[itemp].size())
	  enint[en_ind] =  exp(-ener_arr[en_ind] / t_val) 
	    * e_mlux[itemp][face];
	else
	  enint[en_ind] = 0.;
      }
      davint_(ener_arr.data(), enint.data(), ener_arr.size(), 
	      ener_arr[0], ener_arr[ener_arr.size()-1], 
	      mc_rate(t_ind, face), itemp);
      if(itemp != 1)
	cout << "main: davint integration error" << endl;
      mc_rate(t_ind, face) *= 4. * M_PI / temper_fac;
    } // facet cycle
  }// t_cycle
  

  std::cout << "microcannonical integration done\n";

  // EJ-resolved integration
  Array<double> amom_arr(amom_grid.size());
  for (int am_ind = 0; am_ind < amom_grid.size(); ++am_ind)
    amom_arr[am_ind] = amom_grid[am_ind];

  Array<double> amint(amom_grid.size());

  // angular momentum integration
  for(int i = 0; i < e_flux.size(); ++i)
    e_flux[i].resize(ds_size_max);
  for(int face = 0; face < ds_size_max; ++face) 
    for(int en_ind = 0; en_ind < ener_grid.size(); ++en_ind) {
      for(int am_ind = 0; am_ind < amom_grid.size(); ++am_ind)
	if(face < ej_mlux(en_ind, am_ind).size())
	  amint[am_ind] = amom_arr[am_ind] * amom_arr[am_ind] 
	    * ej_mlux(en_ind, am_ind)[face];
	else
	  amint[am_ind] = 0.;
      davint_(amom_arr.data(), amint.data(), amom_arr.size(), 
	      amom_arr[0], amom_arr[amom_arr.size()-1], 
	      e_flux[en_ind][face], itemp);
      if (itemp != 1)
	cout << "main: davint integration error" << endl;
    }

  // energy integration
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind)  {//t_cycle
    double t_val = temperature[t_ind];
    double temper_fac = power(t_val, dof_num/2); 
    if(dof_num%2)
      temper_fac *= sqrt(t_val); 
    for(int face = 0; face < ds_size_max; ++face) {// facet cycle
      for (int en_ind = 1; en_ind < ener_arr.size(); ++en_ind)
	enint[en_ind] =  exp(-ener_arr[en_ind] / t_val) 
	  * e_flux[en_ind + en_start - 1][face];
      davint_(ener_arr.data(), enint.data(), ener_arr.size(), 
	      ener_arr[0], ener_arr[ener_arr.size()-1], 
	      ej_rate(t_ind, face), itemp);
      if (itemp != 1)
	cout << "main: davint integration error" << endl;
      ej_rate(t_ind, face) *= 4. * M_PI / temper_fac;
    }// facet cycle
  }//t_cycle

  
  std::cout << "EJ-resolved integration done\n";

  // output
  ofstream main_out ("tst.out");
  main_out.precision(4);

  const double conv_fac = 612.6;
  // EJ-resolved rate constant
  main_out << "EJ-resolved rate constant (10^11 cm^3/sec):\n\n";
  main_out << std::setw(5) << "T, K" << std::setw(10) << "Total";
  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face)
      main_out << std::setw(10) << face;
  main_out << "\n";

  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setw(5) << tmpr;
      dtemp = 0.;
      for(int face = 0; face < ds_size_max; ++face) {
	dtemp += ej_rate(t_ind, face);
      }
      main_out << std::setw(10) << dtemp * conv_fac;
      if(ds_size_max != 1)
	for(int face = 0; face < ds_size_max; ++face) 
	  main_out << std::setw(10) << ej_rate(t_ind, face) * conv_fac;
      main_out << "\n";
    }
  }
  main_out << "\n";

  // Microcannonical rate constant
  main_out << "Microcannonical rate constant (10^11 cm^3/sec):\n\n";
  main_out << std::setw(5) << "T, K" << std::setw(10) << "Total";
  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face)
      main_out << std::setw(10) << face;
  main_out << "\n";

  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setw(5) << tmpr;
      dtemp = 0.;
      for(int face = 0; face < ds_size_max; ++face) {
	dtemp += mc_rate(t_ind, face);
      }
      main_out << std::setw(10) << dtemp * conv_fac;
      if(ds_size_max != 1)
	for(int face = 0; face < ds_size_max; ++face) 
	  main_out << std::setw(10) << mc_rate(t_ind, face) * conv_fac;
      main_out << "\n";
    }
  }
  main_out << "\n";

  // Cannonical rate constant
  main_out << "Cannonical rate constant (10^11 cm^3/sec):\n\n";
  main_out << std::setw(5) << "T, K" << std::setw(17) << "Total";
  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face)
      main_out << std::setw(10) << face;
  main_out << "\n";

  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setw(5) << tmpr;
      dtemp = dtemp1 = 0.;
      for(int face = 0; face < mlux[t_ind].size(); ++face) {
	dtemp  += mlux[t_ind][face];
	dtemp1 += mluc[t_ind][face];
      }
      main_out << " " << std::setw(9) 
	       << dtemp * conv_fac;
      dtemp1 = 100. * sqrt(dtemp1) / dtemp;
      main_out << "(";
      if(dtemp1 >= 99.)
	main_out << std::setw(3) << "***";
      else if(dtemp1 >= 1.)
	main_out << std::setprecision(2) << std::setw(3) 
		 << dtemp1;
      else
	main_out << std::setw(3) << "< 1";
            
      main_out << "|" << std::setw(2) << min_ds[t_ind] << ")";

      if(ds_size_max != 1)
	for(int face = 0; face < mlux[t_ind].size(); ++face) 
	  main_out << std::setw(10) << mlux[t_ind][face] * conv_fac;
      main_out << "\n";
    }
  }
  main_out << "\n";
}
