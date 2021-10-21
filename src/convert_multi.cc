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
  const InputKey  multi_key = "MultiInputFile"; // Multi input file
  std::string multi_input_file;

  const InputKey weight_key = "Weight"; // electronic state weights
  std::vector<std::vector<double> > weight;

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
    else if(key == weight_key) {
      from >> itemp;
      weight.resize(itemp);
      std::getline(from, comment);
      for(int i = 0; i < weight.size(); ++i) {
	std::getline(from, line);
	std::istringstream ist(line);
	while(ist >> dtemp)
	  weight[i].push_back(dtemp);
      }
      if(!from) {
	std::cout << funame << key << " reading failed\n";
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

  // checking weight dimensions
  for(int i = 0; i < weight.size(); ++i)
    if(weight[i].size() != pes_size) {
      std::cout << funame << "weights for " << i 
		<< "-th interpolation are inconsistent with PES size\n";
      exit(1);
    }

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


  const int wsize = weight.size() ? weight.size() : pes_size;

  // thermal fluxes
  Array_2<std::vector<double> > flux(temperature.size(), pes_size);
  Array_2<std::vector<double> > wlux(temperature.size(), wsize);
  Array_2<std::vector<double> > mlux(temperature.size(), wsize);

  // thermal flux fluctuation
  Array_2<std::vector<double> > fluc(temperature.size(), pes_size);
  Array_2<std::vector<double> > wluc(temperature.size(), wsize);
  Array_2<std::vector<double> > mluc(temperature.size(), wsize);

  // microcannonical fluxes
  Array_2<std::vector<double> > e_flux(ener_grid.size(), pes_size);
  Array_2<std::vector<double> > e_wlux(ener_grid.size(), wsize);
  Array_2<std::vector<double> > e_mlux(ener_grid.size(), wsize);

  // ej-resolved fluxes
  Array_3<std::vector<double> > ej_flux(ener_grid.size(), amom_grid.size(), pes_size);
  Array_3<std::vector<double> > ej_wlux(ener_grid.size(), amom_grid.size(), wsize);
  Array_3<std::vector<double> > ej_mlux(ener_grid.size(), amom_grid.size(), wsize);

  // minimum flux dividing surface
  Array_2<int> min_ds(temperature.size(), wsize);
  
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

      flux_in >> itemp;// sampling numbers
      std::getline(flux_in, token);

      for(int pes = 0; pes < pes_size; ++pes)// minimal energies
	flux_in >> dtemp;

      for(int pes = 0; pes < pes_size; ++pes)// minimal configurations
	flux_in >> temp_dv;


      for(int i = 0; i < flux.size(); ++i) {// thermal flux
	flux_in >> flux[i][face] >> dtemp;
	dtemp *= flux[i][face] / 100.;
	fluc[i][face] = dtemp * dtemp;
      }

      for(int i = 0; i < e_flux.size(); ++i)// microcannonical flux
	flux_in >> e_flux[i][face] >> dtemp;


      for(int i = 0; i < ej_flux.size(); ++i)// ej-resolved flux
	flux_in >> ej_flux[i][face] >> dtemp;
      
    }// face cycle

    if(!flux_in)
      error("main: flux file is corrupted");

    // convert to weighted fluxes
    if(weight.size()) {
      for(int i=0; i < wlux.size(); ++i) {
	wlux[i].resize(ds_size);
	wluc[i].resize(ds_size);
	for(int face = 0; face < ds_size; ++face) 
	  wluc[i][face] = wlux[i][face] = 0.;
	
      }
      for(int i=0; i < e_wlux.size(); ++i) {
	e_wlux[i].resize(ds_size);
	for(int face = 0; face < ds_size; ++face)
	  e_wlux[i][face] = 0.;
      }
      for(int i=0; i < ej_wlux.size(); ++i) {
	ej_wlux[i].resize(ds_size);
	for(int face = 0; face < ds_size; ++face)
	  ej_wlux[i][face] = 0.;
      }

      for(int face = 0; face < ds_size; ++face) {
	for(int w = 0; w < weight.size(); ++w)
	  for(int t = 0; t < temperature.size(); ++t)
	    for(int pes = 0; pes < weight[w].size(); ++pes) {
	      wlux(t,w)[face] += weight[w][pes] * flux(t, pes)[face];
	      wluc(t,w)[face] += weight[w][pes] * weight[w][pes] * fluc(t, pes)[face];
	    }

	for(int w = 0; w < weight.size(); ++w)
	  for(int en = 0; en < ener_grid.size(); ++en)
	    for(int pes = 0; pes < weight[w].size(); ++pes)
	      e_wlux(en,w)[face] += weight[w][pes] 
		* e_flux(en, pes)[face];
    

	for(int w = 0; w < weight.size(); ++w)
	  for(int en = 0; en < ener_grid.size(); ++en)
	    for(int am = 0; am < amom_grid.size(); ++am)
	      for(int pes = 0; pes < weight[w].size(); ++pes)
		ej_wlux(en, am, w)[face] += weight[w][pes] 
		  * ej_flux(en, am, pes)[face]; 
      }
    }
    else {
      wlux = flux;
      wluc = fluc;
      e_wlux = e_flux;
      ej_wlux = ej_flux;
    }

    // find minima
    for(int i = 0; i < wlux.size(); ++i)
      if(is_first_ds) {
	mlux[i] = wlux[i];
	mluc[i] = wluc[i];
        min_ds[i] = ds;
      }
      else {
	dtemp = 0.;
	for(int face = 0; face < wlux[i].size(); ++face)
	  dtemp += wlux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < mlux[i].size(); ++face)
	  dtemp1 += mlux[i][face];
	
	if(dtemp < dtemp1) {
	  mlux[i] = wlux[i];
	  mluc[i] = wluc[i];
	  min_ds[i] = ds;
	}
      }
    
    for (int i = 0; i < e_wlux.size(); ++i)
      if(is_first_ds) 
	e_mlux[i] = e_wlux[i];
      else {
	dtemp = 0.;
	for(int face = 0; face < e_wlux[i].size(); ++face)
	  dtemp += e_wlux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < e_mlux[i].size(); ++face)
	  dtemp1 += e_mlux[i][face];
	if(dtemp < dtemp1)
	  e_mlux[i] = e_wlux[i];
      }

    for (int i = 0; i < ej_wlux.size(); ++i)
      if(is_first_ds) 
	ej_mlux[i] = ej_wlux[i];
      else {
	dtemp = 0.;
	for(int face = 0; face < ej_wlux[i].size(); ++face)
	  dtemp += ej_wlux[i][face];
	dtemp1 = 0.;
	for(int face = 0; face < ej_mlux[i].size(); ++face)
	  dtemp1 += ej_mlux[i][face];
	if(dtemp < dtemp1)
	  ej_mlux[i] = ej_wlux[i];
      }

    is_first_ds = false;
  }//ds cycle

  std::cout << "flux through " << count << " surfaces was read\n";

  //Array_3<double>   cn_rate(temperature.size(), wsize);
  Array_3<double>   mc_rate(temperature.size(), wsize, ds_size_max);
  Array_3<double>   ej_rate(temperature.size(), wsize, ds_size_max);


  int en_start = 0;
  for (; en_start < ener_grid.size(); ++en_start)
    if (Phys_const::hart2kelv(ener_grid[en_start]) > 0.1)
      break;

  Array<double> ener_arr(ener_grid.size() - en_start + 1);
  ener_arr[0] = 0.;
  for (int en_ind = 1; en_ind < ener_arr.size(); ++en_ind)
    ener_arr[en_ind] = ener_grid[en_ind + en_start - 1];

  Array<double> enint(ener_arr.size());
  enint[0] = 0.;

  // microcanonical integration
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {// t_cycle
    double t_val = temperature[t_ind];
    double temper_fac = power(t_val, dof_num/2); 
    if (dof_num%2)
      temper_fac *= sqrt(t_val); 
    for(int face = 0; face < ds_size_max; ++face)
      for(int pes = 0; pes < wsize; ++pes) {//pes_cycle
	for(int en_ind = 1; en_ind < ener_arr.size(); ++en_ind) {
	  itemp = en_ind + en_start - 1;
	  if(face < e_mlux(itemp, pes).size())
	    enint[en_ind] =  exp(-ener_arr[en_ind] / t_val) 
	      * e_mlux(itemp, pes)[face];
	  else
	    enint[en_ind] = 0.;
	}
	davint_(ener_arr.data(), enint.data(), ener_arr.size(), 
		ener_arr[0], ener_arr[ener_arr.size()-1], 
		mc_rate(t_ind, pes, face), itemp);
	if(itemp != 1)
	  cout << "main: davint integration error" << endl;
	mc_rate(t_ind, pes, face) *= 4. * M_PI / temper_fac;
      }//pes_cycle
  }// t_cycle
  

  std::cout << "microcannonical integration done\n";

  // EJ-resolved integration
  Array<double> amom_arr(amom_grid.size());
  for (int am_ind = 0; am_ind < amom_grid.size(); ++am_ind)
    amom_arr[am_ind] = amom_grid[am_ind];

  Array<double> amint(amom_grid.size());

  // angular momentum integration
  for(int i = 0; i < e_wlux.size(); ++i)
    e_wlux[i].resize(ds_size_max);
  for(int face = 0; face < ds_size_max; ++face) 
    for(int pes = 0; pes < wsize; ++pes)
      for(int en_ind = 0; en_ind < ener_grid.size(); ++en_ind) {
	for(int am_ind = 0; am_ind < amom_grid.size(); ++am_ind)
	  if(face < ej_mlux(en_ind, am_ind, pes).size())
	    amint[am_ind] = amom_arr[am_ind] * amom_arr[am_ind] 
	      * ej_mlux(en_ind, am_ind, pes)[face];
	  else
	    amint[am_ind] = 0.;
	davint_(amom_arr.data(), amint.data(), amom_arr.size(), 
		amom_arr[0], amom_arr[amom_arr.size()-1], 
		e_wlux(en_ind, pes)[face], itemp);
	if (itemp != 1)
	  cout << "main: davint integration error" << endl;
      }

  // energy integration
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind)  {//t_cycle
    double t_val = temperature[t_ind];
    double temper_fac = power(t_val, dof_num/2); 
    if(dof_num%2)
      temper_fac *= sqrt(t_val); 
    for(int face = 0; face < ds_size_max; ++face)
      for(int pes = 0; pes < wsize; ++pes) {//pes_cycle
	for (int en_ind = 1; en_ind < ener_arr.size(); ++en_ind)
	  enint[en_ind] =  exp(-ener_arr[en_ind] / t_val) 
	    * e_wlux(en_ind + en_start - 1, pes)[face];
	davint_(ener_arr.data(), enint.data(), ener_arr.size(), 
		ener_arr[0], ener_arr[ener_arr.size()-1], 
		ej_rate(t_ind, pes, face), itemp);
	if (itemp != 1)
	  cout << "main: davint integration error" << endl;
	ej_rate(t_ind, pes, face) *= 4. * M_PI / temper_fac;
      }//pes_cycle
  }//t_cycle

  
  std::cout << "EJ-resolved integration done\n";

  // output
  ofstream main_out ("tst.out");
  main_out.precision(4);

  const double conv_fac = 612.6;
  // EJ-resolved rate constant
  main_out << "EJ-resolved rate constant (10^11 cm^3/sec):\n\n";
  main_out << "TOTAL:\n";
  main_out << std::setw(5) << "T, K";
  for(int pes = 0; pes < wsize; ++pes)
    main_out << std::setw(10) << pes;
  main_out << "\n";
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setw(5) << tmpr;
      for(int pes = 0; pes < wsize; ++pes) {
	dtemp = 0.;
	for(int face = 0; face < ds_size_max; ++face) {
	  dtemp += ej_rate(t_ind, pes, face);
	}
	main_out << std::setw(10) << dtemp * conv_fac;
      }
      main_out << "\n";
    }
  }
  main_out << "\n";

  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face) {
      main_out << "FACE = " << face << "\n";
      main_out << std::setw(5) << "T, K";
      for(int pes = 0; pes < wsize; ++pes)
	main_out << std::setw(10) << pes;
      main_out << "\n";
      for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
	double tmpr = temperature[t_ind] / Phys_const::kelv;
	if(tmpr < 5000.) {
	  main_out << std::setw(5) << tmpr;
	  for(int pes = 0; pes < wsize; ++pes)
	    main_out << std::setw(10) 
		     << ej_rate(t_ind, pes, face) * conv_fac;
	  main_out << "\n";
	}
      }
    }
  main_out << "\n";
  
  
  // Microcannonical rate constant
  main_out << "Microcannonical rate constant (10^11 cm^3/sec):\n\n";
  main_out << "TOTAL:\n";
  main_out << std::setw(5) << "T, K";
  for(int pes = 0; pes < wsize; ++pes)
    main_out << std::setw(10) << pes;
  main_out << "\n";
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setw(5) << tmpr; 
      for(int pes = 0; pes < wsize; ++pes) {
	dtemp = 0.;
	for(int face = 0; face < ds_size_max; ++face) {
	  dtemp += mc_rate(t_ind, pes, face);
	}
	main_out << std::setw(10) << dtemp * conv_fac;
      }
      main_out << "\n";
    }
  }
  main_out << "\n";
  
  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face) {
      main_out << "FACE = " << face << "\n";
      main_out << std::setw(5) << "T, K";
      for(int pes = 0; pes < wsize; ++pes)
	main_out << std::setw(10) << pes;
      main_out << "\n";
      for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
	double tmpr = temperature[t_ind] / Phys_const::kelv;
	if(tmpr < 5000.) {
	  main_out << std::setw(5) << tmpr;
	  for(int pes = 0; pes < wsize; ++pes)
	    main_out << std::setw(10) 
		     << mc_rate(t_ind, pes, face) * conv_fac;
	  main_out << "\n";
	}
      }
    }
  main_out << "\n";

  // Cannonical rate constant
  main_out << "Cannonical rate constant (10^11 cm^3/sec):\n\n";
  main_out << "TOTAL:\n";
  main_out << std::setw(5) << "T, K";
  for(int pes = 0; pes < wsize; ++pes)
    main_out << std::setw(18) << pes;
  main_out << "\n";
  for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
    double tmpr = temperature[t_ind] / Phys_const::kelv;
    if(tmpr < 5000.) {
      main_out << std::setprecision(4) << std::setw(5) << tmpr;
      for(int pes = 0; pes < wsize; ++pes) {
	dtemp1 = dtemp = 0.;
	for(int face = 0; face < ds_size_max; ++face) {
	  dtemp  += mlux(t_ind, pes)[face];
	  dtemp1 += mluc(t_ind, pes)[face];
	}
	dtemp1 = 100. * sqrt(dtemp1) / dtemp;
        main_out << setprecision(4) << std::setw(10) 
                 << dtemp * conv_fac;
        main_out << "(";
        if(dtemp1 >= 99.)
          main_out << std::setw(3) << "***";
        else if(dtemp1 >= 1.)
          main_out << setprecision(2) << std::setw(3) 
                   << dtemp1;
        else
          main_out << std::setw(3) << "< 1";
            
        main_out << "|" << std::setw(2) << min_ds(t_ind, pes) << ")";
      }
      main_out << "\n";
    }
  }
  main_out << "\n" << setprecision(4);
  if(ds_size_max != 1)
    for(int face = 0; face < ds_size_max; ++face) {
      main_out << "FACE = " << face << "\n";
      main_out << std::setw(5) << "T, K";
      for(int pes = 0; pes < wsize; ++pes)
	main_out << std::setw(10) << pes;
      main_out << "\n";
      for(int t_ind = 0; t_ind < temperature.size(); ++t_ind) {
	double tmpr = temperature[t_ind] / Phys_const::kelv;
	if(tmpr < 5000.) {
	  main_out << std::setw(5) << tmpr;
	  for(int pes = 0; pes < wsize; ++pes)
	    main_out << std::setw(10) 
		     << mlux(t_ind, pes)[face] * conv_fac;
	  main_out << "\n";
	}
      }
    }
  main_out << "\n";
}
