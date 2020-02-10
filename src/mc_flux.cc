#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <set>

#include "input.hh"
#include "units.hh"
#include "rotd.hh"
#include "math.hh"
#include "tmatrix.hh"
#include "slatec.hh"
#include "error.hh"
#include "divsur.hh"
#include "interflux.hh"

int main (int argc, char* argv[])
{
  const char funame [] = "mc_flux: ";

  if (argc != 2)
    error("usage: mc_flux inp_file");

  // temporary variables
  int itemp;
  double dtemp, dtemp1;
  std::string token, line, comment, key;

  // input 
  const InputKey multi_key  = "MultiInputFile";// input file for multi
  const InputKey  corr_key = "EnergyCorrection(kcal/mol)";// energy correction
  const InputKey   out_key = "OutputFile";// ej_flux output file
  const InputKey  face_key = "Face"; 
  const InputKey   pes_key = "ElectronicSurface";
  const InputKey egrid_key = "EnergyGrid(1/cm)";// new energy grid

  Grid new_ener_grid;
  std::string multi_input_file;
  std::map<std::pair<int, int>, double> ener_shift;
  std::string out_file;
  std::set<int> ds_face;
  int pes_val = 0;

  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
    exit(1);
  }

  double min_val, step;
  int ds, ds_size, nstep, size;

  while(from >> token) {
    if(multi_key == token) {
      from >> multi_input_file;
    }
    else if(key == corr_key) {
      from >> size;
      if(!from) {
	std::cout << funame << token << "failed to read number of entries\n";
	exit(1);
      }
      std::getline(from, comment);
      for(int i = 0; i < size; ++i) {
	std::getline(from, line);
	std::istringstream ist(line);
	int sur, face;
	ist >> sur >> face >> dtemp;
	if(!ist) {
	  std::cout << funame << token
		    << ": entry format: sur# face# shift\n";
	  exit(1);
	}
	ener_shift[std::make_pair(sur, face)] = dtemp * Phys_const::kcal;
      }
      if(!from) {
	std::cout << funame << token << ": wrong number of entries\n";
	exit(1);
      }
    }
    else if(pes_key == token) {// surface
      from >> pes_val;
      if(!from) {
	std::cout << funame << token << " reading failed\n";
	exit(1);
      }
    }
    else if(out_key == token) {// output file name
      from >> out_file;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
    }
    else if (face_key == token) {// read face
      std::getline(from, line);
      std::istringstream iss(line);      
      while(iss >> itemp)
	ds_face.insert(itemp);
    }
    else if(egrid_key == token) {// energy grid
      from >> min_val >> nstep >> step;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
      new_ener_grid.set(Phys_const::incm*min_val, nstep, Phys_const::incm*step);
    }
    else {
      std::cout << funame << "unknown keyword " << token << ". Known keywords are:\n";
      InputKey::show_all();
      exit(1);
    }
  }
  from.close();
  from.clear();

  // reading multi input
  std::map <string, Read> input_data;
  typedef map<string, Read>::iterator input_data_it;
 
  std::vector<double> temperature;
  input_data ["tmpr_grid"]   = Read(temperature);
  std::vector<double> ener_grid;
  input_data ["ener_grid"]   = Read(ener_grid);
  std::vector<double> amom_grid;
  input_data ["amom_grid"]  = Read(amom_grid);

  std::string flux_file; // flux output file
  input_data ["flux_file"]  = Read(flux_file, "flux.dat");      
  std::string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");

  int pes_size;
  input_data ["pes_size"] = Read(pes_size, 1);

  from.open(multi_input_file.c_str());
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
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
      std::cout << funame << it->first << " is not initialized\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  //std::cout << "Default parameters:\n";
  //for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
  //if(it->second.is_default())
  //std::cout << it->first << "  " << it->second << "\n";
  //std::cout << "\n";

  // read molecular structure from the file  
  mol_init(mol_spec_file);

  Dynvar temp_dv;
  Array_2<double> ej_flux(ener_grid.size(), amom_grid.size());
  Array_2<double> total_ej_flux(ener_grid.size(), amom_grid.size());
  Matrix<double> min_flux(ener_grid.size(), amom_grid.size());
  Matrix<double> total_min_flux(ener_grid.size(), amom_grid.size());
  double vmin;

  for(int en = 0; en < ener_grid.size(); ++en)
    for(int am = 0; am < amom_grid.size(); ++am)
      total_min_flux[en][am] = min_flux[en][am] = -1.;

  std::map<std::pair<int, int>, double>::const_iterator esp;
  // shifted grid and flux
  Array<double> shift_grid(ener_grid.size());
  Array<double> shift_flux(ener_grid.size());

  // read from flux file
  from.open(flux_file.c_str());
  if(!from) {
    std::cout << funame << "cannot open " << flux_file.c_str() << " file\n";
    exit(1);
  }

  while (from >> ds >> ds_size) {//input cycle
    std::cout << "ds = " << ds << " ... ";
    total_ej_flux.init();
    ej_flux.init();

    // read ej-resolved flux
    for(int face = 0; face < ds_size; ++face) {//face cycle      
      // sampling numbers
      from >> itemp;
      std::getline(from, token);

      // minimum energy value
      for(int pes = 0; pes < pes_size; ++pes)
	if(pes != pes_val)
	  from >> dtemp;
	else
	  from >> vmin;

      // minimal configurations
      for(int pes = 0; pes < pes_size; ++pes)
	from >> temp_dv;

      // thermal flux
      for(int pes = 0; pes < pes_size; ++pes)
	for(int i = 0; i < temperature.size(); ++i)
	  from >> dtemp >> dtemp;

      // microcannonical flux
      for(int pes = 0; pes < pes_size; ++pes)
	for(int i = 0; i < ener_grid.size(); ++i)
	  from >> dtemp >> dtemp;

      // ej_resolved flux
      esp = ener_shift.find(std::make_pair(ds, face));
      // energy shift applied
      if(esp != ener_shift.end()) {
	vmin += esp->second;
	for(int en = 0; en < ener_grid.size(); ++en)
	  shift_grid[en] +=  esp->second + ener_grid[en];

	for(int pes = 0; pes < pes_size; ++pes)
	  for(int am = 0; am < amom_grid.size(); ++am)
	    if(pes != pes_val)
	      for(int en = 0; en < ener_grid.size(); ++en)
		from >> dtemp >> dtemp;
	    else {
	      for(int en = 0; en < ener_grid.size(); ++en)
		from >> shift_flux[en] >> dtemp;
	      InterFlux inter_flux(shift_grid, shift_flux, vmin);
	      for(int en = 0; en < ener_grid.size(); ++en) {
		dtemp = inter_flux(ener_grid[en]);
		total_ej_flux(en, am) += dtemp;
		if(ds_face.find(face) != ds_face.end())
		  ej_flux(en, am) += dtemp;
	      }
	    }
      } 
      // no shift
      else {
	for(int pes = 0; pes < pes_size; ++pes)// ej_resolved flux
	  for(int am = 0; am < amom_grid.size(); ++am)
	    if(pes != pes_val)
	      for(int en = 0; en < ener_grid.size(); ++en)
		from >> dtemp >> dtemp;
	    else
	      for(int en = 0; en < ener_grid.size(); ++en) {
		from >> dtemp >> dtemp1;
		total_ej_flux(en, am) += dtemp;
		if(ds_face.find(face) != ds_face.end())
		  ej_flux(en, am) += dtemp;
	      }
      }
    }// face cycle
    std::cout << "done" << std::endl;

    if(!from) {
      std::cout << funame << "flux file " << flux_file << " is corrupted\n";
      exit(1);
    }

    // minimal flux
    for(int am = 0; am < amom_grid.size(); ++am)
      for(int en = 0; en < ener_grid.size(); ++en) {
	dtemp = total_ej_flux(en, am);
	if(total_min_flux[en][am] < 0. || dtemp < total_min_flux[en][am]) {
	  total_min_flux[en][am] = dtemp;
	  if(!ds_face.size())
	    min_flux[en][am] = dtemp;
	  else
	    min_flux[en][am] = ej_flux(en, am);
	}
      }
  }//input cycle
  from.close();
  from.clear();

  // check that all flux values vere initialized
  for(int am = 0; am < amom_grid.size(); ++am)
    for(int en = 0; en < ener_grid.size(); ++en)
      if(min_flux[en][am] < 0.) {
	std::cout << funame << "flux with energy = "
		  << ener_grid[en]/Phys_const::incm	 << "cm^-1 not initialized\n";
	exit(1);
      }

  // normalization factor
  double nfac = 2. * sqrt(2. * M_PI) * mol_array[0]->stat_sum() * mol_array[1]->stat_sum() *
    std::pow(mol_array[0]->mass() * mol_array[1]->mass() / 
	     (mol_array[0]->mass() + mol_array[1]->mass()), 1.5);


  for(int am = 0; am < amom_grid.size(); ++am) {
    dtemp = amom_grid[am] * amom_grid[am];
    for(int en = 0; en < ener_grid.size(); ++en)
      min_flux[en][am] *= dtemp;
  }

  Array<double> amom_array(amom_grid.size());
  for(int am = 0; am < amom_grid.size(); ++am)
    amom_array[am] = amom_grid[am];
  
  // angular momentum integration
  Array<double> ener_array(ener_grid.size());
  Array<double> stat_array(ener_grid.size());
  int ishift = 0;
  if(ener_grid[0] > 0.) {
    ener_array.resize(ener_grid.size() + 1);
    stat_array.resize(ener_grid.size() + 1);
    ener_array[0] = 0.;
    stat_array[0] = 0.;
    ishift = 1;
  }
    
  for(int en = 0; en < ener_grid.size(); ++en) {
    davint_(amom_array.data(), min_flux[en], amom_array.size(), amom_array[0],
	    amom_array[amom_array.size()-1], dtemp, itemp);
    if (itemp != 1)
      cout << "main: davint integration error" << endl;
    itemp = en + ishift;
    ener_array[itemp] = ener_grid[en];
    stat_array[itemp] = dtemp * nfac;
  }

  Spline states(ener_array.data(), stat_array.data(), ener_array.size());
  dtemp = stat_array.back() / states.fit(ener_array.back() * 0.9, 0);
  const double stat_pow = std::log(dtemp) * 10.;
    
  // output
  ofstream ej_out(out_file.c_str());
  if(new_ener_grid.size())
    for(int en = 0; en < new_ener_grid.size(); ++en) {
      double ener = new_ener_grid[en];
      ej_out << std::setw(15) << ener / Phys_const::incm;
      if(ener < ener_array.back())
	ej_out << std::setw(15) << states.fit(ener, 0);
      else
	ej_out << std::setw(15) << stat_array.back() * std::pow(ener/ener_array.back(), stat_pow);
      ej_out << "\n";
    }
  else
    for(int en = 0; en < ener_array.size(); ++en) {
      ej_out << std::setw(15) << ener_array[en] / Phys_const::incm
	     << std::setw(15) << stat_array[en]
	     << "\n";
    }

  return 0;
}
