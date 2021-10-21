#include <fstream>
#include <string>
#include<sys/stat.h>
#include<cstdlib>
#include <sys/time.h>
#include <csignal>
#include <cstring>
#include <dlfcn.h>
#include <cstdio>

#include "force.hh"
#include "rotd.hh"
#include "math.hh"
#include "error.hh"
#include "log.hh"
#include "units.hh"
#include "system.hh"

#include "multipole.hh"
#include "molpro.hh"
#include "gauss.hh"
#include "sjk.hh"

// molpro template files indices
enum {
  MOLPRO_ENER = 0, // single point energy method 
  MOLPRO_OPT1 = 1, // start optimization/forces method
  MOLPRO_OPT2 = 2  // subsequent optimizations
};

//optimize molecular geometry 
double (*qc_method) (const std::vector<D3>& pos, int flags, 
		     std::vector<D3>& force) = 0;
void (*qc_forces) (std::vector<D3>& forces) = 0; 

void molpro_forces(std::vector<D3>& forces)
{
  const char funame [] = "molpro_forces: ";

  if(!Molpro::isinit()) {
    std::cout << funame << "not initialized\n";
    exit(1);
  }

  // set atomic coordinates
  for(int at = 0; at < atoms.size(); ++at)
    for(int i = 0; i < 3; ++i)
      Molpro::atoms[at].pos[i] = atoms[at].lf_pos[i];

  try {
    Molpro::pot(MOLPRO_OPT1, Molpro::FORCE);
  }
  catch(Molpro::Err) {
    Log::out << funame << "potential calculation failed\n"; 
    throw Pot_error();
  }
  for(int at = 0; at < atoms.size(); ++at)
    forces[at] = Molpro::forces[at];
}

double molpro_opt (const std::vector<D3>& pos, int flags, std::vector<D3>& force)
{

  const char funame [] = "molpro_opt: ";

  if(!Molpro::isinit()) {
    std::cout << funame << "not initialized\n";
    exit(1);
  }

  // set atomic coordinates
  for(int at = 0; at < atoms.size(); ++at)
    for(int i = 0; i < 3; ++i)
      Molpro::atoms[at].pos[i] = pos[at][i];

  int method = flags & OPT_READ ? MOLPRO_OPT2 : MOLPRO_OPT1;

  try {
    Molpro::pot(method, Molpro::WF | Molpro::FORCE);
  }
  catch(Molpro::Err) {
    Log::out << funame << "potential failed\n"; 
    throw Pot_error();
  }
  for(int at = 0; at < atoms.size(); ++at)
    force[at] = Molpro::forces[at];

  if(Molpro::energies.size())
    return Molpro::energies[0];
  else {
    Log::out << funame << "no energies found in molpro output\n"; 
    throw Pot_error();
  }
}

double g98_opt (const std::vector<D3>& pos, int flags, std::vector<D3>& force)
{
  const char funame [] = "g98_opt: ";

  double res;
  
  // set atomic coordinates
  //
  for(int at = 0; at < atoms.size(); ++at)
    //
    for(int i = 0; i < 3; ++i)
      //
      Gauss::cluster[at].pos[i] = pos[at][i];

  int g98_flags = Gauss::FORCE | Gauss::DIRECT;
  
  if(flags & OPT_READ)
    //
    g98_flags |= Gauss::READ;
  
  Gauss::opt_method->init();
  
  if(!Gauss::opt_method->apply(g98_flags)){// method failed
    //
    Log::out << funame << "Oops!!! optimization method " 
	     << Gauss::opt_method->id() << " failed\n";
    
    std::cerr << funame << " optimization method " 
	      << Gauss::opt_method->id() << " failed" << std::endl;
    
    throw Pot_error();
  }
  
  res = Gauss::opt_method->rel_energy();
  
  Gauss::opt_method->read_forces(force);

  return res;
}

void g98_forces (std::vector<D3>& forces)
{
  const char funame [] = "g98_forces: ";

  // set atomic coordinates
  //
  for(int at = 0; at < atoms.size(); ++at)
    //
    for(int i = 0; i < 3; ++i)
      //
      Gauss::cluster[at].pos[i] = atoms[at].lf_pos[i];

  int g98_flags = Gauss::FORCE | Gauss::DIRECT;
  
  Gauss::opt_method->init();

  // method failed
  //
  if(!Gauss::opt_method->apply(g98_flags)){
    //
    Log::out << funame << "Oops!!! optimization method " 
	     << Gauss::opt_method->id() << " failed" << std::endl;
    
    std::cerr << funame << "optimization method " 
	      << Gauss::opt_method->id() << " failed" << std::endl;
    
    throw Pot_error();
  }
  
  Gauss::opt_method->read_forces(forces);
}

void set_opt(const string& type, const string& method)
{
  const char funame [] = "set_opt: ";

  if(type == "full") {
    if(method == "molpro")
      qc_method = molpro_opt;
    else if(method == "gauss" || method == "g98")
      qc_method = g98_opt;
    else {
      std::cout << funame << "unknown optimization method: " 
	       << method << "\navailable methods: molpro, gauss\n";
      exit(1);
    }
    
  }
  else if(type == "one-step") {
    if(method == "molpro")
      qc_forces = molpro_forces;
    else if(method == "gauss" || method == "g98")
      qc_forces = g98_forces;
    else {
      std::cout << funame << "unknown optimization method: " 
	       << method << "\navailable methods: molpro, gauss\n";
      exit(1);
    }
  }
  else if(type == "none")
    ;
  else {
    std::cout << funame << "unknown optimization type: " 
	     << method << "\navailable types: full, one-step, none\n";
    exit(1);
  }
}

double relax_corr ()
{
  const char funame [] = "relax_corr: ";
  if(!qc_forces) {
    std::cout << funame << "method is not set, exitting...\n";
    exit(1);
  }
  double dtemp;
  std::vector<D3> forces(atoms.size());
  qc_forces(forces);

  for (int at = 0; at < atoms.size(); ++at)
    for (int i = 0; i < 3; ++i)
      atoms[at].force[i] = forces[at][i];
    
  double corr = 0.;
  for(int frag = 0; frag < 2; ++frag)
    for(int mode = 0; mode < mol_array[frag]->nm_size(); ++mode) {
      dtemp = mol_array[frag]->nm_force(mode);
      corr -= dtemp * dtemp;
    }
  corr /= 2.;
  Log::out << funame << "one-step relaxation correction = " 
	     << corr / Phys_const::kcal << " kcal/mol\n";
  return corr;
}

double optimize(std::vector<D3>& pos)
{
  const char* funame = "optimize: ";

  const int range_num = 8; 
  static const double max_ener_range [range_num] =
    { -10., -8., -6., -4., -2., 2.,   5.,  10.}; // kcal/mol
  static const double min_corr_range [range_num] = 
    { 0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.3,  0.5}; // kcal/mol
  static const int    corr_num_range [range_num] = 
    {  11,  10,   9,   8,   7,  5,    3,    3};

  const int count_max = 4;
  const double step_max = 0.9;

  if(!qc_method) {
    std::cout << funame << "method is not set, exitting...\n";
    exit(1);
  }

  double dtemp;

  NormalModeCoor q, q_new;
  NormalModeCoor q_corr, q_corr_new;
  double opt_ener, opt_ener_new; 
  double tot_corr, tot_corr_new;
  double frag_corr [2];
  double frag_corr_new [2];

  double min_corr = 0.;
  int corr_num = 0;

  double step [2];

  std::vector<D3> force(atoms.size());

  for(int at = 0; at < atoms.size(); ++at)
    for(int i = 0; i < 3; ++i)
      pos[at][i] = atoms[at].lf_pos[i];
  q.init();

  opt_ener = qc_method(pos, 0, force);
  // relaxation correction
  for (int at = 0; at < atoms.size(); ++at)
    for (int i = 0; i < 3; ++i)
      atoms[at].force[i] = force[at][i];

  tot_corr = 0.;
  for(int frag = 0; frag < 2; ++frag) {
    frag_corr[frag] = 0.;
    for(int mode = 0; mode < q.size(frag); ++mode) {
      dtemp = mol_array[frag]->nm_force(mode);
      q_corr(frag, mode) = dtemp;
      frag_corr[frag] += dtemp * dtemp;
    }
    frag_corr[frag] /= 2.;
    tot_corr += frag_corr[frag];
  }

  Log::out << funame << "Configuration:\n";
  Log::out << funame; print_dist(Log::out);
  Log::out << funame; print_geom(Log::out, pos);
  Log::out << funame << "interaction energy = "
	   << opt_ener / Phys_const::kcal << " kcal/mol\n";
  for(int frag = 0; frag < 2; ++frag)
    Log::out << funame  <<"relaxation energy for "
	     << frag << "-th fragment = "
	     <<  frag_corr[frag] / Phys_const::kcal 
	     << " kcal/mol" << endl;

  if(opt_ener >= max_ener_range[range_num-1] * Phys_const::kcal) {
    Log::out << funame << "interaction is strongly repulsive, "
      "no iterations will be performed\n\n";
    return opt_ener;
  }

  //initial step
  for(int frag = 0; frag < 2; ++frag)
    if(frag_corr[frag] > 10. * Phys_const::kcal)
      step[frag] = 0.25;
    else if(frag_corr[frag] > 5. * Phys_const::kcal)
      step[frag] = 0.3;
    else if(frag_corr[frag] > 2. * Phys_const::kcal)
      step[frag] = 0.4;
    else if(frag_corr[frag] > 1. * Phys_const::kcal)
      step[frag] = 0.5;
    else if(frag_corr[frag] > 0.5 * Phys_const::kcal)
      step[frag] = 0.6;
    else if(frag_corr[frag] > 0.2 * Phys_const::kcal)
      step[frag] = 0.7;
    else if(frag_corr[frag] > 0.1 * Phys_const::kcal)
      step[frag] = 0.8;
    else
      step[frag] = step_max;

  int corr_ind = 0;
  while(1) {//iteration cycle

    // adjustable uncertainties
    for(int i = 0; i < range_num; ++i)
      if(opt_ener < max_ener_range[i] * Phys_const::kcal) {
	min_corr = min_corr_range[i] * Phys_const::kcal;
	corr_num = corr_num_range[i];
	break;
      }

    // stop condition
    if(tot_corr < min_corr || corr_ind > corr_num) {
      if(tot_corr < min_corr)
	Log::out << funame << "energy correction calculation "
	  "successfully converged after " << corr_ind << " iterations.\n";
      else
	Log::out << funame << "maximum number of iterations ("  
		 << corr_num << ") has been reached.\n";

      Log::out << funame << "relaxation energy per normal mode (kcal/mol):\n";
      for(int frag = 0; frag < 2; ++frag) {
	Log::out << "   " << frag << "-th fragment:";
	for(int mode = 0; mode < q.size(frag); ++mode)
	  Log::out << " " << q(frag, mode)*q(frag, mode)/2./Phys_const::kcal;
	Log::out << "\n";
      }

      
      // projected minimum configuration
      for(int frag = 0; frag < 2; ++frag)
	for(int mode = 0; mode < q.size(frag); ++mode)
	  q_new(frag, mode) = q(frag, mode) + step[frag] 
	    * q_corr(frag, mode);
      nm2pos(q_new, pos);

      Log::out << funame << "projected configuration:\n";
      Log::out << funame;
      print_dist(Log::out, pos);
      Log::out << funame;
      print_geom(Log::out, pos);

      opt_ener -= tot_corr;
      Log::out << funame << "projected interaction energy = "
	       << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
      return opt_ener;
    }
    
    // iteration index 
    ++corr_ind;
    Log::out << funame << corr_ind << "-th ITERATION:\n";

    // new energy calculation
    try {
      // new cartesian geometry
      for(int frag = 0; frag < 2; ++frag)
	for(int mode = 0; mode < q.size(frag); ++mode)
	  q_new(frag, mode) = q(frag, mode) + step[frag] 
	    * q_corr(frag, mode);
      nm2pos(q_new, pos);
      // potential and forces
      opt_ener_new = qc_method(pos, OPT_READ, force);
      // relaxation correction
      for (int at = 0; at < atoms.size(); ++at)
	for (int i = 0; i < 3; ++i)
	  atoms[at].force[i] = force[at][i];
      tot_corr_new = 0.;
      for(int frag = 0; frag < 2; ++frag) {
	frag_corr_new[frag] = 0.;
	for(int mode = 0; mode < q.size(frag); ++mode) {
	  dtemp  = mol_array[frag]->nm_force(mode);
	  q_corr_new(frag, mode) = dtemp;
	  frag_corr_new[frag] += dtemp * dtemp;
	}
	frag_corr_new[frag] /= 2.;
	tot_corr_new += frag_corr_new[frag];
      }

      Log::out << funame << "preliminary interaction energy = " 
	       << opt_ener_new / Phys_const::kcal << " kcal/mol\n";
      Log::out << funame << "preliminary relaxation energy = "
	       << tot_corr_new / Phys_const::kcal << " kcal/mol\n";
    } catch(Pot_error) {
      Log::out << funame << "optimization method failed\n";
      Log::out << funame << "projected configuration:\n";
      Log::out << funame; print_dist(Log::out, pos);
      Log::out << funame; print_geom(Log::out, pos);
      opt_ener -= tot_corr;
      Log::out << funame << "projected interaction energy = "
	       << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
      return opt_ener;
    }

    // new step calculation
    bool is_step_changed = false;
    for(int frag = 0; frag < 2; ++frag) {
      double self = 0.;
      double cross = 0.;
      for(int mode = 0; mode < q.size(frag); ++mode) {
	self  += q_corr(frag, mode) * q_corr(frag, mode);
	cross += q_corr(frag, mode) * q_corr_new(frag, mode);
      }
      cross /= self;
      if(cross > .66) {
	step[frag] *= 3.;
	if(step[frag] > step_max)
	  step[frag] = step_max;
	is_step_changed = true;
      }
      else if (abs(cross) > 0.3) {
	step[frag] /= 1. - cross;
	if(step[frag] > step_max)
	  step[frag] = step_max;
	is_step_changed = true;
      }
    }

    // energy calculation with a new step
    if(is_step_changed) {// is_step_changed
      Log::out << funame << "adjusting the step:\n";
      for(int frag = 0; frag < 2; ++frag)
	Log::out << funame << "new step for " << frag 
		 << "-th fragment = " << step[frag] 
		 << "\n";
      try{
	for(int frag = 0; frag < 2; ++frag)
	  for(int mode = 0; mode < q.size(frag); ++mode)
	    q_new(frag, mode) = q(frag, mode) + step[frag] 
	      * q_corr(frag, mode);
	nm2pos(q_new, pos);
	opt_ener_new = qc_method(pos, OPT_READ, force);
	// relaxation correction
	for (int at = 0; at < atoms.size(); ++at)
	  for (int i = 0; i < 3; ++i)
	    atoms[at].force[i] = force[at][i];
	tot_corr_new = 0.;
	for(int frag = 0; frag < 2; ++frag) {
	  frag_corr_new[frag] = 0.;
	  for(int mode = 0; mode < q.size(frag); ++mode) {
	    dtemp  = mol_array[frag]->nm_force(mode);
	    q_corr_new(frag, mode) = dtemp;
	    frag_corr_new[frag] += dtemp * dtemp;
	  }
	  frag_corr_new[frag] /= 2.;
	  tot_corr_new += frag_corr_new[frag];
	}
	Log::out << funame << "preliminary interaction energy = " 
		 << opt_ener_new / Phys_const::kcal << " kcal/mol\n";
	Log::out << funame << "preliminary relaxation energy = "
		 << tot_corr_new / Phys_const::kcal << " kcal/mol\n";
      } catch(Pot_error) {
	Log::out << funame << "optimization method failed\n";
	Log::out << funame << "projected configuration:\n";
	Log::out << funame; print_dist(Log::out, pos);
	Log::out << funame; print_geom(Log::out, pos);
	opt_ener -= tot_corr;
	Log::out << funame << "projected interaction energy = "
		 << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
	return opt_ener;
      }
    }// is_step_changed


    int count = 0;
    while(opt_ener_new > opt_ener) {
      if(count++ == count_max) {
	Log::out << funame << "too many steps, bailing out\n\n";
	Log::out << funame << "projected configuration:\n";
	Log::out << funame; print_dist(Log::out, pos);
	Log::out << funame; print_geom(Log::out, pos);
	opt_ener -= tot_corr;
	Log::out << funame << "projected interaction energy = "
		 << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
	return opt_ener;
      }
      Log::out << funame << "preliminary interaction  energy "
	"is higher than the old one: reducing the step size\n";
      for(int frag = 0; frag < 2; ++frag) {
	step[frag] /= 2.;
	Log::out << funame << "new step for " << frag 
		 << "-th fragment is " << step[frag] 
		 << "\n";
      }
      Log::out << "\n";

      try{
	for(int frag = 0; frag < 2; ++frag)
	  for(int mode = 0; mode < q.size(frag); ++mode)
	    q_new(frag, mode) = q(frag, mode) + step[frag] 
	      * q_corr(frag, mode);
	nm2pos(q_new, pos);
	opt_ener_new = qc_method(pos, OPT_READ, force);
	// relaxation correction
	for (int at = 0; at < atoms.size(); ++at)
	  for (int i = 0; i < 3; ++i)
	    atoms[at].force[i] = force[at][i];
	tot_corr_new = 0.;
	for(int frag = 0; frag < 2; ++frag) {
	  frag_corr_new[frag] = 0.;
	  for(int mode = 0; mode < q.size(frag); ++mode) {
	    dtemp  = mol_array[frag]->nm_force(mode);
	    q_corr_new(frag, mode) = dtemp;
	    frag_corr_new[frag] += dtemp * dtemp;
	  }
	  frag_corr_new[frag] /= 2.;
	  tot_corr_new += frag_corr_new[frag];
	}
	Log::out << funame << "preliminary interaction energy = " 
		 << opt_ener_new / Phys_const::kcal << " kcal/mol\n";
	Log::out << funame << "preliminary relaxation energy = "
		 << tot_corr_new / Phys_const::kcal << " kcal/mol\n";
      } catch(Pot_error) {
	Log::out << funame << "optimization method failed\n";
	Log::out << funame << "projected configuration:\n";
	Log::out << funame; print_dist(Log::out, pos);
	Log::out << funame; print_geom(Log::out, pos);
	opt_ener -= tot_corr;
	Log::out << funame << "projected interaction energy = "
		 << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
	return opt_ener;
      }
    }

    opt_ener = opt_ener_new;
    for(int frag = 0; frag < 2; ++frag)
      frag_corr[frag] = frag_corr_new[frag];
    tot_corr = tot_corr_new;
    q = q_new;
    q_corr = q_corr_new;

    Log::out << funame << "new configuration:\n";
    Log::out << funame; print_dist(Log::out, pos);
    Log::out << funame; print_geom(Log::out, pos);

    Log::out << funame << "new interaction energy = " 
	     << opt_ener / Phys_const::kcal << " kcal/mol\n\n";
    for(int frag = 0; frag < 2; ++frag)
      Log::out << funame  <<"new relaxation energy for "
	       << frag << "-th fragment = "
	       <<  frag_corr[frag] / Phys_const::kcal 
	       << " kcal/mol" << endl;

    // correction to the step
    for(int frag = 0; frag < 2; ++frag)
      if(step[frag] < 0.25)
	step[frag] = 0.25;

  }// iteration cycle
}


/****** wrapper for a Klippenstein fortran potential form ******/

void sjk_pot (int flags, Array<double>& result)
{
  const char funame [] = "sjk_pot: ";

  const double incr = 1.0e-4;


  if(!Sjk::isinit()) {
    std::cerr << funame << "not initialized\n";
    std::exit(1);
  }
  
  if(result.size() != 1) {
    std::cerr << funame << "PES size should be 1\n";
    std::exit(1);
  }

  // print configuration
  Log::out << "Configuration:\n\n";
  print_geom(Log::out);
  Log::out << endl;
  print_dist(Log::out);

  Array_2<double> rriv(3, atoms.size()); // atomic coordinates

  for(int a = 0; a < atoms.size(); ++a)
    for (int i = 0; i < 3; ++i)
      rriv(i, a) = atoms[a].lf_pos[i];

  const double* rr = rriv.data();

  if(Sjk::pot(rr, result[0])) {
    Log::out << funame << "potential energy failure" << endl; 
    throw Pot_error();
  }

  // print energy
  Log::out << "Energy = " << result[0] / Phys_const::kcal << " kcal/mol\n\n" << endl;

  // forces
  if (flags & POT_FRC) {
    double vp, vn;

    for(int a = 0; a < atoms.size(); ++a)
      for (int i = 0; i < 3; ++i) {
	
	rriv(i, a) +=  incr;
	if(Sjk::pot(rr, vp)) {
	  Log::out << funame << "potential energy failure" << endl; 
	  throw Pot_error();
	}
	
	rriv(i, a) -= incr + incr;
	if(Sjk::pot(rr, vn)) {
	  Log::out << funame << "potential energy failure" << endl; 
	  throw Pot_error();
	}

	atoms[a].force[i] = (vn - vp) / incr / 2.;
	  
	rriv(i, a) += incr;
      }
  }// forses
}

/*
void sjk_pot (int flags, Array<double>& result) throw (Pot_error)
{
  const char funame [] = "sjk_pot: ";

  if(!Sjk::isinit()) {
    cout << funame << "not initialized\n";
    exit(1);
  }
  
  if(result.size() != 1) {
    std::cerr << funame << "wrong number of potential energy surfaces (" 
	      << result.size() << ") requested\n ";
    exit(1);
  }

  const double incr = 1.0e-4;

  Array_3<double> rriv(2, 100, 3); // atomic coordinates

  // set input parameters
  for (int mol = 0; mol < 2; ++mol)
    for (int at = 0; at < mol_array[mol]->size(); ++at)
      for (int i = 0; i < 3; ++i)
	rriv(mol, at, i) =(mol_array[mol]->begin() + at)->lf_pos[i];

  // calculate forces
  if (flags & POT_FRC) {
    double vp, vn;
    for (int mol = 0; mol < 2; ++mol)
      for (int at = 0; at < mol_array[mol]->size(); ++at) {
	Ater ap = mol_array[mol]->begin() + at;
	for (int i = 0; i < 3; ++i) {
	  rriv(mol, at, i) = ap->lf_pos[i] + incr;
	  vp = Sjk::pot(rriv.data());
	  rriv(mol, at, i) = ap->lf_pos[i] - incr;
	  vn = Sjk::pot(rriv.data());
	  rriv(mol, at, i) = ap->lf_pos[i];

	  ap->force[i] = (vn - vp) / incr / 2.;
	}
      }
  }
  // calculate potential
  if (Log::more(Log::DEBUG)) {
    Log::out << "Configuration:\n\n";
    print_geom(Log::out);
    Log::out << "\n";
    print_dist(Log::out);
  }

  result[0] = Sjk::pot(rriv.data());

  if(!std::isnormal(result[0])) {
    Log::out << funame << "potential energy is not a valid number" << endl; 
    throw Pot_error();
  }

  if (Log::more(Log::DEBUG))
    Log::out << "Energy = " << result[0] / Phys_const::kcal << " kcal/mol\n\n" << endl;
}
*/
/***********************************************************
 ********************** MOLPRO *****************************
 ***********************************************************/

void molpro_pot (int flags, Array<double>& result)
{
  const char funame [] = "molpro_pot: ";
  const double max_opt_ener = 25. * Phys_const::kcal; // kcal/mol

  const int count_max = 1;
  if(!Molpro::isinit()) {
    std::cout << funame << "not initialized\n";
    exit(1);
  }

  // set atomic coordinates
  if(qc_method) {
    std::vector<D3> r(atoms.size());
    double ener = optimize(r);
    if(qc_method == molpro_opt)
      remove(Molpro::wf_name.c_str());
    if(ener > max_opt_ener) {
      Log::out << funame << "optimization energy  = " 
	       << ener / Phys_const::kcal 
	       << " kcal/mol will be used as a final energy\n";
      for(int i = 0; i < result.size(); ++i)
	result[i] = ener;
      return;
    }
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
	Molpro::atoms[at].pos[i] = r[at][i];
  }
  else
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
	Molpro::atoms[at].pos[i] = atoms[at].lf_pos[i];

    //  if(Molpro::debug>=Log::DEBUG) {
    Log::out << funame << "coordinates (angstrom):\n\n";
    Molpro::print_geom(Log::out) << "\n";
    Log::out << funame;
    Molpro::print_dist();
    //}

  int molpro_flags = 0;
  if (flags & POT_FRC)
    molpro_flags |= Molpro::FORCE;

  int count;
  for(count = 0; count < count_max; ++count) {
    try {
      Molpro::pot(MOLPRO_ENER, molpro_flags);
    } 
    catch(Molpro::Err) {
      continue;
    }
    break;
  }
  if(count == count_max) {
    Log::out << funame << "potential failed " 
	     << count << " times\n";
    throw Pot_error();
  }
  if (flags & POT_FRC)
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
      atoms[at].force[i] = Molpro::forces[at][i];

  if(qc_forces) {
    if(2 * Molpro::energies.size() != result.size()) {
      Log::out << funame << "energies number mismatch\n";
      throw Pot_error();
    }
    double corr = relax_corr();
    for(int i = 0; i < result.size(); ++i) {
      result[2*i] = Molpro::energies[i] + corr;
      result[2*i+1] = Molpro::energies[i];
    }
  }
  else {
    if(Molpro::energies.size() != result.size()) {
      Log::out << funame << "energies number mismatch\n";
      throw Pot_error();
    }
    for(int i = 0; i < result.size(); ++i) {
      result[i] = Molpro::energies[i]; 
      Log::out << funame << i << "-th energy = " << result[i]/Phys_const::kcal 
	       << " kcal/mol" << std::endl;
    }
  }
}
void moldyn_pot (int flags, Array<double>& result)
{
  const char funame [] = "moldyn_pot: ";

  if(!Molpro::isinit()) {
    std::cout << funame << "not initialized\n";
    exit(1);
  }

  for (int at = 0; at < atoms.size(); ++at)
    for (int i = 0; i < 3; ++i)
      Molpro::atoms[at].pos[i] = atoms[at].lf_pos[i];

  if(Molpro::debug>=Log::DEBUG) {
    Log::out << funame << "coordinates (angstrom):\n\n";
    Molpro::print_geom(Log::out) << "\n";
    Log::out << funame;
    Molpro::print_dist();
  }

  if(flags & POT_FRC) {
    try {
      Molpro::pot(1, Molpro::FORCE|Molpro::WF);
    } 
    catch(Molpro::Err) {
      Log::out << funame << "potential failed\n"; 
      throw Pot_error();
    }
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
      atoms[at].force[i] = Molpro::forces[at][i];
  }
  else {
    remove(Molpro::wf_name.c_str());
    try {
      Molpro::pot(0, Molpro::WF);
    } 
    catch(Molpro::Err) {
      Log::out << funame << "potential failed\n"; 
      throw Pot_error();
    }
  }

  if(Molpro::energies.size() != result.size()) {
    Log::out << funame << "energies number mismatch\n";
    throw Pot_error();
  }
  for(int i = 0; i < result.size(); ++i)
    result[i] = Molpro::energies[i]; 
  
}
/**************************************************************
 *****************     G98 Potentails     *********************
 **************************************************************/

void g98_min_pot (int g98_flags, Array<double>& result)
{
  const char funame [] = "g98_min_pot: ";

  double dtemp;
  
  // set atomic coordinates
  if(qc_method) {
    std::vector<D3> r(atoms.size());
    optimize(r);
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
	Gauss::cluster[at].pos[i] = r[at][i];
  }
  else
    for(int at = 0; at < atoms.size(); ++at)
      for(int i = 0; i < 3; ++i)
	Gauss::cluster[at].pos[i] = atoms[at].lf_pos[i];

  Log::out.flush();

  if(Gauss::debug >= Log::INFO)
    //
    print_dist(Log::out);

  // set calculate-forces flag
  //
  int flags = 0;
  
  if (g98_flags & POT_FRC)
    flags |= Gauss::FORCE;
  
  try {
    //
    result[0] = Gauss::min_pot(flags)->rel_energy();

    Log::out << "gaussian   energy = " << result[0] / Phys_const::kcal << " kcal/mol\n";
  }
  catch (Gauss::Err) {
    //
    throw Pot_error();
  }

  // set forces
  //
  if (g98_flags & POT_FRC)
    for (int at = 0; at < atoms.size(); ++at)
      for (int i = 0; i < 3; ++i)
	atoms[at].force[i] = Gauss::cluster[at].force[i];

  // potential correction
  //
  if(Gauss::ispcorr()) {
    //
    Array_2<double> coord(3, atoms.size());
    
    for (int at = 0; at < atoms.size(); ++at)
      //
      for (int i = 0; i < 3; ++i)
	//
	coord(i, at) = atoms[at].lf_pos[i];

    Array_2<double> fcorr(3, atoms.size());
    
    double* dp = g98_flags & POT_FRC ? fcorr.data() : 0;

    dtemp = Gauss::pcorr(atoms.size(), coord.data(), dp);
    
    Log::out << "correction energy = " << dtemp / Phys_const::kcal << " kcal/mol\n";
    
    result[0] += dtemp;
    
    Log::out << "total      energy = " << result[0] / Phys_const::kcal << " kcal/mol" << std::endl;

    if(g98_flags & POT_FRC)
      //
      for (int at = 0; at < atoms.size(); ++at)
	//
	for (int i = 0; i < 3; ++i)
	  //
	  atoms[at].force[i] += fcorr(i, at);
  }
  
  // relaxation correction
  if(qc_forces) {
    //
    if(result.size() != 2) {
      //
      Log::out << funame << "energies number mismatch" << std::endl;

      throw Pot_error();
    }
    
    result[1] = result[0];
    
    result[0] += relax_corr();
  }
}

// gaussian potential with one-step relaxation correction
void g98_corr_pot (int g98_flags, Array<double>& result)
{
  const char funame [] = "g98_corr_pot: ";

  if(result.size() != 2) {
    //
    Log::out << funame << "energies number mismatch\n";
    
    throw Pot_error();
  }

  double dtemp;

  // set atomic coordinates
  for(int at = 0; at < atoms.size(); ++at)
    //
    for(int i = 0; i < 3; ++i)
      //
      Gauss::cluster[at].pos[i] = atoms[at].lf_pos[i];

  if(Gauss::debug >= Log::INFO)
    //
    print_dist(Log::out);
  
  Log::out.flush();


  try {
    //
    result[0] = Gauss::min_pot(Gauss::FORCE)->rel_energy();
  }
  catch (Gauss::Err) {
    //
    throw Pot_error();
  }

  // potential correction
  //
  if(Gauss::ispcorr()) {
    //
    Array_2<double> coord(3, atoms.size());

    for (int at = 0; at < atoms.size(); ++at)
      //
      for (int i = 0; i < 3; ++i)
	//
	coord(i, at) = atoms[at].lf_pos[i];
    
    result[0] += Gauss::pcorr(atoms.size(), coord.data(), 0);
  }

  // set forces
  //
  for (int at = 0; at < atoms.size(); ++at)
    for (int i = 0; i < 3; ++i)
      atoms[at].force[i] = Gauss::cluster[at].force[i];

  // correction
  double corr = 0.;
  for(int frag = 0; frag < 2; ++frag)
    for(int mode = 0; mode < mol_array[frag]->nm_size(); ++mode) {
      dtemp = mol_array[frag]->nm_force(mode);
      corr -= dtemp * dtemp;
    }
  corr /= 2.;
  Log::out << funame << "one-step relaxation correction = " 
	     << corr * Phys_const::kcal << "\n";
  result[1] = result[0];
  result[0] += corr;
}

// potential for trajectory calculations
void g98_dyn_pot (int g98_flags, Array<double>& result)
{
  static const char funame [] = "g98_dyn_pot: ";

  static Gauss::Method* g98_dyn_method = 0;

  // set atomic coordinates
  for (int at = 0; at < atoms.size(); ++at)
    for (int i = 0; i < 3; ++i)
      Gauss::cluster[at].pos[i] = atoms[at].lf_pos[i];

  Log::out.flush();

  if (Gauss::debug >= Log::INFO)
    print_dist(Log::out);

  // read checkpoint file
  //
  if (g98_flags & POT_CHK) {
    //
    if(!g98_dyn_method) {
      std::cerr << funame << "method should not be zero, exitting" << std::endl;
      
      exit(1);
    }
    
    g98_dyn_method->init();

    // method failed
    //
    if(!g98_dyn_method->apply(Gauss::FORCE | Gauss::READ | Gauss::DIRECT)) {
      //
      Log::out << funame << "Oops!!! current method " << g98_dyn_method->id() 
	       << " failed, looking for the minimal one" << std::endl;
      
      try {
	//
	g98_dyn_method = Gauss::min_pot(Gauss::FORCE);
      }
      catch (Gauss::Err) {
	//
	Log::out << funame << "Oops!!! no method found, trajectory propagation aborted" << std::endl;
	throw Pot_error();
      }
      
      Log::out << funame << "new method is " 
	       << g98_dyn_method->id() << std::endl;
      //
    }//
    // method succeeded
    //
    else {
      //
      if (Gauss::debug >= Log::INFO)
	//
	g98_dyn_method->print_info(Log::out);
      
      try {
	//
	std::vector<D3> forces(atoms.size());
    
	g98_dyn_method->read_forces(forces);

	for(int at = 0; at < atoms.size(); ++at)
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    Gauss::cluster[at].force[i] = forces[at][i];
      }
      catch(Gauss::Err) {
	//
	Log::out << funame << "Oops!!! cannot read forces, "
	  " trajectory propagation aborted" << std::endl;
	
	throw Pot_error();
      }//
      //
    }//
    //
  }//
  //
  // start new trajectory (no checkpoint file)
  //
  else {
    //
    Log::out << funame << "start new trajectory propagation" << std::endl;
    try {
      //
      g98_dyn_method = Gauss::min_pot(Gauss::FORCE);
    }
    catch (Gauss::Err) {
      //
      Log::out << funame << "Oops!!! no method found to start, "
	" trajectory propagation aborted" << std::endl;
      
      throw Pot_error();
    }
    
    Log::out << funame << "new method is " << g98_dyn_method->id() << std::endl;
  }

  result[0] = g98_dyn_method->rel_energy();

  // set forces
  //
  for (int at = 0; at < atoms.size(); ++at)
    //
    for (int i = 0; i < 3; ++i)
      //
      atoms[at].force[i] = Gauss::cluster[at].force[i];

  // potential correction
  //
  if(Gauss::ispcorr()) {
    //
    Array_2<double> coord(3, atoms.size());
    
    for (int at = 0; at < atoms.size(); ++at)
      //
      for (int i = 0; i < 3; ++i)
	//
	coord(i, at) = atoms[at].lf_pos[i];

    Array_2<double> fcorr(3, atoms.size());
    
    double* dp = g98_flags & POT_FRC ? fcorr.data() : 0;
    
    result[0] += Gauss::pcorr(atoms.size(), coord.data(), dp);

    if(g98_flags & POT_FRC)
      //
      for (int at = 0; at < atoms.size(); ++at)
	//
	for (int i = 0; i < 3; ++i)
	  //
	  atoms[at].force[i] += fcorr(i, at);
  }
}

void g98_lin_pot (int g98_flags, Array<double>& result)
{
  // set atomic coordinates
  //
  for (int at = 0; at < atoms.size(); ++at)
    //
    for (int i = 0; i < 3; ++i)
      //
      Gauss::cluster[at].pos[i] = atoms[at].lf_pos[i];

  Log::out.flush();

  if (Gauss::debug >= Log::INFO)
    //
    print_dist(Log::out);

  // set calculate-forces flag
  //
  int flags = 0;
  
  if (g98_flags & POT_FRC)
    //
    flags |= Gauss::FORCE;

  try {
    //
    result[0] = Gauss::inter_pot(flags);
  }
  catch (Gauss::Err) {
    //
    throw Pot_error();
  }

  if (g98_flags & POT_FRC)
    //
    for (int at = 0; at < atoms.size(); ++at)
      //
      for (int i = 0; i < 3; ++i)
	//
	atoms[at].force[i] = Gauss::cluster[at].force[i];

  // potential correction
  //
  if(Gauss::ispcorr()) {
    //
    Array_2<double> coord(3, atoms.size());
    
    for (int at = 0; at < atoms.size(); ++at)
      //
      for (int i = 0; i < 3; ++i)
	//
	coord(i, at) = atoms[at].lf_pos[i];

    Array_2<double> fcorr(3, atoms.size());
    
    double* dp = g98_flags & POT_FRC ? fcorr.data() : 0;

    result[0] += Gauss::pcorr(atoms.size(), coord.data(), dp);

    if(g98_flags & POT_FRC)
      //
      for (int at = 0; at < atoms.size(); ++at)
	//
	for (int i = 0; i < 3; ++i)
	  //
	  atoms[at].force[i] += fcorr(i, at);
  }
}

/******************************************************************
 *                    Dummy potential
 ******************************************************************/

void zero_pot (int flags, Array<double>& result)
{
  result[0] = 0.;
  if(flags & POT_FRC) {
    for(int at = 0; at < atoms.size(); ++at)
      for(int i = 0; i < 3; ++i)
	atoms[at].force[i] =0.;
  }
}

void dummy_pot(int flags, Array<double>& result)
{


  zero_pot(flags, result);

  // randomly send different signals
  struct timeval tv;
  gettimeofday(&tv, 0);
  if(tv.tv_usec % 100 == 0) {
    throw  Pot_error();
  }
}

/*************************************************************
 ****************** Multipole Potential *********************
 *************************************************************/

void multi_pot (int flags, Array<double>& result)
{
  static double mf_n[2][3];
  static bool is_first = true;
  if(is_first) {
    is_first = false;
    for(int frag = 0; frag < 2; ++frag) {
      for(int i = 0; i < 3; ++i)
	mf_n[frag][i] = 0.;
      if(mol_array[frag]->type() == NONLINEAR) {
	if((mol_array[frag]->iner_mom(1) - mol_array[frag]->iner_mom(0))
	   / mol_array[frag]->iner_mom(2) > 0.01 &&
	   (mol_array[frag]->iner_mom(2) - mol_array[frag]->iner_mom(1))
	   / mol_array[frag]->iner_mom(2) > 0.01)
	  std::cout << "multi_pot: WARNING!: " << frag 
		    << "-th molecule does not look symmetric\n";
	if(mol_array[frag]->iner_mom(1) - mol_array[frag]->iner_mom(0) >
	   mol_array[frag]->iner_mom(2) - mol_array[frag]->iner_mom(1)) {
	  mf_n[frag][0] = 1.;
	  std::cout << "multi_pot: symmetry axis is assumed to be "
	    "directed along the X axis\n";
	} 
	else {
	  mf_n[frag][2] = 1.;
	  std::cout << "multi_pot: symmetry axis is assumed to be "
	    "directed along the Z axis\n";
	}
      }
    }
  }

  double dtemp;
  double lf_n[2][3];
  for(int frag = 0; frag < 2; ++frag)
    switch(mol_array[frag]->type()) {
    case LINEAR:
      for(int i = 0; i < 3; ++i)
	lf_n[frag][i] = mol_array[frag]->read_ang_pos()[i];
      ::normalize(lf_n[frag], 3);
      break;
    case NONLINEAR:
      mol_array[frag]->mf2lf(mf_n[frag], lf_n[frag]);
      break;
    case ATOM:
      for(int i = 1; i < 3; ++i)
        lf_n[frag][i] = 0.;
      lf_n[frag][0] = 1.;
    }

  double cm_n[3];
  for(int i = 0; i < 3; ++i)
    cm_n[i] = mol_array[1]->read_cm_pos()[i] - mol_array[0]->read_cm_pos()[i];
  double r = ::normalize(cm_n, 3);

  // cos_theta1, cos_theta_2, sin_theta_1 * sin_theta_2 * cos_phi
  double cos_theta[2];
  cos_theta[0] = vdot(lf_n[0], cm_n, 3);
  cos_theta[1] = -vdot(lf_n[1], cm_n, 3);
  double phi = vdot(lf_n[0], lf_n[1], 3)  + cos_theta[0] * cos_theta[1];

  result[0] = 0.;
  double v0_fac, ang_fac, r_fac;
  // r^2 terms
  r_fac = r * r;
  if(mol_array[0]->is_charge() && mol_array[1]->is_dipole()) {
    // charge-dipole
    v0_fac = mol_array[0]->charge() * mol_array[1]->dipole();
    ang_fac = cos_theta[1];
    result[0] += v0_fac * ang_fac / r_fac;
  }
  if(mol_array[1]->is_charge() && mol_array[0]->is_dipole()) {
    // charge-dipole
    v0_fac = mol_array[1]->charge() * mol_array[0]->dipole();
    ang_fac = cos_theta[0];
    result[0] += v0_fac * ang_fac / r_fac;
  }

  // r^3 terms
  r_fac = r * r * r;
  if(mol_array[0]->is_dipole() && mol_array[1]->is_dipole()) {
    // dipole-dipole
    v0_fac = mol_array[0]->dipole() * mol_array[1]->dipole();
    ang_fac = phi + 2. * cos_theta[0] * cos_theta[1];
    result[0] += v0_fac * ang_fac / r_fac;
  }
  if(mol_array[0]->is_charge() && mol_array[1]->is_quad()) {
    // charge-quadrupole
    v0_fac = mol_array[0]->charge() * mol_array[1]->quad();
    ang_fac = (3. * cos_theta[1] * cos_theta[1] - 1.)/4.;
    result[0] += v0_fac * ang_fac / r_fac;
  }
  if(mol_array[1]->is_charge() && mol_array[0]->is_quad()) {
    // charge-quadrupole
    v0_fac = mol_array[1]->charge() * mol_array[0]->quad();
    ang_fac = (3. * cos_theta[0] * cos_theta[0] - 1.)/4.;
    result[0] += v0_fac * ang_fac / r_fac;
  }

  // r^4 terms
  r_fac *= r;
  if(mol_array[0]->is_dipole() && mol_array[1]->is_quad()) {
    // dipole-quadrupole
    v0_fac = mol_array[0]->dipole() * mol_array[1]->quad();
    ang_fac = 3. * cos_theta[0] * cos_theta[1] * cos_theta[1]
      + 2. * phi * cos_theta[1] - cos_theta[0];
    result[0] += 0.75 * v0_fac * ang_fac / r_fac;
  }
  if(mol_array[1]->is_dipole() && mol_array[0]->is_quad()) {
    // quadrupole-dipole
    v0_fac = mol_array[1]->dipole() * mol_array[0]->quad();
    ang_fac = 3. * cos_theta[1] * cos_theta[0] * cos_theta[0]
      + 2. * phi * cos_theta[0] - cos_theta[1];
    result[0] += 0.75 * v0_fac * ang_fac / r_fac;
  }
  if(mol_array[0]->is_charge() && mol_array[1]->is_polar()) {
    //charge-induced dipole
    v0_fac = - mol_array[0]->charge() * mol_array[0]->charge() 
      * mol_array[1]->polar() / 2.;
    result[0] += v0_fac / r_fac;
  }
  if(mol_array[1]->is_charge() && mol_array[0]->is_polar()) {
    //charge-induced dipole
    v0_fac = - mol_array[1]->charge() * mol_array[1]->charge() 
      * mol_array[0]->polar() / 2.;
    result[0] += v0_fac / r_fac;
  }

  // r^5 terms
  r_fac *= r;
  if(mol_array[0]->is_quad() && mol_array[1]->is_quad()) {
    //quadrupole-quadrupole
    v0_fac = mol_array[0]->quad() * mol_array[1]->quad();
    dtemp = phi + 4. *  cos_theta[1] * cos_theta[0];
    ang_fac = 1. - 5. * cos_theta[0] * cos_theta[0]
      - 5. * cos_theta[1] * cos_theta[1] 
      - 15. * cos_theta[1] * cos_theta[1] * cos_theta[0] * cos_theta[0]
      + 2. * dtemp * dtemp;
    result[0] += 0.1875 * v0_fac * ang_fac / r_fac;
  }

  // r^6 terms
  r_fac *= r;
  if(mol_array[0]->is_dipole() && mol_array[1]->is_polar()) {
    //dipole-induced dipole
    v0_fac = - mol_array[0]->dipole() * mol_array[0]->dipole() 
      * mol_array[1]->polar();
    ang_fac = 0.5 + 1.5 * cos_theta[0] * cos_theta[0];
    result[0] += v0_fac * ang_fac / r_fac;
  }
  if(mol_array[1]->is_dipole() && mol_array[0]->is_polar()) {
    //induced dipole-dipole
    v0_fac = - mol_array[1]->dipole() * mol_array[1]->dipole() 
      * mol_array[0]->polar();
    ang_fac = 0.5 + 1.5 * cos_theta[1] * cos_theta[1];
    result[0] += v0_fac * ang_fac / r_fac;
  }
  if(mol_array[0]->is_polar() && mol_array[1]->is_polar() &&
     mol_array[0]->is_ip() && mol_array[1]->is_ip()) {
    //dispersion
    v0_fac = - 1.5 * mol_array[0]->polar() * mol_array[1]->polar()
      * mol_array[0]->ip() * mol_array[1]->ip()
      / (mol_array[0]->ip() + mol_array[1]->ip());
    result[0] += v0_fac / r_fac;
  }
}

// charge-dipole potential
//
void cq_pot (int flags, Array<double>& result)
{
  static const char funame[] = "cq_pot: ";

  if(!Multipole::is_init()) {
    cout << funame << "multipole potential is not initialized\n";
    exit(1);
  }

  if(mol_array[0]->type() != ATOM) {
    cout << funame << "wrong fragment type\n";
    exit(1);
  }

  double dtemp;

  double lf_n1[3];// quadrupole direction
  switch(mol_array[1]->type()) {
  case ATOM:
    cout << funame << "wrong fragment type\n";
    exit(1);

  case LINEAR:
    for(int i = 0; i < 3; ++i)
      lf_n1[i] = mol_array[1]->read_ang_pos()[i];
    ::normalize(lf_n1, 3);
    break;

  case NONLINEAR:
    mol_array[1]->mf2lf(Multipole::pos(1), lf_n1);
    break;
  };

  double lf_n12[3];
  for(int i = 0; i < 3; ++i)
    lf_n12[i] = mol_array[1]->read_cm_pos()[i] 
      - mol_array[0]->read_cm_pos()[i];
  double r12 = ::normalize(lf_n12, 3);

  const double cos_1  = vdot(lf_n12, lf_n1, 3);
  const double fac = Multipole::coupling()/r12/r12/r12;
  
  result[0] = -fac * (3. * cos_1 * cos_1 - 1.) / 2.;// potential

  if (flags & POT_FRC) {// forces

    double lf_torque[3];
    vector_product(lf_n12, lf_n1, lf_torque);
    dtemp = 3. * fac * cos_1;
    for(int i = 0; i < 3; ++i)
      lf_torque[i] *= dtemp;

    double lf_force[3];
    dtemp = fac * 3. / 2. / r12;
    for(int i = 0; i < 3; ++i)
      lf_force[i] = dtemp * (cos_1 * (5. * lf_n12[i] - 2. * lf_n1[i])
			     - lf_n12[i]);

    // convert torque and force to atomic forces
    for(int i = 0; i < 3; ++i)
      mol_array[0]->begin()->force[i] = -lf_force[i];

    switch(mol_array[1]->type()) {// fragment type
    case LINEAR:
      for(Ater at = mol_array[1]->begin(); at != mol_array[1]->end(); ++at) {
	vector_product(lf_torque, at->rel_pos, at->force);
	for(int i = 0; i < 3; ++i)
	  at->force[i] *= at->mass() / mol_array[1]->iner_mom(2);
      }
      break;

    case NONLINEAR:

      double mf_x [3];
      double lf_x [3];

      mol_array[1]->lf2mf(lf_torque, mf_x);
      for(int i = 0; i < 3; ++i)
	mf_x[i] /= mol_array[1]->iner_mom(i);
      mol_array[1]->mf2lf(mf_x, lf_x);

      for(Ater at = mol_array[1]->begin(); at != mol_array[1]->end(); ++at) {
	vector_product(lf_x, at->rel_pos, at->force);
	for(int i = 0; i < 3; ++i)
	  at->force[i] *= at->mass();
      }
    }// fragment type

    for(Ater at = mol_array[1]->begin(); at != mol_array[1]->end(); ++at)
      for(int i = 0; i < 3; ++i)
	at->force[i] += at->mass() * lf_force[i] / mol_array[1]->mass();
  }
}

// dipole-dipole potential
//
void dd_pot (int flags, Array<double>& result)
{
  static const char funame[] = "dd_pot: ";

  if(!Multipole::is_init()) {
    cout << funame << "multipole potential is not initialized\n";
    exit(1);
  }

  double lf_n[2][3];// dipole directions
  for(int mol = 0; mol < 2; ++mol)
    switch(mol_array[mol]->type()) {
    case ATOM:
      cout << funame << "wrong fragment type\n";
      exit(1);

    case LINEAR:
      for(int i = 0; i < 3; ++i)
	lf_n[mol][i] = mol_array[mol]->read_ang_pos()[i];
      ::normalize(lf_n[mol], 3);
      break;

    case NONLINEAR:
      mol_array[mol]->mf2lf(Multipole::pos(mol), lf_n[mol]);
      break;
    };

  double lf_n0[3];
  for(int i = 0; i < 3; ++i)
    lf_n0[i] = mol_array[1]->read_cm_pos()[i] - mol_array[0]->read_cm_pos()[i];
  double r12 = ::normalize(lf_n0, 3);

  const double cos_1  = vdot(lf_n0,   lf_n[0], 3);
  const double cos_2  = vdot(lf_n0,   lf_n[1], 3);
  const double cos_12 = vdot(lf_n[0], lf_n[1], 3);

  const double r3 = r12*r12*r12;
  result[0] = Multipole::coupling() * (cos_12 - 3.*cos_1*cos_2) / r3;

  if (flags & POT_FRC) {// forces

    double nvp12[3];
    vector_product(lf_n[0], lf_n[1], nvp12);

    double nvp[2][3];
    for(int mol = 0; mol < 2; ++mol)
      vector_product(lf_n[mol], lf_n0, nvp[mol]);

    // torque
    double lf_torque[2][3];
    for(int i = 0; i < 3; ++i) {
      lf_torque[0][i] = Multipole::coupling() *
	(3.*nvp[0][i]*cos_2 - nvp12[i]) / r3;
      lf_torque[1][i] = Multipole::coupling() *
	(3.*nvp[1][i]*cos_1 + nvp12[i]) / r3;
    }
    
    const double r4 = r3*r12 / 3.;

    // force
    double lf_force [3]; // force exerted on the 2nd fragment
    for(int i = 0; i < 3; ++i)
      lf_force[i] = Multipole::coupling() *
	(lf_n0[i]*(cos_12 - 5.*cos_1*cos_2) +
	 lf_n[0][i]*cos_2 + lf_n[1][i]*cos_1) / r4;

    // convert torque and force to atomic forces
    for(int mol = 0; mol < 2; ++mol) {// fragment cycle

      switch(mol_array[mol]->type()) {// fragment type
      case ATOM:
	for(int i = 0; i < 3; ++i)
	  mol_array[mol]->begin()->force[i] = 0.;
	break;
	
      case LINEAR:
	for(Ater at = mol_array[mol]->begin(); 
	    at != mol_array[mol]->end(); ++at) {
	  vector_product(lf_torque[mol], at->rel_pos, at->force);
	  for(int i = 0; i < 3; ++i)
	    at->force[i] *= at->mass() / mol_array[mol]->iner_mom(i);
	}
	break;

      case NONLINEAR:

	double mf_x [3];
	double lf_x [3];

	mol_array[mol]->lf2mf(lf_torque[mol], mf_x);
	for(int i = 0; i < 3; ++i)
	  mf_x[i] /= mol_array[mol]->iner_mom(i);
	mol_array[mol]->mf2lf(mf_x, lf_x);

	for(Ater at = mol_array[mol]->begin(); 
	    at != mol_array[mol]->end(); ++at) {
	  vector_product(lf_x, at->rel_pos, at->force);
	  for(int i = 0; i < 3; ++i)
	    at->force[i] *= at->mass();
	}
      }// fragment type

      for(Ater at = mol_array[mol]->begin(); at != mol_array[mol]->end(); ++at)
	for(int i = 0; i < 3; ++i)
	  if(mol)
	    at->force[i] += at->mass() * lf_force[i] / mol_array[mol]->mass();
	  else
	    at->force[i] -= at->mass() * lf_force[i] / mol_array[mol]->mass();

    }// fragment cycle
  }// forces
}

pot_f name2pot (const string& name)
{
  const char funame [] = "name2pot: ";

  if (name == "sjk")
    return sjk_pot;
  else if(name == "molpro")
    return molpro_pot;
  else if(name == "moldyn")
    return moldyn_pot;
  else if (name == "g98_min")
    return g98_min_pot;
  else if (name == "g98_dyn")
    return g98_dyn_pot;
  else if (name == "g98_corr")
    return g98_corr_pot;
  else if (name == "g98_lin")
    return g98_lin_pot;
  else if (name == "zero")
    return zero_pot;
  else if (name == "multi")
    return multi_pot;
  else {
    cout << funame << "unknown potential name: "
	 << name << ", exitting\n";
    exit(1);
  }
  return 0;
}
