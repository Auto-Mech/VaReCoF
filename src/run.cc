#include <cmath>
#include <cstdio>
#include <setjmp.h>
#include <ctime>

#include <mpi.h>

#include "slatec.hh"
#include "run.hh"
#include "math.hh"
#include "error.hh"
#include "comm.hh"
#include "units.hh"

namespace Traj {// Trajectory namespace

  double rel_tol;
  double abs_tol; 
  double time_step;

  int dist_out_flag;
  int aux_out_flag;
  int anim_out_flag;
  int dvd_out_flag;
  int ang_out_flag;
  int dyn_out_flag;
  int ener_out_flag;

  Div_surf reac_surf;
  Div_surf tran_surf;
  Div_surf diss_surf;

  int reac_ener_flag;
  vector<double> reac_ener;
  int reac_surf_flag;

  //pot_f pot;
  string base_name;

  int save_time_step;
  int debug = 0;
  string debug_mesg = "DEBUG: ";

  struct Par {// data for derivative function
    double ener;
    jmp_buf jmp;
  };

  extern "C" void update_dvd (const double&, const double*, 
				   double*, void*, void*);

}// Trajectory namespace

/*********************************************************************
 ****************************** Res **********************************
 *********************************************************************/

void Traj::Res::send (int node, int tag) const
{
  static const char funame [] = "Traj::Res::send: ";

  MPI::COMM_WORLD.Send(_lf, 2, MPI::INT, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame
	 << tag << "-th tag sent to " 
	 << node << "-th node" << endl;
}

void Traj::Res::recv (int node, int tag)
{
  static const char funame [] = "Traj::Res::recv: ";

  MPI::COMM_WORLD.Recv(_lf, 2, MPI::INT, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame
	 << tag << "-th tag received from " 
	 << node << "-th node" << endl;
}

void Traj::Res::recv (int* node, int* tag)
{
  static const char funame [] = "Traj::Res::recv(int*, int*): ";

  MPI::Status stat;
  MPI::COMM_WORLD.Recv(&_lf, 2, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);

  *tag  = stat.Get_tag();
  *node = stat.Get_source();

  if(Comm::debug)
    cout << Comm::mesg() << funame
	 << *tag  << "-th tag received from " 
	 << *node << "-th node" << endl;
}

istream& Traj::operator>> (istream& from, Res& r)
{
  const char funame [] = "Traj::operator>> (Res&): ";

  from >> r._lf[0] >> r._lf[1];

  if(!from) {
    cout << Comm::mesg() << funame << "input failed\n";
    throw Form_Err();
  }

  return from;
}

ostream& Traj::operator<< (ostream& to, const Res& r)
{


  to << setw(2) << r._lf[0] << " " << setw(2) << r._lf[1];

  return to;
}


/*********************************************************************
 ***************************** State *********************************
 *********************************************************************/

void Traj::State::send (int node, int tag) const
{
  static const char funame [] = "Traj::State::send: ";

  MPI::COMM_WORLD.Send(&time, 1, MPI::DOUBLE, node, tag);
  static_cast<const Res&>(*this).send(node, tag);
  MPI::COMM_WORLD.Send(data(), size(), MPI::DOUBLE, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame
	 << tag << "-th tag sent to " 
	 << node << "-th node" << endl;
}

void Traj::State::recv (int node, int tag)
{
  static const char funame [] = "Traj::State::recv: ";

  MPI::COMM_WORLD.Recv(&time, 1, MPI::DOUBLE, node, tag);
  static_cast<Res&>(*this).recv(node, tag);
  MPI::COMM_WORLD.Recv(data(), size(), MPI::DOUBLE, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame
	 << tag << "-th tag received from " 
	 << node << "-th node" << endl;
}

istream& Traj::operator>> (istream& from, State& rs)
{
  static const char funame [] = "Traj::operator>> (istream&, State&): ";
  from >> rs.time 
       >> static_cast<Res&>(rs) 
       >> static_cast<Dynvar&>(rs);
  if(!from) {
    cout << funame << "input failed\n";
    throw Form_Err();
  }
  return from;
}

ostream& Traj::operator<< (ostream& to, const State& rs)
{
  to << rs.time << "\t" 
     << static_cast<const Res&>(rs) << "\n"
     << static_cast<const Dynvar&>(rs);
  return to;
}

/*********************************************************************
 *********************** Trajectory Propagator ***********************
 *********************************************************************/

extern "C" void Traj::update_dvd (const double& traj_time, const double* dv,
                                 double* dvd, void* rpar,void* ipar)
{// calculates derivatives of dynamic variables

  Par& dvd_par = *static_cast<Par*>(ipar);
  double* dyn_var = const_cast<double*>(dv);

  // calculate atomic coordinates
  for(int mol = 0; mol < 2; ++mol) {
    mol_array[mol]->set_dv(dyn_var); // set dyn. var. pointers
    mol_array[mol]->update_atoms();
  }

  // calculate forces end energy
  Array<double> pot_ener(1);
  try {
    PES::pot(POT_FRC | POT_CHK, pot_ener);
  } catch (Pot_error perr) {
    longjmp(dvd_par.jmp, perr.stat);
  }

  dvd_par.ener = pot_ener[0] + mol_array[0]->kin_energy()
    + mol_array[1]->kin_energy();

  // calculate derivatives of dynamic variables
  for(int mol = 0; mol < 2; ++mol) {
    mol_array[mol]->set_dvd(dvd);
    mol_array[mol]->update_dvd();
  }
}

void Traj::State::run (int num_id)
{
  static const char funame [] = "Traj::State::run: ";

  // temporary variables
  int itemp;


  // save time
  time_t save_time = std::time(0) 
    + (time_t)(save_time_step * Comm::rank() / Comm::size());

  // trajectory id
  char str_id [10];
  sprintf(str_id, "_%i", num_id / 2);

  // Run the trajectory backward and forward
  int traj_ward = num_id % 2;

  // trajectory name
  string traj_name = base_name + str_id;
  if(traj_ward)
    traj_name += "_dis";
  else
    traj_name += "_ass";

  // output streams
  ofstream aux_out;      // auxiliary output
  if (aux_out_flag)
    aux_out.open((traj_name + ".aux").c_str());
   
  ofstream anim_out;    // animated trajectory output
  if (anim_out_flag)
    anim_out.open((traj_name + ".xyz").c_str());
  
  attach();
  update_atoms();

  Array<double> pot_ener_arr(1);
  try {
    PES::pot(0, pot_ener_arr);
  }
  catch(Pot_error) {
    cout << Comm::mesg() << funame 
    	 << "Oops!!! potential calculation failed, "
 	 << num_id << "-th trajectory"
 	 << endl;
    if(aux_out_flag)
      aux_out << "potential calculation failed\n";
    throw Pot_Err();
  }
  
  double pot_en = pot_ener_arr[0];
  double kin_en = mol_array[0]->kin_energy() 
    + mol_array[1]->kin_energy();
  double tot_en = pot_en + kin_en;

  if(tot_en < 0.) {
    if(reac_surf_flag) {// reactive surface mechanism
      set_face(0);
    }
    if(reac_ener_flag) {// reactive energy mechanism
      set_ener_level(reac_ener.size());
      set_face(-1);
    }
    if(aux_out_flag)
      aux_out << "negative total energy: " 
	      << tot_en / Phys_const::kcal 
	      << " kcal/mol\n";
    return;
  }

  // auxiliary output
  if(aux_out_flag)
    aux_out << "\ntime = " << time << "\n\n";

  if(ener_out_flag) // energy output
    aux_out << "energy (kcal/mol): K = " << kin_en / Phys_const::kcal
	    << "  P = " << pot_en / Phys_const::kcal
	    << " T = " <<  tot_en / Phys_const::kcal
	    << "\n\n";

  if(dist_out_flag) {// distance output
    aux_out << "cm-to-cm distance = "
	    << vdistance(mol_array[0]->read_cm_pos(), 
			 mol_array[1]->read_cm_pos(), 3) 
	    << "\n\n";
    print_dist(aux_out);

    aux_out << "transition surface relative distance = "
	    << get_ref_dist(tran_surf) << "\n";
    if(reac_surf_flag)
      aux_out << "reactive surface relative distance = " 
	      << get_ref_dist(reac_surf) << "\n";
    aux_out << "\n";
  }

  if(ang_out_flag) // angular momentum output
    print_ang_mom_info(aux_out);

  if(dyn_out_flag) // dynamical variables output
    aux_out << "Dynamical information:\n"
	    << static_cast<Dynvar&>(*this)
	    << "\n";

  int istep = 0;
  if(anim_out_flag) {// animated trajectory output
    itemp = anim_out.precision(5);
    anim_out << "\t" << atoms.size() << "\n"
	     << " step # : " << istep++ << "\n";
    for (Ater at = atoms.begin(); at != atoms.end(); ++at) {
      anim_out << at->name();
      for (int i = 0; i < 3; ++i)
	anim_out << " " << setw(12)
		 << Phys_const::bohr2ang(at->lf_pos[i]);
      anim_out << "\n";
    }
    anim_out.precision(itemp);
  }

  // integrator parameters
  int info [15];
  info [0] = 0;     // start new trajectory
  info [1] = 0;     // rel_tol and abs_tol are scalars
  info [2] = 0;     // the solution only at the end point
  info [3] = 0;     // the integration can go beyond final point

  const int lrw = 130 + 21 * size();
  Array<double> rwork(lrw);
  const int liw = 51;
  int iwork [liw];

  Par dvd_par;// parameters to pass to update_dvd
  if(setjmp(dvd_par.jmp)) {// long jump
    cout << Comm::mesg() << funame 
	 << "potential calculation failed, " 
	 << num_id << "-th trajectory" << endl;
    if(aux_out_flag)
      aux_out << "potential calculation failed\n";
    throw Pot_Err();
  }

  while(1) {// integration cycle
      
    int idid;
    double timeout;
    if (traj_ward)
      timeout = time + time_step;  // forward
    else
      timeout = time - time_step;  // backward

    // integration step
    ddeabm_(update_dvd, size(), time, data(), timeout, info, 
	    rel_tol, abs_tol, idid, rwork.data(), lrw, iwork, 
	    liw, 0, static_cast<void*>(&dvd_par));
	  
    if (info[0] != 1 || idid < 0) {
      if (aux_out_flag)
	aux_out << "info[0] = " << info[0] << "\n";
      info [0] = 1; // continue the integration
    }
    if (idid < 0 && aux_out_flag)
      aux_out << " idid = " << idid << "\n";

    attach();
    update_atoms();

    // kinetic and potential energy
    kin_en = mol_array[0]->kin_energy() 
      + mol_array[1]->kin_energy();
    pot_en = dvd_par.ener - kin_en;

    double ts_dist = get_ref_dist(tran_surf);
    int rface = -1;
    double rs_dist = 0.;
    if(reac_surf_flag)
      rs_dist = get_ref_dist(reac_surf, &rface);

    double ds_dist = get_ref_dist(diss_surf);

    // auxiliary output
    if(aux_out_flag)
      aux_out << "\ntime = " << time << "\n\n";

    if(ener_out_flag) // energy output
      aux_out << "energy (kcal/mol): K = " << kin_en / Phys_const::kcal
	      << "  P = " << pot_en / Phys_const::kcal
	      << " T = " << dvd_par.ener / Phys_const::kcal
	      << "\n\n";

    if(dist_out_flag) {// distance output
      aux_out << "cm-to-cm distance = "
	      << vdistance(mol_array[0]->read_cm_pos(), 
			   mol_array[1]->read_cm_pos(), 3) 
	      << "\n\n";
      print_dist(aux_out);

      aux_out << "transition surface relative distance = "
	      << ts_dist << "\n";
      if(reac_surf_flag)
	aux_out << "reactive surface relative distance = " 
		<< rs_dist << "\n";
      aux_out << "\n";
    }

    if(ang_out_flag) // angular momentum output
      print_ang_mom_info(aux_out);

    if(dyn_out_flag) // dynamical variables output
      aux_out << "Dynamical information:\n"
	      << static_cast<Dynvar&>(*this)
	      << "\n";
      
    // Animated trajectory output
    if(anim_out_flag) {
      itemp = anim_out.precision(5);
      anim_out << "\t" << atoms.size() << "\n"
	       << " step # : " << istep++ << "\n";
      for (Ater at = atoms.begin(); at != atoms.end(); ++at) {
	anim_out << at->name();
	for (int i = 0; i < 3; ++i)
	  anim_out << " " << setw(12)
		   << Phys_const::bohr2ang(at->lf_pos[i]);
	anim_out << "\n";
      }
      anim_out.precision(itemp);
    }

    // energy level
    if(reac_ener_flag)
      for(; ener_level() < reac_ener.size(); add_level())
	if(pot_en > reac_ener[ener_level()])
	  break;
	    
    // stop conditions
    if(traj_ward && ds_dist > 1.) {
      set_face(-1);
      if(aux_out_flag)
	aux_out << "fragments dissociated\n";
      return;
    }
    if(!traj_ward && ts_dist > 1.) {
      set_face(-1);
      if(aux_out_flag)
	aux_out << "fragments recrossed transition state\n";
      return;
    }
    if(reac_surf_flag && rs_dist < 1.) {// reactive surface mechanism
      set_face(rface);
      if(aux_out_flag)
	aux_out << "fragments reacted\n";
      return;
    }
    if(reac_ener_flag && ener_level() == reac_ener.size()) {
      // reactive energy mechanism
      set_face(-1);
      if(aux_out_flag)
	aux_out << "fragments reacted\n";
      return;
    }

    // send intermediate results
    time_t curr_time = std::time(0);
    if(curr_time > save_time) {
      save_time = curr_time + (time_t)save_time_step;
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, 0, Comm::RUN_TAG);
      if(Comm::debug)
	cout << Comm::mesg() << funame
	     << Comm::RUN_TAG << "-th tag sent to master" << endl;
      send(0, Comm::RUN_TAG);
    }

    // check if normalization || orthogonalization is needed
    bool is_modified = false;
    for(int mol = 0; mol < 2; ++mol) {
      if(!mol_array[mol]->is_ang_normalized()) {
	mol_array[mol]->normalize();
	if(aux_out_flag)
	  aux_out << "Angular orientation of the " << mol 
		  << "-th molecule normalized\n";
	is_modified = true; 
      } 
      if(!mol_array[mol]->is_ang_vel_orthogonal()) {
	mol_array[mol]->orthogonalize();
	if(aux_out_flag)
	  aux_out << "Angular velocity of the " << mol 
		  << "-th molecule was orthogonalized\n";
	is_modified = true;
      }
    }
    if(is_modified)
      info [0] = 0; // start new integration
  }// integration loop
}

