#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>

#include "flux.hh"
#include "raninit.hh"
#include "integral.hh"
#include "ref.hh"
#include "gauss.hh"
#include "units.hh"

/****************************************************************
 **************************** Sampling **************************
 ****************************************************************/

double Samp::energy (int i) const 
{ 
  return ener[i];
}

void Samp::set_val (double w, const Array<double>& e)
{
  wfac = w; 
  rval = Random::flat();
  for(int i=0; i<PES::size(); ++i)
    ener[i]=e[i];
  need_update = true;
}

double Samp::weight (double temper, int en) const
{
  const char funame [] = "Samp::weight: ";

  double dtemp = ener[en] / temper;
  
  if(dtemp < -Limits::exp_arg_max()) {
    //
    std::cerr << funame << "energy, " << ener[en] / Phys_const::kcal
	      << " kcal/mol, is too negative, trancating" << std::endl;
    
    dtemp = -Limits::exp_arg_max();
  }
  else if(dtemp < Limits::exp_arg_max()) {
    //
    _weight = wfac * exp(-dtemp);
  }
  else
    //
    _weight = 0.;

  need_update = false;
  
  return _weight;
}

double Samp::weight () const
{ 
  static const char funame [] = "Samp::weight: ";

  if(need_update) {
    cout << funame << "should be initialized first\n";
    throw Error::Init_Err();
  }

  return _weight;
}

void Samp::send (int node, int tag) const
{
  const char funame [] = "Samp::send: ";

  MPI::COMM_WORLD.Send(&wfac, 1, MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Send(ener.data(), PES::size(), MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Send(&rval, 1, MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Send(data(), size(), MPI::DOUBLE, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 << tag << "-th tag sent to "
	 << node << "-th node" << endl;
}

void Samp::recv (int node, int tag)
{
  const char funame [] = "Samp::recv: ";

  need_update = true;

  MPI::COMM_WORLD.Recv(&wfac, 1, MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Recv(ener.data(), PES::size(), MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Recv(&rval, 1, MPI::DOUBLE, node, tag);
  MPI::COMM_WORLD.Recv(data(), size(), MPI::DOUBLE, node, tag);
  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 << tag << "-th tag received from "
	 << node << "-th node" << endl;
}

ostream& operator<< (ostream& to, const Samp& smp)
{
  for(int i=0; i<PES::size(); ++i)
    to << smp.ener[i] << "  ";
  to << "\n" << smp.wfac <<  "   " << smp.rval << "\n";
  to << static_cast<const Dynvar&>(smp);

  return to;
}

istream& operator>> (istream& from, Samp& smp)
{
  smp.need_update = true;

  for(int i=0; i<PES::size(); ++i)
    from >> smp.ener[i];

  from >> smp.wfac >> smp.rval;
  from >> static_cast<Dynvar&>(smp);

  return from;
}


/*****************************************************************************
 ******************************* Flux Base  **********************************
 *****************************************************************************/

bool FluxBase::_is_stat_init = false;
bool FluxBase::is_stat_init() { return _is_stat_init;}

int FluxBase::t_size = 0;
int FluxBase::e_size = 0;
int FluxBase::j_size = 0;

double* FluxBase::tmpr = 0;
double* FluxBase::ener = 0;
double* FluxBase::amom = 0;

int FluxBase::smp_out_flag = 0; // raw sampling output flag

 
void FluxBase::stat_init (const std::vector<double>& td, 
			  const std::vector<double>& ed, 
			  const std::vector<double>& jd)
{ 
  const char funame [] = "FluxBase::stat_init: ";

  _is_stat_init = true; 

  t_size = td.size();
  tmpr = new double[td.size()];
  for(int i = 0; i < td.size(); ++i)
    tmpr[i] = td[i];

  e_size = ed.size();
  ener = new double[ed.size()];
  for(int i = 0; i < ed.size(); ++i)
    ener[i] = ed[i];
  

  j_size = jd.size();
  amom = new double[jd.size()];
  for(int i = 0; i < jd.size(); ++i)
    amom[i] = jd[i];
}

// sample active surface area
FluxBase::FluxBase(const Div_surf& ds, int n) 
  : _acct_num(0), _fail_num(0), _fake_num(0), _face_num(0), _dist_num(0)
{
  const char funame [] = "FluxBase::FluxBase: ";

  static Dynvar dv;

  for (int i = 0; i < n; ++i) {
    switch(rand_pos(ds, dv, 0))	{
    case 0:
      ++_fake_num;
      break;
    case SAMP_ATOMS_CLOSE:
      ++_dist_num;
      break;
    case SAMP_FACE_OUT:
      ++_face_num;
      break;
    default:
      cout << Comm::mesg() << funame 
	   << "rand_pos's status unknown, exitting\n";
      exit(1);
    }
  }
}

FluxBase& FluxBase::operator+= (const FluxBase& fl)
{
  _acct_num += fl._acct_num;
  _fail_num += fl._fail_num;
  _fake_num += fl._fake_num;
  _face_num += fl._face_num;
  _dist_num += fl._dist_num;

  return *this;
}

bool FluxBase::is_pot () const
{
  const char funame [] = "FluxBase::is_pot (): ";

  if(!tot_smp()) {
    cout << funame << "no samplings" << endl;
    throw Error::Run_Err();
  }

  if(!pot_smp())
    return false;

  if(!acct_smp()) {
    cout << funame <<  "no accepted samplings" << endl;
    throw Error::Run_Err();
  }

  return true;
}

void FluxBase::recv (int src, Comm::tag_t tag)
{


  int smp_data [5];

  MPI::COMM_WORLD.Recv(smp_data, 5, MPI::INT, src, tag);

  _acct_num = smp_data[0];
  _fail_num = smp_data[1];
  _face_num = smp_data[2];
  _dist_num = smp_data[3];
  _fake_num = smp_data[4];
}
    
void FluxBase::send (Comm::tag_t tag) const
{


  int smp_data [5] = {_acct_num, _fail_num, _face_num, _dist_num, _fake_num};

  MPI::COMM_WORLD.Send(smp_data, 5, MPI::INT, 0, tag);
}

ostream& operator<< (ostream& to, const FluxBase& fl)
{
  to << fl._acct_num  << "   " 
     << fl._fail_num  << "   "
     << fl._face_num  << "   " 
     << fl._dist_num  << "   "
     << fl._fake_num  << "\n"; 
  return to;
}

istream& operator>> (istream& from, FluxBase& fl)
{
  const char funame [] = "operator>>(istream&, FluxBase&): "; 

  from >> fl._acct_num >> fl._fail_num >> fl._face_num >> fl._dist_num;
  if(!from) {
    std::cout << funame << "format error\n";
    throw Error::Form_Err();
  }

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);

  if(!(iss >> fl._fake_num))
    fl._fake_num = 0;

  return from;
}

/*****************************************************************************
 *****************************   Thermal Flux   ******************************
 *****************************************************************************/

ThermalFlux::ThermalFlux() : t_sum(tm_size()), t_var(tm_size())
{ 
  if (!is_stat_init()) 
    error("ThermalFlux::ThermalFlux: static variables are not initialized");

  min_en = 0.;
  min_dv.init();

  t_sum.init();
  t_var.init(); 
}

ThermalFlux& ThermalFlux::operator+= (const ThermalFlux& fl)
{
  if(fl.acct_smp() && (!acct_smp() || fl.min_en < min_en)) {
      min_en = fl.min_en;
      min_dv = fl.min_dv;
  }

  FluxBase::operator+=(fl);

  t_sum += fl.t_sum;
  t_var += fl.t_var;

  return *this;
}

ThermalFlux ThermalFlux::operator+ (const ThermalFlux& fl) const
{
  ThermalFlux res(*this);
  return res += fl;
}

void ThermalFlux::normalize ()
{
  static const double min_flux = 1.e-99;

  if(!acct_smp())
    return;

  double dtemp;

  double smp_fac = (double)pot_smp() / (double)tot_smp()
    / (double)acct_smp();

  double fl_fac = 1./acct_smp() + 1./tot_smp() - 1./pot_smp();

  for(int i = 0; i < tm_size(); ++i) {
    double& fl = t_sum[i];
    double& df = t_var[i];

    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }
}
  
double ThermalFlux::average (int t_ind) const
{
  const char funame [] = "ThermalFlux::average: ";

  if(t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  return t_sum[t_ind] / (double)acct_smp()
    * (double)pot_smp() / (double)tot_smp();
}

double ThermalFlux::pot_var (int t_ind) const
{
  const char funame [] = "ThermalFlux::pot_var: ";
  if(t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  double dtemp = t_sum[t_ind] / (double)acct_smp();
  return sqrt(t_var[t_ind] / (double)acct_smp() - dtemp * dtemp)
    * (double) pot_smp() / (double) tot_smp();
}

double ThermalFlux::vol_var (int t_ind) const
{
  const char funame [] = "ThermalFlux::vol_var: ";

  if(t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  return t_sum[t_ind] / (double)acct_smp()
    * sqrt((double) vol_smp() * (double) pot_smp())
    / (double) tot_smp();
}

ostream& operator<< (ostream& to, const ThermalFlux& fl)
{
  to << static_cast<const FluxBase&>(fl)
     << fl.min_en << "\n"
     << fl.min_dv << "\n";

  for(int i = 0; i <  ThermalFlux::tm_size(); ++i)
    to <<  fl.cn_flux(i) << "   " <<  fl.cn_dflx(i) << "\n";
  to << "\n";

  return to;
}

istream& operator>> (istream& from, ThermalFlux& fl)
{
  from >> static_cast<FluxBase&>(fl)
       >> fl.min_en
       >> fl.min_dv;

  for(int i = 0; i <  ThermalFlux::tm_size(); ++i)
    from >>  fl.t_sum[i] >> fl.t_var[i];

  return from;
}

void ThermalFlux::master_recv (int src)
{
  FluxBase::recv(src, Comm::FLUX_TAG);

  MPI::COMM_WORLD.Recv(&min_en, 1, MPI::DOUBLE, src, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Recv(min_dv.data(), min_dv.size(), MPI::DOUBLE, src, 
			 Comm::FLUX_TAG);

  MPI::COMM_WORLD.Recv(t_sum.data(), tm_size(), 
		       MPI::DOUBLE, src, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Recv(t_var.data(), tm_size(), 
		       MPI::DOUBLE, src, Comm::FLUX_TAG);

}

void ThermalFlux::slave_send () const
{


  FluxBase::send(Comm::FLUX_TAG);
  MPI::COMM_WORLD.Send(&min_en, 1, MPI::DOUBLE, 0, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Send(min_dv.data(), min_dv.size(), MPI::DOUBLE, 0, 
		       Comm::FLUX_TAG);

  MPI::COMM_WORLD.Send(t_sum.data(), tm_size(), MPI::DOUBLE, 
		       0, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Send(t_var.data(), tm_size(), MPI::DOUBLE, 
		       0, Comm::FLUX_TAG);

}

/**********************************************************************
 ************************   Thermal MultiFlux   ***********************
 **********************************************************************/

ThermalMultiFlux::ThermalMultiFlux() : t_sum(tm_size(), PES::size()), 
				       t_var(tm_size(), PES::size()),
				       min_en(PES::size()), min_dv(PES::size())
{ 
  if (!is_stat_init()) 
    error("ThermalMultiFlux::ThermalMultiFlux: "
	  "static variables are not initialized");

  min_en.init();
  for(int i = 0; i < PES::size(); ++i)
    min_dv[i].init();

  t_sum.init();
  t_var.init(); 
}

ThermalMultiFlux& ThermalMultiFlux::operator+= (const ThermalMultiFlux& fl)
{
  for(int i = 0; i < PES::size(); ++i)
    if(fl.acct_smp() && (!acct_smp() || fl.min_en[i] < min_en[i])) {
      min_en[i] = fl.min_en[i];
      min_dv[i] = fl.min_dv[i];
    }

  FluxBase::operator+=(fl);

  t_sum += fl.t_sum;
  t_var += fl.t_var;

  return *this;
}

ThermalMultiFlux ThermalMultiFlux::operator+ (const ThermalMultiFlux& fl) const
{
  ThermalMultiFlux res(*this);
  return res += fl;
}

void ThermalMultiFlux::normalize ()
{
  static const double min_flux = 1.e-99;

  if(!acct_smp())
    return;

  double dtemp;

  double smp_fac = (double)pot_smp() / (double)tot_smp()
    / (double)acct_smp();

  double fl_fac = 1./acct_smp() + 1./tot_smp() - 1./pot_smp();

  for(int i = 0; i < t_sum.size(); ++i) {
    double& fl = t_sum[i];
    double& df = t_var[i];
      
    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }
}
  
double ThermalMultiFlux::average (int t_ind, int pes) const
{
  const char funame [] = "ThermalMultiFlux::average: ";

  if (t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  return t_sum(t_ind, pes) / (double)acct_smp()
    * (double)pot_smp() / (double)tot_smp();
}

double ThermalMultiFlux::pot_var (int t_ind, int pes) const
{
  const char funame [] = "ThermalMultiFlux::pot_var: ";
  if(t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  double dtemp = t_sum(t_ind, pes) / (double)acct_smp();
  return sqrt(t_var(t_ind, pes) / (double)acct_smp() - dtemp * dtemp)
    * (double)pot_smp() / (double)tot_smp();
}

double ThermalMultiFlux::vol_var (int t_ind, int pes) const
{
  const char funame [] = "ThermalMultiFlux::vol_var: ";

  if(t_ind < 0 || t_ind >= tm_size()) {
    cout << funame << "index out of range" << endl;
    throw Error::Range_Err();
  }

  if(!is_pot())
    return 0.;

  return t_sum(t_ind, pes) / (double)acct_smp()
    * sqrt((double) vol_smp() * (double) pot_smp())
    / (double) tot_smp();
}

ostream& operator<< (ostream& to, const ThermalMultiFlux& fl)
{
  to << static_cast<const FluxBase&>(fl)
     << fl.min_en << "\n"
     << fl.min_dv << "\n";

  for(int i = 0; i <  fl.t_sum.size(); ++i)
    to <<  fl.t_sum[i] << "   " <<  fl.t_var[i] << "\n";
  to << "\n";

  return to;
}

istream& operator>> (istream& from, ThermalMultiFlux& fl)
{
  from >> static_cast<FluxBase&>(fl)
       >> fl.min_en
       >> fl.min_dv;

  for(int i = 0; i <  fl.t_sum.size(); ++i)
    from >>  fl.t_sum[i] >> fl.t_var[i];

  return from;
}

void ThermalMultiFlux::master_recv (int src)
{


  FluxBase::recv(src, Comm::FLUX_TAG);

  MPI::COMM_WORLD.Recv(min_en.data(), min_en.size(), MPI::DOUBLE, 
		       src, Comm::FLUX_TAG);
  for(int i = 0; i < min_dv.size(); ++i)
    MPI::COMM_WORLD.Recv(min_dv[i].data(), min_dv[i].size(), 
			 MPI::DOUBLE, src, Comm::FLUX_TAG);

  MPI::COMM_WORLD.Recv(t_sum.data(), t_sum.size(), 
		       MPI::DOUBLE, src, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Recv(t_var.data(), t_var.size(), 
		       MPI::DOUBLE, src, Comm::FLUX_TAG);
}

void ThermalMultiFlux::slave_send () const
{


  FluxBase::send(Comm::FLUX_TAG);
  MPI::COMM_WORLD.Send(min_en.data(), min_en.size(), MPI::DOUBLE, 0, 
		       Comm::FLUX_TAG);
  for(int i = 0; i < min_dv.size(); ++i)
    MPI::COMM_WORLD.Send(min_dv[i].data(), min_dv[i].size(), 
			 MPI::DOUBLE, 0, Comm::FLUX_TAG);

  MPI::COMM_WORLD.Send(t_sum.data(), t_sum.size(), MPI::DOUBLE, 
		       0, Comm::FLUX_TAG);
  MPI::COMM_WORLD.Send(t_var.data(), t_var.size(), MPI::DOUBLE, 
		       0, Comm::FLUX_TAG);
}

/*****************************************************************************
 ************************** Microcannonical Flux * ***************************
 *****************************************************************************/

Flux::Flux() : e_sum(en_size()), e_var(en_size()), 
	       ej_sum(en_size(), am_size()), 
	       ej_var(en_size(), am_size())
{ 
  if (!is_stat_init()) 
    error("Flux::Flux: static variables are not initialized");

  e_sum.init(); 
  e_var.init(); 

  ej_sum.init(); 
  ej_var.init();
}

Flux& Flux::operator+= (const Flux& fl)
{
  ThermalFlux::operator+=(fl);

  e_sum += fl.e_sum;
  e_var += fl.e_var;

  ej_sum += fl.ej_sum;
  ej_var += fl.ej_var;

  return *this;
}

Flux Flux::operator+ (const Flux& fl) const
{
  Flux res(*this);
  return res += fl;
}

Flux::Flux (const Div_surf& ds, int smp_num, std::vector<Samp>& smp_array)
  : e_sum(en_size()), e_var(en_size()), ej_sum(en_size(), am_size()), ej_var(en_size(), am_size())
{
  static const char funame [] = "Flux::Flux: ";

  if (!is_stat_init()) 
    error("Flux::Flux: static variables are not initialized");

  e_sum.init(); 
  e_var.init(); 

  ej_sum.init(); 
  ej_var.init();

  if(smp_array.size() != smp_num)
    smp_array.resize(smp_num);

  // temporary
  double dtemp;
  int itemp;

  // number of degrees of freedom
  static int dof_num = 3;
  static bool is_dof = false;
  if(!is_dof) {
    is_dof = true;
    for (int i = 0; i < 2; ++i) {
      switch(mol_array[i]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dof_num += 2;
	break;
      case NONLINEAR:
	dof_num += 3;
      }
    }
  }

  static double cn_fac = 2.0 * sqrt(2.0 * M_PI);
  static double mc_fac = M_2_SQRTPI * M_SQRT1_2 / 2. / gamma_2(dof_num + 1);
  static double ej_fac = M_1_PI / gamma_2(dof_num - 2);

  static Samp samp;
  double weight;         // sampling statistical weight
  double tim [3];        // total inertial moments
  Array<double> pot_ener(PES::size());       // potential energy
  pot_ener.init();

  int as_num = 0, fs_num = 0, cs_num = 0, ds_num = 0;
  int pot_num, vol_num;
  while(as_num < smp_num) {// sampling cycle

    pot_num = as_num + fs_num;
    vol_num = cs_num + ds_num;

    if (fs_num > fail_num_max && !as_num) {
      cout << "node " << Comm::rank() << ": " << funame
	   << "Oops!!! all potential calculations failed\n";
      add_acct_smp(as_num);
      add_fail_smp(fs_num);
      add_face_smp(cs_num);
      add_dist_smp(ds_num);
      throw Error::Run_Err();
    }

    if (!pot_num && vol_num > vol_num_max)
      break;

    // random orient
    switch(rand_pos(ds, samp, &weight)) {
    case 0:
      break;
    case SAMP_ATOMS_CLOSE:
      ++ds_num;
      continue;
    case SAMP_FACE_OUT:
      ++cs_num;
      continue;
    default:
      cout << Comm::mesg() << funame
	   << "rand_pos's status unknown, exitting\n";
      exit(1);
    }

    if(get_tot_iner_mom(tim)) {
      std::cout << Comm::mesg() << funame
		<< "Oops!!! inertia moments calculation failed\n";
      ++fs_num;
      continue;
    }
    if(tim[0] <= 0. || tim[1] <= 0. || tim[2] <= 0.) {
      std::cout << Comm::mesg() << funame
		<< "Oops!!! some of inertia moments are negative or zeros\n";
      ++fs_num;
      continue;
    }
    
    const double tim_sqrt = sqrt(tim[0]*tim[1]*tim[2]);

    try {
      PES::pot(0, pot_ener);
    }
    catch (Pot_error perr) {

      cout << Comm::mesg() << funame
	   << "Oops!!! Caught exception from potential function\n";

      ++fs_num;
      continue;

    }// exception handler

    // sampling succeeded
    // energy minimum
    double pot_ener_min;
    for(int pes = 0; pes < PES::size(); ++pes)
      if(!pes || pot_ener[pes] < pot_ener_min)
	pot_ener_min = pot_ener[pes];

    if(!(as_num + acct_smp()) || pot_ener_min < min_en) {
      min_en = pot_ener_min;
      min_dv = samp;
    }
	      
    if(smp_out_flag) {// send raw sampling to the master
      samp.set_val(weight, pot_ener_min);
      smp_array[as_num]=samp;
    }

    ++as_num;

    //PES cycle
    //
    for(int pes = 0; pes < PES::size(); ++pes) {
      //
      // cannonical flux
      //
      for (int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind)  {
	//
	double t_val = ThermalFlux::tm_grid(t_ind);
	
	dtemp = pot_ener[pes] / t_val;
	
	if (dtemp < -Limits::exp_arg_max()) {
	  //
	  std::cerr << Comm::mesg() << funame  << "potential energy, " << pot_ener[pes] / Phys_const::kcal 
		    << " kcal/mol, is too negative, trancating" << std::endl;
	  
	  dtemp = -Limits::exp_arg_max();
	}
	
	if (dtemp < Limits::exp_arg_max()) {
	  //
	  // thermal flux contribution
	  //
	  double fval = cn_fac * weight * exp(-dtemp) * sqrt(t_val);

	  t_sum[t_ind] += fval;
	  
	  t_var[t_ind] += fval * fval;
	}
      }

      // microcanonical flux
      //
      for(int en_ind = 0; en_ind < en_size(); ++en_ind) {
	//
	double ken = en_grid(en_ind) - pot_ener[pes];
	
	if (ken <= 0.)
	  //
	  continue;
	
	dtemp = mc_fac * weight * std::pow(ken, double(dof_num - 1) / 2.);
	
	e_sum[en_ind] += dtemp;
	
	e_var[en_ind] += dtemp * dtemp;

	// E,J-resolved flux

	// angular momentum cycle
	//
	for (int am_ind = 0; am_ind < am_size(); ++am_ind) {
	  //
	  dtemp = mc_stat_weight(ken, am_grid(am_ind), tim, dof_num) * ej_fac * weight / tim_sqrt;
	  
	  ej_sum(en_ind, am_ind) += dtemp;
	  
	  ej_var(en_ind, am_ind) += dtemp * dtemp;
	  //
	  //
	}// angular momentum cycle
	//
      }// energy cycle
      //
    }// potential energy surface cycle
    //
  }// sampling cycle		    

  add_acct_smp(as_num);
  add_fail_smp(fs_num);
  add_face_smp(cs_num);
  add_dist_smp(ds_num);
}

void Flux::normalize ()
{
  static const double min_flux = 1.e-99;

  if(!acct_smp())
    return;

  double dtemp;

  ThermalFlux::normalize();

  double smp_fac = (double)pot_smp() / (double)tot_smp()
    / (double)acct_smp();

  double fl_fac = 1./acct_smp() + 1./tot_smp() - 1./pot_smp();

  for (int i = 0; i < en_size(); ++i) {
    double& fl = e_sum[i];
    double& df = e_var[i];

    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }

  for(int i = 0; i < ej_sum.size(); ++i) {
    double& fl = ej_sum[i];
    double& df = ej_var[i];

    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }
}
  
ostream& operator<< (ostream& to, const Flux& fl)
{
  to << static_cast<const ThermalFlux&>(fl);

  for (int i = 0; i <  FluxBase::en_size(); ++i)
    to <<  fl.e_sum[i] << "   " <<  fl.e_var[i] << "\n";
  to << "\n";
  


  for (int i = 0; i <  fl.ej_sum.size(); ++i)
    to <<  fl.ej_sum[i] << "   " <<  fl.ej_var[i] << "\n";
  to << "\n";
  
  return to;
}

istream& operator>> (istream& from, Flux& fl)
{
  from >> static_cast<ThermalFlux&>(fl);

  for (int i = 0; i <  FluxBase::en_size(); ++i)
    from >> fl.e_sum[i] >> fl.e_var[i];

  for (int i = 0; i < fl.ej_sum.size(); ++i)
    from >> fl.ej_sum[i] >> fl.ej_var[i];
  return from;
}

void Flux::slave_send () const
{
  const char funame [] = "Flux::slave_send: ";

  try {

    ThermalFlux::slave_send();

    MPI::COMM_WORLD.Send(e_sum.data(), en_size(), MPI::DOUBLE, 
			 0, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Send(e_var.data(), en_size(), MPI::DOUBLE,
			 0, Comm::FLUX_TAG);
    

    MPI::COMM_WORLD.Send(ej_sum.data(), ej_sum.size(),
			 MPI::DOUBLE, 0, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Send(ej_var.data(), ej_var.size(), 
			 MPI::DOUBLE, 0, Comm::FLUX_TAG);      
  }
  catch (MPI::Exception) {
    cout << "node " << Comm::rank() 
	 << ": error in sending flux to master" << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 <<  Comm::FLUX_TAG 
	 << "-th tag sent to master" << endl;
}

void Flux::master_recv (int src)
{
  const char funame [] = "Flux::master_recv: ";

  try {
    ThermalFlux::master_recv(src);

    MPI::COMM_WORLD.Recv(e_sum.data(), en_size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Recv(e_var.data(), en_size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);


    MPI::COMM_WORLD.Recv(ej_sum.data(), ej_sum.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Recv(ej_var.data(), ej_var.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);

  }
  catch (MPI::Exception) {
    cout << "master: error in receiving flux from "
         << src << "-th node" << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame
         <<  Comm::FLUX_TAG
         << "-th tag received from "
         << src << "-th node" << endl;
}

/**********************************************
 **************** MultiFlux *******************
 **********************************************/

MultiFlux::MultiFlux (const Div_surf& ds, int smp_num, std::vector<Samp>& smp_array)
  : e_sum(en_size(), PES::size()), e_var(en_size(), PES::size()), 
    ej_sum(en_size(), am_size(), PES::size()), ej_var(en_size(), am_size(), PES::size())
{
  static const char funame [] = "MultiFlux::MultiFlux: ";

  if (!is_stat_init()) 
    error("MultiFlux::MultiFlux: static variables are not initialized");

  e_sum.init(); 
  e_var.init(); 

  ej_sum.init(); 
  ej_var.init();

  // temporary
  double dtemp;
  int itemp;

  // number of degrees of freedom
  static int dof_num = 3;
  static bool is_dof = false;
  if(!is_dof) {
    is_dof = true;
    for (int i = 0; i < 2; ++i) {
      switch(mol_array[i]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dof_num += 2;
	break;
      case NONLINEAR:
	dof_num += 3;
      }
    }
  }

  static double cn_fac = 2.0 * sqrt(2.0 * M_PI);
  static double mc_fac = M_2_SQRTPI * M_SQRT1_2 / 2. / gamma_2(dof_num + 1);
  static double ej_fac = M_1_PI / gamma_2(dof_num - 2);

  if(smp_num != smp_array.size())
    smp_array.resize(smp_num);

  static Samp samp;
  static Array<double> pot_ener(PES::size());       // potential energy
  pot_ener.init();

  double weight;         // sampling statistical weight
  
  double tim [3];        // total inertial moments

  int as_num = 0, fs_num = 0, cs_num = 0, ds_num = 0;
  
  int pot_num, vol_num;
  while(as_num < smp_num) {// sampling cycle

    pot_num = as_num + fs_num;
    vol_num = cs_num + ds_num;

    if(fs_num > fail_num_max && !as_num) {
      cout << "node " << Comm::rank() << ": " << funame
	   << "Oops!!! all potential calculations failed\n";
      add_acct_smp(as_num);
      add_fail_smp(fs_num);
      add_face_smp(cs_num);
      add_dist_smp(ds_num);
      throw Error::Run_Err();
    }

    if (!pot_num && vol_num > vol_num_max)
      break;

    // random orient
    //
    switch(rand_pos(ds, samp, &weight)) {
      //
    case 0:
      //
      break;
    case SAMP_ATOMS_CLOSE:
      //
      ++ds_num;
      
      continue;
    case SAMP_FACE_OUT:
      //
      ++cs_num;
      
      continue;
    default:
      //
      std::cerr << Comm::mesg() << funame
		<< "rand_pos's status unknown, exitting" << std::endl;
      
      exit(1);
    }

    if(get_tot_iner_mom(tim)) {
      std::cerr << Comm::mesg() << funame
		<< "Oops!!! inertia moments calculation failed" << std::endl;
      
      ++fs_num;
      
      continue;
    }
    
    if(tim[0] <= 0. || tim[1] <= 0. || tim[2] <= 0.) {
      //
      std::cerr << Comm::mesg() << funame
		<< "Oops!!! some of inertia moments are negative or zeros"
		<< std::endl;
      
      ++fs_num;
      continue;
    }
    
    const double tim_sqrt = sqrt(tim[0] * tim[1] * tim[2]);
    
    try {
      //
      PES::pot(0, pot_ener);
    }
    catch (Pot_error perr) {
      //
      std::cerr << Comm::mesg() << funame
		<< "Oops!!! Caught exception from potential function" << std::endl;

      ++fs_num;
      continue;
    }

    // sampling succeeded

    
    for(int pes = 0; pes < PES::size(); ++pes) {
      //
      // energy minimum
      //
      if(!(as_num + acct_smp()) || pot_ener[pes] < min_en[pes]) {
	//
	min_en[pes] = pot_ener[pes];
	
	min_dv[pes] = samp;
      }
    }

    // send raw sampling to the master
    //
    if(smp_out_flag) {
      //
      samp.set_val(weight, pot_ener);

      smp_array[as_num] = samp;
    }

    ++as_num;

    // potential energy cycle
    //
    for(int pes = 0; pes < PES::size(); ++pes) {
      //
      // canonical flux
      //
      for (int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind)  {
	//
	double t_val = ThermalFlux::tm_grid(t_ind);
	
	dtemp = pot_ener[pes] / t_val;
	
	if(dtemp < -Limits::exp_arg_max()) {
	  //
	  std::cerr << Comm::mesg() << funame  << "potential energy, " << pot_ener[pes] / Phys_const::kcal 
		    << " kcal/mol, is too negative, trancating" << std::endl;
	  
	  dtemp = -Limits::exp_arg_max();
	}
	
	if(dtemp < Limits::exp_arg_max()) {
	  //
	  // thermal flux contribution
	  //
	  double fval = cn_fac * weight * sqrt(t_val) * std::exp(-dtemp);

	  t_sum(t_ind, pes) += fval;
	  
	  t_var(t_ind, pes) +=  fval * fval;
	}
      }

      // energy cycle
      //
      for (int en_ind = 0; en_ind < en_size(); ++en_ind) {
	//
	// kinetic energy
	//
	double kinen = en_grid(en_ind) - pot_ener[pes];
	
	if (kinen <= 0.)
	  //
	  continue;
	
	// microcanonical flux
	//
	dtemp = mc_fac * weight * std::pow(kinen, double(dof_num - 1) / 2.);
	
	e_sum(en_ind, pes) += dtemp;
	
	e_var(en_ind, pes) += dtemp * dtemp;

	// E,J-resolved flux
	//
	// angular momentum cycle
	//
	for (int am_ind = 0; am_ind < am_size(); ++am_ind) {
	  //
	  dtemp = mc_stat_weight(kinen, am_grid(am_ind), tim, dof_num) * ej_fac * weight / tim_sqrt;

	  ej_sum(en_ind, am_ind, pes) += dtemp;
	  
	  ej_var(en_ind, am_ind, pes) += dtemp * dtemp;
	  //
	  //
	}// angular momentum cycle
	//
      }// energy cycle
      //
    }// potential energy surface cycle
    //
  }// sampling cycle		    

  add_acct_smp(as_num);
  add_fail_smp(fs_num);
  add_face_smp(cs_num);
  add_dist_smp(ds_num);
}

MultiFlux::MultiFlux() : e_sum(en_size(), PES::size()), 
			 e_var(en_size(), PES::size()), 
			 ej_sum(en_size(), am_size(), PES::size()), 
			 ej_var(en_size(), am_size(), PES::size())
{ 
  if (!is_stat_init()) 
    error("MultiFlux::MultiFlux: static variables are not initialized");

  e_sum.init(); 
  e_var.init(); 

  ej_sum.init(); 
  ej_var.init();
}

MultiFlux& MultiFlux::operator+= (const MultiFlux& fl)
{
  ThermalMultiFlux::operator+=(fl);

  e_sum += fl.e_sum;
  e_var += fl.e_var;

  ej_sum += fl.ej_sum;
  ej_var += fl.ej_var;

  return *this;
}

MultiFlux MultiFlux::operator+ (const MultiFlux& fl) const
{
  MultiFlux res(*this);
  return res += fl;
}

MultiFlux& MultiFlux::operator*= (double mf) {
  t_sum  *= mf;
  e_sum  *= mf;
  ej_sum *= mf;
  t_var  *= mf*mf;
  e_var  *= mf*mf;
  ej_var *= mf*mf;

  return *this;
}

void MultiFlux::normalize ()
{
  static const double min_flux = 1.e-99;

  if(!acct_smp())
    return;

  double dtemp;

  ThermalMultiFlux::normalize();

  double smp_fac = (double)pot_smp() / (double)tot_smp()
    / (double)acct_smp();

  double fl_fac = 1./acct_smp() + 1./tot_smp() - 1./pot_smp();

  for (int i = 0; i < e_sum.size(); ++i) {
    double& fl = e_sum[i];
    double& df = e_var[i];

    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }

  for(int i = 0; i < ej_sum.size(); ++i) {
    double& fl = ej_sum[i];
    double& df = ej_var[i];

    fl *= smp_fac;
    if(fl < min_flux)
      df = 0.;
    else {
      dtemp = df * smp_fac * smp_fac - fl * fl * fl_fac;
      if(dtemp <= 0.)
	df = 0.;
      else
	df = sqrt(dtemp)/ fl * 100.;
    }
  }
}
  
ostream& operator<< (ostream& to, const MultiFlux& fl)
{
  to << static_cast<const ThermalMultiFlux&>(fl);

  for (int i = 0; i <  fl.e_sum.size(); ++i)
    to <<  fl.e_sum[i] << "   " <<  fl.e_var[i] << "\n";
  to << "\n";
  

  for (int i = 0; i <  fl.ej_sum.size(); ++i)
    to <<  fl.ej_sum[i] << "   " <<  fl.ej_var[i] << "\n";
  to << "\n";
  
  return to;
}

istream& operator>> (istream& from, MultiFlux& fl)
{
  from >> static_cast<ThermalMultiFlux&>(fl);

  for (int i = 0; i <  fl.e_sum.size(); ++i)
    from >> fl.e_sum[i] >> fl.e_var[i];

  for (int i = 0; i < fl.ej_sum.size(); ++i)
    from >> fl.ej_sum[i] >> fl.ej_var[i];
  return from;
}

void MultiFlux::slave_send () const
{
  const char funame [] = "MultiFlux::slave_send: ";

  try {
    ThermalMultiFlux::slave_send();


    MPI::COMM_WORLD.Send(e_sum.data(), e_sum.size(), MPI::DOUBLE, 
			 0, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Send(e_var.data(), e_var.size(), MPI::DOUBLE,
			 0, Comm::FLUX_TAG);



    MPI::COMM_WORLD.Send(ej_sum.data(), ej_sum.size(),
			 MPI::DOUBLE, 0, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Send(ej_var.data(), ej_var.size(), 
			 MPI::DOUBLE, 0, Comm::FLUX_TAG);

  }
  catch (MPI::Exception) {
    cout << "node " << Comm::rank() 
	 << ": error in sending flux to master" << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 <<  Comm::FLUX_TAG 
	 << "-th tag sent to master" << endl;
}

void MultiFlux::master_recv (int src)
{
  const char funame [] = "MultiFlux::master_recv: ";

  try {
    ThermalMultiFlux::master_recv(src);

    MPI::COMM_WORLD.Recv(e_sum.data(), e_sum.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Recv(e_var.data(), e_var.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Recv(ej_sum.data(), ej_sum.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);
    MPI::COMM_WORLD.Recv(ej_var.data(), ej_var.size(),
			 MPI::DOUBLE, src, Comm::FLUX_TAG);

  }
  catch (MPI::Exception) {
    cout << "master: error in receiving flux from "
         << src << "-th node" << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame
         <<  Comm::FLUX_TAG
         << "-th tag received from "
         << src << "-th node" << endl;
}

/**********************************************************************
 ************************ MultiFlux Array *****************************
 **********************************************************************/

void MultiFluxArray::normalize()
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

ostream& operator<< (ostream& to, const MultiFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    to << f[i];
  return to;
}
istream& operator>> (istream& from,  MultiFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    from >> f[i];
  return from;
}

/*************************************************************
 **************** Distributed MultiFlux Array ****************
 *************************************************************/

void DistMultiFluxArray::normalize()
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

int DistMultiFluxArray::node_size () const
{ 
  int res = 0; 
  for (int i = 0; i < _array.size(); ++i)
    res += _array[i].node_size();
  return res;
}

ostream& operator<< (ostream& to, const DistMultiFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    to << f[i];
  return to;
}

istream& operator>> (istream& from,  DistMultiFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    from >> f[i];
  return from;
}

/********************************************************************
 ******************* Thermal MultiFlux Array ************************
 ********************************************************************/

ThermalMultiFluxArray::ThermalMultiFluxArray(const MultiFluxArray& sf) 
  : _array(sf.size()) 
{
  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
}

ThermalMultiFluxArray::ThermalMultiFluxArray(const DistMultiFluxArray& sf) 
  : _array(sf.size())
{
  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
}

ThermalMultiFluxArray& ThermalMultiFluxArray::operator= (const MultiFluxArray& sf)
{
  const char funame [] = "ThermalMultiFluxArray::operator=: ";

  if(sf.size() != _array.size()) {
    std::cout << funame << "different dimensions\n";
    throw Error::Range_Err();
  }

  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
  return *this;
}

ThermalMultiFluxArray& ThermalMultiFluxArray::operator= (const DistMultiFluxArray& sf)
{
  const char funame [] = "ThermalMultiFluxArray::operator=: ";

  if(sf.size() != _array.size()) {
    std::cout << funame << "different dimensions\n";
    throw Error::Range_Err();
  }

  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
  return *this;
}

void ThermalMultiFluxArray::normalize() 
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

/*****************************************************************
 ****************** MultiFlux Calculation ************************
 *****************************************************************/

int DistMultiFluxArrayCalc::check_state(int& face)
{
  double dtemp;
  int itemp;

  // wait for sampling results
  if(iswait())
    return WAIT;

  // initial sampling
  bool is_init_smp_done = true;
  for (face = 0; face < size(); ++face) {
    int pot_smp  = (*this)[face].pot_smp();
    int tot_smp  = (*this)[face].tot_smp();
    int acct_smp = (*this)[face].acct_smp();
    int proj_smp = (*this)[face].node_size() * pot_len() + acct_smp;

    if (tot_smp >= tot_min() && !pot_smp || acct_smp >= pot_min())
      continue; // initial sampling done

    is_init_smp_done = false;
    if (proj_smp <  pot_min())
      return FLUX;
    else
      continue;
  }

  if (!is_init_smp_done)
    return WAIT; // initial sampling stage

  // estimate an error
  double min_val; int min_ind;
  for(int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind) {
    double val = 0.;
    for(int i = 0; i < size(); ++i)
      val += (*this)[i].average(t_ind, 0);
      
    if(!t_ind) {
      min_val = val;
      min_ind = t_ind;
      continue;
    }      
    if(val < min_val) {
      min_val = val;
      min_ind = t_ind;
    }      
  }
  
  double tot_pot_var = 0.;
  for (int i = 0; i < size(); ++i)
    tot_pot_var += (*this)[i].pot_var(min_ind, 0);

  // projected number of samplings
  double proj_smp_num;
  if(min_val < 1.e-99)
    proj_smp_num = 0.;
  else {
    proj_smp_num = 10000. * tot_pot_var * tot_pot_var 
      / min_val / min_val / tol() / tol();
    proj_smp_num =  proj_smp_num < (double)pot_max() ? proj_smp_num 
      : (double)pot_max();
  }

  if(proj_smp_num > 1.) {
    double max_smp = -1.;
    for (int i = 0; i < size(); ++i) {
      dtemp = proj_smp_num * (*this)[i].pot_var(min_ind, 0) / tot_pot_var;
      if(dtemp > 1.)
	dtemp = 1. - double((*this)[i].node_size() * pot_len() + 
			    (*this)[i].acct_smp()) / dtemp;
      else
	dtemp = -1.;

      if(dtemp > max_smp) {
	max_smp = dtemp;
	face = i;
      }
    }
    if(max_smp > 0.)
      return FLUX;
  }

  if(node_size())// wait for running nodes to finish
    return WAIT;

  // SURFACE AREA SAMPLING

  // estimate an error
  double tot_vol_var = 0.;
  for (int i = 0; i < size(); ++i)
    tot_vol_var += (*this)[i].vol_var(min_ind, 0);

  // projected number of  samplings
  if (min_val < 1.e-99)
    proj_smp_num = 0.;
  else {
    proj_smp_num = 100000. * tot_vol_var * tot_vol_var 
      / min_val / min_val / tol() / tol();
    proj_smp_num =  proj_smp_num < (double)tot_max() ? proj_smp_num 
      : (double)tot_max();
  }

  if (proj_smp_num <= 1.)
    return STOP;
  
  bool is_surf_smp = false;
  for(int i = 0; i < size(); ++i) {
    itemp  = int(proj_smp_num * (*this)[i].vol_var(min_ind, 0) 
		       / tot_vol_var)	- (*this)[i].tot_smp();

    if(itemp > 0) {
      is_surf_smp = true;
      _smp_num[i] = itemp;
    }
    else
      _smp_num[i] = 0;
  }

  // need surface sampling
  if (is_surf_smp)
    return SURF;
  // done
  return STOP;

}

int DistMultiFluxCalc::check_state(int& surf_smp_num)
{
  // wait for sampling results
  if (iswait())
    return WAIT;

  if(!pot_smp() && tot_smp() >= tot_min())
    if(node_size())
      return WAIT;
    else
      return STOP;

  // initial flux sampling
  if(acct_smp() < pot_min())
    if (node_size() * pot_len() + acct_smp() < pot_min())
      return FLUX;
    else
      return WAIT;

  // flux sampling
  double dtemp;
  double min_val; int min_ind;
  for(int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind) {
    dtemp = average(t_ind, 0);
    if(!t_ind) {
      min_val = dtemp;
      min_ind = t_ind;
      continue;
    }      
    if(dtemp < min_val) {
      min_val = dtemp;
      min_ind = t_ind;
    }      
  }
  
  int proj_smp_num;
  if(min_val < 1.e-99)
    proj_smp_num = 0;
  else {
    dtemp = pot_var(min_ind, 0);
    proj_smp_num = int(10000. * dtemp * dtemp 
		       / min_val / min_val / tol() / tol());
    proj_smp_num =  proj_smp_num < pot_max() ? proj_smp_num : pot_max();
  }

  if(proj_smp_num - node_size() * pot_len() - acct_smp() > 0)
    return FLUX;

  if(node_size())
    return WAIT;

  // surface sampling
  if(min_val < 1.e-99)
    proj_smp_num = 0;
  else {
    dtemp = vol_var(min_ind, 0);
    proj_smp_num = int(100000. * dtemp * dtemp 
		       / min_val / min_val / tol()/ tol());
    proj_smp_num =  proj_smp_num < tot_max() ? proj_smp_num : tot_max();
  }

  surf_smp_num = proj_smp_num - tot_smp();
  if(surf_smp_num > 0)
    return SURF;
  else
    return STOP;
}

/*****************************************************************
 ************************ Flux Array *****************************
 *****************************************************************/

void FluxArray::normalize()
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

ostream& operator<< (ostream& to, const FluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    to << f[i];
  return to;
}
istream& operator>> (istream& from,  FluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    from >> f[i];
  return from;
}

/***********************************************************
 ***************** Distributed Calculation *****************
 ***********************************************************/

void DistBase::add_node(int n)
{
  const char funame [] = "DistBase::add_node: ";
  if(!_nodes.insert(n).second) {
    std::cout << funame << "node " << n << " is allready in the list\n";
    throw Error::Run_Err();
  }
}

void DistBase::remove_node(int n)
{
  const char funame [] = "DistBase::remove_node: ";
  
  if(!_nodes.erase(n)) {
    std::cout << funame << "node " << n << " is not found in the list\n";
    throw Error::Run_Err();
  }
}

int DistBase::node (int i) const
{ 
  std::set<int>::const_iterator it = _nodes.begin();
  std::advance(it, i);
  return *it;
}

/*****************************************************************************
 ************************ Distributed Flux Array *****************************
 *****************************************************************************/

void DistFluxArray::normalize()
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

int DistFluxArray::node_size () const
{ 
  int res = 0; 
  for (int i = 0; i < _array.size(); ++i)
    res += _array[i].node_size();
  return res;
}
ostream& operator<< (ostream& to, const DistFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    to << f[i];
  return to;
}
istream& operator>> (istream& from,  DistFluxArray& f)
{
  for(int i = 0; i < f.size(); ++i)
    from >> f[i];
  return from;
}

/***************************************************************
 ******************* Thermal Flux Array ************************
 ***************************************************************/

ThermalFluxArray::ThermalFluxArray(const FluxArray& sf) : _array(sf.size()) 
{
  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
}

ThermalFluxArray::ThermalFluxArray(const DistFluxArray& sf) : _array(sf.size())
{
  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
}

ThermalFluxArray& ThermalFluxArray::operator= (const FluxArray& sf) 
{
  const char funame [] = "ThermalFluxArray::operator=: ";

  if(sf.size() != _array.size()) {
    std::cout << funame << "different dimensions\n";
    throw Error::Range_Err();
  }

  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
  return *this;
}

ThermalFluxArray& ThermalFluxArray::operator= (const DistFluxArray& sf) 
{
  const char funame [] = "ThermalFluxArray::operator=: ";

  if(sf.size() != _array.size()) {
    std::cout << funame << "different dimensions\n";
    throw Error::Range_Err();
  }

  for(int i = 0; i < sf.size(); ++i)
    _array[i] = sf[i];
  return *this;
}

void ThermalFluxArray::normalize() 
{
  for(int i = 0; i < _array.size(); ++i)
    _array[i].normalize();
}

/***************************************************************
 ******************* Flux Array Calculation ********************
 ***************************************************************/

double Calc::_tol = 0.; // required accuracy
int Calc::_pmx = 0;  // maximal possible number of potential samplings
int Calc::_pmn = 0;  // minimal number of potential samplings
int Calc::_tmx = 0;  // maximal possible number of total samplings
int Calc::_tmn = 0;  // minimal number of total samplings
int Calc::_pln = 0;  // number of potential samplings in one calculation
bool Calc::_isinit = false;  // are static vars initialized?
  
void Calc::stat_init(double tl, int px, int pn, int tx, int tn, int pl)
{
  const char funame [] = "Calc::stat_init: ";
  if(_isinit) {
    std::cout << funame << "Calculation parameters allready initialized\n";
    throw Error::Init_Err();
  }
  _isinit = true;

  _tol = tl;
  _pmx = px;
  _pmn = pn;
  _tmx = tx;
  _tmn = tn;
  _pln = pl;
}

int DistFluxArrayCalc::check_state(int& face)
{
  double dtemp;
  int itemp;

  // wait for sampling results
  if(iswait())
    return WAIT;

  // initial sampling
  bool is_init_smp_done = true;
  for (face = 0; face < size(); ++face) {
    int pot_smp  = (*this)[face].pot_smp();
    int tot_smp  = (*this)[face].tot_smp();
    int acct_smp = (*this)[face].acct_smp();
    int proj_smp = (*this)[face].node_size() * pot_len() + acct_smp;

    if (tot_smp >= tot_min() && !pot_smp || acct_smp >= pot_min())
      continue; // initial sampling done

    is_init_smp_done = false;
    if (proj_smp <  pot_min())
      return FLUX;
    else
      continue;
  }

  if (!is_init_smp_done)
    return WAIT; // initial sampling stage

  // estimate an error
  double min_flux_val;
  int    min_flux_ind;
  for(int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind) {
    double val = 0.;
    for(int i = 0; i < size(); ++i)
      val += (*this)[i].average(t_ind);
      
    if(!t_ind) {
      min_flux_val = val;
      min_flux_ind = t_ind;
      continue;
    }      
    if(val < min_flux_val) {
      min_flux_val = val;
      min_flux_ind = t_ind;
    }      
  }
  
  double tot_pot_var = 0.;
  for (int i = 0; i < size(); ++i)
    tot_pot_var += (*this)[i].pot_var(min_flux_ind);

  // maximum number of samplings
  double smp_max;
  if(min_flux_val < 1.e-99)
    smp_max = 0.;
  else {
    smp_max = 10000. * tot_pot_var * tot_pot_var 
      / min_flux_val / min_flux_val / tol() / tol();
    smp_max =  smp_max < (double)pot_max() ? smp_max 
      : (double)pot_max();
  }

  // Flux sampling
  face = -1;
  double max_samp_diff = 0.;
  for (int i = 0; i < size(); ++i) {
    dtemp = smp_max * (*this)[i].pot_var(min_flux_ind) / tot_pot_var
      - double((*this)[i].node_size() * pot_len() + (*this)[i].acct_smp());
    if(max_samp_diff < dtemp) {
      max_samp_diff = dtemp;
      face = i;
    }
  }
  if(face >= 0)
    return FLUX;
  
  // wait for running nodes to finish
  if(node_size())
    return WAIT;

  // SURFACE AREA SAMPLING

  // estimate an error
  double tot_vol_var = 0.;
  for(int i = 0; i < size(); ++i)
    tot_vol_var += (*this)[i].vol_var(min_flux_ind);

  // projected number of  samplings
  if(min_flux_val < 1.e-99)
    smp_max = 0.;
  else {
    smp_max = 100000. * tot_vol_var * tot_vol_var 
      / min_flux_val / min_flux_val / tol() / tol();
    smp_max =  smp_max < (double)tot_max() ? smp_max 
      : (double)tot_max();
  }

  if(smp_max <= 1.)
    return STOP;
  
  bool is_surf_smp = false;
  for(int i = 0; i < size(); ++i) {
    itemp  = int(smp_max * (*this)[i].vol_var(min_flux_ind)/tot_vol_var)
      - (*this)[i].tot_smp();

    if(itemp > 0) {
      is_surf_smp = true;
      _smp_num[i] = itemp;
    }
    else
      _smp_num[i] = 0;
  }

  // need surface sampling
  if (is_surf_smp)
    return SURF;
  // done
  return STOP;

}

int DistFluxCalc::check_state(int& surf_smp_num)
{
  // wait for sampling results
  if (iswait())
    return WAIT;

  if(!pot_smp() && tot_smp() >= tot_min())
    if(node_size())
      return WAIT;
    else
      return STOP;

  // initial flux sampling
  if(acct_smp() < pot_min())
    if (node_size() * pot_len() + acct_smp() < pot_min())
      return FLUX;
    else
      return WAIT;

  // flux sampling
  double dtemp;
  double min_flux_val; 
  int    min_flux_ind;
  for(int t_ind = 0; t_ind < FluxBase::tm_size(); ++t_ind) {
    dtemp = average(t_ind);
    if(!t_ind) {
      min_flux_val = dtemp;
      min_flux_ind = t_ind;
      continue;
    }      
    if(dtemp < min_flux_val) {
      min_flux_val = dtemp;
      min_flux_ind = t_ind;
    }      
  }
  
  int proj_smp_num;
  if(min_flux_val < 1.e-99)
    proj_smp_num = 0;
  else {
    dtemp = pot_var(min_flux_ind);
    proj_smp_num = int(10000. * dtemp * dtemp 
		       / min_flux_val / min_flux_val / tol() / tol());
    proj_smp_num =  proj_smp_num < pot_max() ? proj_smp_num : pot_max();
  }

  if(proj_smp_num - node_size() * pot_len() - acct_smp() > 0)
    return FLUX;

  if(node_size())
    return WAIT;

  // surface sampling
  if(min_flux_val < 1.e-99)
    proj_smp_num = 0;
  else {
    dtemp = vol_var(min_flux_ind);
    proj_smp_num = int(100000. * dtemp * dtemp 
		       / min_flux_val / min_flux_val / tol()/ tol());
    proj_smp_num =  proj_smp_num < tot_max() ? proj_smp_num : tot_max();
  }

  surf_smp_num = proj_smp_num - tot_smp();
  if(surf_smp_num > 0)
    return SURF;
  else
    return STOP;
}

/******************************************************
 ******************** Surface Id **********************
 ******************************************************/

void Sid::send (int dest, int tag) const
{
  const char funame [] = "Sid::send: ";

  try {
    MPI::COMM_WORLD.Send(sid_, 2, MPI::INT, dest, tag);
  } 
  catch (MPI::Exception) {
    cout << "node " << Comm::rank() 
	 << ": error in sending sid to " << dest << "-th node with tag " 
	 << tag << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 <<  tag
	 << "-th tag sent to "
	 << dest << "-th node" << endl;
}

int Sid::slave_recv () // slave blocking receive
{
  const char funame [] = "Sid::slave_recv: ";

  if (!Comm::rank())
    error("Sid::master_recv: master node");

  MPI::Status stat;
  try {
    MPI::COMM_WORLD.Recv(sid_, 2, MPI::INT, 0, MPI::ANY_TAG, stat);
  } 
  catch (MPI::Exception) {
    cout << "node " << Comm::rank() 
	 << ": error in receiving sid" << endl;
    throw;
  }

  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 <<  stat.Get_tag()
	 << "-th tag received from master" << endl;

  return stat.Get_tag();
}

MPI::Status Sid::master_recv () // master nonblocking receive
{
  const char funame [] = "Sid::master_recv: ";

  if (Comm::rank())
    error("Sid::master_recv: not a master");

  MPI::Status stat;
  try {
    MPI::COMM_WORLD.Recv(sid_, 2, MPI::INT, MPI::ANY_SOURCE, 
			 MPI::ANY_TAG, stat);
  }
  catch (MPI::Exception) {
    cout << "master: error in receiving sid" << endl;
    throw;
  }
    
  if(Comm::debug)
    cout << Comm::mesg() << funame 
	 <<  stat.Get_tag()
	 << "-th tag received from "
	 << stat.Get_source() << "-th node" << endl;

  return stat;
}

ostream& operator<< (ostream& to, const Sid& s)
{
  to << s.sid_[0] << "   " << s.sid_[1] << "\n";
  return to;
}

istream& operator>> (istream& from, Sid& s)
{
  from >> s.sid_[0] >> s.sid_[1];
  //s.check();

  return from;
}
