#include "raninit.hh"

#include "math.hh"
#include "tmatrix.hh"
#include "random.hh"
#include "error.hh"

#include <cmath>
#include <iostream>

int rand_pos_0  (const Div_surf&, Dynvar&, double* =0);
int rand_pos_plane  (const Div_surf&, Dynvar&, double* =0);
int rand_pos_1a (const Div_surf&, Dynvar&, double* =0);
int rand_pos_2a (const Div_surf&, Dynvar&, double* =0);
int rand_pos_1b (const Div_surf&, Dynvar&, double* =0);
int rand_pos_2b (const Div_surf&, Dynvar&, double* =0);

double MIN_ATOM_DIST;

int (*rand_pos) (const Div_surf&, Dynvar&, double*) = rand_pos_0;

void name2ranp (const string& name)
{
  if (name == "multifacet")
    rand_pos = rand_pos_0;
  else if (name == "plane")
    rand_pos = rand_pos_plane;
  else if (name == "plane_1")
    rand_pos = rand_pos_1a;
  else if (name == "plane_2")
    rand_pos = rand_pos_2a;
  else if (name == "plane_1b")
    rand_pos = rand_pos_1b;
  else
    error("name2ranp: unknown name");
}

int rand_pos_plane (const Div_surf& ds, Dynvar& dynvar, double* weight)
{
  /********************************************************************
   * Sample spherical sector
   * Sphere - 0th face
   ********************************************************************/

  if (ds.face() != 0)
    error("rand_pos_plane: wrong face");

  if (ds.dist() <= 0.) // negative radius
    return SAMP_FACE_OUT;

  dynvar.attach();

  static double rrmass =   1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 

  // temporary variables
  double dtemp;


  double temp_pos [3];
  double pos [3];
  double cm_pos [3];    // cm-to-cm vector
  
  // randomly orient molecules
  for (int frag = 0; frag < 2; ++frag) {
    if (mol_array[frag]->type() != ATOM)
      Random::orient(mol_array[frag]->write_ang_pos(), 
		     mol_array[frag]->ang_pos_size());
    for (int j = 0; j < 3; ++j)
      mol_array[frag]->write_cm_pos()[j] = 0.0;
    mol_array[frag]->update_atoms();
  }

  // randomly orient pivot-point-to-pivot-point vector
  double pos_12 [3];
  Random::orient(pos_12, 3);
  for(int i = 0; i < 3; ++i)
    pos_12[i] *= ds.dist();

  // working part of the sphere
  for (int frag = 0; frag < 2; ++frag) 
    switch (mol_array[frag]->type()) {
    case ATOM:

      break;

    case LINEAR:

      switch(frag) {
      case 0:

	for(int i = 0; i < 3; ++i)
	  pos[i] = pos_12[i];
	break;

      case 1:

	for(int i = 0; i < 3; ++i)
	  pos[i] = - pos_12[i];
	break;
      }

      for(int i = 0; i < 3; ++i)
	pos[i] += ds.ref_pos(frag)[0] * mol_array[frag]->read_ang_pos()[i];

      for(int i = 1; i < ds.ref_num(frag); ++i) {
	if(frag)
	  dtemp = ds.dist(0,i);
	else
	  dtemp = ds.dist(i,0);
	if(vdot(pos, mol_array[frag]->read_ang_pos(), 3) * ds.ref_pos(frag, i)[0] < dtemp)
	  return SAMP_FACE_OUT;
      }
      break;

    case NONLINEAR:

      mol_array[frag]->lf2mf(pos_12, pos);
      if(frag)
	for(int i = 0; i < 3; ++i)
	  pos[i] = -pos[i];

      for(int i = 0; i < 3; ++i)
	pos[i] += ds.ref_pos(frag)[i];

      for(int i = 1; i < ds.ref_num(frag); ++i) {
	if(frag)
	  dtemp = ds.dist(0,i);
	else
	  dtemp = ds.dist(i,0);
	if(vdot(pos, ds.ref_pos(frag,i), 3) < dtemp)
	  return SAMP_FACE_OUT;
      }
      break;
    }

  // center-of-mass position
  for (int i = 0; i < 3; ++i)
    cm_pos[i] = pos_12[i];

  for (int frag = 0; frag < 2; ++frag) {
    switch (mol_array[frag]->type()) {
    case ATOM:
      for (int j = 0; j < 3; ++j)
	pos[j] = 0.0;
      break;
    case LINEAR:
      for (int j = 0; j < 3; ++j)
	pos[j] = mol_array[frag]->read_ang_pos()[j] * ds.ref_pos(frag)[0];
      break;
    case NONLINEAR:
      mol_array[frag]->mf2lf(ds.ref_pos(frag), pos);
    }
    for(int i = 0; i < 3; ++i)
      if(frag)
	cm_pos[i] -= pos[i];
      else
	cm_pos[i] += pos[i];
  }

  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = cm_pos[i];
  for(int frag = 0; frag < 2; ++frag)
    mol_array[frag]->update_atoms();

  // minimum interatomic distance
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;


  // statistical weight
  if(weight) {
    double d2 = ds.dist() * ds.dist();
    // cm motion
    *weight = rrmass * d2;

    double mf_pos_12 [3];
    for(int frag = 0; frag < 2; ++frag)
      switch(mol_array[frag]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dtemp = vdot(pos_12, mol_array[frag]->read_ang_pos(), 3);
	*weight += ds.ref_pos(frag)[0] * ds.ref_pos(frag)[0] 
	  * (d2 - dtemp * dtemp) / mol_array[frag]->iner_mom(2);
	break;
      case NONLINEAR:
	mol_array[frag]->lf2mf(pos_12, mf_pos_12);
	vector_product(mf_pos_12, ds.ref_pos(frag), temp_pos);
	for (int i = 0; i < 3; ++i)
	  *weight += temp_pos[i] * temp_pos[i] / mol_array[frag]->iner_mom(i);
      }
    *weight = sqrt(*weight * d2);
  }
  return 0;
}

int rand_pos_0 (const Div_surf& ds, Dynvar& dynvar, double* weight)
{


  dynvar.attach();

  if (ds.dist() <= 0.)
    return SAMP_FACE_OUT;

  if (weight)
    *weight = 0.;

  // temporary variables
  Molecule* mp_temp;
  double dtemp;
  int itemp;
  double pos_temp[3];

  // initialization of different parameters
  static double m_factor [2]; // mass factors 
  static double iner_mom2 [2] [3]; // square roots of inertia moments
  static bool first_run = true;
  if (first_run)
    {
      if (mol_array.size() != 2)
	error("init_thermal: wrong number of molecules");

      first_run = false;
      m_factor [0] =   mol_array[1]->mass() 
	/ (mol_array[0]->mass() + mol_array[1]->mass());
      m_factor [1] = - mol_array[0]->mass() 
	/ (mol_array[0]->mass() + mol_array[1]->mass());

      for (int i = 0; i < 2; ++i)
	{
	  mp_temp = mol_array[i];
	  switch(mp_temp->type())
	    {
	    case ATOM:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = 0.0;
	      break;
	    case LINEAR:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = sqrt(mp_temp->iner_mom(j));
	      break;
	    case NONLINEAR:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = sqrt(mp_temp->iner_mom(j));
	      break;
	    }
	}
    }
  // effective mass
  static const double emass = m_factor[0] * mol_array[0]->mass();
  static const double emass2 = sqrt(emass);
      
      
  // randomly orient a vector between reference points
  double lf_cos_12 [3];
  Random::orient(lf_cos_12, 3);

  // randomly orient molecules
  for (int frag = 0; frag < 2; ++frag) {
    mp_temp = mol_array[frag];
    if (mp_temp->type() != ATOM)
      Random::orient(mp_temp->write_ang_pos(), mp_temp->ang_pos_size());
    for (int j = 0; j < 3; ++j)
      mp_temp->write_cm_pos()[j] = 0.0;  //cm positions are not needed for now
    mp_temp->update_atoms();
  }

  //  reference points in laboratory frame
  double lfactor = 1.;
  double lf_ref_pos[2][3];
  if(mol_array[0]->type() == LINEAR && mol_array[1]->type() == LINEAR)
    for(int frag = 0; frag < 2; ++frag)
      for (int j = 0; j < 3; ++j)
	lf_ref_pos[frag][j] = mol_array[frag]->read_ang_pos()[j] * 
	  ds.ref_pos(frag)[0];
  else
    for (int frag = 0; frag < 2; ++frag) {
      mp_temp = mol_array[frag];
      switch (mp_temp->type()) {
      case ATOM:
	for (int j = 0; j < 3; ++j)
	  lf_ref_pos [frag] [j] = 0.0;
	break;

      case LINEAR:
	if(ds.ref_pos(frag)[1] > 0.) {// toroidal surface
	  for (int j = 0; j < 3; ++j)
	    pos_temp[j] = lf_cos_12[j];
	
	  orthogonalize(pos_temp, mp_temp->read_ang_pos(), 3);
	  dtemp = normalize(pos_temp, 3);

	  if(Random::flat() > 0.5)
	    itemp = 1;// internal part of the torus
	  else
	    itemp = 0;// external part of the torus
	
	  double u;
	  if(itemp)
	    u = ds.ref_pos(frag)[1] - dtemp*ds.dist();
	  else
	    u = ds.ref_pos(frag)[1] + dtemp*ds.dist();

	  if(u <= 0.)
	    return SAMP_FACE_OUT;

	  if(dtemp > 1.e-14)
	    lfactor = 2. * u / dtemp / ds.dist();
	  else 
	    lfactor = 0.;

	  for (int j = 0; j < 3; ++j) {
	    lf_ref_pos [frag] [j] = mp_temp->read_ang_pos()[j] * 
	      ds.ref_pos(frag)[0];
	    if((itemp+frag)%2)
	      lf_ref_pos [frag] [j] += pos_temp[j] * ds.ref_pos(frag)[1];
	    else
	      lf_ref_pos [frag] [j] -= pos_temp[j] * ds.ref_pos(frag)[1];
	  }
	}
	else // spherical surface
	  for (int j = 0; j < 3; ++j)
	    lf_ref_pos[frag][j] = mol_array[frag]->read_ang_pos()[j] * 
	      ds.ref_pos(frag)[0];
	break;

      case NONLINEAR:
	mp_temp->mf2lf(ds.ref_pos(frag), lf_ref_pos[frag]);
	break;

      default:
	error ("rand_pos_0: wrong type");
      }
    }

  // set cm to cm vector ...
  double pos_0 [3];
  for (int i = 0; i < 3; ++i)
    pos_0[i] = lf_cos_12[i] * ds.dist() - lf_ref_pos[0][i] + lf_ref_pos[1][i];

  // and update atoms coordinates
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 3; ++i)
      mol_array[j]->write_cm_pos()[i] = pos_0[i] * m_factor [j];
  for (int j = 0; j < 2; ++j)
    mol_array[j]->update_atoms();

  // check if the face is an external one
  for (int face = 0; face < ds.end(); ++face) {// external face
    if ( (ds.face() - face) % ds.size() == 0)
      continue;

    // vector between to ref. points in LF
    double tmp_lf_pos_12 [3];
    for (int i = 0; i < 3; ++i)
      tmp_lf_pos_12[i] = pos_0[i];

    if(mol_array[0]->type() == LINEAR && mol_array[1]->type() == LINEAR)
      // both linear fragments
      for(int i = 0; i < 3; ++i)
	tmp_lf_pos_12[i] +=  ds.ref_pos(0, ds.ref_ind(0, face))[0] * 
	  mol_array[0]->read_ang_pos()[i] -
	  ds.ref_pos(1, ds.ref_ind(1, face))[0] * 
	  mol_array[1]->read_ang_pos()[i];

    else {
      for (int frag = 0; frag < 2; ++frag)
	if(mol_array[frag]->type() == NONLINEAR) {
	  double tmp_lf_ref_pos [3];
	  mol_array[frag]->mf2lf(ds.ref_pos(frag, ds.ref_ind(frag, face)), 
				 tmp_lf_ref_pos);
	  if (frag)
	    for (int i = 0; i < 3; ++i)
	      tmp_lf_pos_12 [i] -= tmp_lf_ref_pos[i];
	  else
	    for (int i = 0; i < 3; ++i)
	      tmp_lf_pos_12 [i] += tmp_lf_ref_pos[i];
	}

      for (int frag = 0; frag < 2; ++frag)
	if(mol_array[frag]->type() == LINEAR)
	  if(ds.ref_pos(frag, ds.ref_ind(frag, face))[1] > 0.) {// torus
	    for(int i = 0; i < 3; ++i)
	      pos_temp[i] = tmp_lf_pos_12[i];

	    double p0 = vdot(pos_temp, mol_array[frag]->read_ang_pos(),3);
	    if(frag)
	      p0 -= ds.ref_pos(frag, ds.ref_ind(frag, face))[0];
	    else
	      p0 += ds.ref_pos(frag, ds.ref_ind(frag, face))[0];

	    orthogonalize(pos_temp, mol_array[frag]->read_ang_pos(),3);
	    double p1 = normalize(pos_temp,3);
	    p1 -= ds.ref_pos(frag, ds.ref_ind(frag, face))[1];

	    for(int i = 0; i < 3; ++i)
	      tmp_lf_pos_12[i] = p0 * mol_array[frag]->read_ang_pos()[i] +
		p1 * pos_temp[i];
	  }// torus
	  else if(frag) // sphere
	    for(int i = 0; i < 3; ++i)
	      tmp_lf_pos_12[i] -= mol_array[frag]->read_ang_pos()[i] *
		ds.ref_pos(frag, ds.ref_ind(frag, face))[0];
	  else
	    for(int i = 0; i < 3; ++i)
	      tmp_lf_pos_12[i] += mol_array[frag]->read_ang_pos()[i] *
		ds.ref_pos(frag, ds.ref_ind(frag, face))[0];
    }
    // the distance between two reference points is less than the minimum
    if (vlength(tmp_lf_pos_12, 3) < ds.dist(face))
      return SAMP_FACE_OUT;
  }// external face

  // minimum interatomic distance
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;
  
  if (weight) {// statistical weight

    // transform a vector between ref. points to molecular frame
    double rc_vec [9];    // reaction coordinate
    double* rcp = rc_vec;

    for (int i = 0; i < 3; ++i)// orbital coordinates
      *rcp++ = lf_cos_12[i] / emass2;

    double mf_cos_12 [3];
    for (int frag = 0; frag < 2; ++frag) {// internal coordinates
      mp_temp = mol_array[frag];
      switch (mp_temp->type()) {
      case ATOM:
	break;

      case LINEAR:
	if (frag == 0)
	  vector_product(lf_ref_pos [frag], lf_cos_12, rcp);
	else
	  vector_product(lf_cos_12, lf_ref_pos [frag], rcp);

	for (int j = 0; j < 3; ++j)
	  *rcp++ /= iner_mom2[frag][j];

	break;

      case NONLINEAR:
	mp_temp->lf2mf(lf_cos_12, mf_cos_12);
	if (frag == 0)
	  vector_product(ds.ref_pos(frag), mf_cos_12, rcp);
	else
	  vector_product(mf_cos_12, ds.ref_pos(frag), rcp);

	for (int j = 0; j < 3; ++j)
	  *rcp++ /= iner_mom2[frag][j];

	break;
      }
    }// internal coordinates

    // normalize reaction coordinate vector
    const int rc_dim = rcp - rc_vec;
    *weight = lfactor * ds.dist() * ds.dist() * normalize(rc_vec, rc_dim);

    // generalized velocity
    double gv_vec [9];
    for (int i = 0; i < rc_dim; ++i)
      gv_vec[i] = Random::norm();
    // adjust reactive component
    orthogonalize(gv_vec, rc_vec, rc_dim);
    dtemp = Random::exp();
    for (int i = 0; i < rc_dim; ++i)
      gv_vec [i] += dtemp * rc_vec[i];

    // setting random velocities

    // orbital
    const double* gvp = gv_vec;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 2; ++j)
	mol_array[j]->write_cm_vel()[i] = m_factor[j] * gvp[i] 
	  / emass2;
    gvp += 3;
    // internal
    for (int i = 0; i < 2; ++i) {
      mp_temp = mol_array[i];
      switch (mp_temp->type()) {
      case ATOM:
	break;
      case LINEAR:
	for (int j = 0; j < 3; ++j)
	  mp_temp->write_ang_vel()[j] = gvp[j] / iner_mom2[i][j];
	mp_temp->orthogonalize();
	gvp += 3;
	break;
      case NONLINEAR:
	for (int j = 0; j < 3; ++j)
	  mp_temp->write_ang_vel()[j] = gvp[j] / iner_mom2[i][j];
	gvp += 3;
	break;
      }
    }

  }// statistical weight

  return 0;
}

int rand_pos_1a (const Div_surf& ds, Dynvar& dynvar, double* weight)
{
  dynvar.attach();

  /********************************************************************
   * Sample surface which consists of two semispheres and a plane.
   * Indices: 0,1 - spheres; 2 - plane
   * Faces:   0,1 - spherical elements; 2 - plane element
   ********************************************************************/

  if (ds.ref_num(0) != 3 || ds.ref_num(1) != 1)
    error("rand_pos_1a: wrong dividing surface");

  if (mol_array[0]->type() != NONLINEAR)
    error("rand_pos_1a: wrong molecular specification");

  double n_0 [3];
  for(int i = 0; i < 3; ++i)
    n_0[i] = ds.ref_pos(0,2)[i];
  normalize(n_0, 3);

  const double disp_0 = ds.dist(2);

  static double rrmass =   1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 

  // temporary variables
  double dtemp, dtemp1;


  double temp_pos [3];
  double pos [3];
  double cm_pos [3];

  //standard orient the first fragment
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();

  // randomly orient second fragment
  if (mol_array[1]->type() != ATOM)
    Random::orient(mol_array[1]->write_ang_pos(), 
		   mol_array[1]->ang_pos_size());
  for (int j = 0; j < 3; ++j)
    mol_array[1] ->write_cm_pos()[j] = 0.0;
  mol_array[1]->update_atoms();

  //  reference point of the second fragment in laboratory frame
  double lf_ref_pos [3];
  switch (mol_array[1]->type()) {
  case ATOM:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = 0.0;
    break;
  case LINEAR:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = mol_array[1]->read_ang_pos()[j] * ds.ref_pos(1)[0];
    break;
  case NONLINEAR:
    mol_array[1]->mf2lf(ds.ref_pos(1), lf_ref_pos);
  }

  if (ds.face() == 0 || ds.face() == 1) {// sperical surface
    
    if (ds.dist() <= 0.)
      return SAMP_FACE_OUT;

    // random pivot-point-to-pivot-point direction
    double pos_12 [3]; // 0 -> 1
    Random::orient(pos_12, 3);

    for (int i = 0; i < 3; ++i)
      pos[i]    = ds.ref_pos(0)[i] + ds.dist() * pos_12[i];

    // working part of the sphere
    switch(ds.face()) {
    case 0:
      if (vdot(pos, n_0, 3) < disp_0)
	return SAMP_FACE_OUT;
      break;
    case 1:
      if (vdot(pos, n_0, 3) > disp_0)
	return SAMP_FACE_OUT;
      break;
    }

    for (int i = 0; i < 3; ++i)
      cm_pos[i] = pos[i] - lf_ref_pos[i];

    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = cm_pos[i];
    mol_array[1]->update_atoms();

    // minimum interatomic distance
    if(are_atoms_close(MIN_ATOM_DIST))
      return SAMP_ATOMS_CLOSE;

    // statistical weight
    if (weight) {
      // cm motion
      *weight = rrmass;

      // 1st fragment
      vector_product(pos_12, ds.ref_pos(0), temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);

      // 2nd fragment
      double mf_pos_12 [3];
      switch(mol_array[1]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dtemp = vdot(pos_12, mol_array[1]->read_ang_pos(), 3);
	*weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0] 
	  * (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
	break;
      case NONLINEAR:
	mol_array[1]->lf2mf(pos_12, mf_pos_12);
	vector_product(mf_pos_12, ds.ref_pos(1), temp_pos);
	for (int i = 0; i < 3; ++i)
	  *weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      }

      *weight = sqrt(*weight) * ds.dist() * ds.dist();
    }

    return 0;

  }// spherical surface

  /*********************** plain surface **************************/

  if (ds.face() != 2)
    error("rand_pos_1a: wrong face");

  double circle_pos [2] [3]; // center of the circle
  double circle_rad [2];     // circle radius
  
  for (int i1 = 0; i1 < 2; ++i1) {

    dtemp = vdot(n_0, ds.ref_pos(0, i1), 3) - disp_0;
    for (int i = 0; i < 3; ++i)
      circle_pos [i1] [i] = ds.ref_pos(0, i1)[i] - dtemp * n_0[i];

    dtemp1 = ds.dist(i1);
    if (dtemp1 > 0) {
      circle_rad[i1] = dtemp1 * dtemp1 - dtemp * dtemp;
      if (circle_rad[i1] > 0.)
	circle_rad[i1] = sqrt(circle_rad[i1]);
      else
	circle_rad[i1] = -1.;
    }
    else
      circle_rad[i1] = -1.;

  }

  // circle with bigest radius
  int major = circle_rad[0] < circle_rad[1] ? 1 : 0;
  int minor = 1 - major;

  if (circle_rad[major] <= 0.)
    return SAMP_FACE_OUT;

  double dd [3];
  for (int i = 0; i < 3; ++i)
    dd[i] = circle_pos[1][i] - circle_pos[0][i];
  double ddl = normalize(dd, 3);

  if (ddl < 1.e-4 && circle_rad[major] - circle_rad[minor] < 1.e-4 )
    return SAMP_FACE_OUT;
    
  /************************** One circle **************************/

  if (circle_rad[minor] <= 0. || 
      circle_rad[major] >= ddl + circle_rad[minor]) {// one circle

    // first orth
    double orth1 [3] = {0., 0., 1.};
    orthogonalize(orth1, n_0, 3);
    normalize(orth1, 3);

    // second ort
    double orth2 [3];
    vector_product(n_0, orth1, orth2);

    // random position
    do {
      double c1 = 2. * Random::flat() - 1.;
      double c2 = 2. * Random::flat() - 1.;
      for (int i = 0; i < 3; ++i)
	pos[i] = circle_pos[major] [i] + circle_rad[major] 
	  * (c1 * orth1[i] + c2 * orth2[i]);

      dtemp = 1.;
      if (circle_rad[minor] > 0.) {
	dtemp = 0.;
	for (int i = 0; i < 3; ++i) {
	  dtemp1 = pos[i] - circle_pos[minor][i];
	  dtemp += dtemp1 * dtemp1;
	}
	dtemp -= circle_rad[minor] * circle_rad[minor];
      }

      dtemp1 = c1 * c1 + c2 * c2 - 1.;

    } while (dtemp1 > 0. || dtemp < 0.);

    // cm-to-cm vector
    for (int i = 0; i < 3; ++i)
      cm_pos[i] = pos[i] - lf_ref_pos[i];

    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = cm_pos[i];
    mol_array[1]->update_atoms();

    mol_array[0]->init_dv();
    mol_array[0]->update_atoms();

    // minimum interatomic distance
    if(are_atoms_close(MIN_ATOM_DIST))
      return SAMP_ATOMS_CLOSE;

    // statistical weight
    if(weight) {
      *weight = rrmass;
      // 1st fragment
      vector_product(pos, n_0, temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
      // 2nd fragment
      double mf_n [3];
      switch(mol_array[1]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dtemp = vdot(n_0, mol_array[1]->read_ang_pos(), 3);
	*weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0]
	  * (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
	break;
      case NONLINEAR:
	mol_array[1]->lf2mf(n_0, mf_n);
	vector_product(ds.ref_pos(1), mf_n, temp_pos);
	for (int i = 0; i < 3; ++i)
	  *weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
	break;
      }

      *weight = sqrt(*weight);

      dtemp = circle_rad[major] * circle_rad[major];
      if (circle_rad[minor] > 0.)
	dtemp -= circle_rad[minor] * circle_rad[minor];
      *weight *= dtemp / 4.; // 4 due to 1/(4 pi)
    }

    return 0;
    
  }// one circle

  /********************** two circles ************************/

  // third ort
  double orth2 [3];
  vector_product(dd, n_0, orth2);

  // sampling circle
  double smp_rad = (circle_rad[0] + circle_rad[1] + ddl) / 2.;
  double smp_pos [3];
  for (int i = 0; i < 3; ++i)
    smp_pos[i] = circle_pos[0][i] + (circle_rad[1] + ddl - circle_rad[0])
      /  2. * dd[i];

  // random position in sampling circle
  double c1, c2;
  do {
    c1 = 2. * Random::flat() - 1.;
    c2 = 2. * Random::flat() - 1.;
  } while (c1 * c1 + c2 * c2 >= 1.);
  double rel_pos[3];
  for (int i = 0; i < 3; ++i)
    rel_pos[i] = smp_rad * (c1 * dd[i] + c2 * orth2[i]);

  // distances from cirlces centers
  double r [2];
  for (int i = 0; i < 2; ++i) {
    r[i] = 0.;
    for (int j = 0; j < 3; ++j) {
      if (!i)
	dtemp = rel_pos[j] + (circle_rad[1] + ddl - circle_rad[0]) 
	  / 2. * dd[j];
      else
	dtemp = rel_pos[j] + (circle_rad[1] - ddl - circle_rad[0]) 
	  / 2. * dd[j];
      r[i] += dtemp * dtemp;
    }
    r[i] = sqrt(r[i]);
  }

  // is inside or outside of both circles?
  if ( r[0] > circle_rad[0] && r[1] > circle_rad[1] ||
       r[0] < circle_rad[0] && r[1] < circle_rad[1] )
    return SAMP_FACE_OUT;

  // atom position
  for (int i = 0; i < 3; ++i)
    pos[i] = smp_pos[i] + rel_pos[i];

  // cm-to-cm vector
  for (int i = 0; i < 3; ++i)
    cm_pos[i] = pos[i] - lf_ref_pos[i];
  
  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = cm_pos[i];
  mol_array[1]->update_atoms();
  
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();

  // minimum interatomic distance
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;

  // statistical weight
  if(weight) {
    *weight = rrmass;
    // 1st fragment
    vector_product(pos, n_0, temp_pos);
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
    // 2nd fragment
    double mf_n [3];
    switch(mol_array[1]->type()) {
    case ATOM:
      break;
    case LINEAR:
      dtemp = vdot(n_0, mol_array[1]->read_ang_pos(), 3);
      *weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0]
	* (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
      break;
    case NONLINEAR:
      mol_array[1]->lf2mf(n_0, mf_n);
      vector_product(ds.ref_pos(1), mf_n, temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      break;
    }

    *weight = sqrt(*weight);

    *weight *= smp_rad * smp_rad / 4.;
  }

  return 0;
}

int rand_pos_1b (const Div_surf& ds, Dynvar& dynvar, double* weight)
{
  dynvar.attach();

  /********************************************************************
   * Sample surface which consists of two semispheres and several planes.
   * 0 - outer sphere, 1 - inner sphere, 2, ... - planes
   ********************************************************************/

  if (ds.ref_num(0) < 3 || ds.ref_num(1) != 1)
    error("rand_pos_1b: wrong dividing surface");

  if (mol_array[0]->type() != NONLINEAR)
    error("rand_pos_1b: wrong molecular specification");

  static const double rrmass =   
    1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 

  // temporary variables
  double dtemp, dtemp1;

  double temp_pos [3];
  double pos [3];
  double cm_pos [3];

  //standard orient the first fragment
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();

  // randomly orient second fragment
  if (mol_array[1]->type() != ATOM)
    Random::orient(mol_array[1]->write_ang_pos(), 
		   mol_array[1]->ang_pos_size());
  for (int j = 0; j < 3; ++j)
    mol_array[1] ->write_cm_pos()[j] = 0.0;
  mol_array[1]->update_atoms();

  //  reference point of the second fragment in laboratory frame
  double lf_ref_pos [3];
  switch (mol_array[1]->type()) {
  case ATOM:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = 0.0;
    break;
  case LINEAR:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = mol_array[1]->read_ang_pos()[j] * ds.ref_pos(1)[0];
    break;
  case NONLINEAR:
    mol_array[1]->mf2lf(ds.ref_pos(1), lf_ref_pos);
  }

  /*********************** Spherical Surface **********************/

  if (ds.face() == 0 || ds.face() == 1) {// sperical surface
    
    if (ds.dist() <= 0.)
      return SAMP_FACE_OUT;

    // random pivot-point-to-pivot-point direction
    double pos_12 [3]; // 0 -> 1
    Random::orient(pos_12, 3);

    for (int i = 0; i < 3; ++i)
      pos[i]    = ds.ref_pos(0)[i] + ds.dist() * pos_12[i];

    // working part of the sphere
    for(int i = 2; i < ds.ref_num(0); ++i)
      if (vdot(pos, ds.ref_pos(0,i), 3) < ds.dist(i))
	return SAMP_FACE_OUT;

    for (int i = 0; i < 3; ++i)
      cm_pos[i] = pos[i] - lf_ref_pos[i];

    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = cm_pos[i];
    mol_array[1]->update_atoms();

    // minimum interatomic distance
    if(are_atoms_close(MIN_ATOM_DIST))
      return SAMP_ATOMS_CLOSE;

    // statistical weight
    if (weight) {
      // cm motion
      *weight = rrmass;

      // 1st fragment
      vector_product(pos_12, ds.ref_pos(0), temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);

      // 2nd fragment
      double mf_pos_12 [3];
      switch(mol_array[1]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dtemp = vdot(pos_12, mol_array[1]->read_ang_pos(), 3);
	*weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0] 
	  * (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
	break;
      case NONLINEAR:
	mol_array[1]->lf2mf(pos_12, mf_pos_12);
	vector_product(mf_pos_12, ds.ref_pos(1), temp_pos);
	for (int i = 0; i < 3; ++i)
	  *weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      }

      *weight = sqrt(*weight) * ds.dist() * ds.dist();
    }

    return 0;

  }// spherical surface


  /*********************** plain surface **************************/

  double circle_pos [2] [3]; // center of the circle
  double circle_rad [2];     // circle radius
  const double norm = vdot(ds.ref_pos(0), ds.ref_pos(0), 3);

  for (int i1 = 0; i1 < 2; ++i1) {
    dtemp = (vdot(ds.ref_pos(0), ds.ref_pos(0, i1), 3) - ds.dist()) /
      norm;
    circle_rad[i1] = ds.dist(i1) * ds.dist(i1) - dtemp * dtemp * norm;
    if (circle_rad[i1] > 0.)
      circle_rad[i1] = sqrt(circle_rad[i1]);
    else
      circle_rad[i1] = -1.;

    for(int i = 0; i < 3; ++i)
      circle_pos [i1] [i] = ds.ref_pos(0, i1)[i] - dtemp * ds.ref_pos(0)[i];
  }

  if (circle_rad[0] <= 0.) {
    //std::cout << "negative radius\n";
    return SAMP_FACE_OUT;
  }

  // orts
  int imax = 0;
  dtemp = ds.ref_pos(0)[0] > 0. ? ds.ref_pos(0)[0] : -ds.ref_pos(0)[0];
  for(int i = 1; i < 3; ++i) {
    dtemp1 = ds.ref_pos(0)[i] > 0. ? ds.ref_pos(0)[i] : -ds.ref_pos(0)[i];
    if(dtemp1 > dtemp) {
      imax = i;
      dtemp = dtemp1;
    }
  }
  double oy[3];
  dtemp = 0.;
  for(int i = 0; i < 3; ++i)
    if(i != imax) {
      oy[i] = 1.;
      dtemp += ds.ref_pos(0)[i];
    }
  oy[imax] = -dtemp/ds.ref_pos(0)[imax];
  normalize(oy, 3);

  double oz[3];
  vector_product(ds.ref_pos(0), oy, oz);
  normalize(oz, 3);

  // random position in outer circle
  double c1, c2;
  do {
    c1 = 2. * Random::flat() - 1.;
    c2 = 2. * Random::flat() - 1.;
  } while (c1 * c1 + c2 * c2 >= 1.);

  dtemp = 0.;
  for (int i = 0; i < 3; ++i) { 
    pos[i] = circle_rad[0] * (c1 * oy[i] + c2 * oz[i]) + circle_pos[0][i];
    dtemp1 = pos[i] - circle_pos[1][i];
    dtemp += dtemp1 * dtemp1;
  }
  // is inside of inner circle?
  if(circle_rad[1] > 0. && dtemp < circle_rad[1] * circle_rad[1]) {
    //std::cout << "inner circle\n";
    return SAMP_FACE_OUT;
  }

  // working part of the plane
  for(int i = 2; i < ds.ref_num(0); ++i)
    if(i != ds.face() && vdot(pos, ds.ref_pos(0,i), 3) < ds.dist(i)) {
      //std::cout << "wrong plane side\n";
      return SAMP_FACE_OUT;
    }

  // cm-to-cm vector
  for (int i = 0; i < 3; ++i)
    cm_pos[i] = pos[i] - lf_ref_pos[i];
  
  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = cm_pos[i];
  mol_array[1]->update_atoms();
  
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();

  // minimum interatomic distance
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;

  // statistical weight
  if(weight) {
    *weight = rrmass;
    // 1st fragment
    vector_product(pos, ds.ref_pos(0), temp_pos);
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i]/norm/mol_array[0]->iner_mom(i);
    // 2nd fragment
    double mf_n [3];
    switch(mol_array[1]->type()) {
    case ATOM:
      break;
    case LINEAR:
      dtemp = vdot(ds.ref_pos(0), mol_array[1]->read_ang_pos(), 3);
      *weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0]
	* (1. - dtemp*dtemp/norm) / mol_array[1]->iner_mom(2);
      break;
    case NONLINEAR:
      mol_array[1]->lf2mf(ds.ref_pos(0), mf_n);
      vector_product(ds.ref_pos(1), mf_n, temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i]*temp_pos[i]/norm/mol_array[1]->iner_mom(i);
      break;
    }
    *weight = sqrt(*weight);

    *weight *= circle_rad[0] * circle_rad[0] / 4.;
  }

  return 0;
}

int rand_pos_2a (const Div_surf& ds_0, Dynvar& dynvar, double* weight)
{
  /********************************************************************
   * Sample surface which consists of two semispheres and two plane.
   * Indices: 0,1 - spheres; 2,3 - plane
   * Faces:   0,1 - spherical elements; 2,3 - plane element
   ********************************************************************/

  dynvar.attach();

  Div_surf ds = ds_0;

  if (ds.ref_num(0) != 4 || ds.ref_num(1) != 1)
    error("rand_pos_2a: wrong dividing surface");

  if (mol_array[0]->type() != NONLINEAR)
    error("rand_pos_2a: wrong molecular specification");

  const double* n_0 = ds.ref_pos(0, 2);
  double disp_0 = ds.dist(2);
  const double* n_1 = ds.ref_pos(0, 3);
  double disp_1 = ds.dist(3);

  normalize(ds.write_ref_pos(0, 2), 3);
  normalize(ds.write_ref_pos(0, 3), 3);

  static double rrmass =   1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 

  // temporary variables
  double dtemp, dtemp1;


  double temp_pos [3];
  double pos [3];
  double cm_pos [3];    // cm-to-cm vector

  // standard orientation of the first fragment
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();
 
  // randomly orient second fragment
  if (mol_array[1]->type() != ATOM)
    Random::orient(mol_array[1]->write_ang_pos(), 
		   mol_array[1]->ang_pos_size());
  for (int j = 0; j < 3; ++j)
    mol_array[1] ->write_cm_pos()[j] = 0.0;
  mol_array[1]->update_atoms();

  //  reference point of the second fragment in laboratory frame
  double lf_ref_pos [3];
  switch (mol_array[1]->type()) {
  case ATOM:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = 0.0;
    break;
  case LINEAR:
    for (int j = 0; j < 3; ++j)
      lf_ref_pos[j] = mol_array[1]->read_ang_pos()[j] * ds.ref_pos(1)[0];
    break;
  case NONLINEAR:
    mol_array[1]->mf2lf(ds.ref_pos(1), lf_ref_pos);
  }

  /*************************** Sperical surfaces *****************************/

  if (ds.face() == 0 || ds.face() == 1) {// sperical surface
    
    if (ds.dist() <= 0.) // negative radius
      return SAMP_FACE_OUT;

    // randomly choose postion on a sphere
    double pos_12 [3];
    Random::orient(pos_12, 3);

    for (int i = 0; i < 3; ++i)
      pos[i]    = ds.ref_pos(0)[i] + ds.dist() * pos_12[i];

    // working part of the sphere
    switch(ds.face()) {
    case 0:
      if (vdot(pos, n_0, 3) < disp_0 || vdot(pos, n_1, 3) < disp_1)
	return SAMP_FACE_OUT;
      break;
    case 1:
      if (vdot(pos, n_0, 3) > disp_0 && vdot(pos, n_1, 3) > disp_1)
	return SAMP_FACE_OUT;
      break;
    }

    // cm-to-cm vector, 0 -> 1
    for (int i = 0; i < 3; ++i)
      cm_pos[i] = pos[i] - lf_ref_pos[i];
    
    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = cm_pos[i];
    mol_array[1]->update_atoms();

    // minimum interatomic distance
    if(are_atoms_close(MIN_ATOM_DIST))
      return SAMP_ATOMS_CLOSE;

    // statistical weight
    if(weight) {
      // cm motion
      *weight = rrmass;

      // 1st fragment
      vector_product(pos_12, ds.ref_pos(0), temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);

      // 2nd fragment
      double mf_pos_12 [3];
      switch(mol_array[1]->type()) {
      case ATOM:
	break;
      case LINEAR:
	dtemp = vdot(pos_12, mol_array[1]->read_ang_pos(), 3);
	*weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0] 
	  * (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
	break;
      case NONLINEAR:
	mol_array[1]->lf2mf(pos_12, mf_pos_12);
	vector_product(mf_pos_12, ds.ref_pos(1), temp_pos);
	for (int i = 0; i < 3; ++i)
	  *weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      }

      *weight = sqrt(*weight) * ds.dist() * ds.dist();
    }

    return 0;

  }// spherical surface

  /*********************** plain surface **************************/

  switch (ds.face()) {

  case 2:
    break;

  case 3:
    n_0 = ds.ref_pos(0, 3);
    disp_0 = ds.dist(3);
    n_1 = ds.ref_pos(0, 2);
    disp_1 = ds.dist(2);
    break;

  default:
    error("rand_pos_2a: wrong face");
  }

  double circle_pos [2] [3]; // center of the circle
  double circle_rad [2];     // circle radius
  
  for (int i1 = 0; i1 < 2; ++i1) {

    dtemp = vdot(n_0, ds.ref_pos(0, i1), 3) - disp_0;
    for (int i = 0; i < 3; ++i)
      circle_pos [i1] [i] = ds.ref_pos(0, i1)[i] - dtemp * n_0[i];

    dtemp1 = ds.dist(i1);
    if (dtemp1 > 0) {
      circle_rad[i1] = dtemp1 * dtemp1 - dtemp * dtemp;
      if (circle_rad[i1] > 0.)
	circle_rad[i1] = sqrt(circle_rad[i1]);
      else
	circle_rad[i1] = -1.;
    }
    else
      circle_rad[i1] = -1.;

  }

  // circle with bigest radius
  int major = circle_rad[0] < circle_rad[1] ? 1 : 0;
  int minor = 1 - major;

  if (circle_rad[major] <= 0.)
    return SAMP_FACE_OUT;

  double dd [3];
  for (int i = 0; i < 3; ++i)
    dd[i] = circle_pos[1][i] - circle_pos[0][i];
  double ddl = normalize(dd, 3);

  if (ddl < 1.e-4 && circle_rad[major] - circle_rad[minor] < 1.e-4 )
    return SAMP_FACE_OUT;
    
  /************************** One circle **************************/

  if (circle_rad[minor] <= 0. || 
      circle_rad[major] >= ddl + circle_rad[minor]) {// one circle

    // first orth
    double orth1 [3] = {0., 0., 1.};
    orthogonalize(orth1, n_0, 3);
    normalize(orth1, 3);

    // second ort
    double orth2 [3];
    vector_product(n_0, orth1, orth2);

    // random position
    do {
      double c1 = 2. * Random::flat() - 1.;
      double c2 = 2. * Random::flat() - 1.;
      for (int i = 0; i < 3; ++i)
	pos[i] = circle_pos[major] [i] + circle_rad[major] 
	  * (c1 * orth1[i] + c2 * orth2[i]);

      dtemp = 1.;
      if (circle_rad[minor] > 0.) {
	dtemp = 0.;
	for (int i = 0; i < 3; ++i) {
	  dtemp1 = pos[i] - circle_pos[minor][i];
	  dtemp += dtemp1 * dtemp1;
	}
	dtemp -= circle_rad[minor] * circle_rad[minor];
      }

      dtemp1 = c1 * c1 + c2 * c2 - 1.;

    } while (dtemp1 > 0. || dtemp < 0.);

    // is right part of the plane?
    if ( vdot(pos, n_1, 3) < disp_1 )
      return SAMP_FACE_OUT;

    // cm-to-cm vector, 0 -> 1
    for (int i = 0; i < 3; ++i)
      cm_pos[i] = pos[i] - lf_ref_pos[i];
    
    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = cm_pos[i];
    mol_array[1]->update_atoms();

    mol_array[0]->init_dv();
    mol_array[0]->update_atoms();
 
    // minimum interatomic distance
    if(are_atoms_close(MIN_ATOM_DIST))
      return SAMP_ATOMS_CLOSE;

    // statistical weight
    if(weight) {
    *weight = rrmass;
    // 1st fragment
    vector_product(pos, n_0, temp_pos);
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
    // 2nd fragment
    double mf_n [3];
    switch(mol_array[1]->type()) {
    case ATOM:
      break;
    case LINEAR:
      dtemp = vdot(n_0, mol_array[1]->read_ang_pos(), 3);
      *weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0]
	* (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
      break;
    case NONLINEAR:
      mol_array[1]->lf2mf(n_0, mf_n);
      vector_product(ds.ref_pos(1), mf_n, temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      break;
    }

    *weight = sqrt(*weight);

    dtemp = circle_rad[major] * circle_rad[major];
    if (circle_rad[minor] > 0.)
      dtemp -= circle_rad[minor] * circle_rad[minor];
    *weight *= dtemp / 4.; // 4 due to 1/(4 pi)
    }

    return 0;
    
  }// one circle

  /********************** two circles ************************/

  // third ort
  double orth2 [3];
  vector_product(dd, n_0, orth2);

  // sampling circle
  double smp_rad = (circle_rad[0] + circle_rad[1] + ddl) / 2.;
  double smp_pos [3];
  for (int i = 0; i < 3; ++i)
    smp_pos[i] = circle_pos[0][i] + (circle_rad[1] + ddl - circle_rad[0])
      /  2. * dd[i];

  // random position in sampling circle
  double c1, c2;
  do {
    c1 = 2. * Random::flat() - 1.;
    c2 = 2. * Random::flat() - 1.;
    } while (c1 * c1 + c2 * c2 >= 1.);
  double rel_pos[3];
  for (int i = 0; i < 3; ++i)
    rel_pos[i] = smp_rad * (c1 * dd[i] + c2 * orth2[i]);

  // distances from cirlces centers
  double r [2];
  for (int i = 0; i < 2; ++i) {
    r[i] = 0.;
    for (int j = 0; j < 3; ++j) {
      if (!i)
	dtemp = rel_pos[j] + (circle_rad[1] + ddl - circle_rad[0]) 
	  / 2. * dd[j];
      else
	dtemp = rel_pos[j] + (circle_rad[1] - ddl - circle_rad[0]) 
	  / 2. * dd[j];
      r[i] += dtemp * dtemp;
    }
    r[i] = sqrt(r[i]);
  }

  // is inside or outside of both circles?
  if ( r[0] > circle_rad[0] && r[1] > circle_rad[1] ||
       r[0] < circle_rad[0] && r[1] < circle_rad[1] )
    return SAMP_FACE_OUT;

  // atom position
  for (int i = 0; i < 3; ++i)
    pos[i] = smp_pos[i] + rel_pos[i];

  // is right part of the plane?
  if ( vdot(pos, n_1, 3) < disp_1 )
    return SAMP_FACE_OUT;

  // cm-to-cm vector, 0 -> 1
  for (int i = 0; i < 3; ++i)
    cm_pos[i] = pos[i] - lf_ref_pos[i];
    
  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = cm_pos[i];
  mol_array[1]->update_atoms();
    
  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();
 
  // minimum interatomic distance
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;

  // statistical weight
  if(weight) {

    vector_product(pos, n_0, temp_pos);
    //relative motion
    *weight = rrmass;
    
    // 1st fragment
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);

    // 2nd fragment
    double mf_n [3];
    switch(mol_array[1]->type()) {
    case ATOM:
      break;
    case LINEAR:
      dtemp = vdot(n_0, mol_array[1]->read_ang_pos(), 3);
      *weight += ds.ref_pos(1)[0] * ds.ref_pos(1)[0]
	* (1. - dtemp * dtemp) / mol_array[1]->iner_mom(2);
      break;
    case NONLINEAR:
      mol_array[1]->lf2mf(n_0, mf_n);
      vector_product(ds.ref_pos(1), mf_n, temp_pos);
      for (int i = 0; i < 3; ++i)
	*weight += temp_pos[i] * temp_pos[i] / mol_array[1]->iner_mom(i);
      break;
    }

    *weight = sqrt(*weight);

    *weight *= smp_rad * smp_rad / 4.;
  }

  return 0;
}

int _rand_pos_1b (const Div_surf& ds, Dynvar& dynvar, double* weight)
{
  /********************************************************************
   * Sample surface which consists of two semispheres and a plane.
   * Indices: 0,1 - spheres; 2 - plane
   * Faces:   0,1 - spherical elements; 2,3 - plane elements
   ********************************************************************/

  dynvar.attach();

  if (ds.ref_num(0) != 4 || ds.ref_num(1) != 1)
    error("rand_pos_1b: wrong dividing surface");

  if (mol_array[0]->type() != NONLINEAR || mol_array[1]->type() != ATOM)
    error("rand_pos_1b: wrong molecular specification");

  const double* const n_0 = ds.ref_pos(0, 2);
  const double disp_0 = ds.dist(2);

  if (fabs(vdot(n_0, n_0, 3) - 1.0) > 1.e-10)
    error("rand_pos_1b: surface orth is not normalized");

  static const double rrmass =   1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 
 
  // temporary variables
  double dtemp, dtemp1;


  double temp_pos [3];
  double pos [3];

  if (ds.face() == 0 || ds.face() == 1) {// sperical surfaces
    
    if (ds.dist() <= 0.)
      return SAMP_FACE_OUT;

    // randomly choose postion on a sphere
    double pos_12 [3];
    Random::orient(pos_12, 3);
    for (int i = 0; i < 3; ++i) {
      pos[i] = ds.ref_pos(0)[i] + ds.dist() * pos_12[i];
    }

    // working part of the sphere
    switch(ds.face())
      {
      case 0:
	if (vdot(pos, n_0, 3) < disp_0)
	  return SAMP_FACE_OUT;
	break;
      case 1:
	if (vdot(pos, n_0, 3) > disp_0)
	  return SAMP_FACE_OUT;
	break;
      }

    mol_array[0]->init_dv();
    mol_array[0]->update_atoms();
 
    // minimum interatomic distance
    for (Ater at0 = mol_array[0]->begin(); at0 != mol_array[0]->end(); ++at0)
      if(vdistance(at0->mf_pos, pos, 3) < MIN_ATOM_DIST)
	return SAMP_ATOMS_CLOSE;

    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = pos[i];
    mol_array[1]->update_atoms();

    // statistical weight
    if (weight) {
    vector_product(pos_12, ds.ref_pos(0), temp_pos);

    *weight = rrmass;
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
    *weight = sqrt(*weight);

    *weight *= ds.dist() * ds.dist();
    }

    return 0;

  }// spherical surface

  /*********************** plane surface **************************/

  if (ds.face() != 2 && ds.face() != 3)
    error("rand_pos_1b: wrong face");

  // major is internal circle
  int major = ds.face() == 2 ? 0 : 1;
  int minor = 1 - major;

  double circle_pos [2] [3]; // center of the circle
  double circle_rad [2];     // circle radius
  
  for (int i1 = 0; i1 < 2; ++i1) {

    dtemp = vdot(n_0, ds.ref_pos(0, i1), 3) - disp_0;
    for (int i = 0; i < 3; ++i)
      circle_pos [i1] [i] = ds.ref_pos(0, i1)[i] - dtemp * n_0[i];

    dtemp1 = ds.dist(i1);
    if (dtemp1 > 0) {
      circle_rad[i1] = dtemp1 * dtemp1 - dtemp * dtemp;
      if (circle_rad[i1] > 0.)
	circle_rad[i1] = sqrt(circle_rad[i1]);
    }
    else
      circle_rad[i1] = -1.;

  }

  if (circle_rad[major] <= 0.)
    return SAMP_FACE_OUT;

  // first orth
  double orth1 [3] = {0., 0., 1.};
  orthogonalize(orth1, n_0, 3);
  normalize(orth1, 3);

  // second ort
  double orth2 [3];
  vector_product(n_0, orth1, orth2);

  // random position
  double c1, c2;
  do {
    c1 = 2. * Random::flat() - 1.;
    c2 = 2. * Random::flat() - 1.;
  } while (c1 * c1 + c2 * c2 > 1.0);

  for (int i = 0; i < 3; ++i)
    pos[i] = circle_pos[major] [i] + circle_rad[major] 
      * (c1 * orth1[i] + c2 * orth2[i]);

  // check if outside minor circle
  if (circle_rad[minor] > 0.) {
    dtemp = 0.;
    for (int i = 0; i < 3; ++i) {
      dtemp1 = pos[i] - circle_pos[minor][i];
      dtemp += dtemp1 * dtemp1;
    }
    if (dtemp < circle_rad[minor] * circle_rad[minor])
      return SAMP_FACE_OUT;
  }

  // minimum interatomic distance
  for (Ater at0 = mol_array[0]->begin(); at0 != mol_array[0]->end(); ++at0)
    if(vdistance(at0->mf_pos, pos, 3) < MIN_ATOM_DIST)
      return SAMP_ATOMS_CLOSE;

  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = pos[i];
  mol_array[1]->update_atoms();

  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();
 
  // statistical weight
  if (weight) {
  vector_product(pos, n_0, temp_pos);

  *weight = rrmass;
  for (int i = 0; i < 3; ++i)
    *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
  *weight = sqrt(*weight);

  *weight *= circle_rad[major] * circle_rad[major] / 4.; // 4 due to 1/(4 pi)
  }

  return 0;
    
}

int rand_pos_2b (const Div_surf& ds, Dynvar& dynvar, double* weight)
{
  /********************************************************************
   * Sample surface which consists of two semispheres and two planes.
   * Indices: 0,1 - spheres; 2,3 - planes
   * Faces:   0,1 - spherical elements; 2,3,4,5 - plane elements
   ********************************************************************/

  dynvar.attach();

  if (ds.ref_num(0) != 6 || ds.ref_num(1) != 1)
    error("rand_pos_2b: wrong dividing surface");

  if (mol_array[0]->type() != NONLINEAR || mol_array[1]->type() != ATOM)
    error("rand_pos_2b: wrong molecular specification");

  const double *n_0 = ds.ref_pos(0, 2); // dividing surface plane normal
  const double *n_1 = ds.ref_pos(0, 3); // other plane normal
  double disp_0 = ds.dist(2), disp_1 = ds.dist(3); // plane shifts

  if (fabs(vdot(n_0, n_0, 3) - 1.0) > 1.e-10 ||
      fabs(vdot(n_1, n_1, 3) - 1.0) > 1.e-10)
    error("rand_pos_2b: surface orth is not normalized");

  static const double rrmass =   1./mol_array[0]->mass() + 1./mol_array[1]->mass(); 
 
  // temporary variables
  double dtemp, dtemp1;


  double temp_pos [3];
  double pos [3];

  if (ds.face() == 0 || ds.face() == 1) {// sperical surface
    
    if (ds.dist() <= 0.)
      return SAMP_FACE_OUT;

    // randomly choose postion on a sphere
    double pos_12 [3];
    Random::orient(pos_12, 3);
    for (int i = 0; i < 3; ++i) {
      pos[i] = ds.ref_pos(0)[i] + ds.dist() * pos_12[i];
    }

    // working part of the sphere
    switch(ds.face())
      {
      case 0: // spherical element inside the angle formed by two planes

	if (vdot(pos, n_0, 3) < disp_0 || vdot(pos, n_1, 3) < disp_1)
	  return SAMP_FACE_OUT;
	break;

      case 1: // spherical element outside the angle between two planes

	if (vdot(pos, n_0, 3) > disp_0 && vdot(pos, n_1, 3) > disp_1)
	  return SAMP_FACE_OUT;
	break;
      }

    // minimum interatomic distance
    for (Ater at0 = mol_array[0]->begin(); at0 != mol_array[0]->end(); ++at0)
      if(vdistance(at0->mf_pos, pos, 3) < MIN_ATOM_DIST)
	return SAMP_ATOMS_CLOSE;

    // update atom coordinates
    for (int i = 0; i < 3; ++i)
      mol_array[1]->write_cm_pos()[i] = pos[i];
    mol_array[1]->update_atoms();

    mol_array[0]->init_dv();
    mol_array[0]->update_atoms();
 
    // statistical weight
    if (weight) {
    vector_product(pos_12, ds.ref_pos(0), temp_pos);

    *weight = rrmass;
    for (int i = 0; i < 3; ++i)
      *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
    *weight = sqrt(*weight);

    *weight *= ds.dist() * ds.dist();
    }

    return 0;

  }// spherical surface

  /*********************** plane surfaces **************************/

  int major = 0, minor = 1;             // major internal circle

  switch (ds.face()) {

  case 2:
    break;

  case 3:
    major = 1;
    minor = 0;
    break;

  case 4:
    n_0 = ds.ref_pos(0, 3);
    n_1 = ds.ref_pos(0, 2);
    disp_0 = ds.dist(3);
    disp_1 = ds.dist(2);
    break;

  case 5:
    major = 1;
    minor = 0;
    n_0 = ds.ref_pos(0, 3);
    n_1 = ds.ref_pos(0, 2);
    disp_0 = ds.dist(3);
    disp_1 = ds.dist(2);
    break;

  default:
    error("rand_pos_2b: wrong face");
  }

  double circle_pos [2] [3]; // center of the circle
  double circle_rad [2];     // circle radius
  
  for (int i1 = 0; i1 < 2; ++i1) {

    dtemp = vdot(n_0, ds.ref_pos(0, i1), 3) - disp_0;
    for (int i = 0; i < 3; ++i)
      circle_pos [i1] [i] = ds.ref_pos(0, i1)[i] - dtemp * n_0[i];

    dtemp1 = ds.dist(i1);
    if (dtemp1 > 0) {
      circle_rad[i1] = dtemp1 * dtemp1 - dtemp * dtemp;
      if (circle_rad[i1] > 0.)
	circle_rad[i1] = sqrt(circle_rad[i1]);
    }
    else
      circle_rad[i1] = -1.;

  }

  if (circle_rad[major] <= 0.)
    return SAMP_FACE_OUT;

  // first orth
  double orth1 [3] = {0., 0., 1.};
  orthogonalize(orth1, n_0, 3);
  normalize(orth1, 3);

  // second ort
  double orth2 [3];
  vector_product(n_0, orth1, orth2);

  // random position
  double c1, c2;
  do {
    c1 = 2. * Random::flat() - 1.;
    c2 = 2. * Random::flat() - 1.;
  } while (c1 * c1 + c2 * c2 > 1.0);

  for (int i = 0; i < 3; ++i)
    pos[i] = circle_pos[major] [i] + circle_rad[major] 
      * (c1 * orth1[i] + c2 * orth2[i]);

  // check if outside minor circle
  if (circle_rad[minor] > 0.) {
    dtemp = 0.;
    for (int i = 0; i < 3; ++i) {
      dtemp1 = pos[i] - circle_pos[minor][i];
      dtemp += dtemp1 * dtemp1;
    }
    if (dtemp < circle_rad[minor] * circle_rad[minor])
      return SAMP_FACE_OUT;
  }

  // check if on the right part of the plane
  if(vdot(n_1, pos, 3) < disp_1)
    return SAMP_FACE_OUT;

  // minimum interatomic distance
  for (Ater at0 = mol_array[0]->begin(); at0 != mol_array[0]->end(); ++at0)
    if(vdistance(at0->mf_pos, pos, 3) < MIN_ATOM_DIST)
      return SAMP_ATOMS_CLOSE;

  // update atom coordinates
  for (int i = 0; i < 3; ++i)
    mol_array[1]->write_cm_pos()[i] = pos[i];
  mol_array[1]->update_atoms();

  mol_array[0]->init_dv();
  mol_array[0]->update_atoms();
 
  // statistical weight
  if (weight) {
  vector_product(pos, n_0, temp_pos);

  *weight = rrmass;
  for (int i = 0; i < 3; ++i)
    *weight += temp_pos[i] * temp_pos[i] / mol_array[0]->iner_mom(i);
  *weight = sqrt(*weight);

  *weight *= circle_rad[major] * circle_rad[major] / 4.; // 4 due to 1/(4 pi)
  }

  return 0;
    
}

int rand_vel (const Div_surf& ds, Dynvar& dynvar, 
	      double temper2, double* weight)
{

  dynvar.attach();

  if (ds.dist() <= 0.)
    return SAMP_FACE_OUT;

  *weight = 0.;

  // temporary variables
  Molecule* mp_temp;



  // initialization of different parameters
  static double m_factor [2]; // mass factors 
  static double iner_mom2 [2] [3]; // square roots of inertia moments
  static bool first_run = true;
  if (first_run)
    {
      if (mol_array.size() != 2)
	error("init_thermal: wrong number of molecules");

      first_run = false;
      m_factor [0] =   mol_array[1]->mass() / (mol_array[0]->mass() + mol_array[1]->mass());
      m_factor [1] = - mol_array[0]->mass() / (mol_array[0]->mass() + mol_array[1]->mass());

      for (int i = 0; i < 2; ++i)
	{
	  mp_temp = mol_array[i];
	  switch(mp_temp->type())
	    {
	    case ATOM:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = 0.0;
	      break;
	    case LINEAR:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = sqrt(mp_temp->iner_mom(j));
	      break;
	    case NONLINEAR:
	      for (int j = 0; j < 3; ++j)
		iner_mom2[i][j] = sqrt(mp_temp->iner_mom(j));
	      break;
	    }
	}
    }
  // effective mass
  static const double emass = m_factor[0] * mol_array[0]->mass();
  static const double emass2 = sqrt(emass);
      
  // randomly orient molecules
  // (cm positions are not needed for now)
  double lf_ref_pos [2] [3];  // reference points in laboratory frame
  for (int i = 0; i < 2; ++i)
    {
      mp_temp = mol_array[i];
      if (mp_temp->type() != ATOM)
	Random::orient(mp_temp->write_ang_pos(), mp_temp->ang_pos_size());
      for (int j = 0; j < 3; ++j)
	mp_temp->write_cm_pos()[j] = 0.0;
      mp_temp->update_atoms();
      // transform reference points to laboratory frame
      switch (mp_temp->type())
	{
	case ATOM:
	  for (int j = 0; j < 3; ++j)
	    lf_ref_pos [i] [j] = 0.0;
	  break;
	case LINEAR:
	  for (int j = 0; j < 3; ++j)
	    lf_ref_pos [i] [j] = mp_temp->read_ang_pos()[j] * 
	      ds.ref_pos(i)[0];
	  break;
	case NONLINEAR:
	  mp_temp->mf2lf(ds.ref_pos(i), lf_ref_pos [i]);
	  break;
	default:
	  error ("init_thermal: wrong type");
	}
    }

  // randomly orient a vector between reference points
  double lf_cos_12 [3];
  Random::orient(lf_cos_12, 3);

  // set cm to cm vector ...
  double pos_0 [3];
  for (int i = 0; i < 3; ++i)
    pos_0[i] = lf_cos_12[i] * ds.dist() - lf_ref_pos[0][i] + lf_ref_pos[1][i];

  // and update atoms coordinates
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 3; ++i)
      mol_array[j]->write_cm_pos()[i] = pos_0[i] * m_factor [j];
  for (int j = 0; j < 2; ++j)
    mol_array[j]->update_atoms();

  // check if the face is an external one
  for (int face = 0; face < ds.end(); ++face)
    {
      if ( (ds.face() - face) % ds.size() == 0)
	continue;

      // vector between to ref. points in LF
      double tmp_lf_pos_12 [3];
      for (int i = 0; i < 3; ++i)
	tmp_lf_pos_12[i] = pos_0[i];

      double tmp_lf_ref_pos [3];
      for (int frag = 0; frag < 2; ++frag)
	switch (mol_array[frag]->type())
	  {
	  case ATOM:
	    break;
	  case LINEAR:
	    for (int i = 0; i < 3; ++i)
	      if (frag)
		tmp_lf_pos_12 [i] -= mol_array[frag]->read_ang_pos()[i] * 
		  ds.ref_pos(frag, ds.ref_ind(frag, face))[0];
	      else
		tmp_lf_pos_12 [i] += mol_array[frag]->read_ang_pos()[i] * 
		  ds.ref_pos(frag, ds.ref_ind(frag, face))[0];
	    break;
	  case NONLINEAR:
	    mol_array[frag]->mf2lf(ds.ref_pos(frag, ds.ref_ind(frag, face)), 
				   tmp_lf_ref_pos);
	    for (int i = 0; i < 3; ++i)
	      if (frag)
		tmp_lf_pos_12 [i] -= tmp_lf_ref_pos[i];
	      else
		tmp_lf_pos_12 [i] += tmp_lf_ref_pos[i];
	    break;
	  default:
	    error ("init_thermal: wrong type");
	  }

      // the distance between two reference points is less than the minimum
      if (vlength(tmp_lf_pos_12, 3) < ds.dist(face))
	return SAMP_FACE_OUT;
    }

  // minimum interatomic distance calculation
  if(are_atoms_close(MIN_ATOM_DIST))
    return SAMP_ATOMS_CLOSE;

  // transform a vector between ref. points to molecular frame
  double mf_cos_12 [2] [3];
  for (int i = 0; i < 2; ++i)
    mol_array[i]->lf2mf(lf_cos_12, mf_cos_12 [i]);

  // reaction coordinate
  double rc_vec [9];
  double* rcp = rc_vec;
  // orbital coordinates
  for (int i = 0; i < 3; ++i)
    rcp[i] = lf_cos_12[i] / emass2;
  rcp += 3;
  // internal coordinates
  for (int i = 0; i < 2; ++i)
    {
      mp_temp = mol_array[i];
      switch (mp_temp->type())
	{
	case ATOM:
	  break;
	case LINEAR:
	  if (i == 0)
	    vector_product(lf_ref_pos [i], lf_cos_12, rcp);
	  else
	    vector_product(lf_cos_12, lf_ref_pos [i], rcp);
	  for (int j = 0; j < 3; ++j)
	    rcp[j] /= iner_mom2[i][j];
	  rcp += 3;
	  break;
	case NONLINEAR:
	  if (i == 0)
	    vector_product(ds.ref_pos(i), mf_cos_12 [i], rcp);
	  else
	    vector_product(mf_cos_12 [i], ds.ref_pos(i), rcp);
	  for (int j = 0; j < 3; ++j)
	    rcp[j] /= iner_mom2[i][j];
	  rcp += 3;
	  break;
	}
    }
  static const int rc_dim = rcp - rc_vec;

  // normalize reaction coordinate vector
  *weight = ds.dist() * ds.dist() * normalize (rc_vec, rc_dim);

  // generalized velocity
  double gv_vec [9];
  for (int i = 0; i < rc_dim; ++i)
    gv_vec[i] = Random::norm();
  // adjust reactive component
  orthogonalize (gv_vec, rc_vec, rc_dim);
  const double rcc = Random::exp();
  for (int i = 0; i < rc_dim; ++i)
    gv_vec [i] += rcc * rc_vec[i];

  // setting random velocities

  //    orbital
  const double* gvp = gv_vec;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 2; ++j)
      mol_array[j]->write_cm_vel()[i] = temper2 * m_factor[j] * gvp[i] 
	/ emass2;
  gvp += 3;
  //    internal
  for (int i = 0; i < 2; ++i)
    {
      mp_temp = mol_array[i];
      switch (mp_temp->type())
	{
	case ATOM:
	  break;
	case LINEAR:
	  for (int j = 0; j < 3; ++j)
	    mp_temp->write_ang_vel()[j] = temper2 * gvp[j] / iner_mom2[i][j];
	  mp_temp->orthogonalize();
	  gvp += 3;
	  break;
	case NONLINEAR:
	  for (int j = 0; j < 3; ++j)
	    mp_temp->write_ang_vel()[j] = temper2 * gvp[j] / iner_mom2[i][j];
	  gvp += 3;
	  break;
	}
    }

  return 0;
}

