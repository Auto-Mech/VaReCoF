#include "surface.hh"
#include "error.hh"
//#include "random.hh"
#include "tmatrix.hh"
#include "math.hh"

#include <iostream>

/********************** Dividing surface **********************/

void Div_surf::destroy ()
{
  for(int i = 0; i < 2; ++i)
    delete[] _ref_pos[i];

  delete[] rr_dist;
}

void Div_surf::create (int n0, int n1)
{
  static const char funame [] = "Div_surf::create: ";

  if(n0 < 1 || n1 < 1)  {
    cout << funame << "wrong number of pivot points\n";
    throw File_input_error();
  }
  _ref_num[0] = n0;
  _ref_num[1] = n1;
  
  for(int i = 0; i < 2; ++i)
    _ref_pos[i] = new double[3 * _ref_num[i]];

  rr_num = n0 * n1;
  rr_dist = new double[rr_num]; 
}

Div_surf::Div_surf (const Div_surf& ds)
  :_face(ds.face())
{
  create(ds.ref_num(0), ds.ref_num(1));

  for (int frag = 0; frag < 2; ++frag)
      for (int ri = 0; ri < ds.ref_num(frag); ++ri)
	for (int i = 0; i < 3; ++i)
	  write_ref_pos(frag, ri)[i] = ds.ref_pos(frag, ri)[i]; 
  
  for (int f = 0; f < rr_num; ++f)
    rr_dist[f] = ds.dist(f);
}

Div_surf& Div_surf::operator= (const Div_surf& ds)
{
  _face = ds.face();

  destroy();
  create(ds.ref_num(0), ds.ref_num(1));

  for (int frag = 0; frag < 2; ++frag)
    for (int ri = 0; ri < ds.ref_num(frag); ++ri)
      for (int i = 0; i < 3; ++i)
	write_ref_pos(frag, ri)[i] = ds.ref_pos(frag, ri)[i]; 

  for (int f = 0; f < rr_num; ++f)
    rr_dist[f] = ds.dist(f);
  
  return *this;
}

Div_surf::Div_surf(istream& from)
  :_face(0)
{
  static const char funame [] = "Div_surf::Div_surf (istream&): ";

  int n0, n1;
  from >> n0 >> n1;
  if(!from) {
    cout << funame << "input stream is corrupted at number of points\n";
    throw File_input_error();
  }
  create(n0, n1);

  for(int frag = 0; frag < 2; ++frag)
    for(int rp_ind = 0; rp_ind < ref_num(frag); ++rp_ind)
      for(int i = 0; i < 3; ++i)
	from >> write_ref_pos(frag, rp_ind)[i];

  for (int face = 0; face < size(); ++face)
    from >> write_dist(face);
  
  if(!from) {
    cout << funame << "input stream is corrupted at data\n";
    destroy();
    throw File_input_error();
  }
}

Div_surf::Div_surf(int n0, int n1)
  :_face(0)
{
  create(n0, n1);
}

Div_surf::~Div_surf()
{
  destroy();
}

int Div_surf::ref_ind (int frag, int f) const
{
  f = f % rr_num;
  if (f < 0)
    f += rr_num;

  switch(frag)
    {
    case 0:
      return f % _ref_num[0];
    case 1:
      return f / _ref_num[0];
    default:
      error("Div_surf::ref_ind: wrong fragment index");
    }
  error("Div_surf::ref_ind: should not reach this point");
  return 0;
}
  
int Div_surf::ref_ind (int frag) const
{
  int f = _face % rr_num;
  if (f < 0)
    f += rr_num;

  switch(frag)
    {
    case 0:
      return f % _ref_num[0];
    case 1:
      return f / _ref_num[0];
    default:
      error("Div_surf::ref_ind: wrong fragment index");
    }
  error("Div_surf::ref_ind: should not reach this point");
  return 0;
}

double Div_surf::dist () const
{ 
  int f = _face % rr_num;
  if (f < 0)
    f += rr_num;

  return rr_dist[f]; 
}

double& Div_surf::write_dist () 
{ 
  int f = _face % rr_num;
  if (f < 0)
    f += rr_num;

  return rr_dist[f];
}


double Div_surf::dist (int f) const 
{ 
  f = f % rr_num;
  if (f < 0)
    f += rr_num;

  return rr_dist[f]; 
}

double& Div_surf::write_dist (int f) 
{ 
  f = f % rr_num;
  if (f < 0)
    f += rr_num;

  return rr_dist[f]; 
}

double Div_surf::dist (int i0, int i1) const
{
  i0 %= _ref_num[0];
  if (i0 < 0)
    i0 += _ref_num[0];

  i1 %= _ref_num[1];
  if (i1 < 0)
    i1 += _ref_num[1];

  return rr_dist[i0 + i1 * _ref_num[0]];
}

double& Div_surf::write_dist (int i0, int i1)
{
  i0 %= _ref_num[0];
  if (i0 < 0)
    i0 += _ref_num[0];

  i1 %= _ref_num[1];
  if (i1 < 0)
    i1 += _ref_num[1];

  return rr_dist[i0 + i1 * _ref_num[0]];
}

bool Div_surf::operator== (const Div_surf& ds) const
{

  const double tol = 1.e-7;
  double dtemp;

  if (ref_num(0) != ds.ref_num(0) || ref_num(1) != ds.ref_num(1))
    return false;

  for (int frag = 0; frag < 2; ++frag)
    for (int ref = 0; ref < ref_num(frag); ++ref)
      for (int coor = 0; coor < 3; ++coor)
	{
	  dtemp = ref_pos(frag, ref)[coor] - ds.ref_pos(frag, ref)[coor];
	  if (dtemp > tol || dtemp < -tol)
	    return false;
	}
  
  for (int face = 0; face < size(); ++face)
    {
      dtemp = dist(face) - ds.dist(face);
      if (dtemp > tol || dtemp < -tol)
	return false;
    }
  
  return true;

}

ostream& operator<< (ostream& to, const Div_surf& ds)
{
  to << ds.ref_num(0) << "   " << ds.ref_num(1) << "\n";
  for (int frag = 0; frag < 2; ++frag)
    for (int rp_ind = 0; rp_ind < ds.ref_num(frag); ++rp_ind)
      {
	for (int i = 0; i < 3; ++i)
	  to << ds.ref_pos(frag, rp_ind)[i] << "  ";
	to << "\n";
      }

  for (int face = 0; face < ds.size(); ++face)
    to << ds.dist(face) << "  ";
  to << endl;
  return to;
}

istream& operator>> (istream& from, Div_surf& ds)
{
  Div_surf ds_temp(from);
  ds = ds_temp;
  return from;
}

/*
int Div_surf::sample (int smp_num) const
{
  if (size() == 1)
    return 0;

  int skip_num = 0;

  double pos [3];   // vector between two chosen pivot points
  double quat [4];  // quaternion

  static TMatrix<double,3> rot_mat_0; // rotational matrix of the first fragment
  static TMatrix<double,3> rot_mat_1; // rotational matrix of the second fragment

  double mf_pos_0 [3]; // position of the first pivot point  (mol. frame)
  double mf_pos_1 [3]; // position of the second pivot point (mol. frame)
  double lf_pos_0 [3]; // position of the first pivot point  (lab. frame)
  double lf_pos_1 [3]; // position of the second pivot point (lab. frame)
  double pos_01   [3]; // vector between two pivot points


  for (int smp_ind = 0; smp_ind < smp_num; ++smp_ind) {

    Random::orient(pos, 3);
    for (int i = 0; i < 3; ++i)
      pos[i] *= dist();

    Random::orient(quat, 4);
    quat2mat (quat, rot_mat_0); 
    
    Random::orient(quat, 4);
    quat2mat (quat, rot_mat_1);

    bool is_skip = false;
    for (int ff = 0; ff < size(); ++ff) {

      if ((ff - _face)/size() == 0)                        
	continue;

      int i0 = ref_ind(0, ff);
      int i1 = ref_ind(1, ff);

      for (int i = 0; i < 3; ++i)
	pos_01[i] = pos [i]; 

      if (i0 != ref_ind(0)) {

	for (int i = 0; i < 3; ++i)
	  mf_pos_0[i] = ref_pos(0)[i] - ref_pos(0, i0)[i];

	vector_matrix_product (mf_pos_0, rot_mat_0, lf_pos_0);

	for (int i = 0; i < 3; ++i)
	  pos_01[i] += lf_pos_0[i]; 
      }

      if (i1 != ref_ind(1)) {

	for (int i = 0; i < 3; ++i)
	  mf_pos_1[i] = ref_pos(1)[i] - ref_pos(1, i1)[i];

	vector_matrix_product (mf_pos_1, rot_mat_1, lf_pos_1);

	for (int i = 0; i < 3; ++i)
	  pos_01[i] -= lf_pos_1[i]; 
      }

      if (vlength(pos_01, 3) < dist(ff)) {
	is_skip = true;
	break;
      }

    }

    if (is_skip)
      ++skip_num; 
  }  
  return skip_num;
}
*/
