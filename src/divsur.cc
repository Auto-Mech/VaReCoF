#include "rotd.hh"
#include "divsur.hh"
#include "math.hh"

vector<Div_surf*> ds_array;

void surf_init (const char* in_name, const char* out_name)
{

  static const char funame [] = "surf_init: ";

  static const string ref_key = "Frame";
  static const string pnt_key = "PivotPoints";
  static const string dst_key = "Distances";
  static const string var_key = "Cycles";
  static const string cmp_key = "Conditions";
  static const string end_key = "EndSurface";

  using namespace Expression;

  ds_array.clear();

  string stemp;
  char next;
  int itemp;
  double dtemp;

  ifstream from(in_name);
  
  ofstream to;
  if(out_name)
    to.open(out_name);
  else
    to.open("/dev/null");
  to.precision(2);
  to.setf(ios::fixed, ios::floatfield);

  int frag = 0;

  // reference frame rotation
  //
  bool need_rotation[2];
  TMatrix<double, 3> rot_mat[2];

  // reference frame shift
  //
  bool need_shift[2];
  double ref_pos[2][3];

  vector<pair<string, Expr> > pnt_xpr[2]; // pivot points position of 2 frags

  vector<pair<string, Expr> > dst_xpr;    // interpoint distances
  vector<Cmp> cmp_xpr;                    // conditions
  vector<pair<Vec<Ident>,Vec<Vec<double> > > > var_spec; // variable arrays 


  string token;
  int surface_index = 0;
  skip_space(from);
  while(from >> token) {

    if(token == end_key) {// make surface
      //
      // print atomic coordinates in the dividing surface frame 
      //
      to << "Molecular geometry in the frame associated with the dividing surface(s) definition\n";

      for(int frag = 0; frag < 2; ++frag) {
	//
	to << "Fragment " << frag + 1 << ":\n";

	for(int a = 0; a < mol_array[frag]->size(); ++a) {
	  //
	  to << setw(2) << mol_array[frag]->atom(a).name();

	  // standard frame
	  //
	  if(!need_rotation[frag] && !need_shift[frag]) {
	    //
	    for(int i = 0; i < 3; ++i)
	      //
	      to << setw(13) <<  mol_array[frag]->atom(a).mf_pos[i];
	  }
	  // non-standard frame
	  //
	  else {
	    //
	    double tmp_pos[3];
	
	    for(int i = 0; i < 3; ++i)
	      //
	      tmp_pos[i] =  mol_array[frag]->atom(a).mf_pos[i];

	    if(need_shift[frag])
	      //
	      for(int i = 0; i < 3; ++i)
	        //
	        tmp_pos[i] -=  ref_pos[frag][i];
	
	    if(need_rotation[frag]) {
	      //
	      double pos[3];

	      vector_matrix_product(tmp_pos, rot_mat[frag], pos);

	      for(int i = 0; i < 3; ++i)
	        //
	        to << setw(13) <<  pos[i];
	    }
	    else {
	      //
	      for(int i = 0; i < 3; ++i)
	        //
	        to << setw(13) <<  tmp_pos[i];
	    }
	  }

	  to << "\n";
	}

	to << "\n";
      }
  for(int frag = 0; frag < 2; ++frag) {
    for(int at = 0; at < mol_array[frag]->size(); ++at) {
      cout << "\n";
    }

      }

      // check dimensions
      for(int i = 0; i < var_spec.size(); ++i)
	for(int j = 0; j < var_spec[i].second.size(); ++j)
	  if(var_spec[i].first.size() != var_spec[i].second[j].size()) {
	    cout << funame << "dimensions are inconsistent\n";
	    throw Form_Err();
	  }

      // check number of distances
      if(dst_xpr.size() * 9 != pnt_xpr[0].size() * pnt_xpr[1].size()) {
	cout << funame << "numbers of distances ("
	     << dst_xpr.size() << ") and pivot points (" 
	     << pnt_xpr[0].size() << ", " 
	     << pnt_xpr[1].size() << ") mismatch\n";
	throw Form_Err();
      }

      // check number of pivot points specs
      if(frag != 2) {
	cout << funame << "there should be EXACTLY TWO pivot point specs\n";
	throw Form_Err();
      }

      // check number of pivot points exprs
      if(pnt_xpr[0].size()%3 || pnt_xpr[1].size()%3) {
	cout << funame << "number of coordinates should be 3-proportional\n";
	throw Form_Err();
      }

      // surface specification
      int tot_num = 1;
      for(int i = 0; i < var_spec.size(); ++i)
	tot_num *= var_spec[i].second.size();

      vector<int> cyc_ind(var_spec.size());

      for(int tot_ind = 0; tot_ind < tot_num; ++tot_ind) {// tot_ind cycle

	// find indices of individual cycles
	itemp = tot_ind;
	for(int i = 0; i < var_spec.size(); ++i) {
	  cyc_ind[i] = itemp % var_spec[i].second.size();
	  itemp /= var_spec[i].second.size();
	}
	// set independent variables
	Varlist vars;
	for(int i = 0; i < var_spec.size(); ++i)
	  for(int j = 0; j < var_spec[i].first.size(); ++j)
	    vars.insert(var_spec[i].first[j](), 
			var_spec[i].second[cyc_ind[i]][j]);

	// set all variables;
	Varlist all_vars = vars;
	for(int j = 0; j < 2; ++j)
	  for(int i = 0; i < pnt_xpr[j].size(); ++i) {
	    stemp = pnt_xpr[j][i].first;
	    dtemp = pnt_xpr[j][i].second(vars);
	    all_vars.insert(stemp, dtemp);
	  }
	for(int i = 0; i < dst_xpr.size(); ++i)
	  all_vars.insert(dst_xpr[i].first, dst_xpr[i].second(vars));
	
	// check the conditions;
	bool is_cmp_valid = true;
	for(int i = 0; i < cmp_xpr.size(); ++i)
	  if(!cmp_xpr[i](all_vars)) {
	    is_cmp_valid = false;
	    break;
	  }    
	if(!is_cmp_valid) // do not make a surface
	  continue;
	
	// make surface
	Div_surf* dsp = new Div_surf(pnt_xpr[0].size()/3, pnt_xpr[1].size()/3);
	ds_array.push_back(dsp);

	to << "SURFACE INDEX = " << surface_index++ << "\n\n";
	for(frag = 0; frag < 2; ++frag) {
	  static const char* coor [3] = {"x", " y", " z"}; 
	  to << "pivot points for " << frag+1 <<"-th fragment:\n\n";
	  int pnt_num = pnt_xpr[frag].size()/3;
	  double tmp_pos [3];

	  for(int pnt_ind = 0; pnt_ind < pnt_num; ++pnt_ind) {
	    for(int spc_ind = 0; spc_ind < 3; ++spc_ind) {
	      stemp = pnt_xpr[frag][pnt_ind*3 + spc_ind].first;
	      tmp_pos[spc_ind] = all_vars[stemp];
	    }

	    matrix_vector_product(rot_mat[frag], tmp_pos, 
				  dsp->write_ref_pos(frag, pnt_ind));

	    for(int spc_ind = 0; spc_ind < 3; ++spc_ind) {
	      dsp->write_ref_pos(frag, pnt_ind)[spc_ind] += 
		ref_pos[frag][spc_ind];
	      to << coor[spc_ind] << " = " << setw(5) 
		 << dsp->ref_pos(frag, pnt_ind)[spc_ind];
	    }
	    to << "\n";
	    /*
	    if(mol_array[frag]->type() == LINEAR && 
	       dsp->ref_pos(frag, pnt_ind)[1] < -1.e-9) {
	      cout << funame << "torus radius is negative\n";
	      throw Form_Err();
	    }
	    */
	  }
	  to << "\n";
	}

	to << "distances between pivot points:\n";
	for(int i = 0; i < dst_xpr.size(); ++i) {
	  itemp = pnt_xpr[1].size() / 3;
	  int ri1 = i % itemp;
	  int ri0 = i / itemp;
	  stemp = dst_xpr[i].first;
	  dsp->write_dist(ri0, ri1) = all_vars[stemp];
	  if(!ri1) to << "\n";
	  to << "r[" << ri0+1 << "," << ri1+1 << "] = "
	     << setw(5) << all_vars[stemp]
	     << "   ";
	}
	to << "\n\nindependent variables:\n";
	for(int i = 0; i < vars.size(); ++i) {
	  if(!(i % 7))
	    to << "\n";
	  to << vars.name(i) << " = " 
	     << vars.value(i) << " ";
	} 
	to << "\n\n\n";

      }// tot_ind cycle
     

      // clean up
      frag = 0;
      for(int i = 0; i < 2; ++i)
	pnt_xpr[i].clear();
      dst_xpr.clear();
      cmp_xpr.clear();
      var_spec.clear();

      try {
	skip_space(from);
      } catch(EOF_Err) {}
	 
    }// make surface

    else if (token == pnt_key) { // pivot point specs

      if(frag > 1) {
	cout << funame << "there should be no more than 2 pivot point specs\n";
	throw Form_Err();
      }
      if(!(from >> itemp)) { // pivot points number
	cout << funame << "cannot read pivot points number\n";
	throw Form_Err();
      }
      int pnt_xpr_num = 3*itemp;
      
      skip_space(from);

      // reference frame
      from >> stemp;
      if(stemp != ref_key) {
	cout << funame << "keyword " << stemp 
	     << "does not match " <<ref_key << "\n";
	throw Form_Err();
      }

      // reference frame shift
      //
      from >> itemp; 
      if(itemp <= 0 || mol_array[frag]->type() == ATOM) {
	//
	need_shift[frag] = false;

	for(int i = 0; i < 3; ++i)
	  ref_pos[frag] [i] = 0.;
      }
      else if(itemp > mol_array[frag]->size()) {
	cout << funame << "reference atom "<< itemp << " does not exist\n";
	throw Form_Err();
      }
      else {
	//
	need_shift[frag] = true;

	for(int i = 0; i < 3; ++i)
	  ref_pos[frag] [i] = (mol_array[frag]->begin()+itemp-1)->mf_pos[i];
      }

      // orientation
      //
      int ref_at_num [3];


      if(mol_array[frag]->type() == NONLINEAR) {
	//
        need_rotation[frag] = true;

	for(int i = 0; i < 3; ++i) {
	  //
	  from >> ref_at_num[i];

	  if(ref_at_num[i] <= 0) {
	    //
	    need_rotation[frag] = false;

	    break;
	  }
	  else if(ref_at_num[i] > mol_array[frag]->size()) {
	    cout << funame << "reference atom " << ref_at_num[i] 
		 << " does not exist\n";
	    throw Form_Err();
	  }
	}
	
	getline(from, stemp);
      }
      else {
	//
	need_rotation[frag] = false;

	getline(from, stemp);
      }

      if(need_rotation[frag]) {

	for(int i = 0; i < 3; ++i) {
	  int i1 = (i+1) % 3;
	  int i2 = (i+2) % 3;
	  if(ref_at_num[i1] == ref_at_num[i2]){
	    cout << funame << "reference atoms should not be the same\n";
	    throw Form_Err();
	  }
	}
	
	double mat [3][3];
	for(int i = 0; i < 3; ++i)
	  for(int j = 0; j < 3; ++j)
	    mat[i][j] = (mol_array[frag]->begin()+ref_at_num[i]-1)->mf_pos[j];
      
	for(int i = 0; i < 2; ++i)
	  for(int j = 0; j < 3; ++j)
	    mat[i][j] -= mat[2][j];

	normalize(mat[0], 3);
	vector_product(mat[0], mat[1], mat[2]);
	normalize(mat[2], 3);
	vector_product(mat[2], mat[0], mat[1]);

	for(int i = 0; i < 3; ++i)
	  for(int j = 0; j < 3; ++j)
	    rot_mat[frag] (i, j) = mat[j][i];

      }
      else
	rot_mat[frag] = 1.;

      // pivot points specs
      for (int i = 0; i < pnt_xpr_num; ++i) {
	Ident name;
	from >> name;
	skip_space(from);
	from.get(next);
	if(next != '=') {
	  cout << funame << "invalid separator " << next 
	       << " (should be =)\n";
	  throw Form_Err();
	}
	pnt_xpr[frag].push_back(pair<string,Expr>(name(),Expr(from)));
      }
      
      // update pivot point index
      ++frag;
    }

    else if (token == dst_key) {// distances between pivot points
      if(frag < 2) {
	cout << funame << "distances section should go AFTER "
	  "pivot points specs for BOTH fragments\n";
	throw Form_Err();
      }
      int dst_num = pnt_xpr[0].size() * pnt_xpr[1].size() / 9;
      for(int i = 0; i < dst_num; ++i) {
	Ident name;
	from >> name;
	skip_space(from);
	from.get(next);
	if(next != '=') {
	  cout << funame << "unknown separator " << next << "\n";
	  throw Form_Err();
	}
	dst_xpr.push_back(pair<string,Expr>(name(), Expr(from)));
      }
    }

    else if (token == cmp_key) {// conditions
      int cmp_num;
      if(!(from >> cmp_num)) {
	cout << funame << "cannot read number of conditions\n";
	throw Form_Err();
      }
      for(int cmp_ind = 0; cmp_ind < cmp_num; ++cmp_ind)
	cmp_xpr.push_back(Cmp(from));
    }

    else if (token == var_key) {// variables
      int cyc_num;
      if(!(from >> cyc_num)) {
	cout << funame << "cannot read number of cycles\n";
	throw Form_Err();
      }
      var_spec.resize(cyc_num);
      for(int i = 0; i < cyc_num; ++i) {
	from >> var_spec[i].first;
	skip_space(from);
	from.get(next);
	if(next != '=') {
	  cout << funame << "unexpected separator " << next 
	       << " (should be =)\n";
	  throw Form_Err();
	}
	from >> var_spec[i].second;
      }
    }



    else {
      cout << funame << "unknown key-word " << token << "\n";
      throw Form_Err();
    }
  }
}
