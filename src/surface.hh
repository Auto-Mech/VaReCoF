#ifndef __SURFACE__
#define __SURFACE__

#include <fstream>

using namespace std;

/********************** Dividing surface **********************/

class Div_surf {

public:

  class File_input_error {};

private:

  int _ref_num[2];     // number of the ref. points of 0-th molecule
  double* _ref_pos[2]; // reference points of the i-th molecule
  int rr_num;          // number of reference point pairs
  double* rr_dist;     // distance between ref. points

  int _face;

  void create (int, int);
  void destroy ();

public:
  

  Div_surf(istream&);
  Div_surf(int =1, int =1);

  ~Div_surf();
  Div_surf (const Div_surf&);
  Div_surf& operator= (const Div_surf&);
 
  int size () const { return rr_num; }
  int end () const { return rr_num; }
  void begin () { _face = 0; }

  int face () const { return _face; }
  void set_face (int i) {  _face = i; }

  int operator++ () { return ++_face; }  
  int operator++ (int) { return _face++; }

  int ref_ind (int frag) const;
  int ref_ind (int frag, int f) const;

  int ref_num (int frag) const { return _ref_num[frag]; }

  const double* ref_pos (int frag) const 
  { return _ref_pos[frag] + 3 * ref_ind(frag); }
  const double* ref_pos (int frag, int rp_ind) const 
  { return _ref_pos[frag] + 3 * rp_ind; }
  double* write_ref_pos (int frag, int rp_ind)
  { return _ref_pos[frag] + 3 * rp_ind; }

  double dist () const;
  double dist (int) const;
  double dist (int, int) const;
  double& write_dist ();
  double& write_dist (int);
  double& write_dist (int, int);

  bool operator== (const Div_surf&) const;
  bool operator!= (const Div_surf& ds) const
  { return !operator==(ds); }

//  int sample (int smp_num) const;
  friend istream& operator>> (istream&, Div_surf&);

};

ostream& operator<< (ostream&, const Div_surf&);
istream& operator>> (istream&, Div_surf&);

#endif
