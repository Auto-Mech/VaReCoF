#ifndef __RANINIT__
#define __RANINIT__

#include "rotd.hh"

extern double MIN_ATOM_DIST;

enum {SAMP_OUT = 1, SAMP_AUX = 2}; //output flags
enum {SAMP_ATOMS_CLOSE = 1, SAMP_FACE_OUT}; // result flags

extern int (*rand_pos) (const Div_surf&, Dynvar&, double*);

void name2ranp (const string&);

int rand_vel (const Div_surf&, Dynvar&, double,  double*);

#endif
