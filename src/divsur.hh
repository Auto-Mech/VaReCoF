#ifndef __DIVSUR__
#define __DIVSUR__

#include <vector>
#include <string>

#include "surface.hh"
#include "expression.hh"

void surf_init (const char*, const char* =0);

extern vector<Div_surf*> ds_array;

#endif
