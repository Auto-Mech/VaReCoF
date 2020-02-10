#ifndef ROTD_FORCE_HH
#define ROTD_FORCE_HH

#include <string>
#include "pes.hh"

enum /* potential energy function flags */
{
  POT_FRC = 1,   /* evaluate forces */
  POT_CHK = 2,   /* check file is available */
  POT_CNV = 4    /* use slow but well converging method */
};

enum /* potential result status */
{
  G98_EGEN = 1,     /* g98 generic failure */
  G98_ECNV          /* g98 convergence failure */
};

enum // potential type
{
  G98_TYPE = 1,  // g98
  ANAL_TYPE, // analytic
  MODEL_TYPE // model potential
};

struct Pot_error
{
  int stat;
  Pot_error(int i) { stat = i; }
  Pot_error() { stat = G98_EGEN; }
};

// finds potential function by name
pot_f name2pot (const std::string&);

// geometry relaxation
enum {
  OPT_READ = 1 // read electronic stracture from the previous calculation
};

void set_opt(const std::string& type, const std::string& method);

#endif
