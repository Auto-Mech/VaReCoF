#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

int main (int argc, char* argv []) {

  const char   en_key [] = "EnergyGrid";
  const char   am_key [] = "AngularMomentumGrid";
  const char  inp_key [] = "InputFluxFile";
  const char  out_key [] = "OutputFluxFile";
  const char freq_key [] = "Frequencies";

  double dtemp;
  int itemp;

  if (argc != 2) {
    std::cout << "usage: convolute inp_file\n";
    exit(1);
  }

  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << "cannot open " << argv[1] << " file\n";
    exit(1);
  }

  std::vector<int> freq_array;// vibrational frequencies, cm^-1
  double emin, estep; int esize = 0;// energy grid, cm^-1
  double amin, astep; int asize = 0;// angular momentum grid, a.u.
  std::string out_flux, inp_flux;

  // input
  std::string token;
  while(from >> token) {
    if(token == en_key) {
      from >> emin >> esize >> estep;
      if(!from) {
	std::cout << "cannot read " << token << " data\n";
	exit(1);
      }
    }
    else if(token == am_key) {
      from >> amin >> asize >> astep;
      if(!from) {
	std::cout << "cannot read " << token << " data\n";
	exit(1);
      }
    }
    else if(token == out_key) {
      from >> out_flux;
      if(!from) {
	std::cout << "cannot read " << token << " data\n";
	exit(1);
      }
    }
    else if(token == inp_key) {
      from >> inp_flux;
      if(!from) {
	std::cout << "cannot read " << token << " data\n";
	exit(1);
      }
    }
    else if(token == freq_key) {
      int freq_size;
      from >> freq_size;
      for(int i = 0; i < freq_size; ++i) {
        from >> dtemp;
	itemp = lrint(dtemp);
        if(itemp <= 0) {
          std::cout << "negative or zero " << i+1 << "-th frequency: "
                    << itemp << "\n";
          exit(1);
        }
        freq_array.push_back(itemp);
      }
      if(!from || !freq_array.size()) {
        std::cout << "cannot read " << freq_key << " data\n";
        exit(1);
      }
      //std::cout << freq_array.size() << " frequencies were read\n";
    }
    else {
      std::cout << "unknown key: " << token << "\n";
      exit(1);
    }
  }
  from.close();
  from.clear();

  if(!esize) {
    std::cout << en_key << " is not initialized\n";
    exit(1);
  }
  if(!asize) {
    std::cout << am_key << " is not initialized\n";
    exit(1);
  }
  if(!out_flux.size()) {
    std::cout << out_key << " is not initialized\n";
    exit(1);
  }
  if(!inp_flux.size()) {
    std::cout << inp_key << " is not initialized\n";
    exit(1);
  }
  if(!freq_array.size()) {
    std::cout << freq_key << " is not initialized\n";
    exit(1);
  }

  std::cout << "parameters reading done\n";

  
  // vibrational density of states
  if(emin > 0.)
    itemp = int(emin + estep * double(esize - 1));
  else
    itemp = int(estep * double(esize));
  std::vector<double> vib_dos(itemp, 0.);
  vib_dos[0] = 1.;
  for(int i = 0; i < freq_array.size(); ++i)
    for(int en = freq_array[i]; en < vib_dos.size(); ++en)
      vib_dos[en] += vib_dos[en - freq_array[i]];

  from.open(inp_flux.c_str());
  if(!from) {
    std::cout << "cannot open input flux file " << inp_flux << "\n";
    exit(1);
  }
  std::ofstream to(out_flux.c_str());
  std::vector<double> tm_nos(esize);
  double ener, amom, nos;
  double tse, pos;
  int ipos, ve_max;
  for(int am = 0; am < asize; ++am) {
    for(int en = 0; en < esize; ++en)
      from >> ener >> amom >> tm_nos[en];
    if(!from) {
      std::cout << "input flux file " << inp_flux << " is corrupted\n";
      exit(1);
    }

    ener = emin;
    for(int en = 0; en < esize; ++en) {
      if(emin > 0.)
	ve_max = int(ener);
      else
	ve_max = int(ener - emin + estep);
      ve_max = ve_max > vib_dos.size() ? vib_dos.size() : ve_max;

      nos = tm_nos[en];
      for(int ve = 1; ve < ve_max; ++ve) {
	tse = ener - double(ve);
	if(tse >= emin) {
	  pos = (tse-emin)/estep;
	  ipos = (int)floor(pos);
	  dtemp = tm_nos[ipos] + (tm_nos[ipos + 1] - tm_nos[ipos]) * 
	    (pos - double(ipos));
	  nos += vib_dos[ve] * dtemp;
	}
	else if(emin > 0.)
	  nos += vib_dos[ve] * tm_nos[0] * tse/emin;
	else
	  nos += vib_dos[ve] * tm_nos[0] * (estep + tse - emin)/estep;
      }
      to << ener << " " << amom << " " << nos << "\n";
      ener += estep;
    }    
  }
  return 0;
}


