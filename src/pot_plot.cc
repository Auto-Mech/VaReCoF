#include <cmath>
#include <fstream>
#include <iomanip>

extern "C" void spl_c2h6__(const double&, const double&, 
			   const double&, double&);

int main (int argc, char* argv[]) {

  const double r2g = 180./M_PI;
  const double au2kcal = 627.503;

  double xmin, xmax, xstep, ymin, ymax, ystep, z;

  if (argc != 2) {
    cout << "main: usage: prog input_file\n";
    return 1;
  }
  ifstream from(argv[1]);
  if(!from) {
    cout << "main: input file " << argv[1] << " is not found\n";
    return 1;
  }
  from >> xmin >> xmax >> xstep
       >> ymin >> ymax >> ystep
       >> z;
  from.close();

  double r, theta, phi;
  double vtot;

  for(double x = xmin; x <= xmax; x+= xstep) {
    if(x == 0.)
      continue;
    for(double y = ymin; y <= ymax; y += ystep) {
      r = sqrt(x*x+y*y+z*z);
      if(r == 0.) 
	continue;
      theta = r2g*acos(z/r);
      phi = r2g*atan(y/x);
      if(x < 0.) 
	phi += 180.;
      spl_c2h6__(r,theta,phi,vtot);
      vtot *= au2kcal;

      cout << setw(13) << x 
	   << setw(13) << y 
	   << setw(13) << vtot << "\n";
    }
    cout << "\n";
  }
  return 0;
}
