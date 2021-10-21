#ifndef __ROTD_RANDOM__
#define __ROTD_RANDOM__

namespace Random {

  void   init      ();
  void   send_seed (int);
  double flat      ();             // rand_flat ();
  double norm      ();             // rand_norm ();
  double exp       ();             // rand_exp ();
  void   orient    (double*, int); // rand_orient (double*, int);

}


#endif
