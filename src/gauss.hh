#ifndef MY_GAUSS
#define MY_GAUSS

#include <string>
#include <vector>
#include <fstream>
#include <set>
#include <map>
#include <list>
#include "math.hh"

namespace Gauss {

  //using namespace std;

  // Error base class
  //
  class Err {};

  // File handling error
  //
  class File_Err: public Err {};

  // File opening error
  //
  class Open_Err: public File_Err{};

  // Input error
  //
  class Form_Err: public File_Err{};

  // Math errors
  //
  class Math_Err: public Err {};
  
  // Indexing errors
  //
  class Range_Err: public Err {};
  
  // Search errors
  //
  class Find_Err: public Err {};

  // Logic errors
  //
  class Logic_Err: public Err {};
  
  // Initialization errors
  //
  class Init_Err: public Err {};
  
  // Run errors
  //
  class Run_Err: public Err {}; 

  extern int debug;

  // reactive system specification

  struct Atom {
    int num;                 // atomic number
    double pos [3];          // atom's position
    double force [3];        // force exerted on the atom
  };

  std::ostream& operator<< (std::ostream&, const Atom&);

  extern std::vector<Atom> cluster;

  // operational flags
  //
  enum {
    FORCE  = 4, // calculate forces
    READ   = 2, // guess = read
    DIRECT = 1  // apply method literally (do not use failsafe methods); 
  };

  class Method;

  // method id
  //
  class MetId: public std::string {
    //
  public:
    //
    MetId () {}

    template<typename T>
    MetId (const T& t) : std::string(t) {}

    Method* operator-> () const;
  };

  inline std::istream& operator>> (std::istream& from, MetId& mid) { return from >> (std::string&)mid; }

  inline std::ostream& operator<< (std::ostream& to, const MetId& mid) { return to << (const std::string&)mid; }

  typedef std::list<MetId> mid_t;
  
  // search pattern base
  //
  class Pattern: private std::list<std::string> {
    //
  public:
    //
    Pattern (std::istream&);

    bool scan (std::istream&) const;

    void print (std::ostream&) const;
  };

  inline std::ostream& operator<< (std::ostream& to, const Pattern& patt) { patt.print(to); return to; }

  // energy search pattern with the corresponding weight
  //
  class EnergyPattern: public Pattern {
    //
    double _weight;

  public:
    //
    EnergyPattern (std::istream&);

    double weight () const { return _weight; }
  };

  // search pattern with the message
  //
  class MessagePattern: public Pattern {

    std::string _message;

  public:
    //
    MessagePattern (std::istream&);
    
    const std::string& message () const { return _message; }

  };

  // failure specific search pattern with the corresponding list of failsafe methods
  //
  class FailurePattern: public MessagePattern, public mid_t {

  public:
    //
    FailurePattern (std::istream&);
  };

  // calculation method
  //
  class Method {
    //
    std::string _id; // method identification string
    
    std::string _template; // gaussian input template
    
    int    _geom_pos; // geometry template position
    
    // methods to calculate the initial guess
    //
    mid_t guess_methods;

    // is initial guess is a minimum of all init guesses or just the first that works
    //
    bool guess_min;

    // is initial guess required or one can start without checkpoint file
    //
    bool guess_required;


    // energy search patterns
    //
    typedef std::list<EnergyPattern> epat_t;

    epat_t ener_patt;

    // failure patterns and corresponding failsafe methods
    //
    typedef std::list<FailurePattern> fpat_t;

    fpat_t fail_patt;

    // rerun patterns
    //
    typedef std::list<MessagePattern> rpat_t;
    
    rpat_t run_patt;

    // non-specific failsafe methods
    //
    mid_t default_failsafe_methods;
    
    // reference energy
    //
    double _ref_ener;

    // calculated energy
    //
    double _abs_ener;
    
    // saved energy
    //
    double _save_ener;

    // method application status
    //
    int _status;

    void _check_status () const; 
    
    // failsafe method
    //
    Method* _failsafe;

    // all available methods
    //
    typedef std::map<std::string, Method*> _mpool_t;
 
    static _mpool_t _method_pool;
    
    // check for infinite loops
    //
    void _check_for_loops (std::set<std::string>&) const; 

    // check status
    //
    mutable bool _checked;

  public:
    //
    // status flags
    //
    enum {
      INIT,    // method was not applied yet
      SUCCESS, // calculation succeeded
      FAIL  // calculation failed
    };

    Method (std::istream&);
    
    Method* apply (int =0);
    
    const std::string& id () const { return _id; }
    
    void read_forces (std::vector<D3>&) const;
    
    void save_log (int) const; // save gauss com, log, etc. files
    
    void print_info (std::ostream&) const; // print calculated energies, etc.

    void backup  ();
    
    void restore ();
    
    int status () const { return _status; }
    
    void init () { _status = INIT; }

    double rel_energy  () const {  _check_status(); return _abs_ener - _ref_ener; }
    
    double abs_energy  () const {  _check_status(); return _abs_ener;             }

    static Method* find_method (const std::string&);

    static void print_available_methods (std::ostream&);

    static void init_all ();

    static int pool_size() { return _method_pool.size(); }

    // check all for infinite loops
    //
    static void check_all ();

    // print geometry, calculated energies, etc. from all methods
    //
    static void print_all_info (std::ostream&); 
  };

  // linear interpolation of several methods
  //
  class Interpolation: private std::list<std::pair<MetId, double> > {
    //
  public:
    //
    Interpolation () {}
    
    void init (std::ifstream&);
    
    bool isinit () const { return size(); }
    
    double apply (int) const;
    
    void   check ()    const;
  };

  // initialize Gaussian from the file
  //
  void init (const std::string&);

  // minimization potential
  //
  Method* min_pot (int =0);

  // interpolation potential
  //
  double  inter_pot (int =0);

  // correction to the potential
  //
  double pcorr(const int&, const double [], double []);
  
  bool  ispcorr ();
  
  void print_geom ();

  // optimization
  //
  extern MetId opt_method;
}

#endif
