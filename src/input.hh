#ifndef READ_HH
#define READ_HH

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "error.hh"

class ReadBase {
public:
    virtual void operator() (std::istream&)        =0;
    virtual void operator() (std::ostream&) const  =0;
    virtual ~ReadBase () {}
};

class ReadInt : public ReadBase {
  int& _data;
public:
  explicit ReadInt (int& v) : _data(v) {}
  void operator() (std::istream&);
  void operator() (std::ostream& to) const { to << _data; }
  ~ReadInt () {}
};

class ReadDouble : public ReadBase {
  double& _data;
public:
  explicit ReadDouble (double& v) : _data(v) {}
  void operator() (std::istream&);
  void operator() (std::ostream& to) const { to << _data; }
  ~ReadDouble () {}
};

class ReadString : public ReadBase {
  std::string& _data;
public:
  explicit ReadString (std::string& v) : _data(v) {}
  void operator() (std::istream&);
  void operator() (std::ostream& to) const { to << _data; }
  ~ReadString () {}
};

class ReadIarr : public ReadBase {
  std::vector<int>& _data;
public:
  explicit ReadIarr (std::vector<int>& v) : _data(v) {}
  void operator() (std::istream&);
  void operator() (std::ostream&) const;
  ~ReadIarr () {}
};

inline void ReadIarr::operator () (std::ostream& to) const
{
  for(int i = 0; i < _data.size(); ++i)
    to << _data[i] << " ";
}

class ReadDarr : public ReadBase {
  std::vector<double>& _data;
public:
  explicit ReadDarr (std::vector<double>& v) : _data(v) {}
  void operator() (std::istream&);
  void operator() (std::ostream&) const;
  ~ReadDarr () {}
};

inline void ReadDarr::operator () (std::ostream& to) const
{
  for(int i = 0; i < _data.size(); ++i)
    to << _data[i] << " ";
}

class ReadFun : public ReadBase {
public:
  typedef void (*fun_t) (std::istream&);
private:
    fun_t _data;
    bool  _used;
public:
  explicit ReadFun (fun_t v) : _data(v) , _used(false) {}
  void operator () (std::istream& from) { _data(from); _used == true; }    
  void operator () (std::ostream&) const;
  ~ReadFun () {}
};

inline void ReadFun::operator () (std::ostream& to) const
{
    if(_used)
	to << "has been processed";
    else
	to << "has not been processed";
}


class Read
{
  enum State {NOVAL, DEFAULT, READ};// reading status

  ReadBase* _read;
  State*    _state;
  int*      _count;

  void _delete_ref ();
  void _create_ref (const Read&);

public:
  
  void operator () (std::ostream&) const;
  void operator () (std::istream&);

  bool is_read    () const;
  bool is_init    () const;
  bool is_default () const;

  Read () : _read(0), _state(0),  _count(0) {}
  
  Read (const Read& r) { _create_ref(r); }
  Read& operator= (const Read& r) { _delete_ref(); _create_ref(r); return *this; }
  ~Read () { _delete_ref(); }

  explicit Read(int& var) 
    : _read(new ReadInt(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read (int& var, int val) 
    : _read(new ReadInt(var)), _state(new State(DEFAULT)), _count(new int(1)) { var = val; }
 
  explicit Read(double& var) 
    : _read(new ReadDouble(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read(double& var, double val) 
    : _read(new ReadDouble(var)), _state(new State(DEFAULT)), _count(new int(1)) { var = val; }
  
  explicit Read(std::string& var) 
    : _read(new ReadString(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read (std::string& var, const std::string& val) 
    : _read(new ReadString(var)), _state(new State(DEFAULT)), _count(new int(1)) { var = val; }

  explicit Read(std::vector<double>& var) 
    : _read(new ReadDarr(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read (std::vector<double>& var, const std::vector<double>& val) 
    : _read(new ReadDarr(var)), _state(new State(DEFAULT)), _count(new int(1)) { var = val; }
  
  explicit Read(std::vector<int>& var) 
    : _read(new ReadIarr(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read (std::vector<int>& var, const std::vector<int>& val) 
    : _read(new ReadIarr(var)), _state(new State(DEFAULT)), _count(new int(1)) { var = val; }
  
  explicit Read(ReadFun::fun_t var) 
    : _read(new ReadFun(var)), _state(new State(NOVAL)), _count(new int(1)) {} 
  Read  (ReadFun::fun_t var, int) 
    : _read(new ReadFun(var)), _state(new State(DEFAULT)), _count(new int(1)) {} 

};

inline void Read::_delete_ref ()
{
  if(!_count)
    return;
  
  if(*_count > 1) {
    --(*_count);
    return;
  }

  delete _read;
  delete _state;
  delete _count;
}
inline void Read::_create_ref (const Read& r)
{
  _read  = r._read;
  _state = r._state;
  _count = r._count;
  
  if(_count)
    ++(*_count);
}

inline void Read::operator () (std::istream& from)
{
  const char funame []  = "Read::operator () (std::istream&): ";

  if(!_read) {
    std::cout << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init_Err();
  }

  *_state = READ;
  (*_read)(from);
}

inline void Read::operator () (std::ostream& to) const
{
  const char funame []  = "Read::operator () (std::ostream&) const: ";

  if(!_read) {
    std::cout << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init_Err();
  }
  if(*_state == NOVAL)
    to << "NoN";
  else
    (*_read)(to);
}

inline std::ostream& operator<< (std::ostream& s, const Read& r) { r(s); return s; }
inline std::istream& operator>> (std::istream& s,       Read& r) { r(s); return s; }

inline bool Read::is_read () const
{
  const char funame []  = "Read::is_read() const: ";

  if(!_state) {
    std::cout << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init_Err();
  }


  if(*_state == READ)
    return true;
  else
    return false;
}

inline bool Read::is_init () const
{
  const char funame []  = "Read::is_init() const: ";

  if(!_state) {
    std::cout << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init_Err();
  }

  if(*_state == NOVAL)
    return false;
  else
    return true;
}

inline bool Read::is_default () const
{
  const char funame []  = "Read::is_default() const: ";

  if(!_state) {
    std::cout << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init_Err();
  }

  if(*_state == DEFAULT)
    return true;
  else
    return false;
}

class InputKey : public std::string {

  static std::set<const std::string*> pool;

public:
  InputKey () { pool.insert(this); }
  InputKey (const std::string& s) : std::string(s) { pool.insert(this); }
  InputKey (const char* s) : std::string(s) { pool.insert(this); }
  ~InputKey () {pool.erase(this); }
  static void show_all ();
};

#endif
