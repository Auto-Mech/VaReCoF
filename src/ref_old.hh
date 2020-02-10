#include <map>
#include <iostream>

/*****************************************************************************
 **************************     Reference   **********************************
 *****************************************************************************/

class Ref_error {};

template<class S>
class Ref
{
  S* sp; // object pointer
  int id_;   // object id

  static map<int, int> ref_num; // number of references
  static int find_id ();        // find available id
 
public:
  
  Ref () : sp(0),  id_(0) {}

  template<class T> 
  Ref (const T&);
  Ref (const Ref&);
  ~Ref ();

  void init () throw (Ref_error);

  template<class T> 
  void init (const T&) throw (Ref_error);

  int is_init () const;

  S& operator() () throw (Ref_error);
  S* operator&  () throw (Ref_error);

  Ref& operator= (const Ref&);
  template<class T> 
  Ref& operator= (const T&);

  S& operator+= (const Ref&) throw (Ref_error);
  S& operator-= (const Ref&) throw (Ref_error);

  template<class T> 
  S& operator+= (const T&) throw (Ref_error);
  template<class T> 
  S& operator-= (const T&) throw (Ref_error);

  void operator++ ()  throw (Ref_error);
  void operator++ (int) throw (Ref_error);

  void operator-- ()  throw (Ref_error);
  void operator-- (int) throw (Ref_error);

  friend istream& operator>> <S> (istream&, Ref<S>&) throw (Ref_error);
  friend ostream& operator<< <S> (ostream&, const Ref<S>&) throw (Ref_error);
  
};

template<class S>
S& Ref<S>::operator+= (const Ref& r) throw (Ref_error) 
{ 
  if (!id_)
    {
      cout << "Ref::operator+=: not initialized" << endl;
      throw Ref_error(); 
    }
  *sp += *r.sp;
  return *sp;
} 

template<class S>
S& Ref<S>::operator-= (const Ref& r) throw (Ref_error) 
{
  if (!id_)
    {
      cout << "Ref::operator-=: not initialized" << endl;
      throw Ref_error(); 
    }
  *sp -= *r.sp;
  return *sp;
} 

template<class S>
template<class T> 
S& Ref<S>::operator+= (const T& t) throw (Ref_error) 
{ 
  if (!id_)
    {
      cout << "Ref::operator+=: not initialized" << endl;
      throw Ref_error(); 
    }
  *sp += t;
  return *sp;
} 

template<class S>
template<class T> 
S& Ref<S>::operator-= (const T& t) throw (Ref_error) 
{
  if (!id_)
    {
      cout << "Ref::operator-=: not initialized" << endl;
      throw Ref_error(); 
    }
  *sp -= t;
  return *sp;
} 

template<class S>
void Ref<S>::operator++ ()  throw (Ref_error)   
{
  if (!id_)
    {
      cout << "Ref::operator++  not initialized" << endl;
      throw Ref_error(); 
    }
  ++(*sp); 
}

template<class S>
void Ref<S>::operator++ (int) throw (Ref_error) 
{
  if (!id_)
    {
      cout << "Ref::operator++: not initialized" << endl;
      throw Ref_error(); 
    }
  (*sp)++;
}

template<class S>
void Ref<S>::operator-- ()  throw (Ref_error)
{
  if (!id_)
    {
      cout << "Ref::operator--: not initialized" << endl;
      throw Ref_error(); 
    }
  --(*sp);
}

template<class S>
void Ref<S>::operator-- (int) throw (Ref_error)
{
  if (!id_)
    {
      cout << "Ref::operator--: not initialized" << endl;
      throw Ref_error(); 
    }
  (*sp)--;
}

template<class S>
S& Ref<S>::operator() () throw (Ref_error)
{ 
  if (!id_)
    {
      cout << "Ref::operator(): not initialized" << endl;
      throw Ref_error(); 
    }
  return *sp; 
}

template<class S>
S* Ref<S>::operator&  () throw (Ref_error) 
{ 
  if (!id_)
    {
      cout << "Ref::operator&: not initialized" << endl;
      throw Ref_error(); 
    }
  return  sp; 
}

template<class S>
istream& operator>> (istream& from, Ref<S>& r) throw (Ref_error)
{
  if (!r.id_)
    {
      cout << "Ref::operator>>: not initialized" << endl;
      throw Ref_error(); 
    }
  
  if (!r.sp)
    error("Ref::operator>>: reference database is corrupted");
      
  from >> *r.sp;

  return from;
}
  
template<class S>
ostream& operator<< (ostream& to, const Ref<S>& r) throw (Ref_error)
{
  if (!r.id_)
    {
      cout << "Ref::operator<<: not initialized" << endl;
      throw Ref_error(); 
    }
  
  if (!r.sp)
    error("Ref::operator<<: reference database is corrupted");
      
  to << *r.sp; 
  return to;
}
  
template<class S>
int Ref<S>::find_id () // find next available id
{
  int res = 0;
  while (1)
    if (ref_num.find(++res) == ref_num.end())       
      return res;
}

template<class S>
template<class T>
Ref<S>::Ref (const T& t)
  : sp(new S(t))
{
  id_ = find_id(); 
  ref_num[id_] = 1;
}

template<class S>
void Ref<S>::init () throw (Ref_error)
{
  if(id_)
    {
      cout << "Ref::init: reference is already initialized" << endl;
      throw Ref_error();
    }

  id_ = find_id();
  ref_num[id_] = 1;
  sp = new S();
}

template<class S>
template<class T>
void Ref<S>::init (const T& t) throw (Ref_error)
{
  if(id_)
    {
      cout << "Ref::init: reference is already initialized" << endl;
      throw Ref_error();
    }

  id_ = find_id();
  ref_num[id_] = 1;
  sp = new S(t);
}

template<class S>
int Ref<S>::is_init () const
{
  return id_;
}

template<class S>
Ref<S>::Ref (const Ref& r)
  : id_(r.id_), sp(r.sp)
{  
  if (!id_)
    return;
  
  if (ref_num.find(id_) == ref_num.end())
    error("Ref(const Ref&): object id is not found");

  ++ref_num[id_];
}

template<class S>
Ref<S>& Ref<S>::operator= (const Ref& r)
{  
  if (id_)
    {
      if (ref_num.find(id_) == ref_num.end())
	error("operator=: object 1 id is not found");

      switch (ref_num[id_])
	{
	case 0:
	  error("operator=: wrong number of references");
	case 1:
	  delete sp;
	  ref_num.erase(id_);
	default:
	  --ref_num[id_];
	}
    }

  id_ = r.id_; sp = r.sp;

  if (!id_)
    return *this;

  if (ref_num.find(id_) == ref_num.end())
    error("operator=: object 2 id is not found");

  ++ref_num[id_];

  return *this;
}

template<class S>
template<class T>
Ref<S>& Ref<S>::operator= (const T& t)
{  
  if (id_)
    {
      if (ref_num.find(id_) == ref_num.end())
	error("operator=: object 1 id is not found");

      switch (ref_num[id_])
	{
	case 0:
	  error("operator=: wrong number of references");
	case 1:
	  delete sp;
	  ref_num.erase(id_);
	default:
	  --ref_num[id_];
	}
    }

  sp = new S(t);
  id_ = find_id(); 
  ref_num[id_] = 1;

  return *this;
}

template<class S>
Ref<S>::~Ref ()
{
  if (!id_)
    return;

  if (ref_num.find(id_) == ref_num.end())
    error("~Ref: object id is not found");

  switch (ref_num[id_])
    {
    case 0:
      error("~Ref: wrong number of references");
    case 1:
      delete sp;
      ref_num.erase(id_);
    default:
      --ref_num[id_];
    }
}

