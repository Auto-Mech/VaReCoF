#include <iostream>

/*****************************************************************************
 **************************     Reference   **********************************
 *****************************************************************************/

class Ref_error {};

template<class S> class Ref;
template<class S>
std::istream& operator>> (std::istream&, Ref<S>&);
template<class S>
std::ostream& operator<< (std::ostream&, const Ref<S>&);

template<class S>
class Ref
{
  S* data; // object pointer
  unsigned* count; // number of references
  void clean ();

public:
  
  Ref () : data(0),  count(0) {}
  Ref (const Ref&);
  Ref (const S& s): data(new S(s)), count(new unsigned(1)) {}
  ~Ref () { clean(); }

  unsigned ref_num () const;

  S& value ();
  S* operator&  ();

  Ref& operator= (const Ref&);
  Ref& operator= (const S&);

  friend std::istream& operator>> <> (std::istream&, Ref<S>&);
  friend std::ostream& operator<< <> (std::ostream&, const Ref<S>&);
};

template<class S>
void Ref<S>::clean()
{
  if (!count)
    return;
  if(*count > 1) 
    (*count)--;
  else {
    delete count;
    delete data;
  }
}

template<class S>
Ref<S>::Ref (const Ref& r) : count(r.count), data(r.data)
{  
  if(count) (*count)++;
}

template<class S>
Ref<S>& Ref<S>::operator= (const Ref<S>& r)
{  
  if(data == r.data)
    return *this;

  clean();
  data = r.data;
  count = r.count;
  if(count)
    (*count)++;
  return *this;
}

template<class S>
Ref<S>& Ref<S>::operator= (const S& s)
{
  clean();
  data = new S(s);
  count = new unsigned(1);
  return *this;
}

template<class S>
S& Ref<S>::value ()
{ 
  if(!data) {
    std::cout << "Ref::value: not initialized\n" << std::endl;
    throw Ref_error(); 
  }
  return *data; 
}

template<class S>
S* Ref<S>::operator&  ()
{ 
  if (!data) {
    std::cout << "Ref::operator&: not initialized\n" << std::endl;
    throw Ref_error(); 
  }
  return  data;
}

template<class S>
unsigned Ref<S>::ref_num () const
{
  if(count)
    return *count;
  else
    return 0;
}

template<class S>
std::istream& operator>> (std::istream& from, Ref<S>& r)
{
  if(!r.count) {
    std::cout << "Ref::operator>>: not initialized" << std::endl;
    throw Ref_error(); 
  }
  from >> *r.data;
  return from;
}
  
template<class S>
std::ostream& operator<< (std::ostream& to, const Ref<S>& r)
{
  if(!r.count) {
    std::cout << "Ref::operator<<: not initialized" << std::endl;
    throw Ref_error(); 
  }
  to << *r.data; 
  return to;
}
