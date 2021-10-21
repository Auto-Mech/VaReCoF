#ifndef ROTD_TMATRIX_HH
#define ROTD_TMATRIX_HH

#include "lapack.h"
#include <iomanip>
#include <exception>
#include <iostream>

/************************* ARRAY *******************************/

class Array_error {};

template <typename T>
class Array
{// class Array

  int real_sz;
  int sz;
  T* aa;
  T* next_to_last;

public:
      
  Array (int =0);
  Array (const Array&);
  ~Array () { if(aa) delete[] aa; }

  T* begin () { return aa; }
  const T* begin () const { return aa; }
  T* end () { return next_to_last; }
  const T* end () const { return next_to_last; }

  T* data () { return aa; }
  const T* data () const { return aa; }

  int size () const { return sz; }
  int capacity () const { return real_sz; }
  void init ();                         // initialize to zero
  void resize(int);
  void compact (); // remove excecive elements

  T& operator[] (int);
  const T&  operator[] (int) const;

  T& front ();
  const T& front () const;

  T& back ();
  const T& back () const;

  Array& operator- ();

  Array& operator = (const Array&);
  Array& operator+= (const Array&);
  Array& operator-= (const Array&);
  Array& operator = (const T*);
  Array& operator+= (const T*);
  Array& operator-= (const T*);

  Array& operator*= (const T&); // multiply each element
  Array& operator/= (const T&); // divide each element
};// class Array

template <typename T>
Array<T>::Array (int s)
  : sz(s), real_sz(s)
{
  const char funame [] = "Array<T>::Array: ";

  if(sz < 0) {
    std::cout << funame << "negative size";
    throw Array_error();
  }

  if(sz)
    aa = new T[sz];
  else
    aa = 0;
  next_to_last = aa + sz;
}

template <typename T>
Array<T>::Array (const Array& ar) : sz(ar.sz), real_sz(ar.real_sz)
{
  if(sz)
    aa = new T[real_sz];
  else
    aa = 0;
  next_to_last = aa + sz;

  const T* it1 = ar.aa;
  for (T* it = aa; it != next_to_last; ++it)
    *it = *it1++;
}

template <typename T>
void Array<T>::compact ()
{
  if(real_sz == sz)
    return;

  T* new_aa = new T[sz];
  for(int i = 0; i < sz; ++i)
    new_aa[i] = aa[i];

  delete[] aa;
  aa = new_aa;
  real_sz = sz;
  next_to_last = aa + sz;
}

template <typename T>
void Array<T>::resize (int s)
{
  const char funame [] = "Array<T>::resize: ";
  if(s == sz)
    return;

  if(s < 0) {
    std::cout << funame << "negative array size\n";
    throw Array_error();
  }

  if(!s) {
    delete[] aa;
    aa = 0;
    sz = 0;
    real_sz = 0;
    next_to_last = 0;
    return;
  }

  if(s <= real_sz) {
    sz = s;
    next_to_last = aa + sz;
    return;
  }

  T* aa_new = new T[s];
  for(int i = 0; i < sz; ++i)
    aa_new[i] = aa[i];

  if(aa)
    delete[] aa;
  aa = aa_new;
  sz = s;
  real_sz = s;
  next_to_last = aa + sz;
}

template <typename T>
void Array<T>::init ()
{
  for (T* it = begin(); it != end(); it++)
    *it = T(0);
}

template <typename T>
T& Array<T>::operator[] (int i)
{
  const char funame [] = "Array<T>::operator[]: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if (i >= sz || i < 0) {
    std::cout << funame << "out of range\n";
    throw Array_error();
  }
#endif

  return aa [i];
}

template <typename T>
const T& Array<T>::operator[] (int i) const
{
  const char funame [] = "Array<T>::operator[]: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if (i >= sz || i < 0) {
    std::cout << funame << "out of range\n";
    throw Array_error();
  }
#endif

  return aa [i];
}

template <typename T>
T& Array<T>::front ()
{
  if(!sz) {
    std::cout << "Array<T>::front(): array is empty\n";
    throw Array_error();
  }
  return *aa;
}
template <typename T>
const T& Array<T>::front () const
{
  if(!sz) {
    std::cout << "Array<T>::front(): array is empty\n";
    throw Array_error();
  }
  return *aa;
}

template <typename T>
T& Array<T>::back ()
{
  if(!sz) {
    std::cout << "Array<T>::back(): array is empty\n";
    throw Array_error();
  }
  return aa[sz-1];
}

template <typename T>
const T& Array<T>::back () const
{
  if(!sz) {
    std::cout << "Array<T>::back(): array is empty\n";
    throw Array_error();
  }
  return aa[sz-1];
}

template <typename T>
Array<T>& Array<T>::operator- ()
{
  for (T* it = begin(); it != end(); ++it)
    *it *= -*it;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator= (const Array& a1)
{
  const char funame [] = "Array<T>::operator=: ";

  if (size() != a1.size()) {
    std::cout << funame << "different size\n";
    throw Array_error();
  }

  const T* it1 = a1.begin();
  for (T* it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const Array& a1)
{
  const char funame [] = "Array<T>::operator+=: ";

  if (size() != a1.size()) {
    std::cout << funame << "different size\n";
    throw Array_error();
  }

  const T* it1 = a1.begin();
  for (T* it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const Array& a1)
{
  const char funame [] = "Array<T>::operator-=: ";

  if (size() != a1.size()) {
    std::cout << funame << "different size\n";
    throw Array_error();
  }

  const T* it1 = a1.begin();
  for (T* it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator*= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it *= val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator/= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it /= val;
  return *this;
}

/****************************** Slice Iterator ****************************/

template <typename T>
class Slice_iter
{// class Slice_iter

  T* pointer_value;
  int stride_value;

public:

  Slice_iter (T* i, int s) : pointer_value(i), stride_value(s) {}

  T* pointer () const { return pointer_value; }
  int stride () const { return stride_value; }

  T& operator* ()       const { return *pointer_value; }
  T& operator[] (int i) const { return pointer_value[i*stride_value]; }

  Slice_iter& operator++ () { pointer_value += stride_value; return *this; }
  Slice_iter  operator++ (int i) 
  { Slice_iter temp=*this; pointer_value += stride_value; return temp; }
  Slice_iter& operator-- () { pointer_value -= stride_value; return *this; }
  Slice_iter  operator-- (int i) 
  { Slice_iter temp=*this; pointer_value -= stride_value; return temp; }

  bool operator!= (const Slice_iter&) const;
  bool operator== (const Slice_iter&) const;

};// class Slice_iter

template <typename T>
bool Slice_iter<T>::operator!= (const Slice_iter& s1) const
{
  const char funame [] = "Slice_iter<T>::operator!=: ";

  if(s1.stride() != stride()) {
    std::cout << funame << "different stride\n";
    throw Array_error();
  }

  return s1.pointer() != pointer();
}

template <typename T>
bool Slice_iter<T>::operator== (const Slice_iter& s1) const
{
  const char funame [] = "Slice_iter<T>::operator==: ";

  if(s1.stride() != stride()) {
    std::cout << funame << "different stride\n";
    throw Array_error();
  }

  return s1.pointer() == pointer();
}

/****************************** Slice **********************************/

template <typename T> class Const_slice;
template <typename T> class Const_slice_iter;

template <typename T>
class Slice
{// class Slice
  typedef Slice_iter<T> Iter;
  Iter start;
  int num;
  //auxiliary
  Iter finish;

public:
	   
  Slice (T* f, int n, int s) : start(f, s), num(n) , finish(f + n*s, s) {}
  Slice (T* f, int n) : start(f, 1), num(n), finish(f + n, 1) {}

  T& operator[] (int i) const {return start[i];}
  int size () const {return num;}

  const Iter& begin () const { return start; }
  const Iter& end   () const { return finish; }
  
  void operator  = (const T&) const;
  void operator += (const T&) const;
  void operator -= (const T&) const;
  void operator *= (const T&) const;
  void operator /= (const T&) const;

  T square () const;
  T length () const { return sqrt(square()); }
  T sum    () const;
  T operator * (const Slice&) const;
  T operator * (const Const_slice<T>&) const;
  T operator * (const T*) const;

  void operator  = (const Slice&) const;
  void operator += (const Slice&) const;
  void operator -= (const Slice&) const;
  void operator  = (const Const_slice<T>&) const;
  void operator += (const Const_slice<T>&) const;
  void operator -= (const Const_slice<T>&) const;
  void operator  = (const T*) const;
  void operator += (const T*) const;
  void operator -= (const T*) const;
  
};// class Slice

template <typename T>
void Slice<T>::operator= (const T& val) const
{
  for (Iter it = begin(); it != end(); ++it)
    *it = val;
}

template <typename T>
void Slice<T>::operator+= (const T& val) const
{
  for (Iter it = begin(); it != end(); ++it)
    *it += val;
}

template <typename T>
void Slice<T>::operator-= (const T& val) const
{
  for (Iter it = begin(); it != end(); ++it)
    *it -= val;
}

template <typename T>
void Slice<T>::operator*= (const T& val) const
{
  for (Iter it = begin(); it != end(); ++it)
    *it *= val;
}

template <typename T>
void Slice<T>::operator/= (const T& val) const
{
  for (Iter it = begin(); it != end(); ++it)
    *it /= val;
}


template <typename T>
T Slice<T>::square () const
{
  T prod = 0.0;
  for (Iter it = begin(); it != end(); ++it)
    prod += (*it)*(*it);
  return prod;
}

template <typename T>
T Slice<T>::sum () const
{
  T res = 0.0;
  for (Iter it = begin(); it != end(); ++it)
    res += *it;
  return res;
}

template <typename T>
T Slice<T>::operator* (const Slice& s1) const
{
  const char funame [] = "Slice<T>::operator*: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  T prod = 0.0;
  Iter it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      prod += (*it1)*(*it);
      ++it1;
    }
  return prod;
}

template <typename T>
T Slice<T>::operator * (const T* it1) const
{
  T prod = 0.0;
  for (Iter it = begin(); it != end(); ++it)
      prod += *(it1++)*(*it);
  return prod;
}

template <typename T>
T Slice<T>::operator* (const Const_slice<T>& s1) const
{
  const char funame [] = "Slice<T>::operator*: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  T prod = 0.0;
  Const_slice_iter<T> it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      prod += (*it1)*(*it);
      ++it1;
    }
  return prod;
}

template <typename T>
void Slice<T>::operator = (const T* it1) const
{
  for (Iter it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator += (const T* it1) const
{
  for (Iter it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator -= (const T* it1) const
{
  for (Iter it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator = (const Slice& s1) const
{
  const char funame [] = "Slice<T>::operator =: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Iter it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator += (const Slice& s1) const
{
  const char funame [] = "Slice<T>::operator +=: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Iter it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator -= (const Slice& s1) const
{
  const char funame [] = "Slice<T>::operator -=: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Iter it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator = (const Const_slice<T>& s1) const
{
  const char funame [] = "Slice<T>::operator =: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Const_slice_iter<T> it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator += (const Const_slice<T>& s1) const
{
  const char funame [] = "Slice<T>::operator +=: ";

  if(s1.size() != size()) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Const_slice_iter<T> it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
}

template <typename T>
void Slice<T>::operator -= (const Const_slice<T>& s1) const
{
  const char funame [] = "Slice<T>::operator -=: ";

  if ( s1.size() != size() ) {
    std::cout << funame << "different length\n";
    throw Array_error();
  }

  Const_slice_iter<T> it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
}

/****************************** Const_slice Iterator ****************************/

template <typename T>
class Const_slice_iter
{// class Const_slice_iter
  const T* pointer_value;
  int stride_value;

public:

  Const_slice_iter (const T* p, int s) : pointer_value(p), stride_value(s) {}
  Const_slice_iter (const Slice_iter<T> s) 
     : pointer_value(s.pointer()), stride_value(s.stride()) {}

  const T* pointer () const { return pointer_value; }
  int stride () const { return stride_value; }

  const T& operator * ()       const { return *pointer_value; }
  const T& operator [] (int i) const { return pointer_value[i*stride_value]; }

  Const_slice_iter& operator ++ () {
    pointer_value += stride_value; return *this; }
  Const_slice_iter  operator ++ (int i) { 
    Const_slice_iter temp=*this; pointer_value += stride_value; return temp; }
  Const_slice_iter& operator -- () {
    pointer_value -= stride_value; return *this; }
  Const_slice_iter  operator -- (int i) {
    Const_slice_iter temp=*this; pointer_value -= stride_value; return temp; }

  bool operator != (const Const_slice_iter&) const;
  bool operator == (const Const_slice_iter&) const;
};// class Const_slice_iter

template <typename T>
bool Const_slice_iter<T>::operator!= (const Const_slice_iter& s1) const
{
  const char funame [] = "Const_slice_iter<T>::operator!=: ";

  if ( s1.stride() != stride() ) {
    std::cout << funame << "not comparable slices\n";
    throw Array_error();
  }

  return s1.pointer() != pointer();
}

template <typename T>
bool Const_slice_iter<T>::operator== (const Const_slice_iter& s1) const
{
  const char funame [] = "Const_slice_iter<T>::operator==: ";

  if ( s1.stride() != stride() ) {
    std::cout << funame << "not comparable slices\n";
    throw Array_error();
  }

  return s1.pointer() == pointer();
}

/************************** Constant Slice **************************/

template <typename T>
class Const_slice
{ // class Const_slice
  typedef Const_slice_iter<T> Iter;
  Iter start;
  int num;
  //auxiliary
  Iter finish;

public:
	   
  Const_slice (const T* f, int n, int s) : start(f, s), num(n) , finish(f + n*s, s) {}
  Const_slice (const T* f, int n) : start(f, 1), num(n), finish(f + n, 1) {}
  Const_slice (const Slice<T>& s) : start(s.begin()), num(s.size()), finish(s.end()) {}

  T operator[] (int i) const { return start[i]; }
  int size () const { return num; }

  const Iter& begin () const { return start; }
  const Iter& end   () const { return finish; }
  
  T square () const;
  T sum () const;
  T operator* (const Const_slice&) const;
  T operator* (const Slice<T>&) const;
  T operator* (const T*) const;
}; // class Const_slice

template <typename T>
T Const_slice<T>::square () const
{
  T prod = 0.0;
  for (Iter it = begin(); it != end(); ++it)
    prod += (*it)*(*it);
  return prod;
}

template <typename T>
T Const_slice<T>::sum () const
{
  T res = 0.0;
  for (Iter it = begin(); it != end(); ++it)
    res += *it;
  return res;
}

template <typename T>
T Const_slice<T>::operator* (const Const_slice& s1) const
{
  const char funame [] = "Const_slice<T>::operator*: ";
  if (s1.size() != size()) {
    std::cout << funame << "different size\n";
    throw Array_error();
  }

  T prod = 0.0;
  Iter it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      prod += (*it1)*(*it);
      ++it1;
    }
  return prod;
}

template <typename T>
T Const_slice<T>::operator* (const Slice<T>& s1) const
{
  const char funame [] = "Const_slice<T>::operator*: ";
  if (s1.size() != size()) {
    std::cout << funame << "different size\n";
    throw Array_error();
  }

  T prod = 0.0;
  Slice_iter<T> it1 = s1.begin();
  for (Iter it = begin(); it != end(); ++it)
    {
      prod += (*it1)*(*it);
      ++it1;
    }
  return prod;
}

template <typename T>
T Const_slice<T>::operator* (const T* it1) const
{
  T prod = 0.0;
  for (Iter it = begin(); it != end(); ++it)
    {
      prod += (*it1)*(*it);
      ++it1;
    }
  return prod;
}

/*********************** MATRIX ********************************/

template <typename T, int DIM>
class TVector;

template <typename T, int DIM> 
class TMatrix : public Array<T>
{
public:

  // Constructors
  TMatrix () : Array<T>(DIM*DIM) {}          

  explicit TMatrix (const T& val) : Array<T>(DIM*DIM) 
  { *this = val; }

  int dim () const { return DIM; }// matrix' dimension

  Array<T>& base () { return *this; }// reference to the base
  const Array<T>& base () const { return *this; }

  // Fortran style indexing
  //
  T& operator() (int, int);
  const T&  operator() (int, int) const;

  // Slices
  //
  Slice<T>       column (int);
  Const_slice<T> column (int) const;
  Slice<T>       row (int);
  Const_slice<T> row (int) const;
  Slice<T>       diagonal (int);
  Const_slice<T> diagonal (int) const;

  //arithmetics
  TMatrix& operator = (const TMatrix& m)
  { base()  = m.base(); return *this; }
  TMatrix& operator+= (const TMatrix& m)
  { base() += m.base(); return *this; }
  TMatrix& operator-= (const TMatrix& m)
  { base() -= m.base(); return *this; }

  TMatrix& operator = (const T* m) { base()  = m; return *this; }
  TMatrix& operator+= (const T* m) { base() += m; return *this; }
  TMatrix& operator-= (const T* m) { base() -= m; return *this; }

  TMatrix& operator*= (const T& val) { base() *= val; return *this; }
  TMatrix& operator/= (const T& val) { base() /= val; return *this; }

  TMatrix& operator = (const T& val) { Array<T>::init(); diagonal(0) = val; return *this; }
  TMatrix& operator+= (const T& val) { diagonal(0) += val; return *this; }
  TMatrix& operator-= (const T& val) { diagonal(0) -= val; return *this; }

  // matrix multiplication
  //
  void multiply (const TMatrix&, const TMatrix&);
  

  // operations

  // inverts matrix: 1/A
  //
  void invert ();

  // transpose matrix elements: A(i,j) = A(j,i)
  //
  void transpose ();

  // eigenvalues and eigenvectors
  //
  void eigenv (T*);
  T trace () const { return diagonal(0).sum(); }

};// class TMatrix

template <typename T, int DIM> 
inline Slice<T> TMatrix<T, DIM>::column (int i) 
{ 
  const char funame [] = "TMatrix<T, DIM>::column: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i < 0 || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif

  return Slice<T> (Array<T>::begin() + DIM*i, DIM);
}

template <typename T, int DIM> 
inline Const_slice<T> TMatrix<T, DIM>::column (int i) const 
{ 
  const char funame [] = "TMatrix<T, DIM>::column: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i < 0 || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif

  return Const_slice<T> (Array<T>::begin() + DIM*i, DIM); 
}

template <typename T, int DIM> 
inline Slice<T> TMatrix<T, DIM>::row (int i) 
{ 
  const char funame [] = "TMatrix<T, DIM>::row: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i < 0 || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif

  return Slice<T> (Array<T>::begin() + i, DIM, DIM); 
}

template <typename T, int DIM> 
inline Const_slice<T> TMatrix<T, DIM>::row (int i) const 
{ 
  const char funame [] = "TMatrix<T, DIM>::row: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i < 0 || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif

  return Const_slice<T> (Array<T>::begin() + i, DIM, DIM); 
}

template <typename T, int DIM> 
inline Slice<T> TMatrix<T, DIM>::diagonal (int i) 
{
  const char funame [] = "TMatrix<T, DIM>::diagonal: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i <= -DIM || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif
  if (i < 0)
    return Slice<T>(Array<T>::begin() - i, DIM - i, DIM + 1);
  else
    return Slice<T>(Array<T>::begin() + i*DIM, DIM - i, DIM + 1);
}

template <typename T, int DIM> 
inline Const_slice<T> TMatrix<T, DIM>::diagonal (int i) const 
{
  const char funame [] = "TMatrix<T, DIM>::diagonal: ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if ( i <= -DIM || i >= DIM ) {
    std::cout << funame << "index i = " << i << " is out of range\n";
    throw Array_error();
  }
#endif
  if (i < 0)
    return Const_slice<T>(Array<T>::begin() - i, DIM + i, DIM + 1);
  else
    return Const_slice<T>(Array<T>::begin() + i*DIM, DIM - i, DIM + 1);
}

template <typename T, int DIM> 
T& TMatrix<T, DIM>::operator() (int i, int j)
{
  const char funame [] = "TMatrix<T, DIM>::operator(): ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if (i >= DIM || j >= DIM || i < 0 || j < 0) {
    std::cout << funame << "indices (i=" << i << ", j=" << j 
	      << ") out of range\n";
    throw Array_error();
  }
#endif
  return base()[i+j*DIM];
}

template <typename T, int DIM> 
const T&  TMatrix<T, DIM>::operator() (int i, int j) const 
{
  const char funame [] = "TMatrix<T, DIM>::operator(): ";

#ifndef NO_ARRAY_SUBSCRIPT_CHECKING
  if (i >= DIM || j >= DIM || i < 0 || j < 0) {
    std::cout << funame << "indices (i=" << i << ", j=" << j 
	      << ") out of range\n";
    throw Array_error();
  }
#endif
  return base()[i+j*DIM];
}

template <typename T, int DIM> 
void TMatrix<T, DIM>::multiply (const TMatrix& m1, const TMatrix& m2)
{
  for (int i = 0; i < dim(); ++i)
    for (int j = 0; j < dim(); ++j)    
      (*this) (i, j) = m1.row(i)*m2.column(j);
}

template <typename T, int DIM>
void TMatrix<T, DIM>::invert()
{
  const char funame [] = "TMatrix<T, DIM>::invert: ";
  
  Array<integer> ipiv (DIM);
  TMatrix res(1.0);
  integer info;
  integer dim = DIM;
  // Lapack subroutine dgesv
  dgesv_ (&dim, &dim, Array<T>::data(), &dim, ipiv.data(), res.data(), 
	  &dim, &info);
  if (info) {
    std::cout << funame << "lapack error: dgesv info = " << info << "\n";
    throw Array_error();
  }
  *this = res;
}

template <typename T, int DIM>
void TMatrix<T, DIM>::transpose()
{
  for (int i1 = 0; i1 < dim (); i1++)
    for (int i2 = i1+1; i2 < dim (); i2++)
      std::swap((*this)(i1, i2), (*this)(i2, i1));
}

template <typename T, int DIM>
void TMatrix<T, DIM>::eigenv (T* eval)
{
  const char funame [] = "TMatrix<T, DIM>::eigenv: ";
 
  integer lwork = DIM*DIM+5;
  static Array<double> work (lwork);
  integer info = 0;
  // Lapack subroutine dsyev
  const char* job = "V";
  const char* uplo = "U";
  integer dim = DIM;
  dsyev_ (job, uplo, &dim, Array<T>::data(), &dim, eval, work.data(), &lwork, &info);
  if (info) {
    std::cout << funame << "lapack error: dsyev info = " << info << "\n";
    throw Array_error();
  }
}

/**************************** Vector **************************************/

template <typename T, int DIM>
class TVector : public Array<T>
{// class TVector

public:

  // Auxiliary
  int dim () const { return DIM; }        // vector's dimension

  Array<T>& base () { return *this; }        // reference to the base
  const Array<T>& base () const { return *this; }

  Slice<T> self () { return Slice<T> (Array<T>::begin(), Array<T>::size()); } // vector as a slice
  Const_slice<T> self () const { return Const_slice<T> (Array<T>::begin(), Array<T>::size()); }

  // Constructors
  TVector () : Array<T>(DIM) {}          
  TVector (const T& val) : Array<T>(DIM) { const_cast<TVector*>(this)->self () = val; }

  T& operator [] ( int i) { return base()[i]; }               // subscribing
  const T&  operator [] (int i) const { return base()[i]; }

  // Arithmetics
  TVector& operator  = (const TVector& v) { base()  = v.base(); return *this; }
  TVector& operator += (const TVector& v) { base() += v.base(); return *this; }
  TVector& operator -= (const TVector& v) { base() -= v.base(); return *this; }

  TVector& operator  = (const T* v) { base()  = v; return *this; }
  TVector& operator += (const T* v) { base() += v; return *this; }
  TVector& operator -= (const T* v) { base() -= v; return *this; }

  TVector& operator  = (const Slice<T>& s) { self()  = s; return *this; }
  TVector& operator += (const Slice<T>& s) { self() += s; return *this; }
  TVector& operator -= (const Slice<T>& s) { self() -= s; return *this; }

  TVector& operator  = (const Const_slice<T>& s) { self()  = s; return *this; }
  TVector& operator += (const Const_slice<T>& s) { self() += s; return *this; }
  TVector& operator -= (const Const_slice<T>& s) { self() -= s; return *this; }

  TVector& operator *= (const T& val) { base() *= val; return *this; }
  TVector& operator /= (const T& val) { base() /= val; return *this; }

  TVector& operator  = (const T& val) { self()  = val; return *this; }
  TVector& operator += (const T& val) { self() += val; return *this; }
  TVector& operator -= (const T& val) { self() -= val; return *this; }

};// class TVector

/*************************  Tensors ***************************************/

template<class C>
class Array_2 : public Array<C>
{// class Array_2

      int dd_1;
      int dd_2;

   public:

      Array_2 (int d1, int d2) : Array<C>(d1*d2), dd_1(d1), dd_2(d2) {}

      C operator () (int i1, int i2) const
      {
	 return Array<C>::data() [i1 + dd_1 * i2];
      }

      C& operator () (int i1, int i2)
      {
	 return Array<C>::data() [i1 + dd_1 * i2];
      }

      int dim_1 () const { return dd_1; }
      int dim_2 () const { return dd_2; }
};

template <class C>
class Array_3 : public Array<C>
{// class Array_3

      int dd_1;
      int dd_2;
      int dd_3;

   public:

      Array_3 (int d1, int d2, int d3)
	: Array<C>(d1*d2*d3), dd_1(d1), dd_2(d2), dd_3(d3)
      {}

      C operator () (int i1, int i2, int i3) const
      {
	 return Array<C>::data() [i1 + dd_1 * i2 + dd_1 *dd_2 * i3];
      }

      C& operator () (int i1, int i2, int i3)
      {
	 return Array<C>::data() [i1 + dd_1 * i2 + dd_1 *dd_2 * i3];
      }

      int dim_1 () const { return dd_1; }
      int dim_2 () const { return dd_2; }
      int dim_3 () const { return dd_3; }

};

template <typename T>
std::istream& operator>> (std::istream& from, Array<T>& a)
{
  for (int i = 0; i < a.size(); ++i)
    from >> a.data()[i];
  return from;
}

template <typename T>
std::ostream& operator<< (std::ostream& to, const Array<T>& a)
{
  const int width = to.precision() + 7;
  const int count_max = 80 / (width + 2);
  int count = 0;
  for (int i = 0; i < a.size(); ++i)
    {
      to << std::setw(width) << a.data()[i] << "  ";
      if (++count == count_max)
	{
	  to << "\n";
	  count = 0;
	}
    }
  to << std::endl;
  return to;
}

/** Seems that template specialization does not work with g++
template <>
std::ostream& operator<< (std::ostream& to, const Array<int>& a)
{
  const int width = 6;
  const int count_max = 80 / (width + 2);
  int count = 0;
  for (int i = 0; i < a.size(); ++i)
    {
      to << std::setw(width) << a.data()[i] << "  ";
      if (++count == count_max)
	{
	  to << "\n";
	  count = 0;
	}
    }
  to << std::endl;
  return to;
} 
**/ 

/**************** matrix by vector multiplication **************************/

template <typename T, int DIM>
TVector<T, DIM> operator * (const TMatrix <T, DIM>& m, const TVector<T, DIM>& v) 
// Matrix x Vector
{
  TVector<T, DIM> temp;
  Const_slice<T> cs = v.self();
  for (int i = 0; i < DIM; ++i)
    temp[i] = cs*m.row(i);
  return temp;
}

template <typename T, int DIM>
TVector<T, DIM> operator * (const TMatrix <T, DIM>& m, const T* v) 
// Matrix x Vector
{
  TVector<T, DIM> temp;
  for (int i = 0; i < DIM; ++i)
    temp[i] = m.row(i)*v;
  return temp;
}

template <typename T, int DIM>
TVector<T, DIM> operator * (const TVector<T, DIM>& v, const TMatrix <T, DIM>& m) 
// Vector x Matrix
{
  TVector<T, DIM> temp;
  Const_slice<T> cs = v.self();
  for (int i = 0; i < DIM; ++i)
    temp[i] = cs*m.column(i);
  return temp;
}

template <typename T, int DIM>
TVector<T, DIM> operator * (const T* v, const TMatrix <T, DIM>& m) 
// Vector x Matrix
{
  TVector<T, DIM> temp;
  for (int i = 0; i < DIM; ++i)
    temp[i] = m.column(i)*v;
  
  return temp;
}

template <typename T, int DIM>
void matrix_vector_product (const TMatrix <T, DIM>& m, const T* v, T* res)
{
  for (int i = 0; i < DIM; ++i)
    res[i] = m.row(i)*v;
}

template <typename T, int DIM>
void vector_matrix_product (const T* v, const TMatrix <T, DIM>& m, T* res)
{
  for (int i = 0; i < DIM; ++i)
    res[i] = m.column(i)*v;
}

template <typename T>
T vector_product (int i, const T* v1, const T* v2)
{
  const char funame [] = "vector_product: ";
  switch(i) {
  case 0:
    return v1 [1] * v2 [2] - v2 [1] * v1 [2];
  case 1:
    return v1 [2] * v2 [0] - v2 [2] * v1 [0];
  case 2:
    return v1 [0] * v2 [1] - v2 [0] * v1 [1];
  default:
    std::cout << funame << "index " << i << " out of range\n";
    throw Array_error();
  }
}

template <typename T>
void vector_product (const T* v1, const T* v2, T* res)
{
  register double 
    x1 = v1 [0], y1 = v1 [1], z1 = v1 [2],
    x2 = v2 [0], y2 = v2 [1], z2 = v2 [2];

  res [0] = y1 * z2 - z1 * y2;
  res [1] = z1 * x2 - x1 * z2;
  res [2] = x1 * y2 - y1 * x2;
}

/************************* Matrix *******************************/

template<typename T>
class Matrix {
  int _size[2];
  T* data;
public:
  Matrix(int d0, int d1) : data(new T[d0*d1]) 
  {_size[0] = d0; _size[1] = d1; }
  Matrix(const Matrix&);
  Matrix& operator= (const Matrix&);
  ~Matrix() { delete[] data; }

  const T* operator[] (int i) const { return data + i * _size[1]; }
  T* operator[] (int i) { return data + i * _size[1]; }
  int size (int dim) const { return _size[dim]; }

  void read_by_line    (std::istream&);
  void read_by_column  (std::istream&);
  void write_by_line   (std::ostream&);
  void write_by_column (std::ostream&);
};

template<typename T>
Matrix<T>::Matrix (const Matrix<T>& m) 
{
  for(int i = 0; i < 2; ++i)
    _size[i] = m._size[i];

  int n = _size[0] * _size[1];
  data = new T[n];
  for(int i = 0; i < n; ++i)
    data[i] = m.data[i];
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
  delete[] data;

  for(int i = 0; i < 2; ++i)
    _size[i] = m._size[i];

  int n = _size[0] * _size[1];
  data = new T[n];
  for(int i = 0; i < n; ++i)
    data[i] = m.data[i];

  return *this;
}
    
template<typename T>
void Matrix<T>::write_by_line(std::ostream& to)
{
  for(int i = 0; i < _size[0]; ++i)
    for(int j = 0; j < _size[1]; ++j) {
      to << data[i*_size[1]+j] << "  ";
    to << "\n";
  }
  to << "\n";
}

template<typename T>
void Matrix<T>::write_by_column(std::ostream& to)
{
  for(int j = 0; j < _size[1]; ++j) {
    for(int i = 0; i < _size[0]; ++i)
      to << data[i*_size[1]+j] << "  ";
    to << "\n";
  }
  to << "\n";
}

template<typename T>
void Matrix<T>::read_by_column(std::istream& from)
{
  for(int j = 0; j < _size[1]; ++j)
    for(int i = 0; i < _size[0]; ++i)
      from >> data[i*_size[1]+j];
}

template<typename T>
void Matrix<T>::read_by_line(std::istream& from)
{
  for(int i = 0; i < _size[0]; ++i)
    for(int j = 0; j < _size[1]; ++j)
      from >> data[i*_size[1]+j];
}

#endif
