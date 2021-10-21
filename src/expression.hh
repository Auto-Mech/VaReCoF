#ifndef __EXPRESSION__
#define __EXPRESSION__

#include <string>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

namespace Expression {

  // error classes
  class Err {};
  class Find_Err: public Err {}; // search error
  class Form_Err: public Err {}; // format error
  class EOF_Err: public Form_Err {}; // end of file error

  // name-to-value map
  class Varlist 
  {
    map<string, double> vl;

    typedef map<string, double>::iterator It;
    typedef map<string, double>::const_iterator Cit;

  public:
    double& operator[] (const string&);
    double  operator[] (const string&) const;
  
    void insert (const string&, double =0.);
    int size () const { return vl.size(); }
    double value(int i) const;
    const string&  name(int i) const;
  };

  // separator
  class Separ {
    char val;
  public:
    friend istream& operator>> (istream&, Separ&);
    operator char () { return val; }
  };

  istream& operator>> (istream&, Separ&);

  // identifier
  class Ident {
    string val;
  public:
    friend istream& operator>> (istream&, Ident&);
    const string& operator() () const { return val; }
  };

  istream& operator>> (istream&, Ident&);

  char skip_space(istream&);

  class Expr;

  // primitive expression
  class Prim {
    enum Prim_type {NUM, VAR, XPR, FUN};
    Prim_type type;
    void* data;
  public:
    Prim (istream&);
    Prim (const Prim&);
    Prim& operator= (const Prim&);
    ~Prim();
    double operator() (const Varlist&) const;
  };

  // term expression
  class Term {
    vector<pair<Prim,bool> > prims;

  public:
    Term(istream&);
    double operator() (const Varlist& vars) const;
  };

  // expression
  class Expr {
    vector<pair<Term,bool> > terms;

  public:
    Expr(istream&);
    double operator() (const Varlist& vars) const;
  };

  // condition
  class Cmp {
    Expr* xp[2];

  public:
    Cmp (istream&);
    Cmp (const Cmp&);
    Cmp& operator= (const Cmp&);
    ~Cmp ();
    bool operator() (const Varlist&) const;
  };    

  template <class C> class Vec;
  template <class C>  istream& operator>> (istream& from, Vec<C>&);

  template <class C>
  class Vec {
    vector<C> data;
  public:
    friend istream& operator>><C> (istream&, Vec&);
    C& operator[] (int i) { return data[i]; }
    const C& operator[] (int i) const { return data[i]; }
    int size() const { return data.size(); }
  };

  template <class C>  
  istream& operator>> (istream& from, Vec<C>& vec)
  {
    static const char funame [] = "expression::operator>>(istream&, Vec&): ";

    vec.data.clear();

    // opening parenthesis
    char next = skip_space(from);

    if(next != '(') {
      vec.data.resize(1);
      from >> vec.data[0];
      return from;
    }

    from.get(next);

    C cval;
    while(from >> cval) {
      vec.data.push_back(cval);
      skip_space(from);
      from.get(next);
      if(next == ')')
	return from;
      if(next != ',') {
	cout << funame << "unexpected separator " << next 
	     << " encountered (should be ,)\n";
	throw Form_Err();
      }
    }

    // should not be here
    cout << funame << "input error\n";
    throw Form_Err();
  }

}
#endif
