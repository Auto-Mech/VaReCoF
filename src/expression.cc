#include "expression.hh"

#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>

/*********************** class Varlist ****************************/

double& Expression::Varlist::operator[] (const string& s)
{
  static const char funame [] = "Expression::Varlist::operator[]: ";

  It it = vl.find(s);
  if (it == vl.end()) {
    cout << funame << "variable " << s << " is not defined\n";
    throw Find_Err();
  }
  else
    return it->second;
}

double Expression::Varlist::operator[] (const string& s) const
{
  static const char funame [] = "Expression::Varlist::operator[]: ";

  Cit it = vl.find(s);
  if (it == vl.end()) {
    cout << funame << "variable " << s << " is not defined\n";
    throw Find_Err();
  }
  else
    return it->second;
}

void Expression::Varlist::insert(const string& s, double d)
{
  static const char funame [] = "Expression::Varlist::insert: ";

  It it = vl.find(s);
  if (it != vl.end()) {
    cout << funame << "variable " << s << " is allready in the list\n";
    throw Find_Err();
  }
  else
    vl[s] = d;
}

double  Expression::Varlist::value(int i) const
{ 
  Cit it = vl.begin();
  advance(it, i);
  return it->second;
}

const string&  Expression::Varlist::name(int i) const
{ 
  Cit it = vl.begin();
  advance(it, i);
  return it->first;
}

/****************** Separators and Identifiers *******************/

char Expression::skip_space(istream& from)
{
  static const char funame [] = "Expression::skip_space: ";

  char next;
  string comment;
  while(from.get(next))
    if(isspace(next))
      continue;
    else if(next == '#')
      getline(from, comment);
    else
      break;

  if(from.eof()) {
    //cout << funame << "EOF encountered\n";
    throw EOF_Err();
  }
  else if(!from) {
    cout << funame << "input stream error\n";
    throw Form_Err();
  }

  from.putback(next);    
  return next;
}

istream& Expression::operator>> (istream& from, Separ& sep)
{
  skip_space(from);
  return from.get(sep.val);
}

istream& Expression::operator>> (istream& from, Ident& id)
{
  static const char funame [] = "Expression::operator>>(istream&, Ident&): ";

  skip_space(from);

  char next;
  from >> next;
  if(!isalpha(next)) {
    cout << funame << "identifier does not begin with letter";
    throw Form_Err();
  }

  id.val = next;
  while(from.get(next) && isalnum(next))
    id.val += next;

  if(!from)
    return from;

  return from.putback(next);
}

/************************* Expr ****************************/

Expression::Expr::Expr (istream& from)
{
  static const char funame [] = "Expression::Expr::Expr: ";

  string comment;
  char next;
  bool is_first_term = true;

  while(from.get(next)) {
    if (isspace(next)) continue;

    switch(next) {
    case '#':
      getline(from, comment);
      continue;
    case '+':
      terms.push_back(pair<Term,bool>(Term(from),true));
      break;
    case '-':
      terms.push_back(pair<Term,bool>(Term(from),false));
      break;
    default:
      if (!is_first_term) {
	from.putback(next);
	return;
      }
      from.putback(next);
      terms.push_back(pair<Term,bool>(Term(from),true));
    }
    is_first_term = false;
  }

  if(!from.eof()) {
    cout << funame << "stream is corrupted\n";
    throw Form_Err();
  }
}

double Expression::Expr::operator() (const Varlist& vars) const
{
  double res = 0.;
  for(int i = 0; i < terms.size(); ++i)
    if(terms[i].second)
      res += terms[i].first(vars);
    else
      res -= terms[i].first(vars);

  return res;
}

/********************* class Term *************************/

Expression::Term::Term (istream& from)
{
  static const char funame [] = "Expression::Term::Term: ";

  prims.push_back(pair<Prim,bool>(Prim(from), true));

  string comment;
  char next;
  while(from.get(next)) {
    if (isspace(next)) continue;
    switch(next) {
    case '#':
      getline(from, comment);
      continue;
    case '*':
      prims.push_back(pair<Prim,bool>(Prim(from), true));
      break;
    case '/':
      prims.push_back(pair<Prim,bool>(Prim(from), false));
      break;
    default:
      from.putback(next);
      return;
    }
  }

  if(!from.eof()) {
    cout << funame << "stream is corrupted\n";
    throw Form_Err();
  }
}

double Expression::Term::operator() (const Varlist& vars) const 
{  
  double res = prims[0].first(vars);
  for(int i = 1; i < prims.size(); ++i)
    if(prims[i].second)
      res *= prims[i].first(vars);
    else
      res /= prims[i].first(vars);
  return res;
}

/********************* class Prim *************************/


Expression::Prim::Prim (istream& from)
{
  static const char funame [] = "Expression::Prim::Prim: ";

  string comment;
  char next = skip_space(from);
  double dtemp;

  // double const
  if (isdigit(next) || next == '.') {
    type = NUM;
    if (!(from >> dtemp)) {
      cout << funame << "cannot read the double\n";
      throw Form_Err();
    }
    data = new double(dtemp);
    return;
  }

  // variable or function
  if (isalpha(next)) {
    Ident name;
    from >> name;
    try {
      next = skip_space(from);
    }
    catch(EOF_Err) {
      type = VAR;
      data = new string(name());
      return;
    }

    // variable
    if (next != '(') {
      type = VAR;
      data = new string(name());
      return;
    }

    // function
    from.get(next);
    Expr xtemp(from);
    skip_space(from);
    from.get(next);
    if(next != ')') {
      cout << funame << "format error in the " << name() <<" function\n";
      throw Form_Err();
    }
    type = FUN;
    data = new pair<string,Expr>(name(), xtemp);
    return;
  }// variable or function

  // subexpression
  if(next == '(') {
    from.get(next);
    Expr xtemp(from);
    skip_space(from);
    from.get(next);
    if(next != ')') {
      cout << funame << "closing parenthesis is expected\n";
      throw Form_Err();
    }
    type = XPR;
    data = new Expr(xtemp);
    return;
  }// subexpression

  cout << funame << "unexpected character " << next << "\n";
  throw Form_Err();
}

Expression::Prim::~Prim ()
{
  switch(type) {
  case NUM:
    delete static_cast<double*>(data);
    break;
  case VAR:
    delete static_cast<string*>(data);
    break;
  case XPR:
    delete static_cast<Expr*>(data);
    break;
  case FUN:
    delete static_cast<pair<string,Expr>*>(data);
    break;
  }
}

Expression::Prim::Prim(const Prim& pe)
{
  static const char funame [] = "Expression::Prim::Prim: ";

  type = pe.type;

  switch(type) {
  case NUM:
    data = new double(*static_cast<double*>(pe.data));
    break;
  case VAR:
    data = new string(*static_cast<string*>(pe.data));
    break;
  case XPR:
    data = new Expr(*static_cast<Expr*>(pe.data));
    break;
  case FUN:
    data = new pair<string,Expr>(*static_cast<pair<string,Expr>*>(pe.data));
    break;
  default:
    cout << funame << "unknown type " << type << "\n";
    throw Form_Err();
  }
}

Expression::Prim& Expression::Prim::operator= (const Prim& pe)
{
  static const char funame [] = "Expression::Prim::operator=: ";

  type = pe.type;

  switch(type) {
  case NUM:
    delete static_cast<double*>(data);
    data = new double(*static_cast<double*>(pe.data));
    break;
  case VAR:
    delete static_cast<string*>(data);
    data = new string(*static_cast<string*>(pe.data));
    break;
  case XPR:
    delete static_cast<Expr*>(data);
    data = new Expr(*static_cast<Expr*>(pe.data));
    break;
  case FUN:
    delete static_cast<pair<string,Expr>*>(data);
    data = new pair<string,Expr>(*static_cast<pair<string,Expr>*>(pe.data));
    break;
  default:
    cout << funame << "unknown type " << type << "\n";
    throw Form_Err();
  }
  return *this;
}

double Expression::Prim::operator() (const Varlist& vars) const
{
  static const char funame [] = "Expression::Prim::operator(): ";

  string stemp;
  double dtemp;

  switch (type) {
  case NUM:
    return *static_cast<double*>(data);
  case VAR:
    return vars[*static_cast<string*>(data)];
  case XPR:
    return (*static_cast<Expr*>(data))(vars);
  case FUN:
    stemp = static_cast<pair<string,Expr>*>(data)->first;
    dtemp = static_cast<pair<string,Expr>*>(data)->second(vars);
    if (stemp == "cos")
      return cos(dtemp*M_PI/180.);
    else if (stemp == "sin")
      return sin(dtemp*M_PI/180.);
    else if (stemp == "tan")
      return tan(dtemp*M_PI/180.);
    else if (stemp == "acos")
      return acos(dtemp)*180./M_PI;
    else if (stemp == "asin")
      return asin(dtemp)*180./M_PI;
    else if (stemp == "atan")
      return atan(dtemp)*180./M_PI;
    else if (stemp == "exp")
      return exp(dtemp);
    else if (stemp == "log")
      return log(dtemp);
    else {
      cout << funame << "unknown function " << stemp << "\n";
      throw Find_Err();
    }
  default:
    cout << funame << "unknown type " << type << "\n";
    throw Find_Err();
  }
}

Expression::Cmp::Cmp(istream& from)
{
  static const char funame [] = "Expression::Cmp::Cmp: ";
  Expr xpr0(from);
  skip_space(from);
  char next = from.get();
  Expr xpr1(from);

  switch(next) {
  case '<':
    xp[0] = new Expr(xpr0);
    xp[1] = new Expr(xpr1);
    break;
  case '>':
    xp[0] = new Expr(xpr1);
    xp[1] = new Expr(xpr0);
    break;
  default:
    cout << funame << "wrong symbol\n";
    throw Form_Err();
  }
}

Expression::Cmp::Cmp(const Cmp& ce)
{
  for(int i = 0; i < 2; ++i)
    xp[i] = new Expr(*ce.xp[i]);
}

Expression::Cmp& Expression::Cmp::operator= (const Cmp& ce)
{
  for(int i = 0; i < 2; ++i) {
    delete xp[i];
    xp[i] = new Expr(*ce.xp[i]);
  }
  return *this;
}

Expression::Cmp::~Cmp()
{
  for(int i = 0; i < 2; ++i)
    delete xp[i];
}

bool Expression::Cmp::operator() (const Varlist& vars) const
{


  if ((*xp[0])(vars) < (*xp[1])(vars))
    return true;
  else
    return false;
}

