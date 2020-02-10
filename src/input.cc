#include <cstdio>
#include <sstream>
#include <cmath>

#include "input.hh"
#include "units.hh"

char skip_space(std::istream& from)
{
  const char funame [] = "skip_space: ";

  char next;
  std::string comment;
  while(from.get(next))
    if(isspace(next))
      continue;
    else if(next == '#')
      getline(from, comment);
    else
      break;

  if(from.eof())
    return EOF;  

  if(from.fail()) {
    std::cout << funame << "input stream error\n";
    throw Error::Form_Err();
  }

  from.putback(next);
  return next;
}

void ReadInt::operator () (std::istream& from)
{
  const char funame [] = "ReadInt::operator () (std::istream&): ";

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cout << funame << "cannot read integer\n";
    throw Error::Form_Err();
  }

  std::string comment;
  iss >> comment;
  if(iss && comment[0] != '#') {
    std::cout << funame << "comment does not start with #\n";
    throw Error::Form_Err();
  }
}

void ReadIarr::operator () (std::istream& from)
{
  const char funame [] = "ReadIarr::operator () (std::istream&): ";
  
  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  int itemp;
  while(iss >> itemp)
    _data.push_back(itemp);

  if(!_data.size()) {
    std::cout << funame << "cannot read integer array\n";
    throw Error::Form_Err();
  }
}

void ReadDouble::operator () (std::istream& from)
{
  const char funame [] = "ReadDouble::operator () (std::istream&): ";

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cout << funame << "cannot read double\n";
    throw Error::Form_Err();
  }

  std::string unit;
  iss >> unit;
  if(iss && unit[0] !='#')
    _data *= Phys_const::str2fac(unit);
}

void ReadString::operator () (std::istream& from)
{
  const char funame [] = "ReadString::operator () (std::istream&): ";

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cout << funame << "cannot read string\n";
    throw Error::Form_Err();
  }
}

void ReadDarr::operator () (std::istream& from)
{
  const char funame [] = "ReadDarr::operator () (std::istream&): ";

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  char next = skip_space(iss);

  if(next == EOF) {
    std::cout << funame << "no array specification is found\n";
    throw Error::Form_Err();
  }

  // progression with variable step (old fashion input)
  if(next == '.' || next == '-' || next == '+' || isdigit(next)) { 
    double d, dd, incr; 
    int n;

    iss >> d >> dd >> incr >> n;

    if(!iss) {
      std::cout << funame << "old style array specification is corrupted\n";
      throw Error::Form_Err();
    }

    std::string unit;
    iss >> unit;
    if(iss && unit[0] != '#') {
      double fac = Phys_const::str2fac(unit);
      d *= fac;
      dd *= fac;
    }

    double val = d;
    double step = dd;
    for(int i = 0; i < n; ++i) {
      _data.push_back(val);
      val += step;
      step *= incr;
    }
    return;
  }
  
  // new style input
  while(next == '@' || next  == '*') {
    bool is_geom = true;;
    if(next == '@')// arithmetic progression
      is_geom = false;

    double term;
    iss >> term;
    if(!iss) {
      std::cout << funame << "cannot read the progression first term\n";
      throw Error::Form_Err();
    }
    
    bool is_step = false;
    bool is_size = false;
    bool is_end  = false;
    double step = 0.;
    double  end = 0.;
    int    size = 0;

    for(int i = 0; i < 2; ++i) {
      next = skip_space(iss);
      if(next == EOF) {
	std::cout << funame << "unexpected end of stream\n";
	throw Error::Form_Err();
      }
      iss.get();
      switch(toupper(next)) {
      case 'N':
	if(is_size) {
	  std::cout << "the progression size has been allready defined\n";
	  throw Error::Form_Err();
	}
	is_size = true;
	iss >> size;
	if(!iss) {
	  std::cout << funame << "cannot read the progression size\n";
	  throw Error::Form_Err();
	}
	break;
      case 'E':
	if(is_end) {
	  std::cout << "the progression last term has been allready defined\n";
	  throw Error::Form_Err();
	}
	is_end = true;
	iss >> end;
	if(!iss) {
	  std::cout << funame << "cannot read the progression last term\n";
	  throw Error::Form_Err();
	}	
	break;
      case 'S':
	if(is_step) {
	  std::cout << "the progression step has been allready defined\n";
	  throw Error::Form_Err();
	}
	is_step = true;
	iss >> step;
	if(!iss) {
	  std::cout << funame << "cannot read the progression step\n";
	  throw Error::Form_Err();
	}
	break;
      default:
	std::cout << funame << "unknown control character: " << next << "\n";
	throw Error::Form_Err();
      }
    }
    
    if(!is_step && size > 1) {
      if(is_geom)
	step = exp(log(end/term)/(size-1));
      else
	step = (end-term)/(size-1);
    }
    if(!is_size) {
      if(is_geom) {
	if(step == 1.) {
	  std::cout << funame << "infinite geometrical progression\n";
	  throw Error::Form_Err();
	}
	size = int(log(end/term) / log(step)) + 1;
      }
      else {
	if(step == 0.) {
	  std::cout << funame << "infinite arithmetical progression\n";
	  throw Error::Form_Err();
	}
	size = int((end-term)/step) + 1;
      }
    }

    for(int i = 0; i < size; ++i) {
      _data.push_back(term);
      if(is_geom)
	term *= step;
      else
	term += step;
    }
  }

  if(next == EOF || next == '#')
    return;
  
  std::string unit;
  iss >> unit;
  double fac = Phys_const::str2fac(unit);
  for(int i = 0; i < _data.size(); ++i)
    _data[i] *= fac;
}

/********************************** InputKey ***********************************/

std::set<const std::string*> InputKey::pool;

void InputKey::show_all ()
{
  for(std::set<const std::string*>::iterator it = pool.begin(); it != pool.end(); ++it)
    std::cout << *(*it) << "\n";
}
