/****************************************************************
 ******************* Random Number Generators *******************
 ****************************************************************/

double Random::norm ()
{// normal distribution RNG

  static bool first = true;
  static double x1, x2;
  double r;
  if (first)
    {
      first = false;
      do
	{
	  x1 = 2.0 * Random::flat() - 1.0;
	  x2 = 2.0 * Random::flat() - 1.0;
	  r = x1 * x1 + x2 * x2;
	}
      while (r > 1.0);
      r = sqrt(-2.0 * log(r)/r);
      x1 *= r;
      x2 *= r;
      return x1;
    }
  else
    {
      first = true;
      return x2;
    }
}

double Random::exp ()
{// RNG for the distribution exp(-u**2/2)u*du
  return sqrt(-2.0 * log(Random::flat()));
}

void Random::orient (double* vec, int dim)
{ // randomly orients vector _vec_ of the dimension _dim_ of unity length
  double norm, temp;
  do {
    norm = 0.0;
    for (int i = 0; i < dim; ++i)
      {
	temp = 2.0 * Random::flat () - 1.0;
	norm += temp * temp;
	vec[i] = temp;
      }
  } while (norm >= 1.0 || norm < 0.0001);

  norm = sqrt(norm);
  for (int i = 0; i < dim; ++i)
    vec[i] /= norm;
}

