namespace slatec {
  extern "C" {
    typedef double (*dgaus8_t) (const double&);
    typedef void   (*ddeabm_t) (const double& x, const double* y, 
				double* dy, void* rp, void* ip);
  }
}

extern "C" 
{

  /***************************************************************
   *                  Grid integrator
   ***************************************************************/
  void davint_ (const double* x, const double* y, const int& n, 
		const double& xlo, const double& xup, double& ans, 
		int& ierr); 
  
  /***************************************************************
   *                  Spline fit
   ***************************************************************/

  void  dbint4_ (const double* x, const double* y, const int& ndata, 
		 const int& ibcl, const int& ibcr, 
		 const double& fbcl, const double& fbcr, 
		 const int& kntopt, double* t,
		 double* bcoef, int& n, int& k, double* w);

  double dbvalu_ (const double* t, const double* a, const int& n, 
		  const int& k, const int& ideriv, const double& x, 
		  int& inbv, double* work);

  /****************************************************************
   *     Differential equations solver by the Adams-Bashforth-Moulton
   *     Predictor-Corrector formulas of orders one through twelve
   ****************************************************************/

  void ddeabm_(slatec::ddeabm_t DF, const int& NEQ, double& T, double* Y, 
	       const double& TOUT, const int* INFO, double& RTOL, 
	       double& ATOL, int& IDID, double* RWORK, const int& LRW, 
	       int* IWORK, const int& LIW, void* RPAR, void* IPAR);

  /****************************************************************
   *     Differential equations solver by Runge-Kutta method
   ****************************************************************/
  void dderkf_(slatec::ddeabm_t df, const int& neq, double& t, double* y, 
	       const double& tout, const int* info, double& rtol, 
	       double& atol, int& idid, double*  rwork, const int& lrw,
	       int* iwork, const int& liw, void* rpar, void* ipar);

  /*****************************************************************
   *         machine specific integer constants
   *****************************************************************/
  int i1mach_(const int&);

  /*****************************************************************
   * machine specific double precision constants (see also dlamch)
   *****************************************************************/
  double d1mach_(const int&);

  /*****************************************************************
   * symmetric system of linear equations solver (together with dsifa)
   ****************************************************************/
  void dsisl_(double* a, const int& lda, const int& n, int* kpvt, 
	      double* b);

  /*****************************************************************
   * symmetric system of linear equations solver (with dsisl)
   *****************************************************************/
  void dsifa_(double* a, const int& lda, const int& n, int* kpvt, 
	      int& info);

  /*****************************************************************
   * general system of linear equations solver (with dsisl)
   *****************************************************************/
  void dgefs_(double* a,const int& lda,const int& n, double* v, 
	      const int& itask, int& ind, double* work,int* iwork);
  /*****************************************************************
   *                     gaus integrator
   *****************************************************************/
  void dgaus8_(slatec::dgaus8_t fun, const double& a, const double& b, 
	       double& err, double& ans, int& ierr);

  /*******************************************************************
   *                     elliptic integral
   *******************************************************************/
  double drf_(const double& X,const double& Y,const double& Z,int& IERR);

  /*******************************************************************
   *             bessel of first kind, zero order
   ******************************************************************/
  double dbesj0_(const double&);

  /*******************************************************************
   *            bessel of first kind, first order
   ******************************************************************/
  double dbesj1_(const double&);
  
  /*******************************************************************
   *       modified bessel of first kind, zero order
   ******************************************************************/
  double dbesi0_(const double&);

  /*******************************************************************
   *   modified bessel of first kind, zero order (predexponent)
   ******************************************************************/
  double dbsi0e_(const double&);

  /*******************************************************************
   *        modified bessel of third kind, zero order
   ******************************************************************/
  double dbesk0_(const double&);

  /*******************************************************************
   *     modified bessel of third kind, zero order (predexponent)
   ******************************************************************/
  double dbsk0e_(const double&);
}
