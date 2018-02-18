SEXP evalf_r,evalg_r,evalh_r,evalc_r,evaljac_r,evalhc_r,evalfc_r,evalgjac_r,
  evalhl_r,evalhlp_r,inip_r,endp_r,param_r,environment_r;

void evalf(int n,double *x,double *f,int *flag);

void evalg(int n,double *x,double *g,int *flag);

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
	   int *flag);

void evalc(int n,double *x,int ind,double *c,int *flag);

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,int *jcnnz,
	     int *flag);

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
	    int *hcnnz,int *flag);

void evalfc(int n,double *x,double *f,int m,double *constr,int *flag);

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
	      double *jcval,int *jcnnz,int *flag);

void evalhl(int n,double *x,int m,double *lambda,double sf,
            double *sc,int *hllin,int *hlcol,double *hlval,int *hlnnz,
            int *flag);

void evalhlp(int n,double *x,int m,double *lambda,double sf,
	     double *sc,double *p,double *hp,int *gothl,int *flag);

void inip(int* n,double** x,double** l, double** u,int* m,
	  double** lambda,int** equatn,int** linear,int* coded,
	  int* checkder);

void endp(int n,double* x,double* l,double* u,int m,double* lambda,
     int* equatn,int* linear);

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp);

FCALLSCSUB4(evalf,EVALF,evalf,INT,DOUBLEV,PDOUBLE,PINT)

FCALLSCSUB4(evalg,EVALG,evalg,INT,DOUBLEV,DOUBLEV,PINT)

FCALLSCSUB7(evalh,EVALH,evalh,INT,DOUBLEV,INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB5(evalc,EVALC,evalc,INT,DOUBLEV,INT,PDOUBLE,PINT)

FCALLSCSUB7(evaljac,EVALJAC,evaljac,INT,DOUBLEV,INT,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB8(evalhc,EVALHC,evalhc,INT,DOUBLEV,INT,INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB6(evalfc,EVALFC,evalfc,INT,DOUBLEV,PDOUBLE,INT,DOUBLEV,PINT)

FCALLSCSUB9(evalgjac,EVALGJAC,evalgjac,INT,DOUBLEV,DOUBLEV,INT,INTV,INTV,\
DOUBLEV,PINT,PINT)

FCALLSCSUB11(evalhl,EVALHL,evalhl,INT,DOUBLEV,INT,DOUBLEV,DOUBLE,DOUBLEV, \
INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB10(evalhlp,EVALHLP,evalhlp,INT,DOUBLEV,INT,DOUBLEV,DOUBLE,DOUBLEV,\
DOUBLEV,DOUBLEV,PLOGICAL,PINT)

PROTOCCALLSFSUB19(ALGENCAN,algencan,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,\
DOUBLEV,DOUBLEV,INT,DOUBLEV,LOGICALV,LOGICALV,LOGICALV,LOGICAL,PDOUBLE,\
PDOUBLE,PDOUBLE,PDOUBLE,PINT)

#define Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,\
equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform) \
CCALLSFSUB19(ALGENCAN,algencan,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,DOUBLEV,\
DOUBLEV,INT,DOUBLEV,LOGICALV,LOGICALV,LOGICALV,LOGICAL,PDOUBLE,\
PDOUBLE,PDOUBLE,PDOUBLE,PINT,\
epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,linear,coded,checkder,\
f,cnorm,snorm,nlpsupn,inform)

SEXP createRRealScalar(double x);

SEXP createRRealVector(int size, double* x);

SEXP createRIntScalar(int x);

SEXP createRIntVector(int size, int* x);

SEXP ralgencan(SEXP evalf_ptr,SEXP evalg_ptr,SEXP evalh_ptr,
	       SEXP evalc_ptr,SEXP evaljac_ptr,SEXP evalhc_ptr,SEXP evalfc_ptr,
	       SEXP evalgjac_ptr, SEXP evalhl_ptr, SEXP evalhlp_ptr,
	       SEXP inip_ptr, SEXP endp_ptr,SEXP param_ptr,
	       SEXP environment_ptr);
