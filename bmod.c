/* file bmod.c */
#include <R.h>
#include <Rmath.h>
static double parms[35];
#define OUx parms[0]
#define OUy parms[1]
#define Wx parms[2]
#define Wy parms[3]
#define midx parms[4]
#define midy parms[5]
#define xstr parms[6]
#define Q parms[7]
#define r parms[8]
#define S parms[9]
#define a parms[10]
#define m parms[11]
#define K parms[12]
#define w parms[13]
// 12 species - parms[14-25] are xopt, and parms[26-37] are yopt
/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=38;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
GetRNGstate();
  double gix; //Growth rates for each species
  double rgain;
  double rtotal;
  int xoptstart = 14, yoptstart = 26;
  int spstart = 3;
  
  if (ip[0] <14) error("nout should be at least 14");
  //Climate variables
  ydot[0] = (midx-y[0])*OUx+rnorm(0, Wx);
  ydot[1] = (midy-y[1])*OUy+rnorm(0, Wy);
    
  //Species and resources
  rtotal = 0; //Sum growth for resource change
  for(unsigned int i = 0; i < 12; i++) { // 12 species
    gix = r*exp(-0.5*pow((abs(parms[xoptstart+i]-y[0])*(xstr)+abs(parms[yoptstart+i]-y[1])*(1-xstr))/w, 2)); //optimum-dependent growth rate

    //Species increase and decrease
    yout[i] = gix*y[spstart+i];
    rgain = yout[i]*y[2]/(y[2]+K); //Increase in species i
    rtotal = rtotal + yout[i]; //Sum growth for resource change
    ydot[spstart+i] = rgain-m*y[spstart+i]; //Total growth in species i
  }
  ydot[2] = a*(S-y[2])-Q*rtotal*y[2]/(y[2]+K); //Change in resource
  
  yout[12] = rtotal;
  yout[13] = gix;
  
PutRNGstate();
}

/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip) {
const double gix = yout[13];
const double gixR = gix*y[2]/(y[2]+K);
const double gixRm = gix*y[2]/(y[2]+K)-m;
const double QgixR = Q*gix*y[2]/(y[2]+K);
const double devR = K/(K + pow(y[2],2));
unsigned int i;

  // col 1, Climate x
  pd[0] = -OUx;
    //for(i = 1; i < (*nrowpd); i++) {
    //  pd[i] = 0.0;
    //}
    
  // col 2, Climate y
    //pd[(*nrowpd)] = 0.0;
  pd[(*nrowpd)+1] = -OUy;
    //for(i = (*norwpd)+2; i < 2*(*nrowpd); i++) {
      //pd[i] = 0.0;
    //}
    
  // col 3, Resource
    //pd[2*(*nrowpd)] = 0.0;
    //pd[2*(*nrowpd)+1] = 0.0;
  pd[2*(*nrowpd)+2] = a - devR*Q*yout[12];
  for(i = 2*(*nrowpd)+3; i < 3*(*nrowpd); i++) {
    pd[i] = -QgixR;
  }
  
  // col 4+, Species
  for(i = 3; i < 15; i++) {
    pd[i*(*nrowpd)+2] = yout[i-3]*devR; //dSi/dR
    pd[i*(*nrowpd)+i] = gixRm; //dSi/dSi
  }
}
/* END file mymod.c */
