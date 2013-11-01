/* file bmod.c */
#include <R.h>
#include <Rmath.h>
static double parms[35];
static double forc[2];

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
#define temp forc[0]
#define precip forc[1]

/* initializer */
void parmsc(void (* odeparms)(int *, double *)) {
  int N=38;
  odeparms(&N, parms);
}

void forcc(void (* odeforcs)(int *, double *)) {
  int N=2;
  odeforcs(&N, forc);
}

/* Derivatives and output variables */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
GetRNGstate();
  double gix; //Growth rates for each species
  double rgain;
  double rtotal;
  const int xoptstart = 14, yoptstart = 26;
  const int spstart = 1;
  const int respos = 0;
  
  if (ip[0] <13) error("nout should be at least 13");

  //Species and resources
  rtotal = 0; //Sum growth for resource change
  //fprintf(stderr, "temp = %.1f\n", forc[0]);
  //fprintf(stderr, "precip = %.1f\n", forc[1]);
  for(unsigned int i = 0; i < 12; i++) { // 12 species
    if(y[spstart+i] > 0.0001) {
    gix = r*exp(-0.5*pow((abs(parms[xoptstart+i]-temp)*(xstr)+abs(parms[yoptstart+i]-precip)*(1-xstr))/w, 2)); //optimum-dependent growth rate

    //Species increase and decrease
    yout[i] = gix*y[spstart+i];
    rgain = yout[i]*y[respos]/(y[respos]+K); //Increase in species i
    rtotal = rtotal + yout[i]; //Sum growth for resource change
    ydot[spstart+i] = rgain-m*y[spstart+i]; //Total growth in species i
    } else {
      y[spstart+i] = 0;
      gix = 0.0; //optimum-dependent growth rate
      //Species increase and decrease
      yout[i] = 0;
      rgain = 0.0; //Increase in species i
      ydot[spstart+i] = 0.0; //Total growth in species i
    }
    
  }
  ydot[respos] = a*(S-y[respos])-Q*rtotal*y[respos]/(y[respos]+K); //Change in resource
  
  yout[12] = rtotal;
PutRNGstate();
}

/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip) {
unsigned int i;
const int spstart = 1;
const int respos = 0;

 // col 0, Resource
  pd[0] = (Q*y[respos]*yout[12])/pow(K + y[respos],2) - (Q*yout[12])/(K + y[respos]) - a;
  unsigned int n = 0;
  for(i = 1; i < (*nrowpd); i++) {
    pd[i] = (Q*yout[n])/(K + y[respos]) - (Q*y[respos]*yout[n])/pow(K + y[respos],2);
    n = n+1;
  }
  // col 1+, Species
  for(i = 0; i < 12; i++) {
    pd[(i+spstart)*(*nrowpd)] = (K*yout[i])/pow(K + y[respos],2); //dSi/dR
    pd[(i+spstart)*(*nrowpd)+i+spstart] = (y[respos]*yout[i]/y[spstart+i])/(K + y[respos]) - m; //dSi/dSi
  }
}
/* END file mymod.c */
