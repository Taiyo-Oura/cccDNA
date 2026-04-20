#include <R.h>
/* a trick to keep up with the parameters */
static double parms[7];
#define x0 parms[0]
#define dx parms[1]
#define lambda parms[2]
#define f parms[3]
#define alpha parms[4]
#define rho parms[5]
#define rr parms[6]

/* initializers */
void initparms( void (* odeparms)(int *, double *) ) {
    int N = 7;
    odeparms(&N,parms);
}
/* names for states and derivatives */
#define x var[0]
#define y var[1]

#define dxdt vardot[0]
#define dydt vardot[1]

void derivs( int *neq, double *t, double *var, double *vardot, double *varout, int *ip ) {

    if( ip[0]<1 ) {
        error("nout should be at least 1");
    }
        // #y show the dynamics without Tet
        dxdt = f*rho*y-dx*x;
        dydt = lambda*(1-exp(-rr*(*t)))+alpha*x-rho*y;

        double lamv = lambda*(1-exp(-rr*(*t)));

        varout[0] = lamv;
    
}
