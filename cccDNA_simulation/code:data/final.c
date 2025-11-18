#include <R.h>
/* a trick to keep up with the parameters */
static double parms[15];
#define x0 parms[0]
#define x40 parms[1]
#define y0 parms[2]
#define dx parms[3]
#define epsilon parms[4]
#define lambda parms[5]
#define f parms[6]
#define alpha parms[7]
#define rho parms[8]
#define rr parms[9]
#define rr1 parms[10]
#define tme30 parms[11]
#define tme45 parms[12]
#define tme9 parms[13]
#define tme18 parms[14]

/* initializers */
void initparms( void (* odeparms)(int *, double *) ) {
    int N = 15;
    odeparms(&N,parms);
}
/* names for states and derivatives */
#define x var[0]
#define y var[1]
#define x1s var[2]
#define y1s var[3]
#define x2s var[4]
#define y2s var[5]
#define x3s var[6]
#define y3s var[7]
#define x4s var[8]
#define y4s var[9]
// #define x6s var[10]
// #define y6s var[11]

#define dxdt vardot[0]
#define dydt vardot[1]
#define dx1dt vardot[2]
#define dy1dt vardot[3]
#define dx2dt vardot[4]
#define dy2dt vardot[5]
#define dx3dt vardot[6]
#define dy3dt vardot[7]
#define dx4dt vardot[8]
#define dy4dt vardot[9]
// #define dx6dt vardot[10]
// #define dy6dt vardot[11]

void derivs( int *neq, double *t, double *var, double *vardot, double *varout, int *ip ) {

    if( ip[0]<1 ) {
        error("nout should be at least 1");
    }
        // #y show the dynamics without Tet
        dxdt = f*rho*y-dx*x;
        dydt = lambda*(1-exp(-rr*(*t)))+alpha*x-rho*y;
        // #y1 show the dynamics with ETV without Tet
        dx1dt = f*rho*y1s-dx*x1s;
        dy1dt = (1-epsilon)*(lambda*(1-exp(-rr*(*t)))+alpha*x1s)-rho*y1s;
        // #This is case of add ETV after 30 days
        // #When time<30days, y2 show without Tet
        // #When time>30days, y2 show with ETV without Tet
        dx2dt = f*rho*y2s-dx*x2s;
    if( (*t)<tme30 ) {
        dy2dt = lambda*(1-exp(-rr*(*t))) + alpha*x2s- rho*y2s;
    } else {
        dy2dt = (1-epsilon)*(lambda*(1-exp(-rr*(*t))) + alpha*x2s)-rho*y2s;
    }
        // #This is case of add Tet after 30 days
        // #When time<30days, y3 show without Tet
        // #When time>30days, y3 show with Tet
        dx3dt = f*rho*y3s-dx*x3s;
    if( (*t)<tme30 ) {
        dy3dt = lambda*(1-exp(-rr*(*t))) + alpha*x3s - rho*y3s;
    } else {
        dy3dt = alpha*x3s - rho*y3s;
    }
        // #This is case of add Tet after 9 days
        // #When time<9days, y4 show without Tet
        // #When time>9days, y4 show with Tet
        dx4dt = f*rho*y4s-dx*x4s;
    if( (*t)<tme9 ) {
        dy4dt = lambda*(1-exp(-rr1*(*t))) + alpha*x4s - rho*y4s;
    } else {
        dy4dt = alpha*x4s - rho*y4s;
    }

}

void derivs1( int *neq, double *t, double *var, double *vardot, double *varout, int *ip ) {

    if( ip[0]<1 ) {
        error("nout should be at least 1");
    }
        // #y show the dynamics without Tet
        dxdt = f*rho*y-dx*x;
        dydt = lambda*(1-exp(-rr*(*t)))+alpha*x-rho*y;
        // #y1 show the dynamics with ETV without Tet
        dx1dt = f*rho*y1s-dx*x1s;
        dy1dt = (1-epsilon)*(lambda*(1-exp(-rr*(*t)))+alpha*x1s)-rho*y1s;
        // #This is case of add ETV after 30 days
        // #When time<30days, y2 show without Tet
        // #When time>30days, y2 show with ETV without Tet
        dx2dt = f*rho*y2s-dx*x2s;
    if( (*t)<tme30 ) {
        dy2dt = lambda*(1-exp(-rr*(*t))) + alpha*x2s- rho*y2s;
    } else {
        dy2dt = (1-epsilon)*(lambda*(1-exp(-rr*(*t))) + alpha*x2s)-rho*y2s;
    }
        // #This is case of add Tet after 30 days
        // #When time<30days, y3 show without Tet
        // #When time>30days, y3 show with Tet
        dx3dt = f*rho*y3s-dx*x3s;
    if( (*t)<tme30 ) {
        dy3dt = lambda*(1-exp(-rr*(*t))) + alpha*x3s - rho*y3s;
    } else {
        dy3dt = alpha*x3s - rho*y3s;
    }
        // #This is case of add Tet after 9 days
        // #When time<9days, y4 show without Tet
        // #When time>9days, y4 show with Tet
        dx4dt = f*rho*y4s-dx*x4s;
    if( (*t)<tme9 ) {
        dy4dt = lambda*(1-exp(-rr1*(*t))) + alpha*x4s - rho*y4s;
    } else {
        dy4dt = alpha*x4s - rho*y4s;
    }

}
