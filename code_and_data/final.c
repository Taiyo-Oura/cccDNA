#include <R.h>

/* parameter vector */
static double parms[10];

#define x0      parms[0]
#define y0      parms[1]
#define dx      parms[2]
#define epsilon parms[3]
#define lambda  parms[4]
#define f       parms[5]
#define alpha   parms[6]
#define rho     parms[7]
#define rr      parms[8]
#define tme30   parms[9]

void initparms(void (* odeparms)(int *, double *)) {
    int N = 10;
    odeparms(&N, parms);
}

/* states for fitting conditions (1-4) */
#define x   var[0]
#define y   var[1]
#define x1s var[2]
#define y1s var[3]
#define x2s var[4]
#define y2s var[5]
#define x3s var[6]
#define y3s var[7]

/* derivatives for fitting conditions (1-4) */
#define dxdt  vardot[0]
#define dydt  vardot[1]
#define dx1dt vardot[2]
#define dy1dt vardot[3]
#define dx2dt vardot[4]
#define dy2dt vardot[5]
#define dx3dt vardot[6]
#define dy3dt vardot[7]

/* Conditions 1-4 for fitting */
void derivs(int *neq, double *t, double *var, double *vardot, double *varout, int *ip) {
    if (ip[0] < 1) {
        error("nout should be at least 1");
    }

    /* condition1: untreated */
    dxdt = f * rho * y - dx * x;
    dydt = lambda * (1.0 - exp(-rr * (*t))) + alpha * x - rho * y;

    /* condition2: ETV from day 0 */
    dx1dt = f * rho * y1s - dx * x1s;
    dy1dt = (1.0 - epsilon) * (lambda * (1.0 - exp(-rr * (*t))) + alpha * x1s) - rho * y1s;

    /* condition3: ETV from day 30 */
    dx2dt = f * rho * y2s - dx * x2s;
    if ((*t) < tme30) {
        dy2dt = lambda * (1.0 - exp(-rr * (*t))) + alpha * x2s - rho * y2s;
    } else {
        dy2dt = (1.0 - epsilon) * (lambda * (1.0 - exp(-rr * (*t))) + alpha * x2s) - rho * y2s;
    }

    /* condition4: Tet from day 30 */
    dx3dt = f * rho * y3s - dx * x3s;
    if ((*t) < tme30) {
        dy3dt = lambda * (1.0 - exp(-rr * (*t))) + alpha * x3s - rho * y3s;
    } else {
        dy3dt = alpha * x3s - rho * y3s;
    }
}

/* states for validation conditions (5-8) */
#define xv5  var[0]
#define yv5  var[1]
#define xv6  var[2]
#define yv6  var[3]
#define xv7  var[4]
#define yv7  var[5]
#define xv8  var[6]
#define yv8  var[7]

/* derivatives for validation conditions (5-8) */
#define dxv5dt vardot[0]
#define dyv5dt vardot[1]
#define dxv6dt vardot[2]
#define dyv6dt vardot[3]
#define dxv7dt vardot[4]
#define dyv7dt vardot[5]
#define dxv8dt vardot[6]
#define dyv8dt vardot[7]

/* Conditions 5-8 for validation */
void derivs1(int *neq, double *t, double *var, double *vardot, double *varout, int *ip) {
    if (ip[0] < 1) {
        error("nout should be at least 1");
    }

    /* condition5: ETV from day 0, Tet from day 30 */
    dxv5dt = f * rho * yv5 - dx * xv5;
    if ((*t) < tme30) {
        dyv5dt = (1.0 - epsilon) * (lambda * (1.0 - exp(-rr * (*t))) + alpha * xv5) - rho * yv5;
    } else {
        dyv5dt = (1.0 - epsilon) * (alpha * xv5) - rho * yv5;
    }

    /* condition6: ETV from day 30, Tet from day 30 */
    dxv6dt = f * rho * yv6 - dx * xv6;
    if ((*t) < tme30) {
        dyv6dt = lambda * (1.0 - exp(-rr * (*t))) + alpha * xv6 - rho * yv6;
    } else {
        dyv6dt = (1.0 - epsilon) * (alpha * xv6) - rho * yv6;
    }

    /* condition7: ETV from day 30, Tet from day 45 */
    dxv7dt = f * rho * yv7 - dx * xv7;
    if ((*t) < tme30) {
        dyv7dt = lambda * (1.0 - exp(-rr * (*t))) + alpha * xv7 - rho * yv7;
    } else if ((*t) < 45.0) {
        dyv7dt = (1.0 - epsilon) * (lambda * (1.0 - exp(-rr * (*t))) + alpha * xv7) - rho * yv7;
    } else {
        dyv7dt = (1.0 - epsilon) * (alpha * xv7) - rho * yv7;
    }

    /* condition8: ETV from day 45, Tet from day 30 */
    dxv8dt = f * rho * yv8 - dx * xv8;
    if ((*t) < tme30) {
        dyv8dt = lambda * (1.0 - exp(-rr * (*t))) + alpha * xv8 - rho * yv8;
    } else if ((*t) < 45.0) {
        dyv8dt = alpha * xv8 - rho * yv8;
    } else {
        dyv8dt = (1.0 - epsilon) * (alpha * xv8) - rho * yv8;
    }
}
