#include <iostream>
#include <math.h>

#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
#define   Pi         3.141592653589793238462643383

using namespace std;

void gauleg(double x1, double x2, double x[], double w[], int n);
void gauss_laguerre(double *, double *, int, double);
double gammln(double);
double function(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    return  1/sqrt(r1*r1+r2*r2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)));
}

int main()
{
    int n = 25;                         //Number of gridpoints
    double *xr = new double [n+1];
    double *xphi = new double [n];
    double *xtheta = new double [n];
    double *wr = new double [n+1];
    double *wphi = new double [n];
    double *wtheta = new double [n];
    double integralsum = 0;
    double alpha = 4;

    gauleg(0, 2*Pi, xphi, wphi, n);
    gauleg(0, Pi, xtheta, wtheta, n);
    gauss_laguerre(xr, wr, n, alpha);


    cout << "calculating...";
    for(int f = 0; f < n; f++)
    {
        for(int g = 0; g < n; g++)
        {
            for(int h = 0; h < n; h++)
            {
                for(int i = 0; i < n; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        for(int k = 0; k < n; k++)
                        {
                            if(xr[f+1]*xr[f+1]+xr[g+1]*xr[g+1] != 2*xr[f+1]*xr[g+1]*(cos(xtheta[h])*cos(xtheta[i])+sin(xtheta[h])*sin(xtheta[i])*cos(xphi[j]-xphi[k])) )  //We do not want to divide by 0
                            {
                                integralsum += wr[f+1]*wr[g+1]*wtheta[h]*wtheta[i]*wphi[j]*wphi[k]*function(xr[f+1], xr[g+1], xtheta[h], xtheta[i], xphi[j], xphi[k]);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << " finished!" << endl << endl;

    cout << "The calculated value for the integral is: " << integralsum << endl << endl
         << "The value of the integral calculated analyticly is: " << (5* Pi*Pi)/(16*16) << endl;



    delete [] xr;
    delete [] xphi;
    delete [] xtheta;
    delete [] wr;
    delete [] wphi;
    delete [] wtheta;

    return 0;
}



void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
// end function gammln

void gauleg(double x1, double x2, double x[], double w[], int n)
{
int         m,j,i;
double      z1,z,xm,xl,pp,p3,p2,p1;
double      const  pi = 3.14159265359;
double      *x_low, *x_high, *w_low, *w_high;

m  = (n + 1)/2;                             // roots are symmetric in the interval
xm = 0.5 * (x2 + x1);
xl = 0.5 * (x2 - x1);

x_low  = x;                                       // pointer initialization
x_high = x + n - 1;
w_low  = w;
w_high = w + n - 1;

for(i = 1; i <= m; i++) {                             // loops over desired roots
z = cos(pi * (i - 0.25)/(n + 0.5));

    /*
** Starting with the above approximation to the ith root
    ** we enter the mani loop of refinement bt Newtons method.
    */

do {
  p1 =1.0;
p2 =0.0;

/*
** loop up recurrence relation to get the
    ** Legendre polynomial evaluated at x
    */

for(j = 1; j <= n; j++) {
 p3 = p2;
 p2 = p1;
 p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
}

/*
** p1 is now the desired Legrendre polynomial. Next compute
    ** ppp its derivative by standard relation involving also p2,
    ** polynomial of one lower order.
    */

pp = n * (z * p1 - p2)/(z * z - 1.0);
z1 = z;
z  = z1 - p1/pp;                   // Newton's method
} while(fabs(z - z1) > ZERO);

   /*
** Scale the root to the desired interval and put in its symmetric
   ** counterpart. Compute the weight and its symmetric counterpart
   */

*(x_low++)  = xm - xl * z;
*(x_high--) = xm + xl * z;
*w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
*(w_high--) = *(w_low++);
}
} // End_ function gauleg()

