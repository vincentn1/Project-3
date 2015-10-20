#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10

using namespace std;

double function_cartesian (double xone, double xtwo, double yone, double ytwo, double zone, double ztwo);
double function_spherical(double, double, double, double, double, double);
void gauleg(double, double, double *, double *, int);
double legendre(int, double, double);
double laguerre_legendre(int);
double gammln(double);
void gauss_laguerre(double *, double *, int, double);
/*
 *
 * Main programm starts HERE
 *
 */
int main()
{
    int num_threads = omp_get_num_procs();
    omp_set_num_threads(num_threads);
    //Initialize borders a & b (double) and number of grid points n (int) by user
    int n;
    double a, b;
    cout << "Please set the number of grid points!" << endl;
    cin >> n;
    cout << "Please set the integration limits" << endl;
    cin >> a;
    cin >> b;
    cout << a << endl << b << endl;
    clock_t start, finish;  //  declare start and final time
    //First, brute force Gauss-Legendre: calculate integral and print
    start = clock();
    double gaulegint=legendre(n, a, b);
    finish = clock();
    cout << "Legendre: " << gaulegint << endl << "time: " << double ( (finish - start)/(double)CLOCKS_PER_SEC ) << endl;
    //Now, mixture of Laguerre (r-part) and Legendre (angles): calculate integral and print
    start = clock();
    double gaulagint=laguerre_legendre(n);
    finish = clock();
    cout << "Laguerre-Legendre: " << gaulagint << endl << "time: " << double ( (finish - start)/(double)CLOCKS_PER_SEC ) << endl;
    return 0;
}
/*
 *
 * END of the main programm
 *
 */

//Calculate the integral by only using Gauss-Legendre
double legendre(int n, double a, double b){
    double * x = new double [n];
    double * w = new double [n];
    gauleg(a, b, x, w, n);
    double integral = 0;
    int i, j, k, l, y, z;
# pragma omp parallel for reduction(+:integral)  private (i, j, k, l, y, z)
    for (i=1; i<= n; i++){
        for (j=1; j<= n; j++){
            for (k=0; k<n; k++){
                for (l=0; l<n; l++){
                    for (y=0; y<n; y++){
                        for (z=0; z<n; z++){
                            integral += (w[i] * w[j] * w[k] * w[l] * w[y]* w[z] * function_cartesian(x[i], x[j], x[k], x[l], x[y], x[z]));
                        }
                    }
                }
            }
        }
    }
    return integral;
    delete x;
    delete w;
}

//Calculating the integral by using mainly Laguerre (for r_1 and r_2) and Legendre (for the angles)
double laguerre_legendre(int n){
    double * xlag = new double[n+1];
    double * wlag = new double[n+1];
    double * xphi = new double[n];
    double * wphi = new double[n];
    double * xtheta = new double[n];
    double * wtheta = new double[n];
    double a, b;
    a=0;
    b=2*M_PI;
    gauleg(a, b, xphi, wphi, n);
    a=0;
    b=M_PI;
    gauleg(a, b, xtheta, wtheta, n);
    //alpha=2 because of r^2 in Jacobian!
    gauss_laguerre(xlag, wlag, n, 2.);
    double integral=0;
    int i, j, k, l, y, z;
# pragma omp parallel for reduction(+:integral)  private (i, j, k, l, y, z)
    for (i=1; i<= n; i++){
        for (j=1; j<= n; j++){
            for (k=0; k<n; k++){
                for (l=0; l<n; l++){
                    for (y=0; y<n; y++){
                        for (z=0; z<n; z++){
                            //integral+=(laguerre weights)*(legendre weights*function value)*(Jacobian)
                            //NOTE: Substituted r_1 and r_2 by r_1'=4*r_1 and r_2'=4*r_2! Same in "function_spherical"! Changes Jacobian!
                            integral+=((wlag[i]*wlag[j])*(wphi[k]*wphi[l]*wtheta[y]*wtheta[z]*function_spherical(xlag[i], xlag[j], xphi[k], xphi[l], xtheta[y], xtheta[z]))*(sin(xtheta[y])*sin(xtheta[z])))/4096;
                        }
                    }
                }
            }
        }
    }
    delete wphi;
    delete wtheta;
    delete wlag;
    delete xphi;
    delete xtheta;
    delete xlag;
    return integral;
}

//This is to calculate the function value at a specific point in carthesian coordinates
double function_cartesian (double xone, double xtwo, double yone, double ytwo, double zone, double ztwo){
    double alpha = 2.0;
    double rone, rtwo;
    rone = sqrt(xone*xone+yone*yone+zone*zone);
    rtwo = sqrt(xtwo*xtwo+ytwo*ytwo+ztwo*ztwo);
    double temp =sqrt(pow((xone-xtwo),2)+pow((yone-ytwo),2)+pow((zone-ztwo),2));
    if(temp < pow(10.,-6.)){
        return 0;
    }else{
        return exp(-2*alpha*(rone+rtwo))/temp;
    }
    return 0;
}

//This is to calculate the function value (=value of the denominator) at a given point in spherical coordinates
//This will be used for the Laguerre version
double function_spherical(double rone, double rtwo, double phione, double phitwo, double thetaone, double thetatwo){
    double beta =(cos(thetaone)*cos(thetatwo))+(sin(thetaone)*sin(thetatwo)*cos(phione - phitwo));
    double root = sqrt(((rone*rone)/16)+((rtwo*rtwo)/16)-((2*rone*rtwo*beta)/16));
    if(root > pow(10.,-6.)){
        return 1/root;
    }else{
        return 0;
    }

}

//This function is used (with Morten's permission) from the example program... It calculates the zeros of Legendre-functions
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

//This function is again used with Morten's permission from his source file. It calculates the weights and integration points for the Laguere-method
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

//Helper for the Gauss-Laguerre function by Morten
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
