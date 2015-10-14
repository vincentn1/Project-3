/* This program calculates the integral of a function with the Gauss-Legendre method.
   The Gauss-Legendre algorithm is copied from the Project3 folder(ComputationalPhysics1 repository on GitHub).
*/

#include <iostream>
#include <math.h>

#define   ZERO       1.0E-10
#define   Pi         3.141592653589793238462643383

using namespace std;

void gauleg(double, double, double *, double *, int);
double function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)))* 1/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

int main()
{
    int n = 35;                         //Number of gridpoints
    double *x = new double [n];
    double *w = new double [n];
    double a, b;                        //Integrationborders
    double integralsum = 0;
    b = 5;
    a = -b;

    gauleg(a, b, x, w, n);

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
                            if(f != i && g != j && h != k)  //We do not want to divide by 0
                            {
                                integralsum += w[f]*w[g]*w[h]*w[i]*w[j]*w[k]*function(x[f], x[g], x[h], x[i], x[j], x[k]);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << " finished!" << endl << endl;

    cout << "The integral calculated with the Gauss-Legendre using " << n << " gridpoints and using as" << b << " and " << a << " borders is: "
         << integralsum << "." << endl
         << "The value of the integral calculated analyticly is: " << (5* Pi*Pi)/(16*16) << endl;




    delete [] x;
    delete [] w;

    return 0;
}


/*
** The function
**              gauleg()
** takes the lower and upper limits of integration x1, x2, calculates
** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
** of length n of the Gauss--Legendre n--point quadrature formulae.
*/

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
