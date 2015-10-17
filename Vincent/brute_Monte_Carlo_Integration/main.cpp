/*This programm computes the the integral by applying the brute force Monte Carlo algorithm*/

#include <iostream>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;

double function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return exp(-4.0*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)))* 1/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

int main()
{
    double integralsum = 0;
    double a = 5;
    int n = 10000000;
    int m = 1;
    double **x = new double *[6];
    for(int i = 0; i < 6; i++)
    {
        x[i] = new double [n];
    }
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution(-a,a);

    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < n; j++)
        {
            x[i][j] = distribution(generator);
        }
    }


    /* for(int i = 0; i < n; i++)  //Test
    {
        cout << x[0][i] << endl;
    }
    */




    for(int i = 0; i < n; i++)
    {

        if(x[0][i] != x[3][i] || x[1][i] != x[4][i] || x[2][i] != x[5][i])  //We do not want to divide by 0
        {
            integralsum += function(x[0][i], x[1][i], x[2][i], x[3][i], x[4][i], x[5][i]);
        }
    }


    //integralsum /= n*n*n*n*n*n;

    while(n != 10)
    {
        m++;
        n /= 10;
    }

    cout << "The calculated value for the integral is: " << integralsum  << "* 10^-" << m << endl;

    for(int i = 0;i < 6; i++)
    {
        delete [] x[i];
    }
    delete [] x;

    return 0;
}



