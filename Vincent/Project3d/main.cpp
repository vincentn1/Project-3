/*This programm computes the the integral by applying the brute force Monte Carlo algorithm*/

#include <iostream>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;

double function(double r1, double theta1, double phi1, double r2, double theta2, double phi2)
{
    return exp(-4*(r1+r2))* 1/sqrt(r1*r1+r2*r2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)));
}

int main()
{
    double integralsum = 0;
    double a = 5;
    int n = 30;
    double **x = new double *[6];
    for(int i = 0; i < 6; i++)
    {
        x[i] = new double [n];
    }
    srand(time(NULL));
    default_random_engine generator;
    uniform_real_distribution<double> distribution(-a,a);
    srand(time(NULL));

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
                                integralsum += function(x[0][f], x[1][g], x[2][h], x[3][i], x[4][j], x[5][k]);

                            }
                        }
                    }
                }
            }
        }
    }

    integralsum /= n*n*n*n*n*n;

    cout << "The calculated value for the integral is: " << integralsum << endl;

    for(int i = 0;i < 6; i++)
    {
        delete [] x[i];
    }
    delete [] x;

    return 0;
}

