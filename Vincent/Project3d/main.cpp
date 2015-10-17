/*This programm computes the the integral by applying the brute force Monte Carlo algorithm*/

#include <iostream>
#include <math.h>
#include <time.h>
#include <random>

#define Pi 3.14159265359

using namespace std;

double function(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    return 1.0/256.0*log(1-4.0*r1)*log(1-4.0*r1)*log(1-4.0*r2)*log(1-4.0*r2)*sin(theta1)*sin(theta2)*1/sqrt(log(1-4.0*r1)*log(1-4.0*r1)+log(1-4.0*r2)*log(1-4.0*r2)-2.0*log(1-4.0*r1)*log(1-4.0*r2)*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)));
}

int main()
{
    double variance, summand = 0, value, standarddeviation;
    double integralsum = 0;
    int n = 10000000;
    double **x = new double *[6];
    for(int i = 0; i < 6; i++)
    {
        x[i] = new double [n];
    }
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution1(0,2*Pi), distribution2(0,Pi), distribution3(0,0.25);

    //random numbers get generated
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < n; j++)
        {
            x[i][j] = distribution3(generator);
        }
    }
    for(int i = 2; i < 4; i++)
    {
        for(int j = 0; j < n; j++)
        {
            x[i][j] = distribution2(generator);
        }
    }
    for(int i = 4; i < 6; i++)
    {
        for(int j = 0; j < n; j++)
        {
            x[i][j] = distribution1(generator);
        }
    }

     /*for(int i = 0; i < n; i++)  //Test
    {
        cout << x[0][i] << endl;
        //cout << x[2][i] << endl;
    }*/



    //expectation value gets calculated
    for(int i = 0; i < n; i++)
    {

        if(x[0][i] != x[1][i] || x[2][i] != x[3][i] || x[4][i] != x[5][i])  //We do not want to divide by 0
        {
            value = function(x[0][i], x[1][i], x[2][i], x[3][i], x[4][i], x[5][i]);
            integralsum += value;
            //summand += value*value;   // <- wrong approach
        }
    }


    integralsum /= n;
 /* variance = summand - integralsum*integralsum;
    standarddeviation = sqrt(variance)/sqrt((double)n);*/   //in work

    cout << "The calculated value for the integral is: " << integralsum /* << " +/- " << standarddeviation  */ << endl;

    for(int i = 0;i < 6; i++)
    {
        delete [] x[i];
    }
    delete [] x;

    return 0;
}

