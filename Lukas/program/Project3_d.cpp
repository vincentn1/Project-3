#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <time.h>

using namespace std;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
float ran0(long *idum)
/*“Minimal” random number generator of Park and Miller. Returns a uniform random deviate
between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)to initialize the sequence; idum must not be altered between calls for successive deviates in a sequence./**/
{
long k;
float ans;
*idum ^= MASK; //XORing with MASK allows use of zero and other simple bit patterns 
k=(*idum)/IQ; //idum.
*idum=IA*(*idum-k*IQ)-IR*k; //Compute idum=(IA*idum) % IM without overif
if(*idum < 0) *idum += IM;// flows by Schrage’s method.
ans=AM*(*idum);// Convert idum to a floating result.
*idum ^= MASK; //Unmask before return.
return ans;
}

int main(){
	long k=pow((time(0)%100),4.3);
	/*for(int i=0;i<20;i++){
		cout << ran0(&k)<<"  "<<k<<endl;
	}*/ 
	
	int n;
	cin >> n;
	double r1 ;
	double r2 ;
	double theta1 ;	
	double theta2 ;	
	double phi1 ;	
	double phi2 ;
	double sum1=0,sum2=0,fofx;
	double Variation;
	
	clock_t start, finish;  //  declare start and final time
    start = clock();
    
	for(int i=0;i<n;i++){
		theta1=ran0(&k)*3.14157;
		theta2=ran0(&k)*3.14157;
		phi1=ran0(&k)*3.14157*2;
		phi2=ran0(&k)*3.14157*2;
		r1=(-1.0/4)*log(1-ran0(&k));
		r2=(-1.0/4)*log(1-ran0(&k));
		fofx=r1*r1*r2*r2*sin(theta1)*sin(theta2)/(sqrt(r1*r1+r2*r2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2))));
		sum1+=fofx/n;
		sum2+=fofx*fofx;
	}
	finish = clock();
	
	sum1=sum1*(pow(M_PI,4)/4); // Normalization and Volume
	sum2=sum2*(pow(M_PI,4)/4)*(pow(M_PI,4)/4); // Volume 
	Variation= sum2/n-sum1*sum1;

	// n // Integral     //standart deviation of the Integral   //  Compilation time
	ofstream Zieldatei1("Data_d.txt", ios::app);
	Zieldatei1 <<endl<<n<<"  "<< sum1<<"  ";
	Zieldatei1 <<sqrt(Variation/n)<<"  ";
	Zieldatei1 <<double ( (finish - start)/(double)CLOCKS_PER_SEC );
	Zieldatei1.close();
	
}
