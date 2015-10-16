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
	double x1 ;
	double x2 ;
	double y1 ;	
	double y2 ;	
	double z1 ;	
	double z2 ;
	double sum1=0,sum2=0;
	double fofx,Variation;
	int factor=4;// from -factor/2 to factor/2
	
	clock_t start, finish;  //  declare start and final time
    start = clock();
	for(int i=0;i<n;i++){
		x1=ran0(&k)*factor-0.5*factor;
		x2=ran0(&k)*factor-0.5*factor;
		y1=ran0(&k)*factor-0.5*factor;
		y2=ran0(&k)*factor-0.5*factor;
		z1=ran0(&k)*factor-0.5*factor;
		z2=ran0(&k)*factor-0.5*factor;
		/*if(x1<100){
			cout << x1 << endl;
		};*/
		fofx=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)))/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
		sum1+=fofx;
		sum2+=fofx*fofx;
	}
	sum1=sum1/n*pow(factor,6);
	Variation= (sum2*pow(factor,6)*pow(factor,6))/n-sum1*sum1;
	
	finish = clock();
	// n // Integral     //standart deviation of the Integral   //  Compilation time
	ofstream Zieldatei1("Data_c.txt", ios::app);
	Zieldatei1 <<endl<<n<<"  "<< sum1<<"  ";
	Zieldatei1 <<sqrt(Variation/n)<<"  ";
	Zieldatei1 <<double ( (finish - start)/(double)CLOCKS_PER_SEC );
	Zieldatei1.close();
	
}
