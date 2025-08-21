#include <iostream>
#include <fstream>
#include <complex>
#include "Timer.hxx"


using namespace std;

int main()
{
 const double max = 2;
 const double min = -2;
 const int N = 2000;       // Number of points along axis
 const double h = (max-min)/(N-1);
 const int iterMax = 256; // Maximum number of iterations

 int* F = new int[N*N];


 Timer<> T;
 T.tick();

 #pragma omp parallel for 
 for(int i=0; i<N; i++){
  for(int j=0; j<N; j++)
   {
    double im = min+i*h;
    double re = min+j*h;
    complex<double> z,  z0 =  complex<double>(re,im);
    int k;
    for(k = 0; k<iterMax; k++)
    {
      z = z0 - (z0*z0*z0 - 1.0)/(3.0*z0*z0);
      if (abs(z-z0) < 1e-6) break;
      z0 = z;
    }
    F[i*N+j] = k;
   }
 }
 T.tock();
 

double tau = T.duration().count();

cout << "tau = " << tau << " ms" << endl;


T.tick();

 ofstream out("out");
 for(int i=0; i<N; i++)
   for(int j=0; j<N; j++)
   {
       out << i << "\t" << j << "\t" << F[i*N+j] << endl;
   }
 
   out.close();
 T.tock();
tau = T.duration().count();
cout << "IO = " << tau << " ms" << endl;

   delete[] F;
   return 0;
}
