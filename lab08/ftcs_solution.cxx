#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u, const double dx, const double xmin,
                const int N);
void step(double* const u1,const double* const u0,
          const double V, const double dt, const double dx, const int N);
//---------------------------------------
int main(){

  const double tEnd = 15;
  const double V = 1;

  const int N  = 256;
  const double xmin = -10;
  const double xmax =  10;
  const double dx = (xmax-xmin)/(N-1);
  double dt = dx/V;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* h;
  stringstream strm;

  initialize(u0,dx, xmin,N);

  writeToFile(u0, "f_0", dx, xmin, N);

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){
      step(u1,u0,V,dt,dx,N);

      h = u0;
      u0 = u1;
      u1 = h;
   }
   strm.str("");
   strm << "f_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }


  delete[] u0;
  delete[] u1;
  return 0;
}
//-----------------------------------------------
void step(double* const u1,const double* const u0,
          const double V,
          const double dt, const double dx, const int N)
{
  const double c = V*dt/(2*dx);
  u1[0]   = -c * (u0[1]-u0[N-1]) + u0[0]; // periodic boundary
  u1[N-1] = -c * (u0[0]-u0[N-1]) + u0[N-1]; // periodic boundary
  for(int i=0; i<N; i++)
  {
     u1[i] = -c * (u0[i+1]-u0[i-1]) + u0[i];
  }

}
//-----------------------------------------------
void initialize(double* const u, const double dx, const double xmin,
                const int N)
{
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     if (fabs(x)<=1.0)
       u[i] = 1;
     else
      u[i] =0;
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}
