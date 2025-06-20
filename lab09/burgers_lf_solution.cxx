#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
#ifndef M_PI
const double M_PI = 3.14159265359;
#endif
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u1, double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N);
void analytic_solution(const double dx, const double tend, const double N,
                                const double xmin);
//---------------------------------------
int main(){

  const double tEnd = 1/(2*M_PI) ;


  const int N  = 64;
  const double xmin = 0;
  const double xmax = 1;
  const double dx = (xmax-xmin)/N ;
  double dt = dx/2;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* u2 = new double[N];
  double* h;
  stringstream strm;

  initialize(u1,u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);

  analytic_solution(dx,tEnd,N,xmin);

  cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){
        step(u2,u1,u0,dt,dx,N);

        h = u0;
        u0 = u1;
        u1 = u2;
        u2 = h;

        t += dt;
   }
   strm.str("");
   strm << "lfu_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }

  cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  delete[] u2;
  return 0;
}
//-----------------------------------------------
void analytic_solution(const double dx, const double tend, const double N,
                      const double xmin)
{
  ofstream out("ana");

  for(int i=0; i<N; i++)
  {
    double x = xmin + i*dx;
    double u = sin(2*M_PI*x);
    double xi = x + u*tend ;
    out << xi << "\t" << u << endl;
  }

  out.close();
}
//-----------------------------------------------

void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N)
{
  const double c = dt/dx;
  u2[0]   =   u0[0] - c *   u1[0]*(  u1[1] - u1[N-1]); // periodic boundary
  u2[N-1] = u0[N-1] - c * u1[N-1]*(  u1[0] - u1[N-2]); // periodic boundary
  for(int i=1; i<N-1; i++)
  {
     u2[i] = u0[i] - dt/dx * u1[i]*( u1[i+1] - u1[i-1]);
  }

}
//-----------------------------------------------
void initialize(double* const u1, double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     u   =         sin(2*M_PI*x);
     ux  =  2*M_PI*cos(2*M_PI*x);
     uxx =       - 4*M_PI*M_PI*u;
     u1[i] = sin(2*M_PI*x);
     u0[i] = u1[i] + dt*u*ux + dt*dt*0.5 * (u * u * uxx + 2*u* ux*ux);
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s);
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}
