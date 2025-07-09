#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <complex>

//---------------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//---------------------------------------
void initialize(cmplx* const psi0, const double eta, const double dx,
                const double dt,const double xmin, const int N);

void Lstep(cmplx* const psi, const double dt, const double dx, const int N);
void NLstep(cmplx* const psi, const double dt, const int N);

void writeToFile(cmplx* const psi, const string s, const double dx,
                 const int Nx, const double xmin);
//---------------------------------------
int main(){

        const int N = 1000;
        const double L = 200;
        const double xmin = -L/2;
        const double Tend = 200;
        const double dx = L / (N - 1);
        const double dt = dx  / 10;
        const int Na = 100;
        int Nk = int(Tend / Na / dt + 0.5);
        double t=0;

        const double eta = 0.2;

        cmplx* psi = new cmplx[N];

        stringstream strm;

        initialize(psi, eta, dx , dt,xmin, N);

        writeToFile(psi, "psi_0", dx, N, xmin);

        cout << "Nk = " << Nk << endl;

        for(int i=1; i<=Na; i++)
        {
            Lstep(psi , dt*0.5, dx,N);
            NLstep(psi, dt,N);

            for(int j=1; j<=Nk-2; j++){
               Lstep(psi, dt, dx,N);
               NLstep(psi, dt,N);

                t +=dt;
            }

            Lstep(psi, dt*0.5, dx, N);


            strm.str("");
            strm << "psi_" << i;
            writeToFile(psi, strm.str(), dx, N, xmin);
        }

        cout << "t = " << t << endl;

        delete[] psi;
        return 0;
}
//-----------------------------------------------
void initialize(cmplx* const psi0, const double eta, const double dx,
                const double dt,const double xmin, const int N)
{
  const double x0 = 0;
	const double psi = sqrt(2) * eta;

	for(int i=0; i<N; i++){
		double x = i*dx + xmin;
    psi0[i] =  psi/cosh( eta * x);
    //psi0[i] =  3 * psi/cosh( eta * x);
	}
}
//-----------------------------------------------
void Lstep(cmplx* const psi, const double dt, const double dx, const int N)
{
  cmplx* d = new cmplx[N];
  const cmplx alpha = cmplx(0.0, dt/dx/dx);

  const cmplx alpha2 = alpha*alpha;

  for(int i=0;i<N;i++) d[i] = 1.0 + 2.0*alpha;

  for(int i=1;i<N;i++){
    d[i]  -=  alpha2/d[i-1];
    psi[i] -= - alpha/d[i-1]*psi[i-1];
  }

  psi[N-1] = psi[N-1]/d[N-1];
  for(int i=N-2;i>0; i--)
    psi[i] = (psi[i] + alpha*psi[i+1])/d[i];


}
//-----------------------------------------------
void NLstep(cmplx* const psi,const double dt, const int N){
    for(int i=0;i<N;i++){
        double phi = norm(psi[i])*dt;
        psi[i] *= cmplx(cos(phi), sin(phi));
    }
}
//-----------------------------------
void writeToFile(cmplx* const psi, const string s, const double dx,
                 const int Nx, const double xmin)
{
  ofstream out(s.c_str());

	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
    out << x << "\t" << norm(psi[i])
             << "\t" << psi[i].real()
             << "\t" << psi[i].imag() <<  endl;
	}
	out.close();
}
