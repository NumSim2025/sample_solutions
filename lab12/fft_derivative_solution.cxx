#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
//-------------------------------
using namespace std;
//-------------------------------

void writeData(const double* const s, const int N, const double dx, 
               const double xmin, const string fname);

void initialize(double* const s, const int N, const double dx, const double xmin);
//-------------------------------
int main(int argc, char** argv){

	const int N = 128;
    const double xmin = -10;
    const double xmax = 10;
    const double L = xmax - xmin;
    const double dx = L/N;
    const double dk = 2*M_PI/L;
	
	// Allocate memory
	fftw_complex* f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
	double* signal  = (double*) malloc(sizeof(double)*N);
    // Equivalent
    //double* signal  = new double[N];

  	// Create plan
	fftw_plan FW  = fftw_plan_dft_r2c_1d(N, signal, f, FFTW_ESTIMATE);
    fftw_plan BW  = fftw_plan_dft_c2r_1d(N, f, signal, FFTW_ESTIMATE);

    initialize(signal, N, dx, xmin);
	
    writeData(signal, N, dx, xmin, "f");

  	// Calculate FFT
  	fftw_execute(FW);

  	// Calculate derivative
    for(int i=0; i<N/2+1; i++){
        double k = i*dk;
        double re = f[i][0];
        f[i][0] = -f[i][1]*k /N;
        f[i][1] =  re*k /N;
    }


    // Calculate FFT
  	fftw_execute(BW);

    writeData(signal, N, dx, xmin, "df");

	// Clean up
	fftw_destroy_plan(FW);
    fftw_destroy_plan(BW);
	fftw_free(f);
	free(signal);

	return 0;
}
//-------------------------------
void initialize(double* const s, const int N, const double dx, const double xmin){
    for(int i=0; i<N; i++){
        double x = xmin + i*dx;
        s[i] = exp(-x*x);
    }
}
//-------------------------------
void writeData(const double* const s, const int N, const double dx, const double xmin, const string fname){
	ofstream out(fname);
    for(int i=0; i<N; i++){
        double x = xmin + i*dx;
		out << x << "\t" << s[i] + 2*x*exp(-x*x) <<  endl;
	}
	out.close();
}
//-------------------------------