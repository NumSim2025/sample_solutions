#include <random>
#include <iostream>
#include <cmath>

using namespace std;

double mean(const double* const sample, const int N){
    double m = 0;
    for(int i = 0; i < N; i++)
        m += sample[i];
    m /= N;    
    return m;
}

double stdev(const double* const sample, const double m, const int N){
    double s;
    for(int i=0; i<N; i++)
        s += (sample[i] - m) * (sample[i] - m);
    s /= (N-1);
    return sqrt(s);
}

int main(int, char**)
{
    
    random_device rd; 
    mt19937 gen(rd()); 

    double mu = 0;
    double sigma = 1;
    normal_distribution<double> d(mu, sigma); 

    int i;

    const int N = 50;  // Number of measurements per experiment
    const int M = 10000; // Number of experiments with N measurements
    double s_mean;
    double s[M];

    for(int n = 2; n <= N; n++){
        double* sample = new double[n];

        for(int k=0; k<M; k++ ){
            for(i = 0; i < n; i++)
                sample[i] = d(gen); 

            s[k] = stdev(sample, mean(sample,n), n);
        }
    
        s_mean = mean(s, M);
        cout << n << "\t" << s_mean << endl;
        delete[] sample;
    }
       
    return 0;
}
