#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void read_file( double*&  p, const int N, const string file_in){
   ifstream in(file_in);
   p = new double[N];
   for(int i=0; i<N; i++){
     in >> p[i];
   }
   in.close();
}

void write_file(const double* const p, const int N, const string file_out){
   ofstream out(file_out);
   for(int i=0; i<N; i++)
      out << p[i] << endl;
   out.close();
}

void filter(double* const p, const int N){
   double temp;
   const double p0   = p[0];
   double prev = p[N-1];
   for(int i=0; i<N-1; i++){
     temp = p[i];
     p[i] += prev + p[i+1];
     p[i] /= 3.0;
     prev = temp;
   }
   p[N-1] += prev + p0;
   p[N-1] /= 3.0;
}

int main(){
  double* p=NULL;
  const int N=237;
  const string file_in = "noisy.txt";
  const string file_out = "filtered.txt";

  cout << "pointer before: " << p << endl;
  read_file(p, N, file_in);
  cout << "pointer after: " << p << endl;
  filter(p,N);
  filter(p,N);
  filter(p,N);

  write_file(p,N,file_out);
  delete[] p;
  return 0;
}
