#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;


int main(){
    const double tmax = 20;
    const double Nfe = 500;
    const double dt = tmax/Nfe;
    const double q0 = 3;
    const double p0 = 0;
    double q_n,q_old,p_n,H,q_halb;
    
    //Forward Euler
    q_n = q0;
    p_n = p0;
    ofstream ofdataFE("dataFE");
    
    for(double i = 0; i < tmax; i+=dt){
        q_old = q_n;
        q_n += dt*(p_n/sqrt(1+p_n*p_n));
        p_n -= dt*q_old;
        H = sqrt(1+p_n*p_n)+0.5*q_n*q_n;
        ofdataFE << i << "\t" << q_n << "\t" << p_n << "\t" << H << endl;
    }
    ofdataFE.close();
    
    //Sympletic Euler
    q_n = q0;
    p_n = p0;    
    ofstream ofdataSE("dataSE");    
    
    for(double i = 0; i < tmax; i+=dt){
        q_n += dt*(p_n/sqrt(1+p_n*p_n));
        p_n -= dt*q_n;
        H = sqrt(1+p_n*p_n)+0.5*q_n*q_n;
        ofdataSE << i << "\t" << q_n << "\t" << p_n << "\t" << H << endl;
    }
    ofdataSE.close();
    
    
    //Störmer-Verlet
    q_n = q0;
    p_n = p0;
    ofstream ofdataSV("dataSV");
    
    for(double i = 0; i < tmax; i+=dt){
        q_halb = q_n + dt/2 * (p_n/sqrt(1+p_n*p_n));
        p_n -= dt*q_halb;
        q_n = q_halb + dt/2 * (p_n/sqrt(1+p_n*p_n));
        H = sqrt(1+p_n*p_n)+0.5*q_n*q_n;
        ofdataSV << i << "\t" << q_n << "\t" << p_n << "\t" << H << endl;
    }
    ofdataSV.close();
    
    return 0;
}