#include "Fraction.hxx"
#include <iostream>

using namespace std;
//---------------------------------------------------------
void Fraction::operator+=(const Fraction& f){
    n = (n * f.d) + (f.n * d);
    d = d * f.d;
    normalize();
}
//---------------------------------------------------------
void Fraction::operator*=(const Fraction& f){
    n *= f.n;
    d *= f.d;
    normalize();
}
//---------------------------------------------------------
void Fraction::operator/=(const Fraction& f){
    n *= f.d;
    d *= f.n;
    normalize();
}
//---------------------------------------------------------
double Fraction::decimal() const{
    return double(n)/d;
}
//---------------------------------------------------------
void Fraction::invert(){
    int temp = n;
    n = d;
    d = temp;
}
//---------------------------------------------------------
int Fraction::gcd() const{
  int h;
  int a = n;
  int b = d;
  if (a == 0) return abs(b);
  if (b == 0) return abs(a);
  
  do {
    h = a % b;
    a = b;
    b = h;
  } while(b != 0);
  return abs(a);
}
//---------------------------------------------------------
void Fraction::normalize(){
  int factor = gcd();
  n /= factor;
  d /= factor;
  if ( (n>0) && (d<0) ) {n*=1; d*=1;}

}
//---------------------------------------------------------
void Fraction::print() const {
    cout << n << "/" << d << endl;
}
//---------------------------------------------------------
Fraction operator*(const Fraction& lhs, const Fraction& rhs){
   Fraction t(lhs); //Make use of implicitely present copy-constructor
   t *= rhs;
   return t;
}
//---------------------------------------------------------
Fraction operator+(const Fraction& lhs, const Fraction& rhs){
   Fraction t(lhs); //Make use of implicitely present copy-constructor
   t += rhs;
   return t;
}
//---------------------------------------------------------
ostream& operator<<(std::ostream& stream, const Fraction& rhs){
  stream << rhs.n << "/" << rhs.d;
  return stream;
};
