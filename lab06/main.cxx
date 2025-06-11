#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Fraction.hxx"

using namespace std;

int main(){
    Fraction f(2,3);
    cout << "This should be 2/3:" << endl;
    f.print();

    Fraction f2(2,-3);
    cout << "This should be -2/3:" << endl;
    f2.print();


    f.invert();
    cout << "This should be 3/2:" << endl;
    f.print();

    cout << "In decimal 3/2 is :" << endl;
    cout << f.decimal() << endl;
    cout << "------------" << endl;

    Fraction b(12,4);
    cout << "This should be 3/1:" << endl;
    b.print();

    Fraction b2(-4,12);
    cout << "This should be -1/3:" << endl;
    b2.print();

    b2 /= b;
    cout << "This should be -1/9:" << endl;
    b2.print();
    cout << "------------" << endl;



    Fraction res(1,1);
    Fraction a(1,2);
    Fraction c(1,4);
    res = a * c;
    cout << "This should be 1/8:" << endl;
    res.print();

    res = a + c;
    cout << "This should be 3/4:" << endl;
    res.print();

    cout << "------------" << endl;
    cout << "This is a test of the ostream << operator: " << res << endl;


    return 0;
}
