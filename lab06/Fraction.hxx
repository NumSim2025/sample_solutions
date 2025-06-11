#ifndef fraction_hxx
#define fraction_hxx

#include <ostream>

class Fraction{
public:
    Fraction(const int n, const int d):n(n),d(d){
        normalize();
    }

    ~Fraction(){};
    void operator+=(const Fraction& f);
    void operator*=(const Fraction& f);
    void operator/=(const Fraction& f);

    friend Fraction operator*(const Fraction& lhs, const Fraction& rhs);
    friend Fraction operator+(const Fraction& lhs, const Fraction& rhs);

    friend std::ostream& operator<<(std::ostream& stream, const Fraction& rhs);

    double decimal() const;
    void invert();
    void print() const;
private:
    int n;
    int d;
    void normalize();
    int gcd() const;
};

#endif
