//
// Created by user on 03/05/2020.
//
#include "solver.hpp"
#include <complex>

using namespace std;
using solver::solve;
using solver::RealVariable;
using solver::ComplexVariable;
using solver::Quadratic_equation;
using solver::Complex_equation;
using solver::Fraction;
using solver::Fraction_complex;

/*
 * Exceptions:
 */

class unsupportedEquationException: public exception
{
    virtual const char* what() const throw()
    {
        return "Equation must be linear or quadratic only!";
    }
} unsupportedEquationException;

class divisionByZeroException: public exception
{
    virtual const char* what() const throw()
    {
        return "Division by zero!";
    }
} divisionByZeroException;

class noSolutionException: public exception
{
    virtual const char* what() const throw()
    {
        return "There is no real solution or no solution at all!";
    }
} noSolutionException;

/*
 * RealVariable class constructors and methods:
 */

RealVariable::RealVariable(){
    this->x = 0.0;
}

void RealVariable::setX(double var){
    this->x = var;
}

double RealVariable::getX(){
    return this->x;
}

/*
 * ComplexVariable class constructors and methods:
 */

ComplexVariable::ComplexVariable(){
    this->x = 0.0+0i;
}

void ComplexVariable::setX(complex<double> var){
    this->x = var;
}

complex<double> ComplexVariable::getX(){
    return this->x;
}

/*
 * Quadratic_equation operators.
 */

Quadratic_equation solver::operator+ (RealVariable r, double d){
    Quadratic_equation q(r);
    q.a = 0;
    q.b = 1;
    q.c = d;
    return q;
}

Quadratic_equation solver::operator+ (double d, RealVariable r){
    return solver::operator+(r,d);
}

Quadratic_equation solver::operator+ (Quadratic_equation r, double d){
    Quadratic_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c + d;
    return q;
}

Quadratic_equation solver::operator+ (double d, Quadratic_equation r){
    return solver::operator+(r,d);
}

Quadratic_equation solver::operator+ (Quadratic_equation a, Quadratic_equation b){
    Quadratic_equation q(a.x);
    q.a = a.a + b.a;
    q.b = a.b + b.b;
    q.c = a.c + b.c;
    return q;
}

Quadratic_equation solver::operator+ (Quadratic_equation r, Fraction f){
    Quadratic_equation q = r * f.divider + f.dividend;
    return q;
}

Quadratic_equation solver::operator+ (Fraction f, Quadratic_equation r){
    return solver::operator+(r,f);
}

Quadratic_equation solver::operator+ (Fraction a, Fraction b){
    Quadratic_equation q;
    q = a.dividend*b.divider + b.dividend*a.divider;
    return q;
}

Quadratic_equation operator- (Fraction a, Fraction b){
    Quadratic_equation q;
    q = a.dividend*b.divider - b.dividend*a.divider;
    return q;
}

Quadratic_equation solver::operator- (RealVariable r, double d){
   Quadratic_equation q(r);
   q.a = 0;
   q.b = 1;
   q.c = -d;
   return q;
}

Quadratic_equation solver::operator- (double  d, RealVariable r){
    Quadratic_equation q(r);
    q.a = 0;
    q.b = -1;
    q.c = d;
    return q;
}

Quadratic_equation solver::operator- (Quadratic_equation r, double d){
    Quadratic_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c-d;
    return q;
}

Quadratic_equation solver::operator- (double  d, Quadratic_equation r){
    Quadratic_equation q(r.x);
    q.a = -1*r.a;
    q.b = -1*r.b;
    q.c = -1*r.c+d;
    return q;
}

Quadratic_equation solver::operator- (Quadratic_equation a, Quadratic_equation b){
    Quadratic_equation q(a.x);
    q.a = a.a - b.a;
    q.b = a.b - b.b;
    q.c = a.c - b.c;
    return q;
}

Quadratic_equation solver::operator- (Quadratic_equation r, RealVariable d){
    Quadratic_equation q(r.x);
    q.a = r.a;
    q.b = r.b-1;
    q.c = r.c;
    return q;
}

Quadratic_equation solver::operator- (RealVariable r, Quadratic_equation d){
    Quadratic_equation q(d.x);
    q.a = -1*d.a;
    q.b = 1-d.b;
    q.c = -1*d.c;
    return q;
}

Quadratic_equation solver::operator* (RealVariable r, double d){
    Quadratic_equation q(r);
    q.a = 0;
    q.b = d;
    q.c = 0;
    return q;
}

Quadratic_equation solver::operator* (double  d, RealVariable r){
    return solver::operator*(r,d);
}

Quadratic_equation solver::operator* (double  d, Quadratic_equation r){
    Quadratic_equation q(r.x);
    q.a = r.a*d;
    q.b = r.b*d;
    q.c = r.c*d;
    return q;
}

Quadratic_equation solver::operator* (Quadratic_equation r, double  d){
    return solver::operator*(d,r);
}

Quadratic_equation solver::operator* (RealVariable a, RealVariable b){
    Quadratic_equation q(a);
    q.a = 1;
    q.b = 0;
    q.c = 0;
    return q;
}

Quadratic_equation solver::operator* (Quadratic_equation a, Quadratic_equation b){
    Quadratic_equation q;
    if((a.a > 0 && b.a > 0) || (a.a > 0 && b.b > 0) || (b.a > 0 && a.b > 0)){
        throw unsupportedEquationException;
    }else{
        q = a.x * b.x * a.b * b.b + a.x * b.c + a.c * b.x + a.c * b.c;
    }
    return q;
}

Fraction solver::operator* (Fraction a, Fraction b){
    Fraction f;
    f.x = a.x;
    f.dividend = a.dividend*b.dividend;
    f.divider = a.divider*b.divider;
    return f;
}

Quadratic_equation solver::operator/ (RealVariable r, double d){
    if(d == 0)
        throw divisionByZeroException;
    Quadratic_equation q(r);
    q.a = 0;
    q.b = 1/d;
    q.c = 0;
    return q;
}

Fraction solver::operator/ (double  d, RealVariable r){
    Fraction f(r,0,0,d,0,1,0);
    return f;
}

Fraction solver::operator/ (double  d, Quadratic_equation q){
    if(q.a == 0 && q.b == 0 && q.c == 0)
        throw divisionByZeroException;
    Fraction f(q.x,0,0,d,q.a,q.b,q.c);
    return f;
}

Quadratic_equation solver::operator/ (Quadratic_equation q, double d){
    if(d == 0)
        throw divisionByZeroException;
    Quadratic_equation equation(q.x);
    equation.a = q.a/d;
    equation.b = q.b/d;
    equation.c = q.c/d;
    return equation;
}

Fraction solver::operator/ (Fraction a, Fraction b){
    if(b.divider.a == 0 && b.divider.b == 0 && b.divider.c == 0)
        throw divisionByZeroException;
    Fraction f;
    f.x = a.x;
    f.dividend = a.dividend*b.divider;
    f.divider = a.divider*b.dividend;
    return f;
}

Fraction solver::operator/ (Quadratic_equation a, Quadratic_equation b){
    if(b.a == 0 && b.b == 0 && b.c == 0)
        throw divisionByZeroException;
    Fraction f;
    f.x = a.x;
    f.dividend = a;
    f.divider = b;
    return f;
}

Quadratic_equation solver::operator== (Quadratic_equation q, double d){
    if((q.a == 0 && q.b == 0 && q.c == 0)||(q.a == 0 && q.b == 0))
        throw noSolutionException;
    Quadratic_equation equation = q - d;
    return equation;
}

Quadratic_equation solver::operator== (double  d, Quadratic_equation q){
    return solver::operator==(q,d);
}

Quadratic_equation solver::operator== (Fraction q, double d){
    Quadratic_equation equation = d * q.divider - q.dividend;
    return equation;
}

Quadratic_equation solver::operator== (double  d, Fraction q){
    return solver::operator==(q,d);
}

Quadratic_equation solver::operator== (Quadratic_equation a, Quadratic_equation b){
    return a-b;
}

Quadratic_equation solver::operator== (RealVariable a, double d){
    return a-d;
}

Quadratic_equation solver::operator== (double d, RealVariable a){
    return solver::operator==(a,d);
}

Quadratic_equation solver::operator^ (RealVariable a, double d){
    if(d == 2) {
        Quadratic_equation q(a);
        q.a = 1;
        q.b = 0;
        q.c = 0;
        return q;
    } else{
        throw unsupportedEquationException;
    }
}

/*
 * Real variable equation solver;
 */
const double solver::solve(Quadratic_equation q){
    if(q.a == 0){ // Linear equation
        double ans = -q.c/q.b;
        return ans;
    }else{ // Quadratic equation
        double a = q.a;
        double b = q.b;
        double c = q.c;
        double delta = pow(b,2) - 4*a*c;
        if(delta < 0)
            throw noSolutionException;
        double ans = (-b+sqrt(delta))/(2*a);
        return ans;
    }
}

int sameComplex(complex<double> a, complex<double> b){
    double diff = __complex_abs(a-b);
    if(diff <= EPSILON){
        return 1;
    } else{
        return -1;
    }
}

//Complex equation operators:
Complex_equation solver::operator+ (ComplexVariable r, double d){
    Complex_equation q(r);
    q.a = 0.0+0i;
    q.b = 1.0+0i;
    q.c = d+0i;
    return q;
}

Complex_equation solver::operator+ (double d, ComplexVariable r){
    return solver::operator+(r,d);
}

Complex_equation solver::operator+ (Complex_equation r, double d){
    Complex_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c + d;
    return q;
}

Complex_equation solver::operator+ (double d, Complex_equation r){
    return solver::operator+(r,d);
}

Complex_equation solver::operator+ (Complex_equation r, complex<double> d){
    Complex_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c + d;
    return q;
}

Complex_equation solver::operator+ (complex<double> d, Complex_equation r){
    return solver::operator+(r,d);
}

Complex_equation solver::operator+ (Complex_equation a, Complex_equation b){
    Complex_equation q(a.x);
    q.a = a.a + b.a;
    q.b = a.b + b.b;
    q.c = a.c + b.c;
    return q;
}

Complex_equation solver::operator+ (Complex_equation r, Fraction_complex f){
    Complex_equation q = r * f.divider + f.dividend;
    return q;
}

Complex_equation operator+ (Fraction_complex f, Complex_equation r){
    return solver::operator+(r,f);
}

Complex_equation operator+ (Fraction_complex a, Fraction_complex b){
    Complex_equation q;
    q = a.dividend*b.divider + b.dividend*a.divider;
    return q;
}

Complex_equation solver::operator+ (ComplexVariable f, complex<double> r){
    Complex_equation q(f);
    q.a = 0.0+0i;
    q.b = 1.0+0i;
    q.c = r;
    return q;
}

Complex_equation solver::operator+ (complex<double> r, ComplexVariable f){
    return solver::operator+(f,r);
}

Complex_equation solver::operator- (Complex_equation r, double d){
    Complex_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c-d;
    return q;
}

Complex_equation solver::operator- (double  d, Complex_equation r){
    Complex_equation q(r.x);
    q.a = -1.0*r.a;
    q.b = -1.0*r.b;
    q.c = -1.0*r.c+d;
    return q;
}

Complex_equation solver::operator- (Complex_equation r, complex<double> d){
    Complex_equation q(r.x);
    q.a = r.a;
    q.b = r.b;
    q.c = r.c-d;
    return q;
}
Complex_equation solver::operator- (complex<double>  d, Complex_equation r){
    Complex_equation q(r.x);
    q.a = -1.0*r.a;
    q.b = -1.0*r.b;
    q.c = -1.0*r.c+d;
    return q;
}

Complex_equation solver::operator- (ComplexVariable r, double d){
    Complex_equation q(r);
    q.a = 0.0+0i;
    q.b = 1.0+0i;
    q.c = -d+0i;
    return q;
}

Complex_equation solver::operator- (double  d, ComplexVariable r){
    Complex_equation q(r);
    q.a = 0.0+0i;
    q.b = -1.0+0i;
    q.c = d+0i;
    return q;
}

Complex_equation solver::operator- (Complex_equation a, Complex_equation b){
    Complex_equation q(a.x);
    q.a = a.a - b.a;
    q.b = a.b - b.b;
    q.c = a.c - b.c;
    return q;
}

Complex_equation solver::operator- (Fraction_complex a, Fraction_complex b){
    Complex_equation q;
    q = a.dividend*b.divider - b.dividend*a.divider;
    return q;
}

Complex_equation solver::operator- (Complex_equation a, ComplexVariable b){
    Complex_equation q(a.x);
    q.a = a.a;
    q.b = a.b-1.0;
    q.c = a.c;
    return q;
}

Complex_equation solver::operator- (ComplexVariable a, Complex_equation b){
    Complex_equation q(b.x);
    q.a = -1.0*b.a;
    q.b = 1.0-b.b;
    q.c = -1.0*b.c;
    return q;
}

Complex_equation solver::operator* (Complex_equation a, Complex_equation b){
    Complex_equation q;
    if((sameComplex(a.a,0.0) == -1 && sameComplex(b.a,0.0) == -1)
    || (sameComplex(a.a,0.0) == -1 && sameComplex(b.b,0.0) == -1)
    || (sameComplex(b.a,0.0) == -1 && sameComplex(a.b,0.0) == -1)){
        throw unsupportedEquationException;
    }else{
        q = a.x * b.x * a.b * b.b + a.x * b.c + a.c * b.x + a.c * b.c;
    }
    return q;
}

Complex_equation solver::operator* (Complex_equation r, complex<double>  d){
    return solver::operator*(d,r);
}

Complex_equation solver::operator* (complex<double>  d, Complex_equation r){
    Complex_equation q(r.x);
    q.a = r.a*d;
    q.b = r.b*d;
    q.c = r.c*d;
    return q;
}

Complex_equation solver::operator* (ComplexVariable r, double d){
    Complex_equation q(r);
    q.a = 0.0+0i;
    q.b = d+0i;
    q.c = 0.0+0i;
    return q;
}

Complex_equation solver::operator* (double  d, ComplexVariable r){
    return solver::operator*(r,d);
}

Complex_equation solver::operator* (double  d, Complex_equation r){
    Complex_equation q(r.x);
    q.a = r.a*d;
    q.b = r.b*d;
    q.c = r.c*d;
    return q;
}

Complex_equation solver::operator* (Complex_equation r, double  d){
    return solver::operator*(d,r);
}

Complex_equation solver::operator* (ComplexVariable a, ComplexVariable b){
    Complex_equation q(a);
    q.a = 1.0 + 0i;
    q.b = 0.0 + 0i;
    q.c = 0.0 + 0i;
    return q;
}

Fraction_complex solver::operator* (Fraction_complex a, Fraction_complex b){
    Fraction_complex f;
    f.x = a.x;
    f.dividend = a.dividend*b.dividend;
    f.divider = a.divider*b.divider;
    return f;
}

Complex_equation solver::operator/ (ComplexVariable r, double d){
    if(d == 0)
        throw divisionByZeroException;
    Complex_equation q(r);
    q.a = 0.0+0i;
    q.b = 1/d+0i;
    q.c = 0.0+0i;
    return q;
}

Fraction_complex solver::operator/ (double  d, ComplexVariable r){
    Fraction_complex f(r,0.0+0i,0.0+0i,d+0i,0.0+0i,1.0+0i,0.0+0i);
    return f;
}

Fraction_complex solver::operator/ (double  d, Complex_equation q){
    if(sameComplex(q.a,0.0+0i) == 1 && sameComplex(q.b,0.0+0i) == 1 && sameComplex(q.c,0.0+0i) == 1)
        throw divisionByZeroException;
    Fraction_complex f(q.x,0.0+0i,0.0+0i,d+0i,q.a,q.b,q.c);
    return f;
}

Fraction_complex solver::operator/ (Fraction_complex a, Fraction_complex b){
    if(sameComplex(b.divider.a,0.0+0i) == 1
    && sameComplex(b.divider.b,0.0+0i) == 1
    && sameComplex(b.divider.c,0.0+0i) == 1)
        throw divisionByZeroException;
    Fraction_complex f;
    f.x = a.x;
    f.dividend = a.dividend*b.divider;
    f.divider = a.divider*b.dividend;
    return f;
}

Fraction_complex solver::operator/ (Complex_equation a, Complex_equation b){
    if(sameComplex(b.a,0.0+0i) == 1
    && sameComplex(b.b,0.0+0i) == 1
    && sameComplex(b.c,0.0+0i) == 1)
        throw divisionByZeroException;
    Fraction_complex f;
    f.x = a.x;
    f.dividend = a;
    f.divider = b;
    return f;
}

Complex_equation solver::operator/ (Complex_equation a, double b){
    if(b == 0)
        throw divisionByZeroException;
    Complex_equation equation(a.x);
    equation.a = a.a/b;
    equation.b = a.b/b;
    equation.c = a.c/b;
    return equation;
}

Complex_equation solver::operator^ (ComplexVariable c, double d){
    if(d == 2) {
        Complex_equation e(c);
        e.a = 1.0 + 0i;
        e.b = 0.0 + 0i;
        e.c = 0.0 + 0i;
        return e;
    } else{
        throw unsupportedEquationException;
    }
}

Complex_equation solver::operator== (Complex_equation c, double d){
    Complex_equation equation = c - d;
    return equation;
}

Complex_equation solver::operator== (double d, Complex_equation c){
    return solver::operator==(c,d);
}

Complex_equation solver::operator== (Complex_equation c, complex<double> d){
    Complex_equation equation = c - d;
    return equation;
}

Complex_equation solver::operator== (complex<double> d, Complex_equation c){
    return solver::operator==(c,d);
}

Complex_equation solver::operator== (Fraction_complex q, double d){
    Complex_equation equation = d * q.divider - q.dividend;
    return equation;
}

Complex_equation solver::operator== (double  d, Fraction_complex q){
    return solver::operator==(q,d);
}

Complex_equation solver::operator== (Complex_equation a, Complex_equation b){
    return a-b;
}

/*
 * Complex variable equation solver;
 */
complex<double> solver::solve(Complex_equation q){
    if(sameComplex(q.a,0.0+0i) == 1){ // Linear equation
        complex<double> ans = -q.c/q.b;
        return ans;
    }else{ // Quadratic equation
        complex<double> a = q.a;
        complex<double> b = q.b;
        complex<double> c = q.c;
        complex<double> ans = (-b+sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
        return ans;
    }
}