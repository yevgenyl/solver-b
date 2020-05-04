//
// Created by user on 03/05/2020.
//

#ifndef SOLVER_2_SOLVER_H
#define SOLVER_2_SOLVER_H

#pragma once
#include <complex>

#define EPSILON 0.0001

using namespace std;

namespace solver{

    class RealVariable{

    private:
        double x;

    public:
        RealVariable();
        void setX(double var);
        double getX();
    };

    class ComplexVariable{

    private:
        complex<double> x;

    public:
        ComplexVariable();
        void setX(complex<double> var);
        complex<double> getX();
    };

    struct Quadratic_equation{ // a*(x^2) + b*x + c
        RealVariable x;
        double a, b, c;
        Quadratic_equation(){};
        Quadratic_equation(RealVariable r){
            x.setX(0.0);
        }
    };

    struct Complex_equation{
        ComplexVariable x;
        complex<double> a, b, c;
        Complex_equation(){};
        Complex_equation(ComplexVariable c){
            x.setX(0.0);
        }
    };

    struct Fraction{
        RealVariable x;
        Quadratic_equation dividend;
        Quadratic_equation divider;
        Fraction(){};
        Fraction(RealVariable real, double a1, double b1, double c1, double a2, double b2, double c2){
            x = real;
            dividend.x = real;
            dividend.a = a1;
            dividend.b = b1;
            dividend.c = c1;
            divider.x = real;
            divider.a = a2;
            divider.b = b2;
            divider.c = c2;
        }
    };

    struct Fraction_complex{
        ComplexVariable x;
        Complex_equation dividend;
        Complex_equation divider;
        Fraction_complex(){};
        Fraction_complex(ComplexVariable real, complex<double> a1, complex<double> b1, complex<double> c1, complex<double> a2, complex<double> b2, complex<double> c2){
            x = real;
            dividend.x = real;
            dividend.a = a1;
            dividend.b = b1;
            dividend.c = c1;
            divider.x = real;
            divider.a = a2;
            divider.b = b2;
            divider.c = c2;
        }
    };
    // Quadratic equation operators:
    Quadratic_equation operator+ (RealVariable r, double d);
    Quadratic_equation operator+ (double d, RealVariable r);
    Quadratic_equation operator+ (Quadratic_equation r, double d);
    Quadratic_equation operator+ (double d, Quadratic_equation r);
    Quadratic_equation operator+ (Quadratic_equation a, Quadratic_equation b);
    Quadratic_equation operator+ (Quadratic_equation r, Fraction f);
    Quadratic_equation operator+ (Fraction f, Quadratic_equation r);
    Quadratic_equation operator+ (Fraction a, Fraction b);
    Quadratic_equation operator- (RealVariable r, double d);
    Quadratic_equation operator- (double  d, RealVariable r);
    Quadratic_equation operator- (Quadratic_equation r, double d);
    Quadratic_equation operator- (double  d, Quadratic_equation r);
    Quadratic_equation operator- (Quadratic_equation a, Quadratic_equation b);
    Quadratic_equation operator- (Fraction a, Fraction b);
    Quadratic_equation operator- (Quadratic_equation r, RealVariable d);
    Quadratic_equation operator- (RealVariable r, Quadratic_equation d);
    Quadratic_equation operator* (RealVariable r, double d);
    Quadratic_equation operator* (double  d, RealVariable r);
    Quadratic_equation operator* (double  d, Quadratic_equation r);
    Quadratic_equation operator* (Quadratic_equation r, double  d);
    Quadratic_equation operator* (RealVariable a, RealVariable b);
    Quadratic_equation operator* (Quadratic_equation a, Quadratic_equation b);
    Fraction operator* (Fraction a, Fraction b);
    Quadratic_equation operator/ (RealVariable r, double d);
    Fraction operator/ (double  d, RealVariable r);
    Fraction operator/ (double  d, Quadratic_equation q);
    Fraction operator/ (Fraction a, Fraction b);
    Fraction operator/ (Quadratic_equation a, Quadratic_equation b);
    Quadratic_equation operator/ (Quadratic_equation q, double d);
    Quadratic_equation operator^ (RealVariable a, double d);
    Quadratic_equation operator== (Quadratic_equation q, double d);
    Quadratic_equation operator== (double  d, Quadratic_equation q);
    Quadratic_equation operator== (Fraction q, double d);
    Quadratic_equation operator== (double  d, Fraction q);
    Quadratic_equation operator== (Quadratic_equation a, Quadratic_equation b);
    Quadratic_equation operator== (RealVariable a, double d);
    Quadratic_equation operator== (double d, RealVariable a);

    //Complex equation operators:
    Complex_equation operator+ (ComplexVariable r, double d);
    Complex_equation operator+ (double d, ComplexVariable r);
    Complex_equation operator+ (Complex_equation r, double d);
    Complex_equation operator+ (double d, Complex_equation r);
    Complex_equation operator+ (Complex_equation r, complex<double> d);
    Complex_equation operator+ (complex<double> d, Complex_equation r);
    Complex_equation operator+ (Complex_equation a, Complex_equation b);
    Complex_equation operator+ (Complex_equation r, Fraction_complex f);
    Complex_equation operator+ (Fraction_complex f, Complex_equation r);
    Complex_equation operator+ (Fraction_complex a, Fraction_complex b);
    Complex_equation operator+ (ComplexVariable f, complex<double> r);
    Complex_equation operator+ (complex<double> r, ComplexVariable f);
    Complex_equation operator- (Complex_equation r, double d);
    Complex_equation operator- (double  d, Complex_equation r);
    Complex_equation operator- (Complex_equation r, complex<double> d);
    Complex_equation operator- (complex<double>  d, Complex_equation r);
    Complex_equation operator- (ComplexVariable r, double d);
    Complex_equation operator- (double  d, ComplexVariable r);
    Complex_equation operator- (Complex_equation a, Complex_equation b);
    Complex_equation operator- (Fraction_complex a, Fraction_complex b);
    Complex_equation operator- (Complex_equation a, ComplexVariable b);
    Complex_equation operator- (ComplexVariable a, Complex_equation b);
    Complex_equation operator* (ComplexVariable r, double d);
    Complex_equation operator* (double  d, ComplexVariable r);
    Complex_equation operator* (double  d, Complex_equation r);
    Complex_equation operator* (Complex_equation r, double  d);
    Complex_equation operator* (ComplexVariable a, ComplexVariable b);
    Fraction_complex operator* (Fraction_complex a, Fraction_complex b);
    Complex_equation operator* (complex<double>  d, Complex_equation r);
    Complex_equation operator* (Complex_equation r, complex<double>  d);
    Complex_equation operator* (Complex_equation a, Complex_equation b);
    Complex_equation operator/ (ComplexVariable r, double d);
    Fraction_complex operator/ (double  d, ComplexVariable r);
    Fraction_complex operator/ (double  d, Complex_equation q);
    Fraction_complex operator/ (Fraction_complex a, Fraction_complex b);
    Fraction_complex operator/ (Complex_equation a, Complex_equation b);
    Complex_equation operator/ (Complex_equation a, double b);
    Complex_equation operator^ (ComplexVariable c, double d);
    Complex_equation operator== (Complex_equation c, double d);
    Complex_equation operator== (double d, Complex_equation c);
    Complex_equation operator== (Complex_equation c, complex<double> d);
    Complex_equation operator== (complex<double> d, Complex_equation c);
    Complex_equation operator== (Fraction_complex q, double d);
    Complex_equation operator== (double  d, Fraction_complex q);
    Complex_equation operator== (Complex_equation a, Complex_equation b);

    const double solve(Quadratic_equation q);
    complex<double> solve(Complex_equation q);
}

#endif //SOLVER_2_SOLVER_H