#pragma once

class Dual {
public:
    double val;
    double der;

    Dual(double val, double der = 0.0);

    Dual operator+(const Dual& other) const;
    Dual operator-(const Dual& other) const;
    Dual operator*(const Dual& other) const;
    Dual operator/(const Dual& other) const;

    Dual operator-() const;

    Dual operator+(double other) const;
    Dual operator-(double other) const;
    Dual operator*(double other) const;
    Dual operator/(double other) const;
};

// Outras operações
Dual exp(const Dual& x);
// se nao for especificado a base, é log na base natural, por isso 2 funções diferentes
Dual log(const Dual& x);
Dual log(const Dual& x, double base);
Dual sqrt(const Dual& x);

// Trigonometricas
Dual sin(const Dual& x);
Dual cos(const Dual& x);
Dual tan(const Dual& x);
Dual csc(const Dual& x);
Dual sec(const Dual& x);
Dual ctg(const Dual& x);

// Inversas trigonométricas
Dual asin(const Dual& x);
Dual acos(const Dual& x);
Dual atan(const Dual& x);


// Operadores Dual e floats

Dual operator+(double lhs, const Dual &x);
Dual operator-(double lhs, const Dual &x);
Dual operator*(double lhs, const Dual &x);
Dual operator/(double lhs, const Dual &x);



